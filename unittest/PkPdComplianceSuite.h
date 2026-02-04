/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2026 Swiss Tropical and Public Health Institute
 * Copyright (C) 2005-2015 Liverpool School Of Tropical Medicine
 * Copyright (C) 2020-2026 University of Basel
 * Copyright (C) 2026 The Kids Research Institute Australia
 *
 * OpenMalaria is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

// Unittest for the LSTM drug model

#ifndef Hmod_PkPdComplianceSuite
#define Hmod_PkPdComplianceSuite

#include <cxxtest/TestSuite.h>
#include "PkPd/LSTMModel.h"
#include "Host/WithinHost/Infection/DummyInfection.h"
#include "UnittestUtil.h"
#include "ExtraAsserts.h"
#include <limits>
#include <cstdio>

using std::pair;
using std::multimap;
using namespace OM;
using namespace OM::PkPd;

// Use one of these to switch verbosity on/off:
#define PCS_VERBOSE( x )
// #define PCS_VERBOSE( x ) x

// Tolerances. We require either abs(a/b-1) < REL_TOL or abs(a-b) < ABS_TOL.
// For drug concentrations, negligible concentrations are defined, below which
// the model is allowed to approximate to zero.
// For drug factors, we don't actually need a huge amount of precision in the
// simulator (models can do at least 1e-3 relative precision, except that
// required accuracy of the integration algorithms has been lowered for speed).
#define PKPD_CONC_REL_TOL 1e-5
#define PKPD_FACT_REL_TOL 5e-3
#define PKPD_FACT_ABS_TOL 1e-20

/** Test outcomes from the PK/PD code in OpenMalaria with LSTM's external
 * model. Numbers should agree (up to rounding errors). */
class PkPdComplianceSuite : public CxxTest::TestSuite
{
public:
    PkPdComplianceSuite() :
            m_rng(0, 0), proxy(0)
    {
        bodymass = 50 /*kg*/;
       
        PCS_VERBOSE(cout << "\n[ Unittest Output Legend: \033[35mDrug "
                "Factor\033[0m, \033[36mDrug Concentration\033[0m ]" << endl;)
    }
    
    void setUp () {
        m_rng.seed(0, 721347520444481703);
        UnittestUtil::initTime(1);
        UnittestUtil::PkPdSuiteSetup();
        proxy = new LSTMModel ();
        inf = createDummyInfection(m_rng, 0, InfectionOrigin::Indigenous);
        schedule.clear();
    }
    
    void tearDown () {
        delete inf;
        delete proxy;
        LSTMDrugType::clear();
    }
    
    void assembleTripleDosageSchedule ( double dosage ){
        schedule.clear();
        schedule.insert(make_pair(0, make_pair(0, dosage)));
        schedule.insert(make_pair(1, make_pair(0, dosage)));
        schedule.insert(make_pair(2, make_pair(0, dosage)));
    }
    void assembleHexDosageSchedule ( double dosage ){
        schedule.clear();
        schedule.insert(make_pair(0, make_pair(0, dosage)));
        schedule.insert(make_pair(0, make_pair(0.5, dosage)));
        schedule.insert(make_pair(1, make_pair(0, dosage)));
        schedule.insert(make_pair(1, make_pair(0.5, dosage)));
        schedule.insert(make_pair(2, make_pair(0, dosage)));
        schedule.insert(make_pair(2, make_pair(0.5, dosage)));
    }
    
    void assembleCQDosageSchedule ( double dosage ){
        // only used for CQ, which needs 10, 10, 5 dosages instead of constant dosages
        schedule.clear();
        schedule.insert(make_pair(0, make_pair(0, dosage)));
        schedule.insert(make_pair(1, make_pair(0, dosage)));
        schedule.insert(make_pair(2, make_pair(0, dosage/2)));
    }
    
    struct debugDrugSimulations {
        double factor, f_error, f_rel_error;
    };
    
    string drugDebugOutputHeader(bool hasSecondDrug, string drugName){
        // title extension
        if ( hasSecondDrug ){
            cout << drugName << endl << "----" << endl;
        }
        
        /*
        *   Verbose output is formatted in markdown, so it can be used in a github-wiki.
        */
        string white = "\033[0m";
        string yellow = "\033[33m";
        string green = "\033[32m";
        string red = "\033[31m";
        
        // Header: |day|conc|rel|abs|factor|rel|abs|
        string fm = white + "|%1$=3d|"
            + white + "%2$-12d|%3%%4$+-9d" + white + "|%5%%6$+-12d" + white + "|"
            + white  + "%7$-12d|%8%%9$+-9d" + white + "|%10%%11$+-12d" + white + "|"
            + white + "%12$=4s|";
        fm += white + "\n";
        // head row
        string type = hasSecondDrug ? "type" : "";
        printf(fm.c_str(), "day", "conc", yellow.c_str(), "rel err %", yellow.c_str(), "abs err", "factor", yellow.c_str(), "rel err %", yellow.c_str(), "abs err", type.c_str());
        // cout << format(fm) % "day" % "conc" % yellow % "rel err %" % yellow % "abs err" % "factor" % yellow % "rel err %" % yellow % "abs err" % type;
        // head row separator
        string sfill = "------------";
        printf(fm.c_str(), "---", "------------", white.c_str(), "---------", white.c_str(), sfill.c_str(), sfill.c_str(), white.c_str(), "---------", white.c_str(), sfill.c_str(), "----");
        // cout << format(fm) % "---" % "------------" % white % "---------" % white % sfill % sfill % white % "---------" % white % sfill % "----";
        return fm;
    }
    
    void drugDebugOutputLine(size_t day,
            double factor, double f_abs_error, double f_rel_error,
            double concentration, double c_abs_error, double c_rel_error,
            double conc_abs_tol, string type, string fm)
    {
        const char* green = "\033[32m";
        const char* red = "\033[31m";
        // rel errors are already percentages
        const char* col_CR = abs(c_rel_error) > PKPD_CONC_REL_TOL * 100 ? red : green;
        const char* col_CA = abs(c_abs_error) > conc_abs_tol ? red : green;
        const char* col_FR = abs(f_rel_error) > PKPD_FACT_REL_TOL * 100 ? red : green;
        const char* col_FA = abs(f_abs_error) > PKPD_FACT_ABS_TOL ? red : green;
        printf(fm.c_str(), day, concentration, col_CR, c_rel_error, col_CA, c_abs_error, factor, col_FR, f_rel_error, col_FA, f_abs_error, type.c_str());
        // cout << format(fm) % day % concentration % col_CR % c_rel_error % col_CA % c_abs_error
        //         % factor % col_FR % f_rel_error % col_FA % f_abs_error % type;
    }
    
    void runDrugSimulations (string drugName, string drug2Name,
                          const double drug_conc[], const double drug2_conc[],
                          const double drug_factors[])
    {
        bool hasSecondDrug = drug2_conc != 0;
        PCS_VERBOSE(
            cout << "\n\033[32mTesting \033[1m" << drugName ;
            if( hasSecondDrug ) {
                cout << " - " << drug2Name << " Conversion";
            }
            cout <<  endl << "====\033[0m" << endl;
        )
        size_t drugIndex = LSTMDrugType::findDrug( drugName );
        size_t drug2Ind = hasSecondDrug ? LSTMDrugType::findDrug( drug2Name ) : 0;
        double conc_abs_tol = LSTMDrugType::get(drugIndex).getNegligibleConcentration();
        double conc_abs_tol2 = hasSecondDrug ? LSTMDrugType::get(drug2Ind).getNegligibleConcentration() : 0;
        const size_t maxDays = 6;
        PCS_VERBOSE(double res_Fac[maxDays];)
        double res_Conc[maxDays];
        double res_Conc2[maxDays];
        double totalFac = 1;
        for( size_t i = 0; i < maxDays; i++){
            // before update (after last step):
            double fac = proxy->getDrugFactor(m_rng, inf, bodymass);
            totalFac *= fac;
            TS_ASSERT_APPROX_TOL (totalFac, drug_factors[i], PKPD_FACT_REL_TOL, PKPD_FACT_ABS_TOL);
            PCS_VERBOSE(res_Fac[i] = totalFac;)
            
            // update (two parts):
            UnittestUtil::incrTime(sim::oneDay());
            proxy->decayDrugs(bodymass);
            
            // after update:
            res_Conc[i] = proxy->getDrugConc(drugIndex);
            TS_ASSERT_APPROX_TOL (res_Conc[i], drug_conc[i], PKPD_CONC_REL_TOL, conc_abs_tol);
            res_Conc2[i] = hasSecondDrug ? proxy->getDrugConc(drug2Ind) : 0.0;
            if( hasSecondDrug ) TS_ASSERT_APPROX_TOL (res_Conc2[i], drug2_conc[i], PKPD_CONC_REL_TOL, conc_abs_tol2);
            
            // medicate (take effect on next update):
            medicate( drugIndex, i );
        }
        PCS_VERBOSE(
            string fmt = drugDebugOutputHeader(hasSecondDrug, drugName);
            for( size_t i = 0; i < maxDays; i++){
                // calculate relative and absolute differences to expected values

                double f_abs_error = res_Fac[i] - drug_factors[i];
                double f_rel_error = floor((res_Fac[i] / drug_factors[i] -1 )*1000000)/10000;
                double c_abs_error = res_Conc[i] - drug_conc[i];
                double c_rel_error = floor((res_Conc[i] / drug_conc[i] - 1 )*1000000)/10000;
                double c2_abs_error = hasSecondDrug ? res_Conc2[i] - drug2_conc[i] : 0.0;
                double c2_rel_error = hasSecondDrug ? floor((res_Conc2[i] / drug2_conc[i] - 1 )*1000000)/10000: 0.0;

                // (parent) drug debug
                const char * type = hasSecondDrug ? "P" : "";
                drugDebugOutputLine(i,
                        res_Fac[i], f_abs_error, f_rel_error, res_Conc[i],
                        c_abs_error, c_rel_error, conc_abs_tol, type, fmt);

                // metabolite debug
                if( hasSecondDrug ) {
                    drugDebugOutputLine(i,
                            res_Fac[i], f_abs_error, f_rel_error,
                            res_Conc2[i], c2_abs_error, c2_rel_error, conc_abs_tol2, "M", fmt);
                }
            }
        )

    }
    
    void runDrugSimulations (string drugName, const double drug_conc[],
                          const double drug_factors[])
    {
        runDrugSimulations(drugName, "", drug_conc, 0, drug_factors);
    }
    
    void medicate ( size_t drugIndex,  size_t i) {
            typedef multimap<size_t,pair<double, double> >::const_iterator iter;
            pair<iter, iter> doses_tmp = schedule.equal_range(i);
            for( iter it = doses_tmp.first; it != doses_tmp.second; it++){
                const double time = it->second.first, qty = it->second.second;
                UnittestUtil::medicate( m_rng, *proxy, drugIndex, qty, time );
            }
    }
    
    void testAR1 () { /* Artemether no conversion */
        const double dose = 1.7 * bodymass;   // 1.7 mg/kg * 50 kg
        assembleHexDosageSchedule(dose);
        const double drug_conc[] = { 0.0, 0.01535201, 0.01564467, 0.01565025, 0.0002983425, 5.687336e-06 };
        const double drug_factors[] = { 1, 1.033933e-12, 1.068873e-24, 1.103296e-36, 1.734223e-42, 1.729046e-42 };
        runDrugSimulations("AR1", drug_conc, drug_factors);
    }
    
    void testAR () { /* Artemether with conversion */
        const double dose = 1.7 * bodymass;   // 1.7 mg/kg * 50 kg
        assembleHexDosageSchedule(dose);
        const double AR_conc[] = { 0, 0.0001825220, 0.0001825231, 0.0001825231, 1.146952e-09, 7.189475e-15 };
        const double DHA_conc[] = { 0, 0.0002013114, 0.0002013126, 0.0002013126, 1.266891e-09, 7.941293e-15 };
        const double drug_factors[] = { 1, 1.695266e-07, 2.838279e-14, 4.740382e-21, 4.751844e-21, 4.751846e-21 };
        runDrugSimulations("AR", "DHA_AR", AR_conc, DHA_conc, drug_factors);
    }
    
    void testAS1 () { /* Artesunate no conversion */
        const double dose = 4 * bodymass;   // 4 mg/kg * 50 kg
        assembleTripleDosageSchedule(dose);
        const double drug_conc[] = { 0, 8.983362e-08, 8.983362e-08, 8.983362e-08, 5.54818e-15, 3.42659e-22 };
        const double drug_factors[] = { 1, 1.204675e-05, 1.451241e-10, 1.748061e-15, 1.748273e-15, 1.748272e-15 };
        runDrugSimulations("AS1", drug_conc, drug_factors);
    }
    
    void testAS () { /* Artesunate with conversion */
        const double dose = 4 * bodymass;   // 4 mg/kg * 50 kg
        assembleTripleDosageSchedule(dose);
        const double AS_conc[] = { 0, 2.301305e-14, 2.301305e-14, 2.301305e-14, 8.245500e-28, 2.954336e-41 };
        const double DHA_conc[] = { 0, 1.142491e-10, 1.142491e-10, 1.142491e-10, 1.067784e-21, 9.940541e-33 };
        // These are the factors produced by Kay et al with a slightly different formula:
        const double drug_factors[] = { 1, 0.0005152782, 2.655117e-07, 1.368124e-10, 1.368124e-10, 1.368124e-10 };
        runDrugSimulations("AS", "DHA_AS", AS_conc, DHA_conc, drug_factors);
    }
    
    void testCQ () {
        const double dose = 10 * bodymass;
        assembleCQDosageSchedule(dose);
        const double drug_conc[] = { 0.0, 0.03257216, 0.06440052, 0.07921600, 0.07740709, 0.07563948 };
        const double drug_factors[] = { 1, 9.259311e-02, 4.623815e-03, 2.057661e-04, 9.262133e-06, 4.218529e-07 };
        runDrugSimulations("CQ", drug_conc, drug_factors);
    }
    
    void testDHA () {
        const double dose = 4 * bodymass;   // 4 mg/kg * 50 kg
        assembleTripleDosageSchedule(dose);
        const double drug_conc[] = { 0, 6.758386e-09, 6.758386e-09, 6.758386e-09, 1.701423e-17, 4.28333e-26 };
        const double drug_factors[] = { 1, 0.0003552336, 1.261909e-07,4.482726e-11, 4.482726e-11, 4.482726e-11 };
        runDrugSimulations("DHA", drug_conc, drug_factors);
    }
    
    void testLF () {
        const double dose = 12 * bodymass;   // 12 mg/kg * 50 kg
        assembleHexDosageSchedule(dose);
        const double drug_conc[] = { 0, 1.014434363, 1.878878305, 2.615508841, 2.228789614, 1.899249226 };
        const double drug_factors[] = { 1, 0.03174632, 0.001007809, 3.199346e-05, 1.015654e-06, 3.224254e-08 };
        runDrugSimulations("LF", drug_conc, drug_factors);
    }
    
    void testMQ () {
        const double dose = 8.3 * bodymass;   // 8.3 mg/kg * 50 kg
        assembleTripleDosageSchedule( dose) ;
        const double drug_conc[] = { 0, 0.378440101, 0.737345129, 1.077723484, 1.022091411, 0.969331065 };
        const double drug_factors[] = { 1, 0.03174581, 0.001007791, 3.199298e-05, 1.015638e-06, 3.224205e-08 };
        runDrugSimulations("MQ", drug_conc, drug_factors);
    }
    
    // PPQ with a 1-compartment model (WinterHasting2011_single; not preferred)
    void testPPQ1C (){
        const double dose = 18 * bodymass;   // 18 mg/kg * 50 kg
        assembleTripleDosageSchedule( dose );
        const double drug_conc[] = { 0, 0.116453464, 0.2294652081, 0.339137, 0.3291139387, 0.3193871518 };
        const double drug_factors[] = { 1, 0.03174892, 0.001007891, 3.199625e-05, 1.015747e-06, 3.224518e-08 };
        runDrugSimulations("PPQ", drug_conc, drug_factors);
    }
    
    // PPQ with a 2-compartment model (Hodel2013)
    void testPPQ_Hodel2013 (){
        const double dose = 18 * bodymass;   // 18 mg/kg * 50 kg
        assembleTripleDosageSchedule( dose );
        const double drug_conc[] = { 0, 0.08022449, 0.1416033, 0.18962337, 0.14792303, 0.11829172 };
        const double drug_factors[] = { 1, 0.03422595, 0.001086594, 3.449438e-05, 1.095144e-06, 3.479034e-08 };
        runDrugSimulations("PPQ2", drug_conc, drug_factors);
    }
    
    // PPQ with a 3-compartment model (Tarning 2012 AAC)
    void testPPQ_Tarning2012AAC (){
        const double dose = 18 * bodymass;   // 18 mg/kg * 50 kg
        assembleTripleDosageSchedule( dose );
        const double drug_conc[] = { 0, 0.075305088, 0.118119866, 0.150210662, 0.100426437, 0.078729041 };
        const double drug_factors[] = { 1, 0.03421756, 0.001086307, 3.448539e-05, 1.094894e-06, 3.478302e-08 };
        runDrugSimulations("PPQ3", drug_conc, drug_factors);
    }
    
private:
    LocalRng m_rng;
    LSTMModel *proxy;
    CommonInfection *inf;
    double bodymass;
    std::multimap<size_t, pair<double, double> > schedule; // < day, pair < part of day, dosage > >
};

#endif
