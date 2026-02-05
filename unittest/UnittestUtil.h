/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2026 Swiss Tropical and Public Health Institute
 * Copyright (C) 2005-2015 Liverpool School Of Tropical Medicine
 * Copyright (C) 2020-2026 University of Basel
 * Copyright (C) 2025-2026 The Kids Research Institute Australia
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

// Utility for unittests, which is granted "friend" access where necessary in the model.

#ifndef Hmod_UnittestUtil
#define Hmod_UnittestUtil

#include "Global.h"
#include "util/ModelOptions.h"

#include "Clinical/ClinicalModel.h"
#include "Host/Human.h"
#include "PkPd/LSTMModel.h"
#include "PkPd/Drug/LSTMDrugType.h"
#include "PkPd/LSTMTreatments.h"
#include "Host/WithinHost/WHInterface.h"
#include "Host/WithinHost/Infection/Infection.h"
#include "Host/WithinHost/WHFalciparum.h"
#include "Host/WithinHost/Infection/MolineauxInfection.h"
#include "Host/WithinHost/Genotypes.h"
#include "mon/management.h"

#include "schema/scenario.h"

using namespace OM;
using namespace WithinHost;

using ::OM::util::ModelOptions;

namespace OM {
    namespace WithinHost {
        extern bool opt_common_whm;
    }
}

namespace xml_helpers{
    /// Construct a helper for setting PK parameters (1-compartment model)
    ///@param Vd Volume of distribution (l/kg)
    ///@param negl_conc Negligible concentration of drug (mg/l)
    ///@param hl Half-life of drug concentration (days)
    struct PK1Chl{
        PK1Chl(double Vd, double negl_conc, double hl):
            Vd(Vd), hl(hl), negl_conc(negl_conc) {}
        double Vd, hl, negl_conc;
    };
    /// Construct a helper for setting PK parameters (1-compartment model)
    ///@param Vd Volume of distribution (l/kg)
    ///@param negl_conc Negligible concentration of drug (mg/l)
    ///@param k Elimination rate of drug (days^-1)
    ///@param m_exponent
    struct PK1C{
        PK1C(double Vd, double negl_conc, double k, double m_exponent):
            Vd(Vd), negl_conc(negl_conc), k(k), me(m_exponent) {}
        double Vd, negl_conc, k, me;
    };
    /// Construct a helper for setting PK parameters (1-compartment plus conversion)
    ///@param Vd Volume of distribution (l/kg)
    ///@param negl_conc Negligible concentration of drug (mg/l)
    ///@param k Direct elimination rate of drug (days^-1)
    ///@param m_exponent
    ///@param k_a Absorbtion rate
    ///@param metabolite Abbreviation of metabolite drug
    ///@param conv Conversion rate of drug (mg/l)
    ///@param mwr Molecular weight ratio (metabolite weight / parent weight)
    struct PKConv{
        PKConv(double Vd, double negl_conc, double k, double m_exponent,
            double k_a, const char *metabolite, double conv, double mwr):
            Vd(Vd), negl_conc(negl_conc), k(k), me(m_exponent), ka(k_a), met(metabolite), conv(conv), mwr(mwr) {}
        double Vd, negl_conc, k, me, ka;
        const char *met;
        double conv, mwr;
    };
    /// Construct a helper for setting PK parameters (2-compartment model)
    ///@param Vd Volume of distribution (l/kg)
    ///@param negl_conc Negligible concentration of drug (mg/l)
    ///@param k Elimination rate of drug (days^-1)
    ///@param m_exponent
    ///@param ka Absorbtion rate from gut
    ///@param k12 Absorbtion rate parameter k12
    ///@param k21 Absorbtion rate parameter k21
    struct PK2C{
        PK2C(double Vd, double negl_conc, double k, double m_exponent,
             double ka, double k12, double k21):
            Vd(Vd), negl_conc(negl_conc), k(k), me(m_exponent), ka(ka), k12(k12), k21(k21) {}
        double Vd, negl_conc, k, me, ka, k12, k21;
    };
    /// Construct a helper for setting PK parameters (2-compartment model)
    ///@param Vd Volume of distribution (l/kg)
    ///@param negl_conc Negligible concentration of drug (mg/l)
    ///@param k Elimination rate of drug (days^-1)
    ///@param m_exponent
    ///@param ka Absorbtion rate from gut
    ///@param k12 Absorbtion rate parameter k12
    ///@param k21 Absorbtion rate parameter k21
    ///@param k13 Absorbtion rate parameter k13
    ///@param k31 Absorbtion rate parameter k31
    struct PK3C{
        PK3C(double Vd, double negl_conc, double k, double m_exponent,
             double ka, double k12, double k21, double k13, double k31):
            Vd(Vd), negl_conc(negl_conc), k(k), me(m_exponent), ka(ka),
            k12(k12), k21(k21), k13(k13), k31(k31) {}
        double Vd, negl_conc, k, me, ka, k12, k21, k13, k31;
    };
    /// Construct a helper for setting PD parameters
    ///@param vmax Max killing rate
    ///@param ic50 IC50
    ///@param slope Slope (n)
    struct PD{
        PD(double vmax, double ic50, double slope):
            vmax(vmax), ic50(ic50), slope(slope) {}
        double vmax, ic50, slope;
    };
    /// Helper for constructing an XML element with parameters for some drug
    ///@param abbrev Drug name (CQ, PPQ, AR, etc)
    ///@param pk Helper struct with pharmaco-kinetic parameters
    ///@param pd Helper struct with pharmaco-dynamic parameters
    scnXml::PKPDDrug drug(const char *abbrev, PK1Chl pk, PD pd){
        scnXml::PK xPK(pk.negl_conc, pk.Vd);
        xPK.setHalf_life(pk.hl);
        scnXml::PD xPD;
        xPD.getPhenotype().push_back(scnXml::Phenotype(pd.vmax, pd.ic50, pd.slope));
        return scnXml::PKPDDrug(xPD, xPK, abbrev);
    }
    /// Helper for constructing an XML element with parameters for some drug
    ///@param abbrev Drug name (CQ, PPQ, AR, etc)
    ///@param pk Helper struct with pharmaco-kinetic parameters
    ///@param pd Helper struct with pharmaco-dynamic parameters
    scnXml::PKPDDrug drug(const char *abbrev, PK1C pk, PD pd){
        scnXml::PK xPK(pk.negl_conc, pk.Vd);
        xPK.setK(scnXml::SampledValueLN(pk.k));
        xPK.setM_exponent(pk.me);
        scnXml::PD xPD;
        xPD.getPhenotype().push_back(scnXml::Phenotype(pd.vmax, pd.ic50, pd.slope));
        return scnXml::PKPDDrug(xPD, xPK, abbrev);
    }
    /// Helper for constructing an XML element with parameters for some drug
    ///@param abbrev Drug name (CQ, PPQ, AR, etc)
    ///@param pk Helper struct with pharmaco-kinetic parameters
    ///@param pd Helper struct with pharmaco-dynamic parameters
    scnXml::PKPDDrug drug(const char *abbrev, PKConv pk, PD pd, double IC50corr){
        scnXml::PK xPK(pk.negl_conc, pk.Vd);
        xPK.setK(scnXml::SampledValueLN(pk.k));
        xPK.setM_exponent(pk.me);
        xPK.setK_a(scnXml::SampledValueLN(pk.ka));
        xPK.setConversion(scnXml::Conversion(pk.met,
            scnXml::SampledValueLN(pk.conv), pk.mwr, IC50corr));
        scnXml::PD xPD;
        xPD.getPhenotype().push_back(scnXml::Phenotype(pd.vmax, pd.ic50, pd.slope));
        return scnXml::PKPDDrug(xPD, xPK, abbrev);
    }
    /// Helper for constructing an XML element with parameters for some drug
    ///@param abbrev Drug name (CQ, PPQ, AR, etc)
    ///@param pk Helper struct with pharmaco-kinetic parameters
    ///@param pd Helper struct with pharmaco-dynamic parameters
    scnXml::PKPDDrug drug(const char *abbrev, PK2C pk, PD pd){
        scnXml::PK xPK(pk.negl_conc, pk.Vd);
        xPK.setK(scnXml::SampledValueLN(pk.k));
        xPK.setCompartment2(scnXml::Compartment2(
            scnXml::SampledValueLN(pk.k12),
            scnXml::SampledValueLN(pk.k21)));
        xPK.setM_exponent(pk.me);
        xPK.setK_a(scnXml::SampledValueLN(pk.ka));
        scnXml::PD xPD;
        xPD.getPhenotype().push_back(scnXml::Phenotype(pd.vmax, pd.ic50, pd.slope));
        return scnXml::PKPDDrug(xPD, xPK, abbrev);
    }
    /// Helper for constructing an XML element with parameters for some drug
    ///@param abbrev Drug name (CQ, PPQ, AR, etc)
    ///@param pk Helper struct with pharmaco-kinetic parameters
    ///@param pd Helper struct with pharmaco-dynamic parameters
    scnXml::PKPDDrug drug(const char *abbrev, PK3C pk, PD pd){
        scnXml::PK xPK(pk.negl_conc, pk.Vd);
        xPK.setK(scnXml::SampledValueLN(pk.k));
        xPK.setCompartment2(scnXml::Compartment2(
            scnXml::SampledValueLN(pk.k12),
            scnXml::SampledValueLN(pk.k21)));
        xPK.setCompartment3(scnXml::Compartment3(
            scnXml::SampledValueLN(pk.k13),
            scnXml::SampledValueLN(pk.k31)));
        xPK.setM_exponent(pk.me);
        xPK.setK_a(scnXml::SampledValueLN(pk.ka));
        scnXml::PD xPD;
        xPD.getPhenotype().push_back(scnXml::Phenotype(pd.vmax, pd.ic50, pd.slope));
        return scnXml::PKPDDrug(xPD, xPK, abbrev);
    }
}

namespace dummyXML{
    scnXml::DemogAgeGroup demAgeGroup(
        numeric_limits<double>::quiet_NaN() /* lower bound */ );
    scnXml::Demography demography( 
        demAgeGroup,
        "dummy" /* name */,
        0 /* pop size */,
        90.0 /*max human age*/ );
    
    scnXml::MonitoringOptions survOpts;
    scnXml::Surveys surveys;
    scnXml::MonAgeGroup monAgeGroup(
        0.0 /* lower bound */ );
    scnXml::Monitoring monitoring(
        survOpts,
        surveys,
        monAgeGroup,
        "dummy" /* name */ );
    
    scnXml::Interventions interventions( "dummy" /* name */ );
    
    scnXml::AgeGroupValues cfr, pSeq;
    scnXml::HealthSystem healthSystem( cfr, pSeq );
    
    scnXml::Entomology entomology( "dummy" /* name */,
        "dummy" /* mode */ );
    
    scnXml::OptionSet modelOpts;
    scnXml::Clinical modelClinical( "dummy" /* HS memory */ );
    scnXml::AgeGroupValues modelHumanAvailMosq;
    scnXml::Human modelHuman( modelHumanAvailMosq );
    scnXml::ComputationParameters computationParams(0 /* iseed */);
    scnXml::Parameters modelParams(0 /* interval*/, "dummy" /* latentP */);
    scnXml::Model model(modelClinical, modelHuman, computationParams);
    
    scnXml::Scenario scenario(
        demography,
        monitoring,
        interventions,
        healthSystem,
        entomology,
        model,
        0 /* schema version */,
        "dummy" /* name */ );
}

class UnittestUtil {
public:
    static void initTime(int daysPerStep){
        dummyXML::model.setModelOptions(dummyXML::modelOpts);
        dummyXML::modelParams.setInterval( daysPerStep );
        dummyXML::model.setParameters(dummyXML::modelParams);
        dummyXML::model.setComputationParameters(dummyXML::computationParams);
        dummyXML::scenario.setModel(dummyXML::model);
        dummyXML::surveys.setDetectionLimit( numeric_limits<double>::quiet_NaN() );
        dummyXML::surveys.getSurveyTime().push_back( scnXml::SurveyTime ( "1t" ) );
        dummyXML::monitoring.setSurveys( dummyXML::surveys );
        dummyXML::scenario.setMonitoring( dummyXML::monitoring );
        util::ModelNameProvider modelNameProvider(dummyXML::scenario.getModel());
        sim::init( dummyXML::scenario, modelNameProvider );
        
        // we could just use zero, but we may spot more errors by using some weird number
        sim::s_t0 = sim::fromYearsN(83.2591);
        sim::s_t1 = sim::s_t0;
#ifndef NDEBUG
        sim::in_update = true;  // may not always be correct but we're more interested in getting around this check than using it in unit tests
#endif
    }
    static void incrTime(SimTime incr){
        //NOTE: for unit tests, we do not differentiate between s_t0 and s_t1
        sim::s_t0 = sim::s_t0 + incr;
        sim::s_t1 = sim::s_t0;
    }
    
    static const scnXml::Parameters& prepareParameters(){
        if( dummyXML::modelParams.getParameter().size() == 0 ){
            dummyXML::modelParams.getParameter().push_back( scnXml::Parameter( 15, 0.177378570987455 ) );
            dummyXML::modelParams.getParameter().push_back( scnXml::Parameter( 34, 4.7601 ) );
            dummyXML::modelParams.getParameter().push_back( scnXml::Parameter( 35, 0.5008 ) );
            dummyXML::modelParams.getParameter().push_back( scnXml::Parameter( 36, 2.2736 ) );
            dummyXML::modelParams.getParameter().push_back( scnXml::Parameter( 37, 0.2315 ) );
        }
        return dummyXML::modelParams;
    }
    
    // Parameterise standard diagnostics
    static void setDiagnostics(){
        // note that this is only ever called after initSurveys(), thus we
        // don't need to call diagnostics::clear() (and shouldn't because it
        // would leave a dangling pointer in Survey::m_diagnostic)
        scnXml::Diagnostic microscopy( "microscopy" );
        scnXml::Stochastic microscopyParams( 20, 0.75 );
        microscopy.setStochastic( microscopyParams );
        scnXml::Diagnostic rdt( "RDT" );
        scnXml::Stochastic rdtParams( 50, 0.942 );
        rdt.setStochastic( rdtParams );
        scnXml::Diagnostics diagsElt;
        diagsElt.getDiagnostic().push_back( microscopy );
        diagsElt.getDiagnostic().push_back( rdt );
        dummyXML::model.setParameters(prepareParameters());
        Parameters parameters( dummyXML::model.getParameters(), util::ModelNameProvider(dummyXML::model) );
        dummyXML::scenario.setDiagnostics(diagsElt);
        dummyXML::surveys.setDetectionLimit(numeric_limits<double>::quiet_NaN());
        dummyXML::monitoring.setSurveys(dummyXML::surveys);
        dummyXML::scenario.setMonitoring(dummyXML::monitoring);
        diagnostics::init( parameters, dummyXML::scenario );
    }
    
    static void PkPdSuiteSetup () {
        ModelOptions::reset();
        WithinHost::Genotypes::initSingle();

        //Note: we fudge this call since it's not so easy to falsely initialize scenario element.
        //PkPdModel::init ();
        using namespace xml_helpers;

        // Drugs:
        scnXml::Drugs drugs;
        
        // Artemether (no conversion model)
        drugs.getDrug().push_back(drug("AR1",
                PK1C(17.4 /*Vd*/, 1e-17 /*negl_conc*/, 3.96 /*k*/, 0.0 /*m_exp*/),
                PD(27.6 /* vmax */, 0.0023 /* IC50 */, 4.0 /* slope */ )));
        // Artemether plus conversion to DHA
        drugs.getDrug().push_back(drug("DHA_AR",
                PK1C(15 /*Vd*/, 1e-17 /*negl_conc*/, 44.15 /*k*/, 0.0 /*m_exp*/),
                PD(27.6 /* vmax */, 0.009 /* IC50 */, 4.0 /* slope */ )));
        drugs.getDrug().push_back(drug("AR",
                PKConv(46.6 /*Vd*/, 1e-17 /*negl_conc*/, 0 /*k*/, 0 /*m_exp*/, 
                       23.98 /*absorption rate*/, "DHA_AR" /*metabolite*/,
                       11.98 /*conversion rate*/, 0.9547587 /*mol. weight ratio*/),
                PD(27.6 /* vmax */, 0.0023 /* IC50 */, 4.0 /* slope */ ),
                1.0 /* IC50 correlation */));
        
        
        // Artesunate (no conversion model)
        drugs.getDrug().push_back(drug("AS1",
                PK1C(2.75 /*Vd*/, 1e-17 /*negl_conc*/, 16.6 /*k*/, 0.0 /*m_exp*/),
                PD(27.6 /* vmax */, 0.0016 /* IC50 */, 4.0 /* slope */ )));
        // Artesunate plus conversion to DHA
        drugs.getDrug().push_back(drug("DHA_AS",
                PK1C(1.49 /*Vd*/, 1e-35 /*negl_conc*/, 25.4 /*k*/, 0.0 /*m_exp*/),
                PD(27.6 /* vmax */, 0.009 /* IC50 */, 4.0 /* slope */ )));
        drugs.getDrug().push_back(drug("AS",
                PKConv(7.1 /*Vd*/, 1e-45 /*negl_conc*/, 0 /*k*/, 0 /*m_exp*/,
                       252 /*absorption rate*/, "DHA_AS" /*metabolite*/,
                       30.96 /*conversion rate*/, 0.741155 /*mol. weight ratio*/),
                PD(27.6 /* vmax */, 0.0016 /* IC50 */, 4.0 /* slope */ ),
                1.0 /* IC50 correlation */));
        
        // Dihydroartemisinin (when not a metabolite)
        drugs.getDrug().push_back(drug("DHA",
                PK1C(1.49 /*Vd*/, 1e-17 /*negl_conc*/, 19.8 /*k*/, 0.0 /*m_exp*/),
                PD(27.6 /* vmax */, 0.009 /* IC50 */, 4.0 /* slope */ )));
        
        drugs.getDrug().push_back(drug("CQ",    // Chloroquine
                PK1Chl(300 /*Vd*/, 0.00036 /*negl_conc*/, 30.006 /*hl*/),
                PD(3.45 /* vmax */, 0.02 /* IC50 */, 1.6 /* slope */ )));
        drugs.getDrug().push_back(drug("LF",    // Lumefantrine
                PK1C(21 /*Vd*/, 0.00032 /*negl_conc*/, 0.16 /*k*/, 0.0 /*m_exp*/),
                PD(3.45 /* vmax */, 0.032 /* IC50 */, 4.0 /* slope */ )));
        drugs.getDrug().push_back(drug("MQ",    // Mefloquine
                PK1Chl(20.8 /*Vd*/, 0.005 /*negl_conc*/, 13.078 /*hl*/),
                PD(3.45 /* vmax */, 0.027 /* IC50 */, 5.0 /* slope */ )));
        
        drugs.getDrug().push_back(drug("PPQ",   // Piperaquine, 1-compartment
                PK1C(150 /*Vd*/, 0.005 /*negl_conc*/, 0.03 /*k*/, 0.0 /*m_exp*/),
                PD(3.45 /* vmax */, 0.020831339 /* IC50 */, 6.0 /* slope */ )));
        drugs.getDrug().push_back(drug("PPQ2",   // Piperaquine, Hodel2013 model
                PK2C(173 /*Vd*/, 0.005 /*negl_conc*/, 0.2452253 /*k*/, 0.25 /*m_exp*/,
                    11.16 /*k_a*/, 0.2014864 /*k12*/, 0.07870968 /*k21*/),
                PD(3.45 /* vmax */, 0.020831339 /* IC50 */, 6.0 /* slope */ )));
        drugs.getDrug().push_back(drug("PPQ3",  // Piperaquine, Tarning 2012 AAC
                PK3C(57.5625 /*Vd*/, 0.005 /*negl_conc*/, 16.314788273615637 /*k*/, 1.0 /*m_exp*/,
                    3.4825 /*k_a*/, 1.854166666666667 /*k12*/, 1.1545945945945946 /*k21*/,
                     0.9027777777777778 /*k13*/, 0.07948639559767655 /*k31*/),
                PD(3.45 /* vmax */, 0.020831339 /* IC50 */, 6.0 /* slope */ )));
        
        PkPd::LSTMDrugType::init (drugs);
        
        // Treatments
        scnXml::PKPDSchedule sched1("sched1");
        sched1.getMedicate().push_back(
            scnXml::PKPDMedication("MQ", 6 /*mg*/,0 /*hour*/ )
        );
        
        scnXml::PKPDSchedule sched2("sched2");
        sched2.getMedicate().push_back(
            scnXml::PKPDMedication("MQ", 2 /*mg*/, 0 /*hour*/));
        sched2.getMedicate().push_back(
            scnXml::PKPDMedication("MQ", 5 /*mg*/, 12 /*hour*/));
        
        // a very basic dosage table, so that we can test it does what's expected
        scnXml::PKPDDosages dosage1("dosage1");
        dosage1.getAge().push_back(scnXml::PKPDDosageRange(0 /*age lb*/,1 /*mult*/));
        dosage1.getAge().push_back(scnXml::PKPDDosageRange(5 /*age lb*/,5 /*mult*/));
        
        scnXml::Treatments treatments;
        treatments.getSchedule().push_back(sched1);
        treatments.getSchedule().push_back(sched2);
        treatments.getDosages().push_back(dosage1);
        PkPd::LSTMTreatments::init(treatments);
    }
    
    // For when infection parameters shouldn't be used; enforce by setting to NaNs.
    // But do set latentP.
    static void Infection_init_latentP_and_NaN () {
	Infection::s_latentP = sim::fromDays(15);
    }
    
    static void DescriptiveInfection_init () {
        ModelOptions::reset();
    }
    
    static void EmpiricalWHM_setup () {
        ModelOptions::reset();
        ModelOptions::set(util::EMPIRICAL_WITHIN_HOST_MODEL);
        OM::WithinHost::opt_common_whm = true;
    }
    
    static void MolineauxWHM_setup( const std::string& mode, bool repl_gamma ){
        ModelOptions::reset();
        ModelOptions::set(util::MOLINEAUX_WITHIN_HOST_MODEL);
        OM::WithinHost::opt_common_whm = true;
        if( mode == "original" ){
        }else if( mode == "1st_max_gamma" ){
            ModelOptions::set(util::FIRST_LOCAL_MAXIMUM_GAMMA);
        }else if( mode == "mean_dur_gamma" ){
            ModelOptions::set(util::MEAN_DURATION_GAMMA);
        }else if( mode == "1st_max_and_mean_dur_gamma" ){
            ModelOptions::set(util::FIRST_LOCAL_MAXIMUM_GAMMA);
            ModelOptions::set(util::MEAN_DURATION_GAMMA);
        }else if( mode == "pairwise" ){
            ModelOptions::set(util::MOLINEAUX_PAIRWISE_SAMPLE);
        }else{
            ETS_ASSERT( false );        // stop this test
        }
        if( repl_gamma ){
            ModelOptions::set(util::PARASITE_REPLICATION_GAMMA);
        }
        
        // Set parameters; all of these were estimated externally from OpenMalaria
        dummyXML::model.setParameters(prepareParameters());
        Parameters parameters( dummyXML::model.getParameters(), util::ModelNameProvider(dummyXML::model) );
        
        // This sets up the model based on parameters and options
        WithinHost::MolineauxInfection::init(parameters);
    }
    
    static void MosqLifeCycle_init() {
        ModelOptions::reset();
        ModelOptions::set(util::VECTOR_LIFE_CYCLE_MODEL);
    }
    
    static double getPrescribedMg( const PkPd::LSTMModel& pkpd ){
        double r = 0.0;
        for( const PkPd::MedicateData& md : pkpd.medicateQueue ){
            r += md.qty;
        }
        return r;
    }
    
    static void medicate(LocalRng& rng, PkPd::LSTMModel& pkpd, size_t typeIndex, double qty,
                         double time){
        pkpd.medicateDrug(rng, typeIndex, qty, time);
    }
    
    static void clearMedicateQueue( PkPd::LSTMModel& pkpd ){
        pkpd.medicateQueue.clear();
    }
    
    static unique_ptr<Host::Human> createHuman(SimTime dateOfBirth){
        return unique_ptr<Host::Human>( new Host::Human(dateOfBirth) );
    }
    // Set the WithinHost model used by the human, and return a pointer to it. Do not delete this!
    static WithinHost::WHInterface* setHumanWH(Host::Human& human, unique_ptr<WithinHost::WHInterface> wh){
        human.withinHostModel = std::move(wh);
        return &*human.withinHostModel;
    }
};

#endif
