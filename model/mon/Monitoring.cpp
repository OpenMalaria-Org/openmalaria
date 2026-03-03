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

#include "mon/Monitoring.h"
#include "mon/AgeGroup.h"
#include "Host/WithinHost/Diagnostic.h"
#include "Host/WithinHost/Genotypes.h"
#include "Clinical/ClinicalModel.h"
#include "Host/Human.h"
#include "interventions/InterventionManager.h"
#include "util/CommandLine.h"
#include "util/errors.h"
#include "util/UnitParse.h"
#include "schema/monitoring.h"
#include "schema/scenario.h"

#include <cctype>
#include <algorithm>
#include <fstream>
#include <gzstream/gzstream.h>
#include <iostream>
#include <numeric>

#include <stdexcept>
namespace OM {
namespace mon {

NamedMeasureMapT namedOutMeasures;
set<Measure> validCondMeasures;
// This method defines output measures accepted by name in the XML (e.g.
// "nHost") and their numeric output identifier (i.e. measure column of
// outputs), type of output (integer or floating point), aggregation, and the
// corresponding internal measure code.
void defineOutMeasures(){
    /* Don't ever make an existing numerical identifier point to a new/different 
    output measure because this would violate users' expectations */
    // Don't reuse old numerical identifiers
    namedOutMeasures.clear();
    validCondMeasures.clear();

    std::set<int> seenIDs;

    struct NamedDef {
        const char* name;
        OutMeasure om;
    };

    static const NamedDef defs[] = {

        /// Total number of humans
        {"nHost", OutMeasure::humanAC( 0, nHost, false )},
        /** The number of human hosts with an infection (patent or not) at the time
        * the survey is taken. */
        {"nInfect", OutMeasure::humanAC( 1, nInfect, false )},
        {"nInfect_Imported", OutMeasure::humanAC( 1001, nInfect_Imported, false )},
        {"nInfect_Introduced", OutMeasure::humanAC( 2001, nInfect_Introduced, false )},
        {"nInfect_Indigenous", OutMeasure::humanAC( 3001, nInfect_Indigenous, false )},
        /** Expected number of infected hosts
        *
        * This is the sum of the probabilities, across all time steps since the
        * last survey, of each host becoming infected on that time step. */
        {"nExpectd", OutMeasure::humanAC( 2, nExpectd, true )},
        /** The number of human hosts whose total (blood-stage) parasite density is
        * above the detection threshold */
        {"nPatent", OutMeasure::humanAC( 3, nPatent, false )},
        {"nPatent_Imported", OutMeasure::humanAC( 1003, nPatent_Imported, false )},
        {"nPatent_Introduced", OutMeasure::humanAC( 2003, nPatent_Introduced, false )},
        {"nPatent_Indigenous", OutMeasure::humanAC( 3003, nPatent_Indigenous, false )},
        /// Sum of log(1 + p) where p is the pyrogenic threshold
        {"sumLogPyrogenThres", OutMeasure::humanAC( 4, sumLogPyrogenThres, true )},
        /** Sum (across hosts) of the natural logarithm of the parasite density of
        * hosts with detectable parasite density (patent according to the
        * monitoring diagnostic). */
        {"sumlogDens", OutMeasure::humanAC( 5, sumlogDens, true )},
        /** The total number of infections in the population: includes both blood
        * and liver stages. Vivax: this is the number of broods. */
        {"totalInfs", OutMeasure::humanACG( 6, totalInfs, false )},
        {"totalInfs_Imported", OutMeasure::humanACG( 1006, totalInfs_Imported, false )},
        {"totalInfs_Introduced", OutMeasure::humanACG( 2006, totalInfs_Introduced, false )},
        {"totalInfs_Indigenous", OutMeasure::humanACG( 3006, totalInfs_Indigenous, false )},
        /** Infectiousness of human population to mosquitoes
        *
        * Number of hosts transmitting to mosquitoes (i.e. proportion of
        * mosquitoes that get infected multiplied by human population size).
        * Single value, not per age-group. */
        {"nTransmit", OutMeasure::value( 7, nTransmit, true )},
        /** The sum of all detectable infections (where blood stage parasite
        * density is above the detection limit) across all human hosts.
        * Vivax: the number of broods with an active blood stage. */
        {"totalPatentInf", OutMeasure::humanACG( 8, totalPatentInf, false )},
        {"totalPatentInf_Imported", OutMeasure::humanACG( 1008, totalPatentInf_Imported, false )},
        {"totalPatentInf_Introduced", OutMeasure::humanACG( 2008, totalPatentInf_Introduced, false )},
        {"totalPatentInf_Indigenous", OutMeasure::humanACG( 3008, totalPatentInf_Indigenous, false )},
        /// Contribuion to immunity functions (removed)
        {"contrib", OutMeasure::obsolete( 9 )},
        /// Sum of the pyrogenic threshold
        {"sumPyrogenThresh", OutMeasure::humanAC( 10, sumPyrogenThresh, true )},
        /// number of blood-stage treatments (1st line)
        {"nTreatments1", OutMeasure::humanAC( 11, nTreatments1, false )},
        /// number of blood-stage treatments (2nd line)
        {"nTreatments2", OutMeasure::humanAC( 12, nTreatments2, false )},
        /// number of blood-stage treatments (inpatient)
        {"nTreatments3", OutMeasure::humanAC( 13, nTreatments3, false )},
        /// number of episodes (uncomplicated)
        {"nUncomp", OutMeasure::humanAC( 14, nUncomp, false )},
        {"nUncomp_Imported", OutMeasure::humanAC( 1014, nUncomp_Imported, false )},
        {"nUncomp_Introduced", OutMeasure::humanAC( 2014, nUncomp_Introduced, false )},
        {"nUncomp_Indigenous", OutMeasure::humanAC( 3014, nUncomp_Indigenous, false )},
        /// Number of severe episodes (severe malaria or malaria + coinfection)
        {"nSevere", OutMeasure::humanAC( 15, nSevere, false )},
        /// cases with sequelae
        {"nSeq", OutMeasure::humanAC( 16, nSeq, false )},
        /// deaths in hospital
        {"nHospitalDeaths", OutMeasure::humanAC( 17, nHospitalDeaths, false )},
        /// Number of deaths indirectly caused by malaria
        {"nIndDeaths", OutMeasure::humanAC( 18, nIndDeaths, false )},
        /// Number of deaths directly caused by malaria
        {"nDirDeaths", OutMeasure::humanAC( 19, nDirDeaths, false )},
        /** Number of vaccine doses given via EPI.
        *
        * Since schema 22, each vaccine type may be deployed independently. To be
        * roughly backwards-compatible, the first type (PEV, BSV or TBV) described
        * (with an "effect" element) will be reported. */
        {"nEPIVaccinations", OutMeasure::humanDeploy( 20, vaccinations, Deploy::CTS )},
        /** All cause infant mortality rate
        *
        * Reports death rate of infants due to all causes (malaria as modelled
        * plus fixed non-malaria attribution). Calculated via Kaplan-Meier method.
        * Units: deaths per thousand births. */
        {"allCauseIMR", OutMeasure::value( 21, allCauseIMR, true )},
        /** Number of vaccine doses given via mass campaign.
        *
        * Since schema 22, each vaccine type may be deployed independently. To be
        * roughly backwards-compatible, the first type (PEV, BSV or TBV) described
        * (with an "effect" element) will be reported. */
        {"nMassVaccinations", OutMeasure::humanDeploy( 22, vaccinations, Deploy::TIMED )},
        /// recoveries in hospital
        {"nHospitalRecovs", OutMeasure::humanAC( 23, nHospitalRecovs, false )},
        /// sequelae in hospital
        {"nHospitalSeqs", OutMeasure::humanAC( 24, nHospitalSeqs, false )},
        /// Number of IPT Doses (removed together with IPT model)
        {"nIPTDoses", OutMeasure::obsolete( 25 )},
        /** Annual Average Kappa
        *
        * Calculated once a year as sum of human infectiousness divided by initial
        * EIR summed over a year. Single value, not per age-group. */
        {"annAvgK", OutMeasure::value( 26, annAvgK, true )},
        /// Number of episodes (non-malaria fever)
        {"nNMFever", OutMeasure::humanAC( 27, nNMFever, false )},
        /// Inoculations per human (all ages) per day of year, over the last year.
        /// (Reporting removed.)
        {"innoculationsPerDayOfYear", OutMeasure::obsolete( 28 )},
        /// Kappa (human infectiousness) weighted by availability per day-of-year for the last year.
        /// (Reporting removed.)
        {"kappaPerDayOfYear", OutMeasure::obsolete( 29 )},
        /** The total number of inoculations, by age group, cohort and parasite
        * genotype, summed over the reporting period. */
        {"innoculationsPerAgeGroup", OutMeasure::humanACG( 30, innoculationsPerAgeGroup, true )},
        /// N_v0: emergence of feeding vectors during the last time step. Units: mosquitoes/day
        {"Vector_Nv0", OutMeasure::species( 31, Vector_Nv0, false )},
        /// N_v: vectors seeking to feed during the last time step. Units: mosquitoes/day
        {"Vector_Nv", OutMeasure::species( 32, Vector_Nv, false )},
        /// N_v: infected vectors seeking to feed during the last time step. Units: mosquitoes/day
        {"Vector_Ov", OutMeasure::species( 33, Vector_Ov, true )},
        /// N_v: infectious vectors seeking to feed during the last time step. Units: mosquitoes/day
        {"Vector_Sv", OutMeasure::species( 34, Vector_Sv, true )},
        /** Input EIR (Expected EIR entered into scenario file)
        *
        * Units: inoculations per adult per time step. */
        {"inputEIR", OutMeasure::value( 35, inputEIR, true )},
        /** Simulated EIR (EIR output by the transmission model)
        *
        * Units: inoculations per adult per time step (children are excluded
        * when measuring). */
        {"simulatedEIR", OutMeasure::value( 36, simulatedEIR, true )},
        {"simulatedEIR_Introduced", OutMeasure::value( 2036, simulatedEIR_Introduced, true )},
        {"simulatedEIR_Indigenous", OutMeasure::value( 3036, simulatedEIR_Indigenous, true )},
        /// Number of Rapid Diagnostic Tests used
        {"Clinical_RDTs", OutMeasure::obsolete( 39 )},
        /* Effective total quanty of each drug used orally, in mg.
        * (Per active ingredient abbreviation.)
        *
        * The quantity is efffective with respect to the cost (see treatment
        * schedule definition).
        *
        * Reporting removed. */
        {"Clinical_DrugUsage", OutMeasure::obsolete( 40 )},
        /// Direct death on first day of CM (before treatment takes effect)
        {"Clinical_FirstDayDeaths", OutMeasure::humanAC( 41, Clinical_FirstDayDeaths, false )},
        /// Direct death on first day of CM (before treatment takes effect); hospital only
        {"Clinical_HospitalFirstDayDeaths", OutMeasure::humanAC( 42, Clinical_HospitalFirstDayDeaths, false )},
        /** The number of actual infections since the last survey. */
        {"nNewInfections", OutMeasure::humanAC( 43, nNewInfections, false )},
        {"nNewInfections_Imported", OutMeasure::humanAC( 1043, nNewInfections_Imported, false )},
        {"nNewInfections_Introduced", OutMeasure::humanAC( 2043, nNewInfections_Introduced, false )},
        {"nNewInfections_Indigenous", OutMeasure::humanAC( 3043, nNewInfections_Indigenous, false )},
        /** The number of ITNs delivered by mass distribution since last survey.
        *
        * These are "modelled ITNs": cover only a single person, cannot be passed
        * to someone else for reuse or used for fishing, etc. */
        {"nMassITNs", OutMeasure::humanDeploy( 44, itn, Deploy::TIMED )},
        /** The number of ITNs delivered through EPI since last survey.
        *
        * Comments from nMassITNs apply. */
        {"nEPI_ITNs", OutMeasure::humanDeploy( 45, itn, Deploy::CTS )},
        /** The number of people newly protected by IRS since last survey.
        *
        * Modelled IRS: affects one person, cannot be plastered over. */
        {"nMassIRS", OutMeasure::humanDeploy( 46, irs, Deploy::TIMED )},
        /** Defunct; was used by "vector availability" intervention (which is now a
        * sub-set of GVI). */
        {"nMassVA", OutMeasure::obsolete( 47 )},
        /// Number of malarial tests via microscopy used
        {"Clinical_Microscopy", OutMeasure::obsolete( 48 )},
        /* As Clinical_DrugUsage, but for quatities of drug delivered via IV. */
        {"Clinical_DrugUsageIV", OutMeasure::obsolete( 49 )},
        /// Number of cohort recruitments removed)
        {"nAddedToCohort", OutMeasure::obsolete( 50 )},
        /// Number of individuals removed from cohort (removed)
        {"nRemovedFromCohort", OutMeasure::obsolete( 51 )},
        /** Number of people (per age group) treated by mass drug administration
        * campaign. (Note that in one day time-step model MDA can be configured
        * as screen-and-treat. This option reports treatments administered not
        * the number of tests used.) */
        {"nMDAs", OutMeasure::humanDeploy( 52, treat, Deploy::TIMED )},
        /// Number of deaths caused by non-malaria fevers
        {"nNmfDeaths", OutMeasure::humanAC( 53, nNmfDeaths, false )},
        /// Number of antibiotic treatments given (disabled — not used)
        {"nAntibioticTreatments", OutMeasure::obsolete( 54 )},
        /** Report the number of screenings used in a mass screen-and-treat
        * operation. */
        {"nMassScreenings", OutMeasure::humanDeploy( 55, screen, Deploy::TIMED )},
        /// Report the number of mass deployments of generic vector interventions.
        {"nMassGVI", OutMeasure::humanDeploy( 56, gvi, Deploy::TIMED )},
        /** Number of IRS deployments via continuous deployment. */
        {"nCtsIRS", OutMeasure::humanDeploy( 57, irs, Deploy::CTS )},
        /** Number of GVI deployments via continuous deployment. */
        {"nCtsGVI", OutMeasure::humanDeploy( 58, gvi, Deploy::CTS )},
        /** Number of "MDA" deployments via continuous deployment.
        *
        * Note: MDA stands for mass drug administration, but the term has come to
        * be used more flexibly by OpenMalaria, including optional screening and
        * deployment through age-based systems. */
        {"nCtsMDA", OutMeasure::humanDeploy( 59, treat, Deploy::CTS )},
        /** Number of diagnostics used by "MDA" distribution through continuous
        * methods. Can be higher than nCtsMDA since drugs are administered only
        * when the diagnostic is positive. Also see nCtsMDA description. */
        {"nCtsScreenings", OutMeasure::humanDeploy( 60, screen, Deploy::CTS )},
        /** Number of removals from a sub-population due to expiry of duration of
        * membership (e.g. intervention too old). */
        {"nSubPopRemovalTooOld", OutMeasure::humanAC( 61, nSubPopRemovalTooOld, false )},
        /** Number of removals from a sub-population due to first
        * infection/bout/treatment (see onFirstBout & co). */
        {"nSubPopRemovalFirstEvent", OutMeasure::humanAC( 62, nSubPopRemovalFirstEvent, false )},
        /** Report the number of liver-stage treatments (likely Primaquine) administered. */
        {"nLiverStageTreatments", OutMeasure::humanAC( 63, nLiverStageTreatments, false )},
        /** Report the number of diagnostics used during treatment.
        *
        * This is not the same as Clinical_RDTs + Clinical_Microscopy: those
        * outputs are used by the "event scheduler" 1-day time step clinical
        * model, whereas this output is used by the 5-day time step model. */
        {"nTreatDiagnostics", OutMeasure::humanAC( 64, nTreatDiagnostics, false )},
        /** Number of "recruitment only" recruitments via timed deployment. */
        {"nMassRecruitOnly", OutMeasure::humanDeploy( 65, recruit, Deploy::TIMED )},
        /** Number of "recruitment only" recruitments via age-based deployment. */
        {"nCtsRecruitOnly", OutMeasure::humanDeploy( 66, recruit, Deploy::CTS )},
        /** Number of deployments (of all intervention components) triggered by
        * treatment (case management). */
        {"nTreatDeployments", OutMeasure::humanDeploy( 67, nTreatDeployments, Deploy::TREAT )},
        /** Report the total age of all humans in this a group (sum across humans,
        * in years). Divide by nHost to get the average age. */
        {"sumAge", OutMeasure::humanAC( 68, sumAge, true )},
        /** The number of human hosts with an infection (patent or not), for each
        * genotype, at the time the survey is taken. */
        {"nInfectByGenotype", OutMeasure::humanACG( 69, nInfectByGenotype, false )},
        /** The number of human hosts whose total (blood-stage) parasite density,
        * for each genotype, is above the detection threshold */
        {"nPatentByGenotype", OutMeasure::humanACG( 70, nPatentByGenotype, false )},
        /** For each infection genotype, sum across humans the natural log of
        * parasite density (like sumlogDens but per genotype). */
        {"logDensByGenotype", OutMeasure::humanACG( 71, logDensByGenotype, true )},
        /** For each drug type in the pharmacology section of the XML, report the
        * number of humans with non-zero concentration of this drug in their
        * blood. */
        {"nHostDrugConcNonZero", OutMeasure::humanACP( 72, nHostDrugConcNonZero, false )},
        /** For each drug type in the pharmacology section of the XML, report the
        * sum of the natural logarithm of the drug concentration in hosts with
        * non-zero concentration. */
        {"sumLogDrugConcNonZero", OutMeasure::humanACP( 73, sumLogDrugConcNonZero, true )},
        /** Expected number of direct malaria deaths, from those with severe
        * disease.
        *
        * This is calculated as the sum over all steps in the reporting period of
        * the sum over humans with severe malaria of the probability of direct
        * death from malaria. */
        {"expectedDirectDeaths", OutMeasure::humanAC( 74, expectedDirectDeaths, true )},
        /** Expected number of direct malaria deaths which occur in hospital.
        *
        * This is the a subset of `expectedDirectDeaths` and the same notes apply.
        */
        {"expectedHospitalDeaths", OutMeasure::humanAC( 75, expectedHospitalDeaths, true )},
        /** Expected number of indirect malaria deaths, from sick humans.
        *
        * This is calculated as the sum over all steps in the reporting period of
        * the sum over humans with a malaria bout (severe or not) of the
        * proability of indirect death due to malaria, assuming that they do not
        * die of another cause in the mean-time.
        *
        * Note that indirect death is only possible in the simulation when the
        * individual is sick, so the expemctation of this event is the same as were
        * it applied to all humans (sick or not).
        *
        * It does not quite tally with reports of indirect death, since the
        * probability of indirect death is calculated ahead of the actual death
        * and death may occur earlier for another reason (direct death,
        * outmigration).
        *
        * Humans already 'doomed' to die as an 'indirect mortality' are excluded
        * from the sum. */
        {"expectedIndirectDeaths", OutMeasure::humanAC( 76, expectedIndirectDeaths, true )},
        /** Expected number of sequelae, from those with severe disease.
        *
        * This is calculated as the sum over all steps in the reporting period of
        * the sum over humans with severe malaria of the probability of sequelae
        * occuring, assuming the human "recovers" from the bout.
        */
        {"expectedSequelae", OutMeasure::humanAC( 77, expectedSequelae, true )},
        /** Expected number of severe bouts of malaria.
        *
        * This is calculated as the sum over all steps in the reporting period of
        * the sum over humans with a malaria bout (severe or not) of the bout
        * becoming severe. For the 5-day time-step this is calculated once per
        * bout (which lasts one time-step). For other time-steps exact behaviour
        * is not yet defined.
        *
        * This includes both "severe malaria" and "complications due to
        * coinfection" (the same as the `nSevere` output).
        *
        * Note that this has the same expectation as the probability of a severe
        * bout when not already given that there will be a malaria bout, but may
        * be more noisy.
        */
        {"expectedSevere", OutMeasure::humanAC( 78, expectedSevere, true )},
        /** The total number of inoculations, by mosquito species, summed over
        * the reporting period. */
        {"innoculationsPerVector", OutMeasure::species( 79, innoculationsPerAgeGroup, false )},
        /** Number of custom intervention reports done */
        {"nCMDTReport", OutMeasure::humanAC( 80, nCMDTReport, false )},
        /// Similar to nSevere. Number of severe episodes WITHOUT coinfection
        {"nSevereWithoutComorbidities", OutMeasure::humanAC( 81, nSevereWithoutComorbidities, false )},
        /** Similar to 'expectedSevere'.
        * Expected number of severe bouts of malaria WITHOUT "complications due
        * to coinfection" (the same as the `nSevereWithoutComorbidities` output). */
        {"expectedSevereWithoutComorbidities", OutMeasure::humanAC( 82, expectedSevereWithoutComorbidities, true )},
    };

    for (const NamedDef& d : defs) {
        auto [it, inserted] = namedOutMeasures.emplace(d.name, d.om);
        if (!inserted) {
            throw std::runtime_error("Duplicate OutMeasure name detected: " + std::string(d.name));
        }
        if (!seenIDs.insert(d.om.outId).second) {
            throw std::runtime_error("Duplicate OutMeasure (outId) detected: " + std::to_string(d.om.outId));
        }
    }

    // Now initialise valid condition measures:
    for( const NamedMeasureMapT::value_type& v : namedOutMeasures ){
        Measure m = v.second.m;
        // Not the following:
        if( m == mon::nSevere ||
            m == mon::nSevereWithoutComorbidities ||
            m == mon::nUncomp ||
            m == mon::nDirDeaths ||
            m == mon::nHospitalDeaths ||
            m == mon::Clinical_FirstDayDeaths ||
            m == mon::Clinical_HospitalFirstDayDeaths ||
            m == mon::nSeq ||
            m == mon::nHospitalSeqs ||
            m == mon::nHospitalRecovs ||
            m == mon::nNMFever ||
            m == mon::nNmfDeaths ||
            m == mon::nSubPopRemovalFirstEvent ||
            m == mon::innoculationsPerAgeGroup ||
            m == mon::inputEIR ||
            m == mon::simulatedEIR ) continue;
        validCondMeasures.insert(m);
    }
}

struct Condition {
    bool value; // whether the condition was satisfied during the last survey
    bool isDouble;
    Measure measure;
    uint8_t method;
    double min, max;
};

namespace impl {
    // Accumulators, variables:
    bool isInit = false;
    size_t surveyIndex = 0;     // index in surveyTimes of next survey
    size_t survNumEvent = NOT_USED, survNumStat = NOT_USED;
    SimTime nextSurveyDate = sim::future();
    
    vector<Condition> conditions;
}

/// One of these is used for every output index, and is specific to a measure
/// and repeated for every survey.
struct MonIndex {
    // Measure number used in output file
    int outMeasure;
    // Number of categories. Must be > 0. If 1, index is set to zero, otherwise
    // indices *should* be less than this.
    // 
    // nAges may include a final, unreported category.
    size_t nAges, nCohorts, nSpecies, nGenotypes, nDrugs;
    // Either Deploy::NA (not tracking deployments) or a binary 'or' of at
    // least one of Deploy::TIMED, Deploy::CTS, Deploy::TREAT.
    uint8_t deployMask;
    
    inline size_t size() const{
        return nAges * nCohorts * nSpecies * nGenotypes * nDrugs;
    }
    // Get the index in the result array to store this data at
    // (age group, cohort, species, genotype, drug).
    size_t index( size_t a, size_t c, size_t sp, size_t g, size_t d ) const{
#ifndef NDEBUG
        if( (nAges > 1 && a >= nAges) ||
            (nCohorts > 1 && c >= nCohorts) ||
            (nSpecies > 1 && sp >= nSpecies) ||
            (nGenotypes > 1 && g >= nGenotypes) ||
            (nDrugs > 1 && d >= nDrugs)
        ){
            cout << "Index out of bounds for age group\t" << a << " of " << nAges
                << "\ncohort set\t" << c << " of " << nCohorts
                << "\nspecies\t" << sp << " of " << nSpecies
                << "\ngenotype\t" << g << " of " << nGenotypes
                << "\ndrug\t" << d << " of " << nDrugs
                << endl;
        }
#endif
        // We use `a % nAges` etc. to enforce `a < nAges` and handle
        // the case `nAges == 1` (i.e. classification is turned off).
        return (d % nDrugs) + nDrugs *
            ((g % nGenotypes) + nGenotypes *
            ((sp % nSpecies) + nSpecies *
            ((c % nCohorts) + nCohorts *
            (a % nAges))));
    }
    
    // Write out some data from results.
    // 
    // @param stream Data sink
    // @param surveyNum Number to write in output (should start from 1 unlike in code)
    // @param results Vector of results
    // @param surveyStart Index in results where data for the current survey starts
    template<typename T>
    void write( ostream& stream, int surveyNum, const OutMeasure& om,
            const vector<T>& results, size_t surveyStart ) const
    {
        assert(results.size() >= surveyStart + size());
        // First age group starts at 1, unless there isn't an age group:
        const int ageGroupAdd = om.byAge ? 1 : 0;
        // Number of *reported* age categories: either no categorisation (1) or there is an extra unreported category
        const size_t nAgeCats = nAges == 1 ? 1 : nAges - 1;

        auto emit = [&](size_t ageGroup, size_t cohortSet, size_t species, size_t genotype, size_t drug, int col2) {
            const T value = results[surveyStart + index(ageGroup, cohortSet, species, genotype, drug)];
            stream << surveyNum << '\t' << col2 << '\t' << om.outId << '\t' << value << lineEnd;
        };

        if( om.bySpecies ){
            assert( nAges == 1 && nCohorts == 1 && nDrugs == 1 );
            for( size_t species = 0; species < nSpecies; ++species ){
            for( size_t genotype = 0; genotype < nGenotypes; ++genotype ){
                const int col2 = species + 1 +
                    1000000 * genotype;
                emit(0, 0, species, genotype, 0, col2);
            } }
            return;
        }

        if( om.byDrug ){
            assert( nSpecies == 1 && nGenotypes == 1 );
            for( size_t cohortSet = 0; cohortSet < nCohorts; ++cohortSet ){
            // Last age category is not reported
            for( size_t ageGroup = 0; ageGroup < nAgeCats; ++ageGroup ){
            for( size_t drug = 0; drug < nDrugs; ++drug ){
                // Yeah, >999 age groups clashes with cohort sets, but unlikely a real issue
                const int col2 = ageGroup + ageGroupAdd +
                    1000 * internal::cohortSetOutputId( cohortSet ) +
                    1000000 * (drug + 1);
                emit(ageGroup, cohortSet, 0, 0, drug, col2);
            } } }
            return;
        }

        assert( nSpecies == 1 && nDrugs == 1 );
        for( size_t cohortSet = 0; cohortSet < nCohorts; ++cohortSet ){
        // Last age category is not reported
        for( size_t ageGroup = 0; ageGroup < nAgeCats; ++ageGroup ){
        for( size_t genotype = 0; genotype < nGenotypes; ++genotype ){
            // Yeah, >999 age groups clashes with cohort sets, but unlikely a real issue
            const int col2 = ageGroup + ageGroupAdd +
                1000 * internal::cohortSetOutputId( cohortSet ) +
                1000000 * genotype;
            emit(ageGroup, cohortSet, 0, genotype, 0, col2);
        } } }
    }
};

template <typename T>
struct State {
    MonIndex layout;
    vector<T> reports;
};

vector<State<int>> intStates;
vector<State<double>> doubleStates;
vector<vector<size_t>> measureToIntStates;
vector<vector<size_t>> measureToDoubleStates;

template <typename State>
void fillStateLayout(State& state, const OutMeasure& om, size_t nSp, size_t nD, bool forceNoCategories)
{
    state.layout.outMeasure = om.outId;
    state.layout.nAges = forceNoCategories ? 1 : (om.byAge ? AgeGroup::numGroups() : 1);
    state.layout.nCohorts = forceNoCategories ? 1 : (om.byCohort ? impl::nCohorts : 1);
    state.layout.nSpecies = forceNoCategories ? 1 : (om.bySpecies ? nSp : 1);
    state.layout.nGenotypes = forceNoCategories ? 1 : (om.byGenotype ? WithinHost::Genotypes::N() : 1);
    state.layout.nDrugs = forceNoCategories ? 1 : (om.byDrug ? nD : 1);
    state.layout.deployMask = om.method;
}

template <typename T>
State<T> makeState(const OutMeasure& om, size_t nSp, size_t nD, bool forceNoCategories)
{
    State<T> state;
    fillStateLayout(state, om, nSp, nD, forceNoCategories);
    state.reports.assign(state.layout.size() * impl::nSurveys, T{});
    return state;
}

void addState(const OutMeasure& om, size_t nSp, size_t nD, bool forceNoCategories)
{
    assert(om.m < MeasureCount);
    if (om.isDouble) {
        size_t idx = doubleStates.size();
        doubleStates.push_back(makeState<double>(om, nSp, nD, forceNoCategories));
        measureToDoubleStates[om.m].push_back(idx);
    } else {
        size_t idx = intStates.size();
        intStates.push_back(makeState<int>(om, nSp, nD, forceNoCategories));
        measureToIntStates[om.m].push_back(idx);
    }
}

void initStates(const vector<OutMeasure>& enabledMeasures, size_t nSp, size_t nD)
{
    intStates.clear();
    doubleStates.clear();
    measureToIntStates.assign(MeasureCount, {});
    measureToDoubleStates.assign(MeasureCount, {});
    for (const OutMeasure& om : enabledMeasures)
    {
        if (om.m >= MeasureCount) continue;
        addState(om, nSp, nD, false);
    }
}

void ensureConditionState(const OutMeasure& om)
{
    assert(om.m < MeasureCount);
    auto hasMethodState = [&](const auto& states, const auto& measureToStates) {
        for (size_t idx : measureToStates[om.m]) {
            if (states[idx].layout.deployMask == om.method) return true;
        }
        return false;
    };
    if ((om.isDouble && hasMethodState(doubleStates, measureToDoubleStates)) ||
        (!om.isDouble && hasMethodState(intStates, measureToIntStates))) return;
    addState(om, 1, 1, true);
}

template <typename T>
inline void addValue(vector<T>& values, size_t index, T val)
{
    assert(index < values.size());
    values[index] += val;
}

template <typename State, typename T>
void recordValue(T val, vector<State>& states, const vector<vector<size_t>>& measureToStates,
                 Measure measure, size_t survey, size_t ageIndex, uint32_t cohortSet,
                 size_t species, size_t genotype, size_t drug, int outId = 0)
{
    if (survey == NOT_USED) return;
    assert(measure < measureToStates.size());
    for (size_t idx : measureToStates[measure])
    {
        State& state = states[idx];
        if (state.layout.deployMask != Deploy::NA) continue;
        if (outId != 0 && state.layout.outMeasure != outId) continue;
        size_t i = survey * state.layout.size() + state.layout.index(ageIndex, cohortSet, species, genotype, drug);
        addValue(state.reports, i, val);
    }
}

void recordDeployValue(int val, Measure measure, size_t survey, size_t ageIndex,
                       uint32_t cohortSet, Deploy::Method method)
{
    if (survey == NOT_USED) return;
    assert(method == Deploy::TIMED || method == Deploy::CTS || method == Deploy::TREAT);
    assert(measure < measureToIntStates.size());
    for (size_t idx : measureToIntStates[measure])
    {
        State<int>& state = intStates[idx];
        if ((state.layout.deployMask & method) == Deploy::NA) continue;
        assert(state.layout.nSpecies == 1 && state.layout.nGenotypes == 1);
        size_t i = survey * state.layout.size() + state.layout.index(ageIndex, cohortSet, 0, 0, 0);
        addValue(state.reports, i, val);
    }
}

template <typename T>
double sumStateMeasureTyped(Measure measure, uint8_t method, size_t survey,
                            const vector<State<T>>& states,
                            const vector<vector<size_t>>& measureToStates)
{
    assert(measure < measureToStates.size());
    for (size_t idx : measureToStates[measure])
    {
        const State<T>& state = states[idx];
        if (state.layout.deployMask != method) continue;
        const size_t off = survey * state.layout.size();
        const size_t end = off + state.layout.size();
        return std::accumulate(state.reports.begin() + off, state.reports.begin() + end, 0.0);
    }
    throw SWITCH_DEFAULT_EXCEPTION;
}

double sumStateMeasure(Measure measure, bool isDouble, uint8_t method, size_t survey)
{
    assert(survey != NOT_USED);
    return isDouble
        ? sumStateMeasureTyped(measure, method, survey, doubleStates, measureToDoubleStates)
        : sumStateMeasureTyped(measure, method, survey, intStates, measureToIntStates);
}

template <typename T>
void writeMeasureStateTyped(ostream& stream, size_t survey, const OutMeasure& om,
                            const vector<State<T>>& states,
                            const vector<vector<size_t>>& measureToStates)
{
    assert(om.m < measureToStates.size());
    for (size_t idx : measureToStates[om.m]) {
        const State<T>& state = states[idx];
        if (state.layout.outMeasure != om.outId) continue;
        state.layout.write(stream, survey + 1, om, state.reports, survey * state.layout.size());
        return;
    }
    assert(false && "measure not found in records");
}

void writeMeasureState(ostream& stream, size_t survey, const OutMeasure& om)
{
    if (om.isDouble)
        writeMeasureStateTyped(stream, survey, om, doubleStates, measureToDoubleStates);
    else
        writeMeasureStateTyped(stream, survey, om, intStates, measureToIntStates);
}

bool isMeasureUsed(Measure measure, bool isDouble)
{
    assert(measure < MeasureCount);
    const auto& measureToStates = isDouble ? measureToDoubleStates : measureToIntStates;
    return !measureToStates[measure].empty();
}

template <typename T>
void checkpointVec(ostream& stream, vector<T>& values, size_t /*expectedSize*/)
{
    values.size() & stream;
    for (T& y : values) y & stream;
}

template <typename T>
void checkpointVec(istream& stream, vector<T>& values, size_t expectedSize)
{
    size_t storedSize = 0;
    storedSize & stream;
    if (storedSize != expectedSize) {
        throw util::checkpoint_error("mon::reports: invalid list size");
    }
    values.resize(storedSize);
    for (T& y : values) y & stream;
}

template <typename Stream>
void checkpointStates(Stream& stream)
{
    for (State<int>& state : intStates)
        checkpointVec(stream, state.reports, state.layout.size() * impl::nSurveys);
    for (State<double>& state : doubleStates)
        checkpointVec(stream, state.reports, state.layout.size() * impl::nSurveys);
}

// Enabled measures:
vector<OutMeasure> reportedMeasures;
int reportIMR = -1; // special output for fitting

struct MeasureByOutId {
    bool operator()(const OutMeasure& i, const OutMeasure& j) { return i.outId < j.outId; }
} measureByOutId;

void initReporting( const scnXml::Scenario& scenario ){
    defineOutMeasures();        // set up namedOutMeasures
    assert(reportedMeasures.empty());
    
    // First we put used measures in this list:
    const scnXml::MonitoringOptions& optsElt = scenario.getMonitoring().getSurveyOptions();
    // This should be an upper bound on the number of options we need:
    reportedMeasures.reserve(optsElt.getOption().size() + namedOutMeasures.size());
    
    auto applyCategory = [](bool& supports, const bool optionalPresent, const bool requested,
                            const string& optionName, const char* label)
    {
        if (!optionalPresent) return;
        if (supports) {
            supports = requested;   // disable or keep
        } else if (requested) {
            throw util::xml_scenario_error("measure " + optionName + " does not support categorisation by " + label);
        }
    };

    set<int> outIds;    // all measure numbers used in output
    for( const scnXml::MonitoringOption& optElt : optsElt.getOption() ){
        if( optElt.getValue() == false ) continue;      // option is disabled
        
        auto it = namedOutMeasures.find( optElt.getName() );
        if( it == namedOutMeasures.end() ){
            throw util::xml_scenario_error("unrecognised survey option: " + string(optElt.getName()));
        }
        OutMeasure om = it->second;     // copy; we may modify below
        
        if( om.m >= MeasureCount ){
            if( om.m == allCauseIMR ){
                if( om.isDouble && !om.byAge && !om.byCohort && !om.bySpecies ){
                    reportIMR = om.outId;
                }else{
                    throw util::xml_scenario_error( "measure allCauseIMR does not support any categorisation" );
                }
            } else if( om.m == obsoleteMeasure ){
                throw util::xml_scenario_error("obsolete survey option: " + string(optElt.getName()));
            } else TRACED_EXCEPTION_DEFAULT("invalid measure code");
        }
        
        if( om.m == sumlogDens || om.m == logDensByGenotype){
            if( WithinHost::diagnostics::monitoringDiagnostic().allowsFalsePositives() ){
                throw util::xml_scenario_error("measure " + string(optElt.getName()) + " may not be used when monitoring diagnostic sensitivity < 1");
            }
        }
        
        // Categorisation can be disabled but not enabled.
        const string optionName(optElt.getName());
        const bool byAgePresent = optElt.getByAge().present();
        const bool byCohortPresent = optElt.getByCohort().present();
        const bool bySpeciesPresent = optElt.getBySpecies().present();
        const bool byGenotypePresent = optElt.getByGenotype().present();
        const bool byDrugPresent = optElt.getByDrugType().present();
        applyCategory(om.byAge, byAgePresent, byAgePresent ? optElt.getByAge().get() : false, optionName, "age group");
        applyCategory(om.byCohort, byCohortPresent, byCohortPresent ? optElt.getByCohort().get() : false, optionName, "cohort");
        applyCategory(om.bySpecies, bySpeciesPresent, bySpeciesPresent ? optElt.getBySpecies().get() : false, optionName, "species");
        applyCategory(om.byGenotype, byGenotypePresent, byGenotypePresent ? optElt.getByGenotype().get() : false, optionName, "genotype");
        applyCategory(om.byDrug, byDrugPresent, byDrugPresent ? optElt.getByDrugType().get() : false, optionName, "drug type");
        
        // Output number may be changed:
        if( optElt.getOutputNumber().present() ) om.outId = optElt.getOutputNumber().get();
        if( outIds.count(om.outId) ){
            throw util::xml_scenario_error("monitoring output number " + to_string(om.outId) + " used more than once");
        }
        outIds.insert( om.outId );
        
        reportedMeasures.push_back( om );
    }
    
    std::sort( reportedMeasures.begin(), reportedMeasures.end(), measureByOutId );
    
    size_t nSpecies = scenario.getEntomology().getVector().present() ?
        scenario.getEntomology().getVector().get().getAnopheles().size() : 1;
    size_t nDrugs = scenario.getPharmacology().present() ?
        scenario.getPharmacology().get().getDrugs().getDrug().size() : 1;
    
    initStates(reportedMeasures, nSpecies, nDrugs);
}

size_t setupCondition( const string& measureName, double minValue,
                       double maxValue, bool initialState )
{
    auto it = namedOutMeasures.find( measureName );
    if( it == namedOutMeasures.end() )
        throw util::xml_scenario_error("unrecognised measure: " + string(measureName));

    OutMeasure om = it->second;         // copy so that we can modify
    // Refuse to use some measures, since these are not reported reliably in
    // "non-reporting" surveys or are reported after the survey is taken:
    if( validCondMeasures.count(om.m) == 0 ){
        throw util::xml_scenario_error("cannot use measure " + string(measureName) + " as condition of deployment");
    }
    ensureConditionState(om);
    
    Condition condition;
    condition.value = initialState;
    condition.isDouble = om.isDouble;
    condition.measure = om.m;
    condition.method = om.method;
    condition.min = minValue;
    condition.max = maxValue;
    impl::conditions.push_back(condition);
    return impl::conditions.size() - 1;
}

void updateConditions() {
    for( Condition& cond : impl::conditions ){
        double val = sumStateMeasure(cond.measure, cond.isDouble, cond.method, impl::survNumStat);
        cond.value = (val >= cond.min && val <= cond.max);
    }
}
bool checkCondition( size_t conditionKey ){
    assert( conditionKey < impl::conditions.size() );
    return impl::conditions[conditionKey].value;
}

void internal::write( ostream& stream ){
    for( size_t survey = 0; survey < impl::nSurveys; ++survey ){
        for( const OutMeasure& om : reportedMeasures ){
            if( om.m >= MeasureCount ){
                // "Special" measures are not reported this way. The only such measure is IMR.
                assert( om.m == allCauseIMR && reportIMR >= 0 );
                continue;
            }
            writeMeasureState(stream, survey, om);
        }
    }
    if( reportIMR >= 0 ){
        // Infant mortality rate is a single number, therefore treated specially.
        // It is calculated across the entire intervention period and used in
        // model fitting.
        stream << 1 << "\t" << 1 << "\t" << reportIMR
            << "\t" << Clinical::InfantMortality::allCause() << lineEnd;
    }
}

void record(Measure measure, size_t survey, size_t age, uint32_t cohort, size_t species, size_t genotype, size_t drug, int val, int outId)
{
    recordValue(val, intStates, measureToIntStates, measure, survey, age, cohort, species, genotype, drug, outId);
}

void record(Measure measure, size_t survey, size_t age, uint32_t cohort, size_t species, size_t genotype, size_t drug, double val)
{
    recordValue(val, doubleStates, measureToDoubleStates, measure, survey, age, cohort, species, genotype, drug);
}

void recordStat(Measure measure, const Host::Human& human, int val, size_t species, size_t genotype, size_t drug, int outId)
{
    record(measure, statSurveyNumber(), human.monitoringAgeGroup.i(), human.getCohortSet(), species, genotype, drug, val, outId);
}

void recordStat(Measure measure, const Host::Human& human, double val, size_t species, size_t genotype, size_t drug)
{
    record(measure, statSurveyNumber(), human.monitoringAgeGroup.i(), human.getCohortSet(), species, genotype, drug, val);
}

void recordEvent(Measure measure, const Host::Human& human, int val)
{
    record(measure, eventSurveyNumber(), human.monitoringAgeGroup.i(), human.getCohortSet(), 0, 0, 0, val);
}

void recordDeploy(Measure measure, const Host::Human& human, Deploy::Method method, int val)
{
    recordDeployValue(val, measure, eventSurveyNumber(), human.monitoringAgeGroup.i(), human.getCohortSet(), method);
    if (measure != nTreatDeployments)
        recordDeployValue(val, nTreatDeployments, eventSurveyNumber(), human.monitoringAgeGroup.i(), human.getCohortSet(), method);
}

bool isUsed(Measure measure) { return isMeasureUsed(measure, false) || isMeasureUsed(measure, true); }

template <typename Stream>
void checkpoint(Stream& stream)
{
    impl::isInit & stream;
    impl::surveyIndex & stream;
    impl::survNumEvent & stream;
    impl::survNumStat & stream;
    impl::nextSurveyDate & stream;

    checkpointStates(stream);
}
template void checkpoint<ostream>(ostream& stream);
template void checkpoint<istream>(istream& stream);

// ———  surveys  ———

struct SurveyDate {
    SimTime date = sim::never();       // date of survey
    size_t num; // if NOT_USED, the survey is not reported; if greater, this is the survey number
    
    /// Construct
    SurveyDate(SimTime date, size_t num) : date(date), num(num) {}
    
    inline bool isReported() const { return num != NOT_USED; }
};

namespace impl{
    // Constants or defined during init:
    size_t nSurveys = 0;        // number of reported surveys
    size_t nCohorts = 1;     // default: just the whole population
    extern size_t surveyIndex;     // index in surveyDates of next survey
    vector<SurveyDate> surveyDates;     // dates of surveys
}

// trim from start (in place)
static inline void ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {
        return !std::isspace(ch);
    }));
}

// trim from end (in place)
static inline void rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
        return !std::isspace(ch);
    }).base(), s.end());
}

// trim from both ends (in place)
static inline void trim(std::string &s) {
    ltrim(s);
    rtrim(s);
}

SimTime readSurveyDates( const scnXml::Monitoring& monitoring ){
    const scnXml::Surveys::SurveyTimeSequence& survs =
        monitoring.getSurveys().getSurveyTime();
    
    map<SimTime, bool> surveys;        // dates of all surveys (from XML) and whether these are reporting
    
    auto addSurvey = [&surveys](SimTime date, bool reporting) {
        if( reporting ) surveys[date] = true;
        else surveys.insert(make_pair(date, false));   // do not override existing reported entries
    };

    for(size_t i = 0; i < survs.size(); ++i) {
        const scnXml::SurveyTime& surv = survs[i];
        try{
            std::string s = surv;
            trim(s);
            SimTime cur = UnitParse::readDate( s, UnitParse::STEPS );
            bool reporting = surv.getReported();
            if( surv.getRepeatStep().present() != surv.getRepeatEnd().present() ){
                throw util::xml_scenario_error( "surveyTime: use of repeatStep or repeatEnd without other" );
            }
            if( surv.getRepeatStep().present() ){
                SimTime step = UnitParse::readDuration( surv.getRepeatStep().get(), UnitParse::NONE );
                if( step < sim::oneTS() ){
                    throw util::xml_scenario_error( "surveyTime: repeatStep must be >= 1" );
                }
                SimTime end = UnitParse::readDate( surv.getRepeatEnd().get(), UnitParse::NONE );
                while(cur < end){
                    addSurvey(cur, reporting);
                    cur = cur + step;
                }
            }else{
                addSurvey(cur, reporting);
            }
        }catch( const util::format_error& e ){
            throw util::xml_scenario_error( string("surveyTime: ").append(e.message()) );
        }
    }
    
    impl::surveyDates.clear();
    impl::surveyDates.reserve(surveys.size());
    size_t n = 0;
    for( auto it = surveys.begin(); it != surveys.end(); ++it ){
        size_t num = NOT_USED;
        if( it->second ){
            num = n;
            n += 1;
        }
        impl::surveyDates.push_back(SurveyDate(it->first, num));
    }
    impl::nSurveys = n;
    
    if( impl::surveyDates.size() == 0 ){
        throw util::xml_scenario_error( "Scenario defines no surveys; at least one is required." );
    }
    if( !impl::surveyDates.back().isReported() ){
        std::cerr << "Warning: the last survey is unreported. Having surveys beyond the last reported survey is pointless." << std::endl;
    }
    
    if( util::CommandLine::option( util::CommandLine::PRINT_SURVEY_TIMES ) ){
        std::cout << "Survey\tsteps\tdate";
        std::cout << std::endl;
        for( size_t i = 0; i < impl::surveyDates.size(); ++i ){
            const SurveyDate& surveyDate = impl::surveyDates[i];
            if( !surveyDate.isReported() ) continue;
            std::cout
                << (surveyDate.num+1) << '\t'
                << sim::inSteps(surveyDate.date - sim::startDate()) << '\t'
                << surveyDate.date << std::endl;
        }
    }
    
    if( monitoring.getCohorts().present() ){
        // this needs to be set early, but we can't set cohortSubPopIds until after InterventionManager is initialised
        impl::nCohorts = static_cast<uint32_t>(1) << monitoring.getCohorts().get().getSubPop().size();
    }
    
    mon::AgeGroup::init( monitoring );
    
    // final survey date:
    return impl::surveyDates[impl::surveyDates.size()-1].date;
}

void updateSurveyNumbers() {
    if( impl::surveyIndex >= impl::surveyDates.size() ){
        impl::survNumEvent = NOT_USED;
        impl::survNumStat = NOT_USED;
        impl::nextSurveyDate = sim::future();
    }else{
        for( size_t i = impl::surveyIndex; i < impl::surveyDates.size(); ++i ){
            impl::survNumEvent = impl::surveyDates[i].num;  // set to survey number or NOT_USED; this happens at least once!
            if( impl::survNumEvent != NOT_USED ) break;        // stop at first reported survey
        }
        const SurveyDate& nextSurvey = impl::surveyDates[impl::surveyIndex];
        impl::survNumStat = nextSurvey.num;     // may be NOT_USED; this is intended
        impl::nextSurveyDate = nextSurvey.date;
    }
}
void initMainSim(){
    impl::surveyIndex = 0;
    impl::isInit = true;
    updateSurveyNumbers();
}
void concludeSurvey(){
    updateConditions();
    impl::surveyIndex += 1;
    updateSurveyNumbers();
}

void writeToStream(ostream& stream) {
    stream.width(0);
    internal::write(stream);
}

void writeSurveyData ()
{
    string filename = util::CommandLine::getOutputName();
    auto mode = std::ios::out | std::ios::binary;
    
    if (util::CommandLine::option( util::CommandLine::COMPRESS_OUTPUT )) {
        filename.append(".gz");
        ogzstream stream(filename.c_str(), mode);
        writeToStream(stream);
    } else {
        ofstream stream(filename, mode);
        writeToStream(stream);
        // Otherwise file may be written after OpenMalaria has returned (Mac OS Xcode 9.4)
        stream.flush();
        stream.close();
    }

    ifstream stream(filename, mode);
    if(stream.is_open() == false || !stream.good())
    {
        cerr << "STREAM BAD" << endl;
    }
}


// ———  AgeGroup  ———

vector<SimTime> AgeGroup::upperBound;

void AgeGroup::init (const scnXml::Monitoring& monitoring) {
    const scnXml::MonAgeGroup::GroupSequence& groups =
        monitoring.getAgeGroup().getGroup();
    if (!(monitoring.getAgeGroup().getLowerbound() <= 0.0))
        throw util::xml_scenario_error ("Expected survey age-group lowerbound of 0");
    
    // The last age group includes individuals too old for reporting
    upperBound.resize( groups.size() + 1 );
    for(size_t i = 0;i < groups.size(); ++i) {
        // convert to SimTime, rounding down to the next time step
        upperBound[i] = sim::fromYearsD( groups[i].getUpperbound() );
    }
    upperBound[groups.size()] = sim::future();
}

void AgeGroup::update (SimTime age) {
    while (age >= upperBound[index]){
        ++index;
    }
}


// ———  Cohort sets  ———

using interventions::ComponentId;
vector<uint32_t> cohortSubPopNumbers;   // value is output number
map<ComponentId,uint32_t> cohortSubPopIds;      // value is internal index (used above)

bool notPowerOfTwo( uint32_t num ){
    return num == 0 || num > (static_cast<uint32_t>(1) << 21) || (num & (num - 1)) != 0;
}
// Init cohort sets. Depends on interventions (initialise those first).
void initCohorts( const scnXml::Monitoring& monitoring )
{
    if( monitoring.getCohorts().present() ){
        const scnXml::Cohorts monCohorts = monitoring.getCohorts().get();
        uint32_t nextId = 0;
        for( auto it = monCohorts.getSubPop().begin(),
            end = monCohorts.getSubPop().end(); it != end; ++it )
        {
            ComponentId compId = interventions::InterventionManager::getComponentId( it->getId() );
            bool inserted = cohortSubPopIds.insert( make_pair(compId,nextId) ).second;
            if( !inserted ){
                throw util::xml_scenario_error(
                    string("cohort specification uses sub-population \"").append(it->getId())
                    .append("\" more than once") );
            }
            if( it->getNumber() < 0 || notPowerOfTwo( it->getNumber() ) ){
                throw util::xml_scenario_error(
                    string( "cohort specification assigns sub-population \"").append(it->getId())
                    .append("\" a number which is not a power of 2 (up to 2^21)") );
            }
            cohortSubPopNumbers.push_back( it->getNumber() );
            nextId += 1;
        }
    }
}

uint32_t updateCohortSet( uint32_t old, ComponentId subPop, bool isMember ){
    auto it = cohortSubPopIds.find( subPop );
    if( it == cohortSubPopIds.end() ) return old;       // sub-pop not used in cohorts
    uint32_t subPopId = static_cast<uint32_t>(1) << it->second;        // 1 bit positive
    return (old & ~subPopId) | (isMember ? subPopId : 0);
}

uint32_t internal::cohortSetOutputId(uint32_t cohortSet){
    uint32_t outNum = 0;
    assert( (cohortSet >> cohortSubPopNumbers.size()) == 0 );
    for( uint32_t i = 0; i < cohortSubPopNumbers.size(); ++i ){
        if( cohortSet & (static_cast<uint32_t>(1) << i) ){
            outNum += cohortSubPopNumbers[i];
        }
    }
    return outNum;
}

} }
