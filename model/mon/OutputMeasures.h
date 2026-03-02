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

namespace OM {
namespace mon {

void findNamedMeasuresUsing( Measure m, ostream& msg ){
    int nMatches = 0;
    for( auto it = namedOutMeasures.begin(); it != namedOutMeasures.end(); ++it ){
        if( it->second.m == m ){
            if( nMatches > 0 ) msg << ", ";
            msg << it->first;
            nMatches += 1;
        }
    }
    if( nMatches == 0 ) msg << "??";
}

// This method defines output measures accepted by name in the XML (e.g.
// "nHost") and their numeric output identifier (i.e. measure column of
// outputs), type of output (integer or floating point), aggregation, and the
// corresponding internal measure code.
void defineOutMeasures(){
    /* Don't ever make an existing numerical identifier point to a new/different 
    output measure because this would violate users' expectations */
    // Don't reuse old numerical identifiers

    /// Total number of humans
    namedOutMeasures["nHost"] = OutMeasure::humanAC( 0, nHost, false );
    /** The number of human hosts with an infection (patent or not) at the time
     * the survey is taken. */
    namedOutMeasures["nInfect"] = OutMeasure::humanAC( 1, nInfect, false );
    namedOutMeasures["nInfect_Imported"] = OutMeasure::humanAC( 1001, nInfect_Imported, false );
    namedOutMeasures["nInfect_Introduced"] = OutMeasure::humanAC( 2001, nInfect_Introduced, false );
    namedOutMeasures["nInfect_Indigenous"] = OutMeasure::humanAC( 3001, nInfect_Indigenous, false );
    /** Expected number of infected hosts
     * 
     * This is the sum of the probabilities, across all time steps since the
     * last survey, of each host becoming infected on that time step. */
    namedOutMeasures["nExpectd"] = OutMeasure::humanAC( 2, nExpectd, true );
    /** The number of human hosts whose total (blood-stage) parasite density is
     * above the detection threshold */
    namedOutMeasures["nPatent"] = OutMeasure::humanAC( 3, nPatent, false );
    namedOutMeasures["nPatent_Imported"] = OutMeasure::humanAC( 1003, nPatent_Imported, false );
    namedOutMeasures["nPatent_Introduced"] = OutMeasure::humanAC( 2003, nPatent_Introduced, false );
    namedOutMeasures["nPatent_Indigenous"] = OutMeasure::humanAC( 3003, nPatent_Indigenous, false );
    /// Sum of log(1 + p) where p is the pyrogenic threshold
    namedOutMeasures["sumLogPyrogenThres"] =
        OutMeasure::humanAC( 4, sumLogPyrogenThres, true );
    /** Sum (across hosts) of the natural logarithm of the parasite density of
     * hosts with detectable parasite density (patent according to the
     * monitoring diagnostic). */
    namedOutMeasures["sumlogDens"] = OutMeasure::humanAC( 5, sumlogDens, true );
    /** The total number of infections in the population: includes both blood
     * and liver stages. Vivax: this is the number of broods. */
    namedOutMeasures["totalInfs"] = OutMeasure::humanACG( 6, totalInfs, false );
        namedOutMeasures["totalInfs_Imported"] = OutMeasure::humanACG( 1006, totalInfs_Imported, false );
        namedOutMeasures["totalInfs_Introduced"] = OutMeasure::humanACG( 2006, totalInfs_Introduced, false );
        namedOutMeasures["totalInfs_Indigenous"] = OutMeasure::humanACG( 3006, totalInfs_Indigenous, false );
    /** Infectiousness of human population to mosquitoes
     *
     * Number of hosts transmitting to mosquitoes (i.e. proportion of
     * mosquitoes that get infected multiplied by human population size).
     * Single value, not per age-group. */
    namedOutMeasures["nTransmit"] = OutMeasure::value( 7, nTransmit, true );
    /** The sum of all detectable infections (where blood stage parasite
     * density is above the detection limit) across all human hosts.
     * Vivax: the number of broods with an active blood stage. */
    namedOutMeasures["totalPatentInf"] = OutMeasure::humanACG( 8, totalPatentInf, false );
        namedOutMeasures["totalPatentInf_Imported"] = OutMeasure::humanACG( 1008, totalPatentInf_Imported, false );
        namedOutMeasures["totalPatentInf_Introduced"] = OutMeasure::humanACG( 2008, totalPatentInf_Introduced, false );
        namedOutMeasures["totalPatentInf_Indigenous"] = OutMeasure::humanACG( 3008, totalPatentInf_Indigenous, false );
    /// Contribuion to immunity functions (removed)
    namedOutMeasures["contrib"] = OutMeasure::obsolete( 9 );
    /// Sum of the pyrogenic threshold
    namedOutMeasures["sumPyrogenThresh"] =
        OutMeasure::humanAC( 10, sumPyrogenThresh, true );
    /// number of blood-stage treatments (1st line)
    namedOutMeasures["nTreatments1"] = OutMeasure::humanAC( 11, nTreatments1, false );
    /// number of blood-stage treatments (2nd line)
    namedOutMeasures["nTreatments2"] = OutMeasure::humanAC( 12, nTreatments2, false );
    /// number of blood-stage treatments (inpatient)
    namedOutMeasures["nTreatments3"] = OutMeasure::humanAC( 13, nTreatments3, false );
    /// number of episodes (uncomplicated)
    namedOutMeasures["nUncomp"] = OutMeasure::humanAC( 14, nUncomp, false );
        namedOutMeasures["nUncomp_Imported"] = OutMeasure::humanAC( 1014, nUncomp_Imported, false );
        namedOutMeasures["nUncomp_Introduced"] = OutMeasure::humanAC( 2014, nUncomp_Introduced, false );
        namedOutMeasures["nUncomp_Indigenous"] = OutMeasure::humanAC( 3014, nUncomp_Indigenous, false );
    /// Number of severe episodes (severe malaria or malaria + coinfection)
    namedOutMeasures["nSevere"] =
        OutMeasure::humanAC( 15, nSevere, false );
    /// cases with sequelae
    namedOutMeasures["nSeq"] = OutMeasure::humanAC( 16, nSeq, false );
    /// deaths in hospital
    namedOutMeasures["nHospitalDeaths"] =
        OutMeasure::humanAC( 17, nHospitalDeaths, false );
    /// Number of deaths indirectly caused by malaria
    namedOutMeasures["nIndDeaths"] =
        OutMeasure::humanAC( 18, nIndDeaths, false );
    /// Number of deaths directly caused by malaria
    namedOutMeasures["nDirDeaths"] =
        OutMeasure::humanAC( 19, nDirDeaths, false );
    /** Number of vaccine doses given via EPI.
     * 
     * Since schema 22, each vaccine type may be deployed independently. To be
     * roughly backwards-compatible, the first type (PEV, BSV or TBV) described
     * (with an "effect" element) will be reported. */
    namedOutMeasures["nEPIVaccinations"] =
        OutMeasure::humanDeploy( 20, vaccinations, Deploy::CTS );
    /** All cause infant mortality rate
     * 
     * Reports death rate of infants due to all causes (malaria as modelled
     * plus fixed non-malaria attribution). Calculated via Kaplan-Meier method.
     * Units: deaths per thousand births. */
    namedOutMeasures["allCauseIMR"] =
        OutMeasure::value( 21, allCauseIMR, true );
    /** Number of vaccine doses given via mass campaign.
     * 
     * Since schema 22, each vaccine type may be deployed independently. To be
     * roughly backwards-compatible, the first type (PEV, BSV or TBV) described
     * (with an "effect" element) will be reported. */
    namedOutMeasures["nMassVaccinations"] =
        OutMeasure::humanDeploy( 22, vaccinations, Deploy::TIMED );
    /// recoveries in hospital
    namedOutMeasures["nHospitalRecovs"] =
        OutMeasure::humanAC( 23, nHospitalRecovs, false );
    /// sequelae in hospital
    namedOutMeasures["nHospitalSeqs"] =
        OutMeasure::humanAC( 24, nHospitalSeqs, false );
    /// Number of IPT Doses (removed together with IPT model)
    namedOutMeasures["nIPTDoses"] = OutMeasure::obsolete( 25 );
    /** Annual Average Kappa
     *
     * Calculated once a year as sum of human infectiousness divided by initial
     * EIR summed over a year. Single value, not per age-group. */
    namedOutMeasures["annAvgK"] = OutMeasure::value( 26, annAvgK, true );
    /// Number of episodes (non-malaria fever)
    namedOutMeasures["nNMFever"] =
        OutMeasure::humanAC( 27, nNMFever, false );
    /// Inoculations per human (all ages) per day of year, over the last year.
    /// (Reporting removed.)
    namedOutMeasures["innoculationsPerDayOfYear"] = OutMeasure::obsolete( 28 );
    /// Kappa (human infectiousness) weighted by availability per day-of-year for the last year.
    /// (Reporting removed.)
    namedOutMeasures["kappaPerDayOfYear"] = OutMeasure::obsolete( 29 );
    /** The total number of inoculations, by age group, cohort and parasite
     * genotype, summed over the reporting period. */
    namedOutMeasures["innoculationsPerAgeGroup"] =
        OutMeasure::humanACG( 30, innoculationsPerAgeGroup, true );
    /// N_v0: emergence of feeding vectors during the last time step. Units: mosquitoes/day
    namedOutMeasures["Vector_Nv0"] = OutMeasure::species( 31, Vector_Nv0, false );
    /// N_v: vectors seeking to feed during the last time step. Units: mosquitoes/day
    namedOutMeasures["Vector_Nv"] = OutMeasure::species( 32, Vector_Nv, false );
    /// N_v: infected vectors seeking to feed during the last time step. Units: mosquitoes/day
    namedOutMeasures["Vector_Ov"] = OutMeasure::species( 33, Vector_Ov, true );
    /// N_v: infectious vectors seeking to feed during the last time step. Units: mosquitoes/day
    namedOutMeasures["Vector_Sv"] = OutMeasure::species( 34, Vector_Sv, true );
    /** Input EIR (Expected EIR entered into scenario file)
     *
     * Units: inoculations per adult per time step. */
    namedOutMeasures["inputEIR"] = OutMeasure::value( 35, inputEIR, true );
    /** Simulated EIR (EIR output by the transmission model)
     *
     * Units: inoculations per adult per time step (children are excluded
     * when measuring). */
    namedOutMeasures["simulatedEIR"] = OutMeasure::value( 36, simulatedEIR, true );
    namedOutMeasures["simulatedEIR_Introduced"] = OutMeasure::value( 2036, simulatedEIR_Introduced, true );
    namedOutMeasures["simulatedEIR_Indigenous"] = OutMeasure::value( 3036, simulatedEIR_Indigenous, true );
    /// Number of Rapid Diagnostic Tests used
    namedOutMeasures["Clinical_RDTs"] = OutMeasure::obsolete( 39 );
    /* Effective total quanty of each drug used orally, in mg.
     * (Per active ingredient abbreviation.)
     * 
     * The quantity is efffective with respect to the cost (see treatment
     * schedule definition).
     * 
     * Reporting removed. */
    namedOutMeasures["Clinical_DrugUsage"] = OutMeasure::obsolete( 40 );
    /// Direct death on first day of CM (before treatment takes effect)
    namedOutMeasures["Clinical_FirstDayDeaths"] =
        OutMeasure::humanAC( 41, Clinical_FirstDayDeaths, false );
    /// Direct death on first day of CM (before treatment takes effect); hospital only
    namedOutMeasures["Clinical_HospitalFirstDayDeaths"] =
        OutMeasure::humanAC( 42, Clinical_HospitalFirstDayDeaths, false );
    /** The number of actual infections since the last survey. */
    namedOutMeasures["nNewInfections"] = OutMeasure::humanAC( 43, nNewInfections, false );
        namedOutMeasures["nNewInfections_Imported"] = OutMeasure::humanAC( 1043, nNewInfections_Imported, false );
        namedOutMeasures["nNewInfections_Introduced"] = OutMeasure::humanAC( 2043, nNewInfections_Introduced, false );
        namedOutMeasures["nNewInfections_Indigenous"] = OutMeasure::humanAC( 3043, nNewInfections_Indigenous, false );
        
    /** The number of ITNs delivered by mass distribution since last survey.
     *
     * These are "modelled ITNs": cover only a single person, cannot be passed
     * to someone else for reuse or used for fishing, etc. */
    namedOutMeasures["nMassITNs"] =
        OutMeasure::humanDeploy( 44, itn, Deploy::TIMED );
    /** The number of ITNs delivered through EPI since last survey.
     *
     * Comments from nMassITNs apply. */
    namedOutMeasures["nEPI_ITNs"] =
        OutMeasure::humanDeploy( 45, itn, Deploy::CTS );
    /** The number of people newly protected by IRS since last survey.
     *
     * Modelled IRS: affects one person, cannot be plastered over. */
    namedOutMeasures["nMassIRS"] =
        OutMeasure::humanDeploy( 46, irs, Deploy::TIMED );
    /** Defunct; was used by "vector availability" intervention (which is now a
     * sub-set of GVI). */
    namedOutMeasures["nMassVA"] = OutMeasure::obsolete( 47 );
    /// Number of malarial tests via microscopy used
    namedOutMeasures["Clinical_Microscopy"] = OutMeasure::obsolete( 48 );
    /* As Clinical_DrugUsage, but for quatities of drug delivered via IV. */
    namedOutMeasures["Clinical_DrugUsageIV"] = OutMeasure::obsolete( 49 );
    /// Number of cohort recruitments removed)
    namedOutMeasures["nAddedToCohort"] = OutMeasure::obsolete( 50 );
    /// Number of individuals removed from cohort (removed)
    namedOutMeasures["nRemovedFromCohort"] = OutMeasure::obsolete( 51 );
    /** Number of people (per age group) treated by mass drug administration
     * campaign. (Note that in one day time-step model MDA can be configured
     * as screen-and-treat. This option reports treatments administered not
     * the number of tests used.) */
    namedOutMeasures["nMDAs"] =
        OutMeasure::humanDeploy( 52, treat, Deploy::TIMED );
    /// Number of deaths caused by non-malaria fevers
    namedOutMeasures["nNmfDeaths"] = OutMeasure::humanAC( 53, nNmfDeaths, false );
    /// Number of antibiotic treatments given (disabled — not used)
    namedOutMeasures["nAntibioticTreatments"] = OutMeasure::obsolete( 54 );
    /** Report the number of screenings used in a mass screen-and-treat
     * operation. */
    namedOutMeasures["nMassScreenings"] =
        OutMeasure::humanDeploy( 55, screen, Deploy::TIMED );
    /// Report the number of mass deployments of generic vector interventions.
    namedOutMeasures["nMassGVI"] =
        OutMeasure::humanDeploy( 56, gvi, Deploy::TIMED );
    /** Number of IRS deployments via continuous deployment. */
    namedOutMeasures["nCtsIRS"] =
        OutMeasure::humanDeploy( 57, irs, Deploy::CTS );
    /** Number of GVI deployments via continuous deployment. */
    namedOutMeasures["nCtsGVI"] =
        OutMeasure::humanDeploy( 58, gvi, Deploy::CTS );
    /** Number of "MDA" deployments via continuous deployment.
     * 
     * Note: MDA stands for mass drug administration, but the term has come to
     * be used more flexibly by OpenMalaria, including optional screening and
     * deployment through age-based systems. */
    namedOutMeasures["nCtsMDA"] =
        OutMeasure::humanDeploy( 59, treat, Deploy::CTS );
    /** Number of diagnostics used by "MDA" distribution through continuous
     * methods. Can be higher than nCtsMDA since drugs are administered only
     * when the diagnostic is positive. Also see nCtsMDA description. */
    namedOutMeasures["nCtsScreenings"] =
        OutMeasure::humanDeploy( 60, screen, Deploy::CTS );
    /** Number of removals from a sub-population due to expiry of duration of
     * membership (e.g. intervention too old). */
    namedOutMeasures["nSubPopRemovalTooOld"] =
        OutMeasure::humanAC( 61, nSubPopRemovalTooOld, false );
    /** Number of removals from a sub-population due to first
     * infection/bout/treatment (see onFirstBout & co). */
    namedOutMeasures["nSubPopRemovalFirstEvent"] =
        OutMeasure::humanAC( 62, nSubPopRemovalFirstEvent, false );
    /** Report the number of liver-stage treatments (likely Primaquine) administered. */
    namedOutMeasures["nLiverStageTreatments"] =
        OutMeasure::humanAC( 63, nLiverStageTreatments, false );
    /** Report the number of diagnostics used during treatment.
     * 
     * This is not the same as Clinical_RDTs + Clinical_Microscopy: those
     * outputs are used by the "event scheduler" 1-day time step clinical
     * model, whereas this output is used by the 5-day time step model. */
    namedOutMeasures["nTreatDiagnostics"] =
        OutMeasure::humanAC( 64, nTreatDiagnostics, false );
    /** Number of "recruitment only" recruitments via timed deployment. */
    namedOutMeasures["nMassRecruitOnly"] =
        OutMeasure::humanDeploy( 65, recruit, Deploy::TIMED );
    /** Number of "recruitment only" recruitments via age-based deployment. */
    namedOutMeasures["nCtsRecruitOnly"] =
        OutMeasure::humanDeploy( 66, recruit, Deploy::CTS );
    /** Number of deployments (of all intervention components) triggered by
     * treatment (case management). */
    namedOutMeasures["nTreatDeployments"] =
        OutMeasure::humanDeploy( 67, nTreatDeployments, Deploy::TREAT );
    /** Report the total age of all humans in this a group (sum across humans,
     * in years). Divide by nHost to get the average age. */
    namedOutMeasures["sumAge"] = OutMeasure::humanAC( 68, sumAge, true );
    /** The number of human hosts with an infection (patent or not), for each
     * genotype, at the time the survey is taken. */
    namedOutMeasures["nInfectByGenotype"] =
        OutMeasure::humanACG( 69, nInfectByGenotype, false );
    /** The number of human hosts whose total (blood-stage) parasite density,
     * for each genotype, is above the detection threshold */
    namedOutMeasures["nPatentByGenotype"] =
        OutMeasure::humanACG( 70, nPatentByGenotype, false );
    /** For each infection genotype, sum across humans the natural log of
     * parasite density (like sumlogDens but per genotype). */
    namedOutMeasures["logDensByGenotype"] =
        OutMeasure::humanACG( 71, logDensByGenotype, true );
    /** For each drug type in the pharmacology section of the XML, report the
     * number of humans with non-zero concentration of this drug in their
     * blood. */
    namedOutMeasures["nHostDrugConcNonZero"] = 
        OutMeasure::humanACP( 72, nHostDrugConcNonZero, false );
    /** For each drug type in the pharmacology section of the XML, report the
     * sum of the natural logarithm of the drug concentration in hosts with
     * non-zero concentration. */
    namedOutMeasures["sumLogDrugConcNonZero"] =
        OutMeasure::humanACP( 73, sumLogDrugConcNonZero, true );
    /** Expected number of direct malaria deaths, from those with severe
     * disease.
     *
     * This is calculated as the sum over all steps in the reporting period of
     * the sum over humans with severe malaria of the probability of direct
     * death from malaria. */
    namedOutMeasures["expectedDirectDeaths"] =
        OutMeasure::humanAC( 74, expectedDirectDeaths, true );
    /** Expected number of direct malaria deaths which occur in hospital.
     * 
     * This is the a subset of `expectedDirectDeaths` and the same notes apply.
     */
    namedOutMeasures["expectedHospitalDeaths"] =
        OutMeasure::humanAC( 75, expectedHospitalDeaths, true );
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
    namedOutMeasures["expectedIndirectDeaths"] =
        OutMeasure::humanAC( 76, expectedIndirectDeaths, true );
    /** Expected number of sequelae, from those with severe disease.
     * 
     * This is calculated as the sum over all steps in the reporting period of
     * the sum over humans with severe malaria of the probability of sequelae
     * occuring, assuming the human "recovers" from the bout.
     */
    namedOutMeasures["expectedSequelae"] =
        OutMeasure::humanAC( 77, expectedSequelae, true );
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
    namedOutMeasures["expectedSevere"] =
        OutMeasure::humanAC( 78, expectedSevere, true );
    /** The total number of inoculations, by mosquito species, summed over
     * the reporting period. */
    namedOutMeasures["innoculationsPerVector"] =
        OutMeasure::species( 79, innoculationsPerAgeGroup, false );

    /** Number of custom intervention reports done */
    namedOutMeasures["nCMDTReport"] =
        OutMeasure::humanAC( 80, nCMDTReport, false );
        
    /// Similar to nSevere. Number of severe episodes WITHOUT coinfection
    namedOutMeasures["nSevereWithoutComorbidities"] =
        OutMeasure::humanAC( 81, nSevereWithoutComorbidities, false );
    /** Similar to 'expectedSevere'.
     * Expected number of severe bouts of malaria WITHOUT "complications due 
     * to coinfection" (the same as the `nSevereWithoutComorbidities` output). */
    namedOutMeasures["expectedSevereWithoutComorbidities"] =
        OutMeasure::humanAC( 82, expectedSevereWithoutComorbidities, true );

    // Check for duplicate outId
    std::set<int> seenIDs;
    for (const auto& pair : namedOutMeasures) {
        int id = pair.second.outId;
        if (!seenIDs.insert(id).second)
            throw std::runtime_error("Duplicate OutMeasure (outId) detected: " + std::to_string(id));
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

}
}
