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

#ifndef H_OM_mon_OutMeasures
#define H_OM_mon_OutMeasures

#include <cstdint>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <string_view>

namespace OM {
namespace mon {

using Measure = uint16_t;
inline constexpr Measure invalidMeasure = std::numeric_limits<Measure>::max() - 2;
inline constexpr Measure obsoleteMeasure = std::numeric_limits<Measure>::max() - 1;
inline constexpr Measure allCauseIMR = std::numeric_limits<Measure>::max();

namespace Deploy {
enum Method {
    NA = 0,
    TIMED = 1 << 0,
    CTS = 1 << 1,
    TREAT = 1 << 2
};
}

enum class Dim : uint8_t {
    None = 0,
    Age = 1 << 0,
    Cohort = 1 << 1,
    Species = 1 << 2,
    Genotype = 1 << 3,
    Drug = 1 << 4
};
constexpr Dim operator|(Dim a, Dim b) { return static_cast<Dim>(static_cast<uint8_t>(a) | static_cast<uint8_t>(b)); }
inline bool hasDim(Dim dims, Dim f) { return (static_cast<uint8_t>(dims) & static_cast<uint8_t>(f)) != 0; }
inline void clearDim(Dim& dims, Dim f) { dims = static_cast<Dim>(static_cast<uint8_t>(dims) & ~static_cast<uint8_t>(f)); }

struct NamedDef {
    const char* name;
    int outId;
    const char* measureName;
    bool isDouble;
    Dim dims;
    uint8_t method = Deploy::NA;
};

struct OutMeasure {
    int outId = -1;
    Measure m = invalidMeasure;
    bool isDouble = false;
    Dim dims = Dim::None;
    uint8_t method = Deploy::NA;
};

using NamedMeasureMapT = std::map<std::string, OutMeasure>;

inline constexpr NamedDef defs[] = {
        {"nHost", 0, "nHost", false, Dim::Age | Dim::Cohort},
        /** The number of human hosts with an infection (patent or not) at the time
        * the survey is taken. */
        {"nInfect", 1, "nInfect", false, Dim::Age | Dim::Cohort},
        {"nInfect_Imported", 1001, "nInfect_Imported", false, Dim::Age | Dim::Cohort},
        {"nInfect_Introduced", 2001, "nInfect_Introduced", false, Dim::Age | Dim::Cohort},
        {"nInfect_Indigenous", 3001, "nInfect_Indigenous", false, Dim::Age | Dim::Cohort},
        /** Expected number of infected hosts
        *
        * This is the sum of the probabilities, across all time steps since the
        * last survey, of each host becoming infected on that time step. */
        {"nExpectd", 2, "nExpectd", true, Dim::Age | Dim::Cohort},
        /** The number of human hosts whose total (blood-stage) parasite density is
        * above the detection threshold */
        {"nPatent", 3, "nPatent", false, Dim::Age | Dim::Cohort},
        {"nPatent_Imported", 1003, "nPatent_Imported", false, Dim::Age | Dim::Cohort},
        {"nPatent_Introduced", 2003, "nPatent_Introduced", false, Dim::Age | Dim::Cohort},
        {"nPatent_Indigenous", 3003, "nPatent_Indigenous", false, Dim::Age | Dim::Cohort},
        /// Sum of log(1 + p) where p is the pyrogenic threshold
        {"sumLogPyrogenThres", 4, "sumLogPyrogenThres", true, Dim::Age | Dim::Cohort},
        /** Sum (across hosts) of the natural logarithm of the parasite density of
        * hosts with detectable parasite density (patent according to the
        * monitoring diagnostic). */
        {"sumlogDens", 5, "sumlogDens", true, Dim::Age | Dim::Cohort},
        /** The total number of infections in the population: includes both blood
        * and liver stages. Vivax: this is the number of broods. */
        {"totalInfs", 6, "totalInfs", false, Dim::Age | Dim::Cohort | Dim::Genotype},
        {"totalInfs_Imported", 1006, "totalInfs_Imported", false, Dim::Age | Dim::Cohort | Dim::Genotype},
        {"totalInfs_Introduced", 2006, "totalInfs_Introduced", false, Dim::Age | Dim::Cohort | Dim::Genotype},
        {"totalInfs_Indigenous", 3006, "totalInfs_Indigenous", false, Dim::Age | Dim::Cohort | Dim::Genotype},
        /** Infectiousness of human population to mosquitoes
        *
        * Number of hosts transmitting to mosquitoes (i.e. proportion of
        * mosquitoes that get infected multiplied by human population size).
        * Single value, not per age-group. */
        {"nTransmit", 7, "nTransmit", true, Dim::None},
        /** The sum of all detectable infections (where blood stage parasite
        * density is above the detection limit) across all human hosts.
        * Vivax: the number of broods with an active blood stage. */
        {"totalPatentInf", 8, "totalPatentInf", false, Dim::Age | Dim::Cohort | Dim::Genotype},
        {"totalPatentInf_Imported", 1008, "totalPatentInf_Imported", false, Dim::Age | Dim::Cohort | Dim::Genotype},
        {"totalPatentInf_Introduced", 2008, "totalPatentInf_Introduced", false, Dim::Age | Dim::Cohort | Dim::Genotype},
        {"totalPatentInf_Indigenous", 3008, "totalPatentInf_Indigenous", false, Dim::Age | Dim::Cohort | Dim::Genotype},
        /// Contribuion to immunity functions (removed)
        {"contrib", 9, nullptr, false, Dim::None},
        /// Sum of the pyrogenic threshold
        {"sumPyrogenThresh", 10, "sumPyrogenThresh", true, Dim::Age | Dim::Cohort},
        /// number of blood-stage treatments (1st line)
        {"nTreatments1", 11, "nTreatments1", false, Dim::Age | Dim::Cohort},
        /// number of blood-stage treatments (2nd line)
        {"nTreatments2", 12, "nTreatments2", false, Dim::Age | Dim::Cohort},
        /// number of blood-stage treatments (inpatient)
        {"nTreatments3", 13, "nTreatments3", false, Dim::Age | Dim::Cohort},
        /// number of episodes (uncomplicated)
        {"nUncomp", 14, "nUncomp", false, Dim::Age | Dim::Cohort},
        {"nUncomp_Imported", 1014, "nUncomp_Imported", false, Dim::Age | Dim::Cohort},
        {"nUncomp_Introduced", 2014, "nUncomp_Introduced", false, Dim::Age | Dim::Cohort},
        {"nUncomp_Indigenous", 3014, "nUncomp_Indigenous", false, Dim::Age | Dim::Cohort},
        /// Number of severe episodes (severe malaria or malaria + coinfection)
        {"nSevere", 15, "nSevere", false, Dim::Age | Dim::Cohort},
        /// cases with sequelae
        {"nSeq", 16, "nSeq", false, Dim::Age | Dim::Cohort},
        /// deaths in hospital
        {"nHospitalDeaths", 17, "nHospitalDeaths", false, Dim::Age | Dim::Cohort},
        /// Number of deaths indirectly caused by malaria
        {"nIndDeaths", 18, "nIndDeaths", false, Dim::Age | Dim::Cohort},
        /// Number of deaths directly caused by malaria
        {"nDirDeaths", 19, "nDirDeaths", false, Dim::Age | Dim::Cohort},
        /** Number of vaccine doses given via EPI.
        *
        * Since schema 22, each vaccine type may be deployed independently. To be
        * roughly backwards-compatible, the first type (PEV, BSV or TBV) described
        * (with an "effect" element) will be reported. */
        {"nEPIVaccinations", 20, "vaccinations", false, Dim::Age | Dim::Cohort, Deploy::CTS},
        /** All cause infant mortality rate
        *
        * Reports death rate of infants due to all causes (malaria as modelled
        * plus fixed non-malaria attribution). Calculated via Kaplan-Meier method.
        * Units: deaths per thousand births. */
        {"allCauseIMR", 21, "allCauseIMR", true, Dim::None},
        /** Number of vaccine doses given via mass campaign.
        *
        * Since schema 22, each vaccine type may be deployed independently. To be
        * roughly backwards-compatible, the first type (PEV, BSV or TBV) described
        * (with an "effect" element) will be reported. */
        {"nMassVaccinations", 22, "vaccinations", false, Dim::Age | Dim::Cohort, Deploy::TIMED},
        /// recoveries in hospital
        {"nHospitalRecovs", 23, "nHospitalRecovs", false, Dim::Age | Dim::Cohort},
        /// sequelae in hospital
        {"nHospitalSeqs", 24, "nHospitalSeqs", false, Dim::Age | Dim::Cohort},
        /// Number of IPT Doses (removed together with IPT model)
        {"nIPTDoses", 25, nullptr, false, Dim::None},
        /** Annual Average Kappa
        *
        * Calculated once a year as sum of human infectiousness divided by initial
        * EIR summed over a year. Single value, not per age-group. */
        {"annAvgK", 26, "annAvgK", true, Dim::None},
        /// Number of episodes (non-malaria fever)
        {"nNMFever", 27, "nNMFever", false, Dim::Age | Dim::Cohort},
        /// Inoculations per human (all ages) per day of year, over the last year.
        /// (Reporting removed.)
        {"innoculationsPerDayOfYear", 28, nullptr, false, Dim::None},
        /// Kappa (human infectiousness) weighted by availability per day-of-year for the last year.
        /// (Reporting removed.)
        {"kappaPerDayOfYear", 29, nullptr, false, Dim::None},
        /** The total number of inoculations, by age group, cohort and parasite
        * genotype, summed over the reporting period. */
        {"innoculationsPerAgeGroup", 30, "innoculationsPerAgeGroup", true, Dim::Age | Dim::Cohort | Dim::Genotype},
        /// N_v0: emergence of feeding vectors during the last time step. Units: mosquitoes/day
        {"Vector_Nv0", 31, "Vector_Nv0", true, Dim::Species},
        /// N_v: vectors seeking to feed during the last time step. Units: mosquitoes/day
        {"Vector_Nv", 32, "Vector_Nv", true, Dim::Species},
        /// N_v: infected vectors seeking to feed during the last time step. Units: mosquitoes/day
        {"Vector_Ov", 33, "Vector_Ov", true, Dim::Species | Dim::Genotype},
        /// N_v: infectious vectors seeking to feed during the last time step. Units: mosquitoes/day
        {"Vector_Sv", 34, "Vector_Sv", true, Dim::Species | Dim::Genotype},
        /** Input EIR (Expected EIR entered into scenario file)
        *
        * Units: inoculations per adult per time step. */
        {"inputEIR", 35, "inputEIR", true, Dim::None},
        /** Simulated EIR (EIR output by the transmission model)
        *
        * Units: inoculations per adult per time step (children are excluded
        * when measuring). */
        {"simulatedEIR", 36, "simulatedEIR", true, Dim::None},
        {"simulatedEIR_Introduced", 2036, "simulatedEIR_Introduced", true, Dim::None},
        {"simulatedEIR_Indigenous", 3036, "simulatedEIR_Indigenous", true, Dim::None},
        /// Number of Rapid Diagnostic Tests used
        {"Clinical_RDTs", 39, nullptr, false, Dim::None},
        /* Effective total quanty of each drug used orally, in mg.
        * (Per active ingredient abbreviation.)
        *
        * The quantity is efffective with respect to the cost (see treatment
        * schedule definition).
        *
        * Reporting removed. */
        {"Clinical_DrugUsage", 40, nullptr, false, Dim::None},
        /// Direct death on first day of CM (before treatment takes effect)
        {"Clinical_FirstDayDeaths", 41, "Clinical_FirstDayDeaths", false, Dim::Age | Dim::Cohort},
        /// Direct death on first day of CM (before treatment takes effect); hospital only
        {"Clinical_HospitalFirstDayDeaths", 42, "Clinical_HospitalFirstDayDeaths", false, Dim::Age | Dim::Cohort},
        /** The number of actual infections since the last survey. */
        {"nNewInfections", 43, "nNewInfections", false, Dim::Age | Dim::Cohort},
        {"nNewInfections_Imported", 1043, "nNewInfections_Imported", false, Dim::Age | Dim::Cohort},
        {"nNewInfections_Introduced", 2043, "nNewInfections_Introduced", false, Dim::Age | Dim::Cohort},
        {"nNewInfections_Indigenous", 3043, "nNewInfections_Indigenous", false, Dim::Age | Dim::Cohort},
        /** The number of ITNs delivered by mass distribution since last survey.
        *
        * These are "modelled ITNs": cover only a single person, cannot be passed
        * to someone else for reuse or used for fishing, etc. */
        {"nMassITNs", 44, "itn", false, Dim::Age | Dim::Cohort, Deploy::TIMED},
        /** The number of ITNs delivered through EPI since last survey.
        *
        * Comments from nMassITNs apply. */
        {"nEPI_ITNs", 45, "itn", false, Dim::Age | Dim::Cohort, Deploy::CTS},
        /** The number of people newly protected by IRS since last survey.
        *
        * Modelled IRS: affects one person, cannot be plastered over. */
        {"nMassIRS", 46, "irs", false, Dim::Age | Dim::Cohort, Deploy::TIMED},
        /** Defunct; was used by "vector availability" intervention (which is now a
        * sub-set of GVI). */
        {"nMassVA", 47, nullptr, false, Dim::None},
        /// Number of malarial tests via microscopy used
        {"Clinical_Microscopy", 48, nullptr, false, Dim::None},
        /* As Clinical_DrugUsage, but for quatities of drug delivered via IV. */
        {"Clinical_DrugUsageIV", 49, nullptr, false, Dim::None},
        /// Number of cohort recruitments removed)
        {"nAddedToCohort", 50, nullptr, false, Dim::None},
        /// Number of individuals removed from cohort (removed)
        {"nRemovedFromCohort", 51, nullptr, false, Dim::None},
        /** Number of people (per age group) treated by mass drug administration
        * campaign. (Note that in one day time-step model MDA can be configured
        * as screen-and-treat. This option reports treatments administered not
        * the number of tests used.) */
        {"nMDAs", 52, "treat", false, Dim::Age | Dim::Cohort, Deploy::TIMED},
        /// Number of deaths caused by non-malaria fevers
        {"nNmfDeaths", 53, "nNmfDeaths", false, Dim::Age | Dim::Cohort},
        /// Number of antibiotic treatments given (disabled — not used)
        {"nAntibioticTreatments", 54, nullptr, false, Dim::None},
        /** Report the number of screenings used in a mass screen-and-treat
        * operation. */
        {"nMassScreenings", 55, "screen", false, Dim::Age | Dim::Cohort, Deploy::TIMED},
        /// Report the number of mass deployments of generic vector interventions.
        {"nMassGVI", 56, "gvi", false, Dim::Age | Dim::Cohort, Deploy::TIMED},
        /** Number of IRS deployments via continuous deployment. */
        {"nCtsIRS", 57, "irs", false, Dim::Age | Dim::Cohort, Deploy::CTS},
        /** Number of GVI deployments via continuous deployment. */
        {"nCtsGVI", 58, "gvi", false, Dim::Age | Dim::Cohort, Deploy::CTS},
        /** Number of "MDA" deployments via continuous deployment.
        *
        * Note: MDA stands for mass drug administration, but the term has come to
        * be used more flexibly by OpenMalaria, including optional screening and
        * deployment through age-based systems. */
        {"nCtsMDA", 59, "treat", false, Dim::Age | Dim::Cohort, Deploy::CTS},
        /** Number of diagnostics used by "MDA" distribution through continuous
        * methods. Can be higher than nCtsMDA since drugs are administered only
        * when the diagnostic is positive. Also see nCtsMDA description. */
        {"nCtsScreenings", 60, "screen", false, Dim::Age | Dim::Cohort, Deploy::CTS},
        /** Number of removals from a sub-population due to expiry of duration of
        * membership (e.g. intervention too old). */
        {"nSubPopRemovalTooOld", 61, "nSubPopRemovalTooOld", false, Dim::Age | Dim::Cohort},
        /** Number of removals from a sub-population due to first
        * infection/bout/treatment (see onFirstBout & co). */
        {"nSubPopRemovalFirstEvent", 62, "nSubPopRemovalFirstEvent", false, Dim::Age | Dim::Cohort},
        /** Report the number of liver-stage treatments (likely Primaquine) administered. */
        {"nLiverStageTreatments", 63, "nLiverStageTreatments", false, Dim::Age | Dim::Cohort},
        /** Report the number of diagnostics used during treatment.
        *
        * This is not the same as Clinical_RDTs + Clinical_Microscopy: those
        * outputs are used by the "event scheduler" 1-day time step clinical
        * model, whereas this output is used by the 5-day time step model. */
        {"nTreatDiagnostics", 64, "nTreatDiagnostics", false, Dim::Age | Dim::Cohort},
        /** Number of "recruitment only" recruitments via timed deployment. */
        {"nMassRecruitOnly", 65, "recruit", false, Dim::Age | Dim::Cohort, Deploy::TIMED},
        /** Number of "recruitment only" recruitments via age-based deployment. */
        {"nCtsRecruitOnly", 66, "recruit", false, Dim::Age | Dim::Cohort, Deploy::CTS},
        /** Number of deployments (of all intervention components) triggered by
        * treatment (case management). */
        {"nTreatDeployments", 67, "nTreatDeployments", false, Dim::Age | Dim::Cohort, Deploy::TREAT},
        /** Report the total age of all humans in this a group (sum across humans,
        * in years). Divide by nHost to get the average age. */
        {"sumAge", 68, "sumAge", true, Dim::Age | Dim::Cohort},
        /** The number of human hosts with an infection (patent or not), for each
        * genotype, at the time the survey is taken. */
        {"nInfectByGenotype", 69, "nInfectByGenotype", false, Dim::Age | Dim::Cohort | Dim::Genotype},
        /** The number of human hosts whose total (blood-stage) parasite density,
        * for each genotype, is above the detection threshold */
        {"nPatentByGenotype", 70, "nPatentByGenotype", false, Dim::Age | Dim::Cohort | Dim::Genotype},
        /** For each infection genotype, sum across humans the natural log of
        * parasite density (like sumlogDens but per genotype). */
        {"logDensByGenotype", 71, "logDensByGenotype", true, Dim::Age | Dim::Cohort | Dim::Genotype},
        /** For each drug type in the pharmacology section of the XML, report the
        * number of humans with non-zero concentration of this drug in their
        * blood. */
        {"nHostDrugConcNonZero", 72, "nHostDrugConcNonZero", false, Dim::Age | Dim::Cohort | Dim::Drug},
        /** For each drug type in the pharmacology section of the XML, report the
        * sum of the natural logarithm of the drug concentration in hosts with
        * non-zero concentration. */
        {"sumLogDrugConcNonZero", 73, "sumLogDrugConcNonZero", true, Dim::Age | Dim::Cohort | Dim::Drug},
        /** Expected number of direct malaria deaths, from those with severe
        * disease.
        *
        * This is calculated as the sum over all steps in the reporting period of
        * the sum over humans with severe malaria of the probability of direct
        * death from malaria. */
        {"expectedDirectDeaths", 74, "expectedDirectDeaths", true, Dim::Age | Dim::Cohort},
        /** Expected number of direct malaria deaths which occur in hospital.
        *
        * This is the a subset of `expectedDirectDeaths` and the same notes apply.
        */
        {"expectedHospitalDeaths", 75, "expectedHospitalDeaths", true, Dim::Age | Dim::Cohort},
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
        {"expectedIndirectDeaths", 76, "expectedIndirectDeaths", true, Dim::Age | Dim::Cohort},
        /** Expected number of sequelae, from those with severe disease.
        *
        * This is calculated as the sum over all steps in the reporting period of
        * the sum over humans with severe malaria of the probability of sequelae
        * occuring, assuming the human "recovers" from the bout.
        */
        {"expectedSequelae", 77, "expectedSequelae", true, Dim::Age | Dim::Cohort},
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
        {"expectedSevere", 78, "expectedSevere", true, Dim::Age | Dim::Cohort},
        /** The total number of inoculations, by mosquito species, summed over
        * the reporting period. */
        {"innoculationsPerVector", 79, "innoculationsPerAgeGroup", true, Dim::Species},
        /** Number of custom intervention reports done */
        {"nCMDTReport", 80, "nCMDTReport", false, Dim::Age | Dim::Cohort},
        /// Similar to nSevere. Number of severe episodes WITHOUT coinfection
        {"nSevereWithoutComorbidities", 81, "nSevereWithoutComorbidities", false, Dim::Age | Dim::Cohort},
        /** Similar to 'expectedSevere'.
        * Expected number of severe bouts of malaria WITHOUT "complications due
        * to coinfection" (the same as the `nSevereWithoutComorbidities` output). */
        {"expectedSevereWithoutComorbidities", 82, "expectedSevereWithoutComorbidities", true, Dim::Age | Dim::Cohort},
    };

constexpr bool cstrEquals(const char* a, const char* b)
{
    if (a == nullptr || b == nullptr) return a == b;
    while (*a != '\0' && *b != '\0') {
        if (*a != *b) return false;
        ++a;
        ++b;
    }
    return *a == *b;
}

constexpr bool cstrEquals(const char* a, std::string_view b)
{
    if (a == nullptr) return false;
    for (size_t i = 0; i < b.size(); ++i) {
        if (a[i] == '\0' || a[i] != b[i]) return false;
    }
    return a[b.size()] == '\0';
}

constexpr bool isBuiltinMeasureName(const char* name)
{
    return name != nullptr && !cstrEquals(name, "allCauseIMR");
}

constexpr bool isFirstBuiltinOccurrence(size_t i)
{
    const char* current = defs[i].measureName;
    if (!isBuiltinMeasureName(current)) return false;
    for (size_t j = 0; j < i; ++j) {
        if (cstrEquals(current, defs[j].measureName)) return false;
    }
    return true;
}

constexpr Measure countBuiltinMeasures()
{
    Measure count = 0;
    for (size_t i = 0; i < sizeof(defs) / sizeof(defs[0]); ++i) {
        if (isFirstBuiltinOccurrence(i)) {
            ++count;
        }
    }
    return count;
}

inline constexpr Measure builtinMeasureCount = countBuiltinMeasures();
inline constexpr Measure MeasureCount = static_cast<Measure>(builtinMeasureCount + 4);

constexpr Measure lookupBuiltinMeasure(std::string_view measureName)
{
    Measure id = 0;
    for (size_t i = 0; i < sizeof(defs) / sizeof(defs[0]); ++i) {
        if (!isFirstBuiltinOccurrence(i)) continue;
        if (cstrEquals(defs[i].measureName, measureName)) {
            return id;
        }
        ++id;
    }
    return invalidMeasure;
}

constexpr Measure measure(std::string_view measureName)
{
    Measure id = lookupBuiltinMeasure(measureName);
    if (id != invalidMeasure) return id;

    if (measureName == "nNMFTreatments") return builtinMeasureCount;
    if (measureName == "pev") return static_cast<Measure>(builtinMeasureCount + 1);
    if (measureName == "bsv") return static_cast<Measure>(builtinMeasureCount + 2);
    if (measureName == "tbv") return static_cast<Measure>(builtinMeasureCount + 3);
    return invalidMeasure;
}

template <size_t N>
constexpr Measure measure(const char (&measureName)[N])
{
    return measure(std::string_view(measureName, N - 1));
}

inline void defineOutMeasures(NamedMeasureMapT& namedOutMeasures, std::set<Measure>& validCondMeasures)
{
    namedOutMeasures.clear();
    validCondMeasures.clear();

    std::set<int> seenIDs;
    for (const NamedDef& d : defs) {
        OutMeasure om;
        om.outId = d.outId;
        if (d.measureName == nullptr) {
            om.m = obsoleteMeasure;
        } else if (cstrEquals(d.measureName, "allCauseIMR")) {
            om.m = allCauseIMR;
        } else {
            om.m = measure(d.measureName);
            if (om.m == invalidMeasure) {
                throw std::runtime_error("Unknown measure in definition: " + std::string(d.measureName));
            }
        }
        om.isDouble = d.isDouble;
        om.dims = d.dims;
        om.method = d.method;

        auto [it, inserted] = namedOutMeasures.emplace(d.name, om);
        if (!inserted) {
            throw std::runtime_error("Duplicate OutMeasure name detected: " + std::string(d.name));
        }
        if (!seenIDs.insert(d.outId).second) {
            throw std::runtime_error("Duplicate OutMeasure (outId) detected: " + std::to_string(d.outId));
        }
    }

    for (const NamedMeasureMapT::value_type& v : namedOutMeasures) {
        Measure m = v.second.m;
        if (m == measure("nSevere") ||
            m == measure("nSevereWithoutComorbidities") ||
            m == measure("nUncomp") ||
            m == measure("nDirDeaths") ||
            m == measure("nHospitalDeaths") ||
            m == measure("Clinical_FirstDayDeaths") ||
            m == measure("Clinical_HospitalFirstDayDeaths") ||
            m == measure("nSeq") ||
            m == measure("nHospitalSeqs") ||
            m == measure("nHospitalRecovs") ||
            m == measure("nNMFever") ||
            m == measure("nNmfDeaths") ||
            m == measure("nSubPopRemovalFirstEvent") ||
            m == measure("innoculationsPerAgeGroup") ||
            m == measure("inputEIR") ||
            m == measure("simulatedEIR")) {
            continue;
        }
        validCondMeasures.insert(m);
    }
}

} // namespace mon
} // namespace OM

#endif
