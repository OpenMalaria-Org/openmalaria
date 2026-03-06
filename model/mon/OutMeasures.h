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

struct OutMeasure {
    const char* userName = nullptr; // defs[] only
    int outId = -1;
    bool isDouble = false;
    Dim dims = Dim::None;
    uint8_t method = Deploy::NA;
    const char* internalName = nullptr; // defs[] only
    bool obsolete = false;       // defs[] only
    Measure m = 0;               // resolved at init-time
};

using NamedMeasureMapT = std::map<std::string, OutMeasure>;

inline constexpr OutMeasure defs[] = {
        {"nHost", 0, false, Dim::Age | Dim::Cohort},
        /** The number of human hosts with an infection (patent or not) at the time
        * the survey is taken. */
        {"nInfect", 1, false, Dim::Age | Dim::Cohort},
        {"nInfect_Imported", 1001, false, Dim::Age | Dim::Cohort},
        {"nInfect_Introduced", 2001, false, Dim::Age | Dim::Cohort},
        {"nInfect_Indigenous", 3001, false, Dim::Age | Dim::Cohort},
        /** Expected number of infected hosts
        *
        * This is the sum of the probabilities, across all time steps since the
        * last survey, of each host becoming infected on that time step. */
        {"nExpectd", 2, true, Dim::Age | Dim::Cohort},
        /** The number of human hosts whose total (blood-stage) parasite density is
        * above the detection threshold */
        {"nPatent", 3, false, Dim::Age | Dim::Cohort},
        {"nPatent_Imported", 1003, false, Dim::Age | Dim::Cohort},
        {"nPatent_Introduced", 2003, false, Dim::Age | Dim::Cohort},
        {"nPatent_Indigenous", 3003, false, Dim::Age | Dim::Cohort},
        /// Sum of log(1 + p) where p is the pyrogenic threshold
        {"sumLogPyrogenThres", 4, true, Dim::Age | Dim::Cohort},
        /** Sum (across hosts) of the natural logarithm of the parasite density of
        * hosts with detectable parasite density (patent according to the
        * monitoring diagnostic). */
        {"sumlogDens", 5, true, Dim::Age | Dim::Cohort},
        /** The total number of infections in the population: includes both blood
        * and liver stages. Vivax: this is the number of broods. */
        {"totalInfs", 6, false, Dim::Age | Dim::Cohort | Dim::Genotype},
        {"totalInfs_Imported", 1006, false, Dim::Age | Dim::Cohort | Dim::Genotype},
        {"totalInfs_Introduced", 2006, false, Dim::Age | Dim::Cohort | Dim::Genotype},
        {"totalInfs_Indigenous", 3006, false, Dim::Age | Dim::Cohort | Dim::Genotype},
        /** Infectiousness of human population to mosquitoes
        *
        * Number of hosts transmitting to mosquitoes (i.e. proportion of
        * mosquitoes that get infected multiplied by human population size).
        * Single value, not per age-group. */
        {"nTransmit", 7, true, Dim::None},
        /** The sum of all detectable infections (where blood stage parasite
        * density is above the detection limit) across all human hosts.
        * Vivax: the number of broods with an active blood stage. */
        {"totalPatentInf", 8, false, Dim::Age | Dim::Cohort | Dim::Genotype},
        {"totalPatentInf_Imported", 1008, false, Dim::Age | Dim::Cohort | Dim::Genotype},
        {"totalPatentInf_Introduced", 2008, false, Dim::Age | Dim::Cohort | Dim::Genotype},
        {"totalPatentInf_Indigenous", 3008, false, Dim::Age | Dim::Cohort | Dim::Genotype},
        /// Contribuion to immunity functions (removed)
        {"contrib", 9, false, Dim::None, Deploy::NA, nullptr, true},
        /// Sum of the pyrogenic threshold
        {"sumPyrogenThresh", 10, true, Dim::Age | Dim::Cohort},
        /// number of blood-stage treatments (1st line)
        {"nTreatments1", 11, false, Dim::Age | Dim::Cohort},
        /// number of blood-stage treatments (2nd line)
        {"nTreatments2", 12, false, Dim::Age | Dim::Cohort},
        /// number of blood-stage treatments (inpatient)
        {"nTreatments3", 13, false, Dim::Age | Dim::Cohort},
        /// number of episodes (uncomplicated)
        {"nUncomp", 14, false, Dim::Age | Dim::Cohort},
        {"nUncomp_Imported", 1014, false, Dim::Age | Dim::Cohort},
        {"nUncomp_Introduced", 2014, false, Dim::Age | Dim::Cohort},
        {"nUncomp_Indigenous", 3014, false, Dim::Age | Dim::Cohort},
        /// Number of severe episodes (severe malaria or malaria + coinfection)
        {"nSevere", 15, false, Dim::Age | Dim::Cohort},
        /// cases with sequelae
        {"nSeq", 16, false, Dim::Age | Dim::Cohort},
        /// deaths in hospital
        {"nHospitalDeaths", 17, false, Dim::Age | Dim::Cohort},
        /// Number of deaths indirectly caused by malaria
        {"nIndDeaths", 18, false, Dim::Age | Dim::Cohort},
        /// Number of deaths directly caused by malaria
        {"nDirDeaths", 19, false, Dim::Age | Dim::Cohort},
        /** Number of vaccine doses given via EPI.
        *
        * Since schema 22, each vaccine type may be deployed independently. To be
        * roughly backwards-compatible, the first type (PEV, BSV or TBV) described
        * (with an "effect" element) will be reported. */
        {"nEPIVaccinations", 20, false, Dim::Age | Dim::Cohort, Deploy::CTS, "vaccinations"},
        /** All cause infant mortality rate
        *
        * Reports death rate of infants due to all causes (malaria as modelled
        * plus fixed non-malaria attribution). Calculated via Kaplan-Meier method.
        * Units: deaths per thousand births. */
        {"allCauseIMR", 21, true, Dim::None},
        /** Number of vaccine doses given via mass campaign.
        *
        * Since schema 22, each vaccine type may be deployed independently. To be
        * roughly backwards-compatible, the first type (PEV, BSV or TBV) described
        * (with an "effect" element) will be reported. */
        {"nMassVaccinations", 22, false, Dim::Age | Dim::Cohort, Deploy::TIMED, "vaccinations"},
        /// recoveries in hospital
        {"nHospitalRecovs", 23, false, Dim::Age | Dim::Cohort},
        /// sequelae in hospital
        {"nHospitalSeqs", 24, false, Dim::Age | Dim::Cohort},
        /// Number of IPT Doses (removed together with IPT model)
        {"nIPTDoses", 25, false, Dim::None, Deploy::NA, nullptr, true},
        /** Annual Average Kappa
        *
        * Calculated once a year as sum of human infectiousness divided by initial
        * EIR summed over a year. Single value, not per age-group. */
        {"annAvgK", 26, true, Dim::None},
        /// Number of episodes (non-malaria fever)
        {"nNMFever", 27, false, Dim::Age | Dim::Cohort},
        /// Inoculations per human (all ages) per day of year, over the last year.
        /// (Reporting removed.)
        {"innoculationsPerDayOfYear", 28, false, Dim::None, Deploy::NA, nullptr, true},
        /// Kappa (human infectiousness) weighted by availability per day-of-year for the last year.
        /// (Reporting removed.)
        {"kappaPerDayOfYear", 29, false, Dim::None, Deploy::NA, nullptr, true},
        /** The total number of inoculations, by age group, cohort and parasite
        * genotype, summed over the reporting period. */
        {"innoculationsPerAgeGroup", 30, true, Dim::Age | Dim::Cohort | Dim::Genotype},
        /// N_v0: emergence of feeding vectors during the last time step. Units: mosquitoes/day
        {"Vector_Nv0", 31, true, Dim::Species},
        /// N_v: vectors seeking to feed during the last time step. Units: mosquitoes/day
        {"Vector_Nv", 32, true, Dim::Species},
        /// N_v: infected vectors seeking to feed during the last time step. Units: mosquitoes/day
        {"Vector_Ov", 33, true, Dim::Species | Dim::Genotype},
        /// N_v: infectious vectors seeking to feed during the last time step. Units: mosquitoes/day
        {"Vector_Sv", 34, true, Dim::Species | Dim::Genotype},
        /** Input EIR (Expected EIR entered into scenario file)
        *
        * Units: inoculations per adult per time step. */
        {"inputEIR", 35, true, Dim::None},
        /** Simulated EIR (EIR output by the transmission model)
        *
        * Units: inoculations per adult per time step (children are excluded
        * when measuring). */
        {"simulatedEIR", 36, true, Dim::None},
        {"simulatedEIR_Introduced", 2036, true, Dim::None},
        {"simulatedEIR_Indigenous", 3036, true, Dim::None},
        /// Number of Rapid Diagnostic Tests used
        {"Clinical_RDTs", 39, false, Dim::None, Deploy::NA, nullptr, true},
        /* Effective total quanty of each drug used orally, in mg.
        * (Per active ingredient abbreviation.)
        *
        * The quantity is efffective with respect to the cost (see treatment
        * schedule definition).
        *
        * Reporting removed. */
        {"Clinical_DrugUsage", 40, false, Dim::None, Deploy::NA, nullptr, true},
        /// Direct death on first day of CM (before treatment takes effect)
        {"Clinical_FirstDayDeaths", 41, false, Dim::Age | Dim::Cohort},
        /// Direct death on first day of CM (before treatment takes effect); hospital only
        {"Clinical_HospitalFirstDayDeaths", 42, false, Dim::Age | Dim::Cohort},
        /** The number of actual infections since the last survey. */
        {"nNewInfections", 43, false, Dim::Age | Dim::Cohort},
        {"nNewInfections_Imported", 1043, false, Dim::Age | Dim::Cohort},
        {"nNewInfections_Introduced", 2043, false, Dim::Age | Dim::Cohort},
        {"nNewInfections_Indigenous", 3043, false, Dim::Age | Dim::Cohort},
        /** The number of ITNs delivered by mass distribution since last survey.
        *
        * These are "modelled ITNs": cover only a single person, cannot be passed
        * to someone else for reuse or used for fishing, etc. */
        {"nMassITNs", 44, false, Dim::Age | Dim::Cohort, Deploy::TIMED, "itn"},
        /** The number of ITNs delivered through EPI since last survey.
        *
        * Comments from nMassITNs apply. */
        {"nEPI_ITNs", 45, false, Dim::Age | Dim::Cohort, Deploy::CTS, "itn"},
        /** The number of people newly protected by IRS since last survey.
        *
        * Modelled IRS: affects one person, cannot be plastered over. */
        {"nMassIRS", 46, false, Dim::Age | Dim::Cohort, Deploy::TIMED, "irs"},
        /** Defunct; was used by "vector availability" intervention (which is now a
        * sub-set of GVI). */
        {"nMassVA", 47, false, Dim::None, Deploy::NA, nullptr, true},
        /// Number of malarial tests via microscopy used
        {"Clinical_Microscopy", 48, false, Dim::None, Deploy::NA, nullptr, true},
        /* As Clinical_DrugUsage, but for quatities of drug delivered via IV. */
        {"Clinical_DrugUsageIV", 49, false, Dim::None, Deploy::NA, nullptr, true},
        /// Number of cohort recruitments removed)
        {"nAddedToCohort", 50, false, Dim::None, Deploy::NA, nullptr, true},
        /// Number of individuals removed from cohort (removed)
        {"nRemovedFromCohort", 51, false, Dim::None, Deploy::NA, nullptr, true},
        /** Number of people (per age group) treated by mass drug administration
        * campaign. (Note that in one day time-step model MDA can be configured
        * as screen-and-treat. This option reports treatments administered not
        * the number of tests used.) */
        {"nMDAs", 52, false, Dim::Age | Dim::Cohort, Deploy::TIMED, "treat"},
        /// Number of deaths caused by non-malaria fevers
        {"nNmfDeaths", 53, false, Dim::Age | Dim::Cohort},
        /// Number of antibiotic treatments given (disabled — not used)
        {"nAntibioticTreatments", 54, false, Dim::None, Deploy::NA, nullptr, true},
        /** Report the number of screenings used in a mass screen-and-treat
        * operation. */
        {"nMassScreenings", 55, false, Dim::Age | Dim::Cohort, Deploy::TIMED, "screen"},
        /// Report the number of mass deployments of generic vector interventions.
        {"nMassGVI", 56, false, Dim::Age | Dim::Cohort, Deploy::TIMED, "gvi"},
        /** Number of IRS deployments via continuous deployment. */
        {"nCtsIRS", 57, false, Dim::Age | Dim::Cohort, Deploy::CTS, "irs"},
        /** Number of GVI deployments via continuous deployment. */
        {"nCtsGVI", 58, false, Dim::Age | Dim::Cohort, Deploy::CTS, "gvi"},
        /** Number of "MDA" deployments via continuous deployment.
        *
        * Note: MDA stands for mass drug administration, but the term has come to
        * be used more flexibly by OpenMalaria, including optional screening and
        * deployment through age-based systems. */
        {"nCtsMDA", 59, false, Dim::Age | Dim::Cohort, Deploy::CTS, "treat"},
        /** Number of diagnostics used by "MDA" distribution through continuous
        * methods. Can be higher than nCtsMDA since drugs are administered only
        * when the diagnostic is positive. Also see nCtsMDA description. */
        {"nCtsScreenings", 60, false, Dim::Age | Dim::Cohort, Deploy::CTS, "screen"},
        /** Number of removals from a sub-population due to expiry of duration of
        * membership (e.g. intervention too old). */
        {"nSubPopRemovalTooOld", 61, false, Dim::Age | Dim::Cohort},
        /** Number of removals from a sub-population due to first
        * infection/bout/treatment (see onFirstBout & co). */
        {"nSubPopRemovalFirstEvent", 62, false, Dim::Age | Dim::Cohort},
        /** Report the number of liver-stage treatments (likely Primaquine) administered. */
        {"nLiverStageTreatments", 63, false, Dim::Age | Dim::Cohort},
        /** Report the number of diagnostics used during treatment.
        *
        * This is not the same as Clinical_RDTs + Clinical_Microscopy: those
        * outputs are used by the "event scheduler" 1-day time step clinical
        * model, whereas this output is used by the 5-day time step model. */
        {"nTreatDiagnostics", 64, false, Dim::Age | Dim::Cohort},
        /** Number of "recruitment only" recruitments via timed deployment. */
        {"nMassRecruitOnly", 65, false, Dim::Age | Dim::Cohort, Deploy::TIMED, "recruit"},
        /** Number of "recruitment only" recruitments via age-based deployment. */
        {"nCtsRecruitOnly", 66, false, Dim::Age | Dim::Cohort, Deploy::CTS, "recruit"},
        /** Number of deployments (of all intervention components) triggered by
        * treatment (case management). */
        {"nTreatDeployments", 67, false, Dim::Age | Dim::Cohort, Deploy::TREAT},
        /** Report the total age of all humans in this a group (sum across humans,
        * in years). Divide by nHost to get the average age. */
        {"sumAge", 68, true, Dim::Age | Dim::Cohort},
        /** The number of human hosts with an infection (patent or not), for each
        * genotype, at the time the survey is taken. */
        {"nInfectByGenotype", 69, false, Dim::Age | Dim::Cohort | Dim::Genotype},
        /** The number of human hosts whose total (blood-stage) parasite density,
        * for each genotype, is above the detection threshold */
        {"nPatentByGenotype", 70, false, Dim::Age | Dim::Cohort | Dim::Genotype},
        /** For each infection genotype, sum across humans the natural log of
        * parasite density (like sumlogDens but per genotype). */
        {"logDensByGenotype", 71, true, Dim::Age | Dim::Cohort | Dim::Genotype},
        /** For each drug type in the pharmacology section of the XML, report the
        * number of humans with non-zero concentration of this drug in their
        * blood. */
        {"nHostDrugConcNonZero", 72, false, Dim::Age | Dim::Cohort | Dim::Drug},
        /** For each drug type in the pharmacology section of the XML, report the
        * sum of the natural logarithm of the drug concentration in hosts with
        * non-zero concentration. */
        {"sumLogDrugConcNonZero", 73, true, Dim::Age | Dim::Cohort | Dim::Drug},
        /** Expected number of direct malaria deaths, from those with severe
        * disease.
        *
        * This is calculated as the sum over all steps in the reporting period of
        * the sum over humans with severe malaria of the probability of direct
        * death from malaria. */
        {"expectedDirectDeaths", 74, true, Dim::Age | Dim::Cohort},
        /** Expected number of direct malaria deaths which occur in hospital.
        *
        * This is the a subset of `expectedDirectDeaths` and the same notes apply.
        */
        {"expectedHospitalDeaths", 75, true, Dim::Age | Dim::Cohort},
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
        {"expectedIndirectDeaths", 76, true, Dim::Age | Dim::Cohort},
        /** Expected number of sequelae, from those with severe disease.
        *
        * This is calculated as the sum over all steps in the reporting period of
        * the sum over humans with severe malaria of the probability of sequelae
        * occuring, assuming the human "recovers" from the bout.
        */
        {"expectedSequelae", 77, true, Dim::Age | Dim::Cohort},
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
        {"expectedSevere", 78, true, Dim::Age | Dim::Cohort},
        /** The total number of inoculations, by mosquito species, summed over
        * the reporting period. */
        {"innoculationsPerVector", 79, true, Dim::Species, Deploy::NA, "innoculationsPerAgeGroup"},
        /** Number of custom intervention reports done */
        {"nCMDTReport", 80, false, Dim::Age | Dim::Cohort},
        /// Similar to nSevere. Number of severe episodes WITHOUT coinfection
        {"nSevereWithoutComorbidities", 81, false, Dim::Age | Dim::Cohort},
        /** Similar to 'expectedSevere'.
        * Expected number of severe bouts of malaria WITHOUT "complications due
        * to coinfection" (the same as the `nSevereWithoutComorbidities` output). */
        {"expectedSevereWithoutComorbidities", 82, true, Dim::Age | Dim::Cohort},
        // Number of treatments for non-malaria infections. Units: treatments (whole courses)
        {"nNMFTreatments", -1, false, Dim::Age | Dim::Cohort},
        // Number of pre-erythrocytic vaccine doses deployed. Units: doses (as above)
        {"pev", -1, false, Dim::Age | Dim::Cohort},
        // Number of blood-stage vaccine doses deployed. Units: doses (as above)
        {"bsv", -1, false, Dim::Age | Dim::Cohort},
        // Number of transmission-blocking vaccine doses deployed. Units: doses (as above)
        {"tbv", -1, false, Dim::Age | Dim::Cohort},
    };

inline constexpr size_t defCount = sizeof(defs) / sizeof(defs[0]);

constexpr bool isAllCauseDef(const OutMeasure& d)
{
    return std::string_view(d.userName) == "allCauseIMR";
}

constexpr std::string_view measureKey(const OutMeasure& d)
{
    return d.internalName != nullptr ? std::string_view(d.internalName) : std::string_view(d.userName);
}

constexpr bool isFirstBuiltinMeasure(size_t i)
{
    if (defs[i].obsolete || isAllCauseDef(defs[i])) return false;
    const std::string_view current = measureKey(defs[i]);
    if (current.empty()) return false;
    for (size_t j = 0; j < i; ++j) {
        if (defs[j].obsolete || isAllCauseDef(defs[j])) continue;
        if (measureKey(defs[j]) == current) return false;
    }
    return true;
}

constexpr Measure countBuiltinMeasures()
{
    Measure count = 0;
    for (size_t i = 0; i < defCount; ++i) {
        if (isFirstBuiltinMeasure(i)) {
            ++count;
        }
    }
    return count;
}

inline constexpr Measure MeasureCount = countBuiltinMeasures();

constexpr Measure measure(std::string_view measureName)
{
    Measure id = 0;
    for (size_t i = 0; i < defCount; ++i) {
        if (!isFirstBuiltinMeasure(i)) continue;
        if (measureKey(defs[i]) == measureName) return id;
        ++id;
    }
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
    for (const OutMeasure& d : defs) {
        OutMeasure om;
        om.outId = d.outId;
        if (d.obsolete) {
            om.m = obsoleteMeasure;
        } else if (isAllCauseDef(d)) {
            om.m = allCauseIMR;
        } else {
            const std::string_view key = measureKey(d);
            om.m = measure(key);
            if (om.m == invalidMeasure) {
                throw std::runtime_error("Unknown measure in definition: " + std::string(key));
            }
        }
        om.isDouble = d.isDouble;
        om.dims = d.dims;
        om.method = d.method;
        if (d.outId < 0) continue;

        auto [it, inserted] = namedOutMeasures.emplace(d.userName, om);
        if (!inserted) {
            throw std::runtime_error("Duplicate OutMeasure name detected: " + std::string(d.userName));
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
