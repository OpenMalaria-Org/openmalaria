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

#include "mon/init.h"

#include "Host/WithinHost/Diagnostic.h"
#include "interventions/InterventionManager.h"
#include "mon/Monitoring.h"
#include "schema/monitoring.h"
#include "schema/scenario.h"
#include "util/CommandLine.h"
#include "util/UnitParse.h"
#include "util/errors.h"

#include <algorithm>
#include <cctype>
#include <iostream>
#include <map>

namespace OM {
namespace mon {

using internal::SurveyDate;
using internal::runtime;

namespace {

bool notPowerOfTwo(uint32_t num)
{
    return num == 0 || num > (static_cast<uint32_t>(1) << 21) || (num & (num - 1)) != 0;
}

}

void initReporting(const scnXml::Scenario& scenario)
{
    defineOutMeasures(runtime.namedOutMeasures, runtime.validCondMeasures);
    assert(runtime.reportedMeasures.empty());
    runtime.reportIMR = -1;
    const scnXml::MonitoringOptions& optsElt = scenario.getMonitoring().getSurveyOptions();
    runtime.reportedMeasures.reserve(optsElt.getOption().size() + runtime.namedOutMeasures.size());
    auto applyCategory = [](Dim& dims, Dim flag, const bool optionalPresent, const bool requested,
                            const std::string& optionName, const char* label)
    {
        if (!optionalPresent) return;
        const bool supports = hasDim(dims, flag);
        if (supports) {
            if (!requested) clearDim(dims, flag);
        } else if (requested) {
            throw util::xml_scenario_error("measure " + optionName + " does not support categorisation by " + label);
        }
    };

    std::set<int> outIds;
    for (const scnXml::MonitoringOption& optElt : optsElt.getOption()) {
        if (!optElt.getValue()) continue;

        auto it = runtime.namedOutMeasures.find(optElt.getName());
        if (it == runtime.namedOutMeasures.end()) {
            throw util::xml_scenario_error("unrecognised survey option: " + std::string(optElt.getName()));
        }
        OutMeasure om = it->second;
        if (om.m >= MeasureCount) {
            if (om.m == obsoleteMeasure) {
                throw util::xml_scenario_error("obsolete survey option: " + std::string(optElt.getName()));
            }
            assert(om.m == allCauseIMR);
            const bool byAge = optElt.getByAge().present() && optElt.getByAge().get();
            const bool byCohort = optElt.getByCohort().present() && optElt.getByCohort().get();
            const bool bySpecies = optElt.getBySpecies().present() && optElt.getBySpecies().get();
            const bool byGenotype = optElt.getByGenotype().present() && optElt.getByGenotype().get();
            const bool byDrug = optElt.getByDrugType().present() && optElt.getByDrugType().get();
            if (!om.isDouble || byAge || byCohort || bySpecies || byGenotype || byDrug) {
                throw util::xml_scenario_error("measure allCauseIMR does not support any categorisation");
            }
            if (optElt.getOutputNumber().present()) om.outId = optElt.getOutputNumber().get();
            if (outIds.count(om.outId)) {
                throw util::xml_scenario_error("monitoring output number " + std::to_string(om.outId) + " used more than once");
            }
            outIds.insert(om.outId);
            runtime.reportIMR = om.outId;
            continue;
        }

        if ((om.m == measure("sumlogDens") || om.m == measure("logDensByGenotype")) &&
            WithinHost::diagnostics::monitoringDiagnostic().allowsFalsePositives())
        {
            throw util::xml_scenario_error("measure " + std::string(optElt.getName()) + " may not be used when monitoring diagnostic sensitivity < 1");
        }
        const std::string optionName(optElt.getName());
        const bool byAgePresent = optElt.getByAge().present();
        const bool byCohortPresent = optElt.getByCohort().present();
        const bool bySpeciesPresent = optElt.getBySpecies().present();
        const bool byGenotypePresent = optElt.getByGenotype().present();
        const bool byDrugPresent = optElt.getByDrugType().present();
        applyCategory(om.dims, Dim::Age, byAgePresent, byAgePresent ? optElt.getByAge().get() : false, optionName, "age group");
        applyCategory(om.dims, Dim::Cohort, byCohortPresent, byCohortPresent ? optElt.getByCohort().get() : false, optionName, "cohort");
        applyCategory(om.dims, Dim::Species, bySpeciesPresent, bySpeciesPresent ? optElt.getBySpecies().get() : false, optionName, "species");
        applyCategory(om.dims, Dim::Genotype, byGenotypePresent, byGenotypePresent ? optElt.getByGenotype().get() : false, optionName, "genotype");
        applyCategory(om.dims, Dim::Drug, byDrugPresent, byDrugPresent ? optElt.getByDrugType().get() : false, optionName, "drug type");

        if (optElt.getOutputNumber().present()) om.outId = optElt.getOutputNumber().get();
        if (outIds.count(om.outId)) {
            throw util::xml_scenario_error("monitoring output number " + std::to_string(om.outId) + " used more than once");
        }
        outIds.insert(om.outId);
        runtime.reportedMeasures.push_back(om);
    }
    std::sort(runtime.reportedMeasures.begin(), runtime.reportedMeasures.end(),
        [](const OutMeasure& lhs, const OutMeasure& rhs) { return lhs.outId < rhs.outId; });
    const size_t nSpecies = scenario.getEntomology().getVector().present()
        ? scenario.getEntomology().getVector().get().getAnopheles().size() : 1;
    const size_t nDrugs = scenario.getPharmacology().present()
        ? scenario.getPharmacology().get().getDrugs().getDrug().size() : 1;
    runtime.surveyStore.init(runtime.reportedMeasures, nSpecies, nDrugs);
}

SimTime readSurveyDates(const scnXml::Monitoring& monitoring)
{
    const scnXml::Surveys::SurveyTimeSequence& survs = monitoring.getSurveys().getSurveyTime();
    std::map<SimTime, bool> surveys;

    auto addSurvey = [&surveys](SimTime date, bool reporting) {
        if (reporting) surveys[date] = true;
        else surveys.insert(std::make_pair(date, false));
    };

    for (size_t i = 0; i < survs.size(); ++i) {
        const scnXml::SurveyTime& surv = survs[i];
        try {
            std::string s = surv;
            s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {
                return !std::isspace(ch);
            }));
            s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
                return !std::isspace(ch);
            }).base(), s.end());
            SimTime cur = UnitParse::readDate(s, UnitParse::STEPS);
            const bool reporting = surv.getReported();
            if (surv.getRepeatStep().present() != surv.getRepeatEnd().present()) {
                throw util::xml_scenario_error("surveyTime: use of repeatStep or repeatEnd without other");
            }
            if (surv.getRepeatStep().present()) {
                const SimTime step = UnitParse::readDuration(surv.getRepeatStep().get(), UnitParse::NONE);
                if (step < sim::oneTS()) {
                    throw util::xml_scenario_error("surveyTime: repeatStep must be >= 1");
                }
                const SimTime end = UnitParse::readDate(surv.getRepeatEnd().get(), UnitParse::NONE);
                while (cur < end) {
                    addSurvey(cur, reporting);
                    cur = cur + step;
                }
            } else {
                addSurvey(cur, reporting);
            }
        } catch (const util::format_error& e) {
            throw util::xml_scenario_error(std::string("surveyTime: ").append(e.message()));
        }
    }

    runtime.surveyDates.clear();
    runtime.surveyDates.reserve(surveys.size());
    runtime.nSurveys = 0;
    for (const auto& [date, reporting] : surveys) {
        const size_t num = reporting ? runtime.nSurveys++ : NOT_USED;
        runtime.surveyDates.push_back({date, num});
    }
    if (runtime.surveyDates.empty()) {
        throw util::xml_scenario_error("Scenario defines no surveys; at least one is required.");
    }
    if (!runtime.surveyDates.back().isReported()) {
        std::cerr << "Warning: the last survey is unreported. Having surveys beyond the last reported survey is pointless." << std::endl;
    }

    if (util::CommandLine::option(util::CommandLine::PRINT_SURVEY_TIMES)) {
        std::cout << "Survey\tsteps\tdate\n";
        for (const SurveyDate& surveyDate : runtime.surveyDates) {
            if (!surveyDate.isReported()) continue;
            std::cout
                << (surveyDate.num + 1) << '\t'
                << sim::inSteps(surveyDate.date - sim::startDate()) << '\t'
                << surveyDate.date << std::endl;
        }
    }

    runtime.nCohorts = 1;
    if (monitoring.getCohorts().present()) {
        runtime.nCohorts = static_cast<uint32_t>(1) << monitoring.getCohorts().get().getSubPop().size();
    }
    initAgeGroups(monitoring);
    return runtime.surveyDates.back().date;
}

void initCohorts(const scnXml::Monitoring& monitoring)
{
    runtime.cohortSubPopNumbers.clear();
    runtime.cohortSubPopIds.clear();
    if (!monitoring.getCohorts().present()) return;

    const scnXml::Cohorts monCohorts = monitoring.getCohorts().get();
    uint32_t nextId = 0;
    for (auto it = monCohorts.getSubPop().begin(), end = monCohorts.getSubPop().end(); it != end; ++it) {
        const interventions::ComponentId compId = interventions::InterventionManager::getComponentId(it->getId());
        const bool inserted = runtime.cohortSubPopIds.insert(std::make_pair(compId, nextId)).second;
        if (!inserted) {
            throw util::xml_scenario_error(
                std::string("cohort specification uses sub-population \"").append(it->getId()).append("\" more than once"));
        }
        if (it->getNumber() < 0 || notPowerOfTwo(it->getNumber())) {
            throw util::xml_scenario_error(
                std::string("cohort specification assigns sub-population \"").append(it->getId())
                    .append("\" a number which is not a power of 2 (up to 2^21)"));
        }
        runtime.cohortSubPopNumbers.push_back(it->getNumber());
        nextId += 1;
    }
}

void initAgeGroups(const scnXml::Monitoring& monitoring)
{
    const scnXml::MonAgeGroup::GroupSequence& groups = monitoring.getAgeGroup().getGroup();
    if (!(monitoring.getAgeGroup().getLowerbound() <= 0.0)) {
        throw util::xml_scenario_error("Expected survey age-group lowerbound of 0");
    }

    runtime.ageGroupUpperBound.resize(groups.size() + 1);
    for (size_t i = 0; i < groups.size(); ++i) {
        runtime.ageGroupUpperBound[i] = sim::fromYearsD(groups[i].getUpperbound());
    }
    runtime.ageGroupUpperBound[groups.size()] = sim::future();
}

} // namespace mon
} // namespace OM
