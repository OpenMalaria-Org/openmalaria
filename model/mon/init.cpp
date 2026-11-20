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
#include <bit>
#include <cctype>
#include <iostream>
#include <map>

namespace OM {
namespace mon {

using internal::SurveyDate;
using internal::runtime;

namespace {

constexpr uint32_t maxCohortNumber = uint32_t{1} << 21;

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

} // namespace

void initReporting(const scnXml::Scenario& scenario)
{
    runtime.reportIMR = -1;
    const scnXml::MonitoringOptions& optsElt = scenario.getMonitoring().getSurveyOptions();
    std::vector<OutMeasure> reportedMeasures;
    reportedMeasures.reserve(optsElt.getOption().size());

    std::set<int> outIds;
    for (const scnXml::MonitoringOption& optElt : optsElt.getOption()) {
        if (!optElt.getValue()) continue;

        OutMeasure om = findOutMeasure(optElt.getName());
        if (om.measure == invalidMeasure) {
            if (isObsoleteMeasure(optElt.getName())) {
                throw util::xml_scenario_error("obsolete survey option: " + std::string(optElt.getName()));
            }
            throw util::xml_scenario_error("unrecognised survey option: " + std::string(optElt.getName()));
        }
        if ((om.measure == measure("sumlogDens") || om.measure == measure("logDensByGenotype")) &&
            WithinHost::diagnostics::monitoringDiagnostic().allowsFalsePositives())
        {
            throw util::xml_scenario_error("measure " + std::string(optElt.getName()) + " may not be used when monitoring diagnostic sensitivity < 1");
        }

        auto applyCategory = [&](Dim flag, const auto& option, const char* label) {
            if (!option.present()) return;
            if (option.get() && !(om.dims & flag)) {
                throw util::xml_scenario_error("measure " + std::string(optElt.getName()) + " does not support categorisation by " + label);
            }
            if (!option.get()) om.dims &= ~flag;
        };
        applyCategory(Dim::Age, optElt.getByAge(), "age group");
        applyCategory(Dim::Cohort, optElt.getByCohort(), "cohort");
        applyCategory(Dim::Species, optElt.getBySpecies(), "species");
        applyCategory(Dim::Genotype, optElt.getByGenotype(), "genotype");
        applyCategory(Dim::Drug, optElt.getByDrugType(), "drug type");

        if (optElt.getOutputNumber().present()) om.outId = optElt.getOutputNumber().get();
        if (om.outId < 0) {
            throw util::xml_scenario_error("monitoring output number must not be negative");
        }
        if (!outIds.insert(om.outId).second) {
            throw util::xml_scenario_error("monitoring output number " + std::to_string(om.outId) + " used more than once");
        }
        if (om.measure == allCauseIMR) {
            runtime.reportIMR = om.outId;
            continue;
        }
        reportedMeasures.push_back(om);
    }
    std::sort(reportedMeasures.begin(), reportedMeasures.end(),
        [](const OutMeasure& lhs, const OutMeasure& rhs) { return lhs.outId < rhs.outId; });
    const size_t nSpecies = scenario.getEntomology().getVector().present()
        ? scenario.getEntomology().getVector().get().getAnopheles().size() : 1;
    const size_t nDrugs = scenario.getPharmacology().present()
        ? scenario.getPharmacology().get().getDrugs().getDrug().size() : 1;
    internal::initStores(reportedMeasures, nSpecies, nDrugs);
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
        const auto number = it->getNumber();
        if (number < 1 || number > maxCohortNumber ||
            !std::has_single_bit(static_cast<uint32_t>(number))) {
            throw util::xml_scenario_error(
                std::string("cohort specification assigns sub-population \"").append(it->getId())
                    .append("\" a number which is not a power of 2 between 1 and ")
                    .append(std::to_string(maxCohortNumber)));
        }
        runtime.cohortSubPopNumbers.push_back(number);
        nextId += 1;
    }
}

} // namespace mon
} // namespace OM
