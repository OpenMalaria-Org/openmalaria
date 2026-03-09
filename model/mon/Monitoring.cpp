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
#include <cmath>
#include <fstream>
#include <gzstream/gzstream.h>
#include <iostream>
#include <numeric>

namespace OM {
namespace mon {

namespace {

using interventions::ComponentId;

struct Condition {
    bool value; // whether the condition was satisfied during the last survey
    Measure measure;
    uint8_t method;
    double min, max;
};

struct SurveyDate {
    SimTime date = sim::never();
    size_t num = NOT_USED;

    bool isReported() const { return num != NOT_USED; }
};

struct RuntimeState {
    bool isInit = false;
    size_t surveyIndex = 0;
    size_t survNumEvent = NOT_USED, survNumStat = NOT_USED;
    SimTime nextSurveyDate = sim::future();
    size_t nSurveys = 0;
    size_t nCohorts = 1;
    NamedMeasureMapT namedOutMeasures;
    set<Measure> validCondMeasures;
    vector<OutMeasure> reportedMeasures;
    int reportIMR = -1;
    vector<Condition> conditions;
    vector<SurveyDate> surveyDates;
    vector<SimTime> ageGroupUpperBound;
    vector<uint32_t> cohortSubPopNumbers;
    map<ComponentId, uint32_t> cohortSubPopIds;
};

RuntimeState runtime;
uint32_t cohortSetOutputId(uint32_t cohortSet){
    uint32_t outNum = 0;
    assert( (cohortSet >> runtime.cohortSubPopNumbers.size()) == 0 );
    for( uint32_t i = 0; i < runtime.cohortSubPopNumbers.size(); ++i ){
        if( cohortSet & (static_cast<uint32_t>(1) << i) ){
            outNum += runtime.cohortSubPopNumbers[i];
        }
    }
    return outNum;
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
    void write( ostream& stream, int surveyNum, const OutMeasure& om,
            const vector<double>& results, size_t surveyStart ) const
    {
        assert(results.size() >= surveyStart + size());
        // First age group starts at 1, unless there isn't an age group:
        const int ageGroupAdd = hasDim(om.dims, Dim::Age) ? 1 : 0;
        // Number of *reported* age categories: either no categorisation (1) or there is an extra unreported category
        const size_t nAgeCats = nAges == 1 ? 1 : nAges - 1;

        auto emit = [&](size_t ageGroup, size_t cohortSet, size_t species, size_t genotype, size_t drug, int col2) {
            const double value = results[surveyStart + index(ageGroup, cohortSet, species, genotype, drug)];
            stream << surveyNum << '\t' << col2 << '\t' << om.outId << '\t';
            if (om.isDouble) {
                stream << value;
            } else {
                assert(std::trunc(value) == value);
                stream << static_cast<long long>(value);
            }
            stream << lineEnd;
        };

        if( hasDim(om.dims, Dim::Species) ){
            assert( nAges == 1 && nCohorts == 1 && nDrugs == 1 );
            for( size_t species = 0; species < nSpecies; ++species ){
            for( size_t genotype = 0; genotype < nGenotypes; ++genotype ){
                const int col2 = species + 1 +
                    1000000 * genotype;
                emit(0, 0, species, genotype, 0, col2);
            } }
            return;
        }

        if( hasDim(om.dims, Dim::Drug) ){
            assert( nSpecies == 1 && nGenotypes == 1 );
            for( size_t cohortSet = 0; cohortSet < nCohorts; ++cohortSet ){
            // Last age category is not reported
            for( size_t ageGroup = 0; ageGroup < nAgeCats; ++ageGroup ){
            for( size_t drug = 0; drug < nDrugs; ++drug ){
                // Yeah, >999 age groups clashes with cohort sets, but unlikely a real issue
                const int col2 = ageGroup + ageGroupAdd +
                    1000 * cohortSetOutputId( cohortSet ) +
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
                1000 * cohortSetOutputId( cohortSet ) +
                1000000 * genotype;
            emit(ageGroup, cohortSet, 0, genotype, 0, col2);
        } } }
    }
};

struct State {
    MonIndex layout;
    vector<double> reports;
};

MonIndex makeLayout(const OutMeasure& om, size_t nSpecies, size_t nDrugs, bool forceNoCategories)
{
    MonIndex layout;
    layout.outMeasure = om.outId;
    layout.nAges = forceNoCategories ? 1 : (hasDim(om.dims, Dim::Age) ? numAgeGroups() : 1);
    layout.nCohorts = forceNoCategories ? 1 : (hasDim(om.dims, Dim::Cohort) ? runtime.nCohorts : 1);
    layout.nSpecies = forceNoCategories ? 1 : (hasDim(om.dims, Dim::Species) ? nSpecies : 1);
    layout.nGenotypes = forceNoCategories ? 1 : (hasDim(om.dims, Dim::Genotype) ? WithinHost::Genotypes::N() : 1);
    layout.nDrugs = forceNoCategories ? 1 : (hasDim(om.dims, Dim::Drug) ? nDrugs : 1);
    layout.deployMask = om.method;
    return layout;
}

void checkpointVec(ostream& stream, vector<double>& values, size_t /*expectedSize*/)
{
    values.size() & stream;
    for (double& y : values) y & stream;
}

void checkpointVec(istream& stream, vector<double>& values, size_t expectedSize)
{
    size_t storedSize = 0;
    storedSize & stream;
    if (storedSize != expectedSize) {
        throw util::checkpoint_error("mon::reports: invalid list size");
    }
    values.resize(storedSize);
    for (double& y : values) y & stream;
}

struct StateRegistry {
    vector<State> states;
    vector<vector<size_t>> measureToStates;

    void init(const vector<OutMeasure>& enabledMeasures, size_t nSpecies, size_t nDrugs)
    {
        states.clear();
        measureToStates.assign(MeasureCount, {});
        for (const OutMeasure& om : enabledMeasures) {
            add(om, nSpecies, nDrugs, false);
        }
    }

    void ensureConditionState(const OutMeasure& om)
    {
        assert(om.m < MeasureCount);
        for (size_t idx : measureToStates[om.m]) {
            if (states[idx].layout.deployMask == om.method) return;
        }
        add(om, 1, 1, true);
    }

    void record(double val, Measure measure, size_t survey, size_t ageIndex, uint32_t cohortSet,
                size_t species, size_t genotype, size_t drug, int outId = 0)
    {
        if (survey == NOT_USED) return;
        assert(measure < measureToStates.size());
        for (size_t idx : measureToStates[measure]) {
            State& state = states[idx];
            if (state.layout.deployMask != Deploy::NA) continue;
            if (outId != 0 && state.layout.outMeasure != outId) continue;
            const size_t offset = survey * state.layout.size();
            const size_t index = offset + state.layout.index(ageIndex, cohortSet, species, genotype, drug);
            assert(index < state.reports.size());
            state.reports[index] += val;
        }
    }

    void recordDeploy(int val, Measure measure, size_t survey, size_t ageIndex,
                      uint32_t cohortSet, Deploy::Method method)
    {
        if (survey == NOT_USED) return;
        assert(method == Deploy::TIMED || method == Deploy::CTS || method == Deploy::TREAT);
        assert(measure < measureToStates.size());
        for (size_t idx : measureToStates[measure]) {
            State& state = states[idx];
            if ((state.layout.deployMask & method) == Deploy::NA) continue;
            assert(state.layout.nSpecies == 1 && state.layout.nGenotypes == 1);
            const size_t offset = survey * state.layout.size();
            const size_t index = offset + state.layout.index(ageIndex, cohortSet, 0, 0, 0);
            assert(index < state.reports.size());
            state.reports[index] += val;
        }
    }

    double sum(Measure measure, uint8_t method, size_t survey) const
    {
        assert(measure < measureToStates.size());
        for (size_t idx : measureToStates[measure]) {
            const State& state = states[idx];
            if (state.layout.deployMask != method) continue;
            const size_t begin = survey * state.layout.size();
            const size_t end = begin + state.layout.size();
            return std::accumulate(state.reports.begin() + begin, state.reports.begin() + end, 0.0);
        }
        throw SWITCH_DEFAULT_EXCEPTION;
    }

    void write(ostream& stream, size_t survey, const OutMeasure& om) const
    {
        assert(om.m < measureToStates.size());
        for (size_t idx : measureToStates[om.m]) {
            const State& state = states[idx];
            if (state.layout.outMeasure != om.outId) continue;
            state.layout.write(stream, survey + 1, om, state.reports, survey * state.layout.size());
            return;
        }
        assert(false && "measure not found in records");
    }

    bool uses(Measure measure) const
    {
        assert(measure < MeasureCount);
        return !measureToStates[measure].empty();
    }

    template <typename Stream>
    void checkpoint(Stream& stream)
    {
        for (State& state : states) {
            checkpointVec(stream, state.reports, state.layout.size() * runtime.nSurveys);
        }
    }

private:
    void add(const OutMeasure& om, size_t nSpecies, size_t nDrugs, bool forceNoCategories)
    {
        assert(om.m < MeasureCount);
        State state;
        state.layout = makeLayout(om, nSpecies, nDrugs, forceNoCategories);
        state.reports.assign(state.layout.size() * runtime.nSurveys, 0.0);
        measureToStates[om.m].push_back(states.size());
        states.push_back(std::move(state));
    }
};

StateRegistry stateRegistry;

void updateConditions()
{
    assert(runtime.survNumStat != NOT_USED);
    for (Condition& cond : runtime.conditions) {
        const double val = stateRegistry.sum(cond.measure, cond.method, runtime.survNumStat);
        cond.value = (val >= cond.min && val <= cond.max);
    }
}

void write(ostream& stream)
{
    for (size_t survey = 0; survey < runtime.nSurveys; ++survey) {
        for (const OutMeasure& om : runtime.reportedMeasures) {
            stateRegistry.write(stream, survey, om);
        }
    }
    if (runtime.reportIMR >= 0) {
        stream << 1 << "\t" << 1 << "\t" << runtime.reportIMR
            << "\t" << Clinical::InfantMortality::allCause() << lineEnd;
    }
}

void updateSurveyNumbers()
{
    if (runtime.surveyIndex >= runtime.surveyDates.size()) {
        runtime.survNumEvent = NOT_USED;
        runtime.survNumStat = NOT_USED;
        runtime.nextSurveyDate = sim::future();
        return;
    }

    for (size_t i = runtime.surveyIndex; i < runtime.surveyDates.size(); ++i) {
        runtime.survNumEvent = runtime.surveyDates[i].num;
        if (runtime.survNumEvent != NOT_USED) break;
    }

    const SurveyDate& nextSurvey = runtime.surveyDates[runtime.surveyIndex];
    runtime.survNumStat = nextSurvey.num;
    runtime.nextSurveyDate = nextSurvey.date;
}

bool notPowerOfTwo(uint32_t num)
{
    return num == 0 || num > (static_cast<uint32_t>(1) << 21) || (num & (num - 1)) != 0;
}

} // namespace

size_t eventSurveyNumber() { return runtime.survNumEvent; }
size_t statSurveyNumber() { return runtime.survNumStat; }
bool isReported() { return !runtime.isInit || runtime.survNumStat != NOT_USED; }
SimTime nextSurveyDate() { return runtime.nextSurveyDate; }
size_t numCohortSets() { return runtime.nCohorts; }

void initReporting(const scnXml::Scenario& scenario)
{
    defineOutMeasures(runtime.namedOutMeasures, runtime.validCondMeasures);
    assert(runtime.reportedMeasures.empty());
    runtime.reportIMR = -1;
    const scnXml::MonitoringOptions& optsElt = scenario.getMonitoring().getSurveyOptions();
    runtime.reportedMeasures.reserve(optsElt.getOption().size() + runtime.namedOutMeasures.size());
    auto applyCategory = [](Dim& dims, Dim flag, const bool optionalPresent, const bool requested,
                            const string& optionName, const char* label)
    {
        if (!optionalPresent) return;
        const bool supports = hasDim(dims, flag);
        if (supports) {
            if (!requested) clearDim(dims, flag);
        } else if (requested) {
            throw util::xml_scenario_error("measure " + optionName + " does not support categorisation by " + label);
        }
    };

    set<int> outIds;
    for (const scnXml::MonitoringOption& optElt : optsElt.getOption()) {
        if (!optElt.getValue()) continue;

        auto it = runtime.namedOutMeasures.find(optElt.getName());
        if (it == runtime.namedOutMeasures.end()) {
            throw util::xml_scenario_error("unrecognised survey option: " + string(optElt.getName()));
        }
        OutMeasure om = it->second;
        if (om.m >= MeasureCount) {
            if (om.m == obsoleteMeasure) {
                throw util::xml_scenario_error("obsolete survey option: " + string(optElt.getName()));
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
                throw util::xml_scenario_error("monitoring output number " + to_string(om.outId) + " used more than once");
            }
            outIds.insert(om.outId);
            runtime.reportIMR = om.outId;
            continue;
        }

        if ((om.m == measure("sumlogDens") || om.m == measure("logDensByGenotype")) &&
            WithinHost::diagnostics::monitoringDiagnostic().allowsFalsePositives())
        {
            throw util::xml_scenario_error("measure " + string(optElt.getName()) + " may not be used when monitoring diagnostic sensitivity < 1");
        }
        const string optionName(optElt.getName());
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
            throw util::xml_scenario_error("monitoring output number " + to_string(om.outId) + " used more than once");
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
    stateRegistry.init(runtime.reportedMeasures, nSpecies, nDrugs);
}

size_t setupCondition(const string& measureName, double minValue, double maxValue, bool initialState)
{
    auto it = runtime.namedOutMeasures.find(measureName);
    if (it == runtime.namedOutMeasures.end()) {
        throw util::xml_scenario_error("unrecognised measure: " + measureName);
    }

    const OutMeasure om = it->second;
    if (runtime.validCondMeasures.count(om.m) == 0) {
        throw util::xml_scenario_error("cannot use measure " + measureName + " as condition of deployment");
    }
    stateRegistry.ensureConditionState(om);

    runtime.conditions.push_back({initialState, om.m, om.method, minValue, maxValue});
    return runtime.conditions.size() - 1;
}

bool checkCondition(size_t conditionKey)
{
    assert(conditionKey < runtime.conditions.size());
    return runtime.conditions[conditionKey].value;
}

void record(Measure measure, size_t survey, size_t age, uint32_t cohort,
            size_t species, size_t genotype, size_t drug, int val, int outId)
{
    stateRegistry.record(val, measure, survey, age, cohort, species, genotype, drug, outId);
}

void record(Measure measure, size_t survey, size_t age, uint32_t cohort,
            size_t species, size_t genotype, size_t drug, double val)
{
    stateRegistry.record(val, measure, survey, age, cohort, species, genotype, drug);
}

void recordStat(Measure measure, const Host::Human& human, int val, size_t species, size_t genotype, size_t drug, int outId)
{
    record(measure, statSurveyNumber(), human.monitoringAgeGroup, human.getCohortSet(), species, genotype, drug, val, outId);
}

void recordStat(Measure measure, const Host::Human& human, double val, size_t species, size_t genotype, size_t drug)
{
    record(measure, statSurveyNumber(), human.monitoringAgeGroup, human.getCohortSet(), species, genotype, drug, val);
}

void recordEvent(Measure measure, const Host::Human& human, int val)
{
    record(measure, eventSurveyNumber(), human.monitoringAgeGroup, human.getCohortSet(), 0, 0, 0, val);
}

void recordDeploy(Measure measure, const Host::Human& human, Deploy::Method method, int val)
{
    stateRegistry.recordDeploy(val, measure, eventSurveyNumber(), human.monitoringAgeGroup, human.getCohortSet(), method);
    const Measure treatDeployments = ::OM::mon::measure("nTreatDeployments");
    if (measure != treatDeployments) {
        stateRegistry.recordDeploy(val, treatDeployments, eventSurveyNumber(), human.monitoringAgeGroup, human.getCohortSet(), method);
    }
}

bool isUsed(Measure measure)
{
    return stateRegistry.uses(measure);
}

template <typename Stream>
void checkpoint(Stream& stream)
{
    runtime.isInit & stream;
    runtime.surveyIndex & stream;
    runtime.survNumEvent & stream;
    runtime.survNumStat & stream;
    runtime.nextSurveyDate & stream;
    stateRegistry.checkpoint(stream);
}
template void checkpoint<ostream>(ostream& stream);
template void checkpoint<istream>(istream& stream);
// ———  surveys  ———

SimTime readSurveyDates(const scnXml::Monitoring& monitoring)
{
    const scnXml::Surveys::SurveyTimeSequence& survs = monitoring.getSurveys().getSurveyTime();
    map<SimTime, bool> surveys;

    auto addSurvey = [&surveys](SimTime date, bool reporting) {
        if (reporting) surveys[date] = true;
        else surveys.insert(make_pair(date, false));
    };

    for (size_t i = 0; i < survs.size(); ++i) {
        const scnXml::SurveyTime& surv = survs[i];
        try {
            string s = surv;
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
            throw util::xml_scenario_error(string("surveyTime: ").append(e.message()));
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

void initMainSim()
{
    runtime.surveyIndex = 0;
    runtime.isInit = true;
    updateSurveyNumbers();
}

void concludeSurvey()
{
    updateConditions();
    runtime.surveyIndex += 1;
    updateSurveyNumbers();
}

void writeSurveyData()
{
    string filename = util::CommandLine::getOutputName();
    const auto mode = std::ios::out | std::ios::binary;

    if (util::CommandLine::option(util::CommandLine::COMPRESS_OUTPUT)) {
        filename.append(".gz");
        ogzstream stream(filename.c_str(), mode);
        write(stream);
    } else {
        ofstream stream(filename, mode);
        write(stream);
    }
}
// ———  Age groups  ———

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

size_t numAgeGroups()
{
    assert(!runtime.ageGroupUpperBound.empty());
    return runtime.ageGroupUpperBound.size();
}

void updateAgeGroup(size_t& index, SimTime age)
{
    while (age >= runtime.ageGroupUpperBound[index]) {
        ++index;
    }
}
// ———  Cohort sets  ———

void initCohorts(const scnXml::Monitoring& monitoring)
{
    runtime.cohortSubPopNumbers.clear();
    runtime.cohortSubPopIds.clear();
    if (!monitoring.getCohorts().present()) return;

    const scnXml::Cohorts monCohorts = monitoring.getCohorts().get();
    uint32_t nextId = 0;
    for (auto it = monCohorts.getSubPop().begin(), end = monCohorts.getSubPop().end(); it != end; ++it) {
        const ComponentId compId = interventions::InterventionManager::getComponentId(it->getId());
        const bool inserted = runtime.cohortSubPopIds.insert(make_pair(compId, nextId)).second;
        if (!inserted) {
            throw util::xml_scenario_error(
                string("cohort specification uses sub-population \"").append(it->getId()).append("\" more than once"));
        }
        if (it->getNumber() < 0 || notPowerOfTwo(it->getNumber())) {
            throw util::xml_scenario_error(
                string("cohort specification assigns sub-population \"").append(it->getId())
                    .append("\" a number which is not a power of 2 (up to 2^21)"));
        }
        runtime.cohortSubPopNumbers.push_back(it->getNumber());
        nextId += 1;
    }
}

uint32_t updateCohortSet(uint32_t old, ComponentId subPop, bool isMember)
{
    auto it = runtime.cohortSubPopIds.find(subPop);
    if (it == runtime.cohortSubPopIds.end()) return old;
    const uint32_t subPopId = static_cast<uint32_t>(1) << it->second;
    return (old & ~subPopId) | (isMember ? subPopId : 0);
}

} // namespace mon
} // namespace OM
