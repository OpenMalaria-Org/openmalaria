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
#include "Host/WithinHost/Genotypes.h"
#include "Clinical/ClinicalModel.h"
#include "Host/Human.h"
#include "util/CommandLine.h"
#include "util/errors.h"
#include "schema/monitoring.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <gzstream/gzstream.h>
#include <iostream>
#include <numeric>

namespace OM {
namespace mon {

using internal::Condition;
using internal::SurveyDate;
using internal::runtime;

namespace {

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

internal::MeasureLayout makeMeasureLayout(const OutMeasure& om, size_t nSpecies, size_t nDrugs, bool forceNoCategories)
{
    internal::MeasureLayout layout;
    layout.outMeasure = om.outId;
    if (!forceNoCategories && hasDim(om.dims, Dim::Age)) {
        assert(!runtime.ageGroupUpperBound.empty());
    }
    layout.nAges = forceNoCategories ? 1 : (hasDim(om.dims, Dim::Age) ? runtime.ageGroupUpperBound.size() : 1);
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

} // namespace

namespace internal {

size_t MeasureLayout::size() const
{
    return nAges * nCohorts * nSpecies * nGenotypes * nDrugs;
}

size_t MeasureLayout::index(size_t a, size_t c, size_t sp, size_t g, size_t d) const
{
#ifndef NDEBUG
    if ((nAges > 1 && a >= nAges) ||
        (nCohorts > 1 && c >= nCohorts) ||
        (nSpecies > 1 && sp >= nSpecies) ||
        (nGenotypes > 1 && g >= nGenotypes) ||
        (nDrugs > 1 && d >= nDrugs))
    {
        cout << "Index out of bounds for age group\t" << a << " of " << nAges
            << "\ncohort set\t" << c << " of " << nCohorts
            << "\nspecies\t" << sp << " of " << nSpecies
            << "\ngenotype\t" << g << " of " << nGenotypes
            << "\ndrug\t" << d << " of " << nDrugs
            << endl;
    }
#endif
    return (d % nDrugs) + nDrugs *
        ((g % nGenotypes) + nGenotypes *
        ((sp % nSpecies) + nSpecies *
        ((c % nCohorts) + nCohorts *
        (a % nAges))));
}

void MeasureLayout::write(ostream& stream, int surveyNum, const OutMeasure& om,
                          const vector<double>& results, size_t surveyStart) const
{
    assert(results.size() >= surveyStart + size());
    const int ageGroupAdd = hasDim(om.dims, Dim::Age) ? 1 : 0;
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

    if (hasDim(om.dims, Dim::Species)) {
        assert(nAges == 1 && nCohorts == 1 && nDrugs == 1);
        for (size_t species = 0; species < nSpecies; ++species) {
        for (size_t genotype = 0; genotype < nGenotypes; ++genotype) {
            const int col2 = species + 1 + 1000000 * genotype;
            emit(0, 0, species, genotype, 0, col2);
        } }
        return;
    }

    if (hasDim(om.dims, Dim::Drug)) {
        assert(nSpecies == 1 && nGenotypes == 1);
        for (size_t cohortSet = 0; cohortSet < nCohorts; ++cohortSet) {
        for (size_t ageGroup = 0; ageGroup < nAgeCats; ++ageGroup) {
        for (size_t drug = 0; drug < nDrugs; ++drug) {
            const int col2 = ageGroup + ageGroupAdd +
                1000 * cohortSetOutputId(cohortSet) +
                1000000 * (drug + 1);
            emit(ageGroup, cohortSet, 0, 0, drug, col2);
        } } }
        return;
    }

    assert(nSpecies == 1 && nDrugs == 1);
    for (size_t cohortSet = 0; cohortSet < nCohorts; ++cohortSet) {
    for (size_t ageGroup = 0; ageGroup < nAgeCats; ++ageGroup) {
    for (size_t genotype = 0; genotype < nGenotypes; ++genotype) {
        const int col2 = ageGroup + ageGroupAdd +
            1000 * cohortSetOutputId(cohortSet) +
            1000000 * genotype;
        emit(ageGroup, cohortSet, 0, genotype, 0, col2);
    } } }
}

void SurveyStore::init(const vector<OutMeasure>& enabledMeasures, size_t nSpecies, size_t nDrugs)
{
    stores.clear();
    measureToStates.assign(MeasureCount, {});
    for (const OutMeasure& om : enabledMeasures) {
        add(om, nSpecies, nDrugs, false);
    }
}

void SurveyStore::ensureConditionState(const OutMeasure& om)
{
    assert(om.m < MeasureCount);
    for (size_t idx : measureToStates[om.m]) {
        if (stores[idx].layout.deployMask == om.method) return;
    }
    add(om, 1, 1, true);
}

void SurveyStore::record(double val, Measure measure, size_t survey, size_t ageIndex, uint32_t cohortSet,
                         size_t species, size_t genotype, size_t drug, int outId)
{
    if (survey == NOT_USED) return;
    assert(measure < measureToStates.size());
    for (size_t idx : measureToStates[measure]) {
        MeasureStore& store = stores[idx];
        if (store.layout.deployMask != Deploy::NA) continue;
        if (outId != 0 && store.layout.outMeasure != outId) continue;
        const size_t offset = survey * store.layout.size();
        const size_t index = offset + store.layout.index(ageIndex, cohortSet, species, genotype, drug);
        assert(index < store.reports.size());
        store.reports[index] += val;
    }
}

void SurveyStore::recordDeploy(int val, Measure measure, size_t survey, size_t ageIndex,
                               uint32_t cohortSet, Deploy::Method method)
{
    if (survey == NOT_USED) return;
    assert(method == Deploy::TIMED || method == Deploy::CTS || method == Deploy::TREAT);
    assert(measure < measureToStates.size());
    for (size_t idx : measureToStates[measure]) {
        MeasureStore& store = stores[idx];
        if ((store.layout.deployMask & method) == Deploy::NA) continue;
        assert(store.layout.nSpecies == 1 && store.layout.nGenotypes == 1);
        const size_t offset = survey * store.layout.size();
        const size_t index = offset + store.layout.index(ageIndex, cohortSet, 0, 0, 0);
        assert(index < store.reports.size());
        store.reports[index] += val;
    }
}

double SurveyStore::sum(Measure measure, uint8_t method, size_t survey) const
{
    assert(measure < measureToStates.size());
    for (size_t idx : measureToStates[measure]) {
        const MeasureStore& store = stores[idx];
        if (store.layout.deployMask != method) continue;
        const size_t begin = survey * store.layout.size();
        const size_t end = begin + store.layout.size();
        return std::accumulate(store.reports.begin() + begin, store.reports.begin() + end, 0.0);
    }
    throw SWITCH_DEFAULT_EXCEPTION;
}

void SurveyStore::write(ostream& stream, size_t survey, const OutMeasure& om) const
{
    assert(om.m < measureToStates.size());
    for (size_t idx : measureToStates[om.m]) {
        const MeasureStore& store = stores[idx];
        if (store.layout.outMeasure != om.outId) continue;
        store.layout.write(stream, survey + 1, om, store.reports, survey * store.layout.size());
        return;
    }
    assert(false && "measure not found in records");
}

bool SurveyStore::uses(Measure measure) const
{
    assert(measure < MeasureCount);
    return !measureToStates[measure].empty();
}

void SurveyStore::checkpoint(ostream& stream)
{
    for (MeasureStore& store : stores) {
        checkpointVec(stream, store.reports, store.layout.size() * runtime.nSurveys);
    }
}

void SurveyStore::checkpoint(istream& stream)
{
    for (MeasureStore& store : stores) {
        checkpointVec(stream, store.reports, store.layout.size() * runtime.nSurveys);
    }
}

void SurveyStore::add(const OutMeasure& om, size_t nSpecies, size_t nDrugs, bool forceNoCategories)
{
    assert(om.m < MeasureCount);
    MeasureStore store;
    store.layout = makeMeasureLayout(om, nSpecies, nDrugs, forceNoCategories);
    store.reports.assign(store.layout.size() * runtime.nSurveys, 0.0);
    measureToStates[om.m].push_back(stores.size());
    stores.push_back(std::move(store));
}

RuntimeState runtime;

} // namespace internal

namespace {

void updateConditions()
{
    assert(runtime.survNumStat != NOT_USED);
    for (Condition& cond : runtime.conditions) {
        const double val = runtime.surveyStore.sum(cond.measure, cond.method, runtime.survNumStat);
        cond.value = (val >= cond.min && val <= cond.max);
    }
}

void write(ostream& stream)
{
    for (size_t survey = 0; survey < runtime.nSurveys; ++survey) {
        for (const OutMeasure& om : runtime.reportedMeasures) {
            runtime.surveyStore.write(stream, survey, om);
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

} // namespace
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
    runtime.surveyStore.ensureConditionState(om);

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
    runtime.surveyStore.record(val, measure, survey, age, cohort, species, genotype, drug, outId);
}

void record(Measure measure, size_t survey, size_t age, uint32_t cohort,
            size_t species, size_t genotype, size_t drug, double val)
{
    runtime.surveyStore.record(val, measure, survey, age, cohort, species, genotype, drug);
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
    runtime.surveyStore.recordDeploy(val, measure, eventSurveyNumber(), human.monitoringAgeGroup, human.getCohortSet(), method);
    const Measure treatDeployments = ::OM::mon::measure("nTreatDeployments");
    if (measure != treatDeployments) {
        runtime.surveyStore.recordDeploy(val, treatDeployments, eventSurveyNumber(), human.monitoringAgeGroup, human.getCohortSet(), method);
    }
}

template <typename Stream>
void checkpoint(Stream& stream)
{
    runtime.isInit & stream;
    runtime.surveyIndex & stream;
    runtime.survNumEvent & stream;
    runtime.survNumStat & stream;
    runtime.nextSurveyDate & stream;
    runtime.surveyStore.checkpoint(stream);
}
template void checkpoint<ostream>(ostream& stream);
template void checkpoint<istream>(istream& stream);
// ———  surveys  ———

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

void updateAgeGroup(size_t& index, SimTime age)
{
    while (age >= runtime.ageGroupUpperBound[index]) {
        ++index;
    }
}

uint32_t updateCohortSet(uint32_t old, interventions::ComponentId subPop, bool isMember)
{
    auto it = runtime.cohortSubPopIds.find(subPop);
    if (it == runtime.cohortSubPopIds.end()) return old;
    const uint32_t subPopId = static_cast<uint32_t>(1) << it->second;
    return (old & ~subPopId) | (isMember ? subPopId : 0);
}

} // namespace mon
} // namespace OM
