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
#include <bit>
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

template <typename T>
void writeBinaryValue(ostream& stream, const T& value)
{
    stream.write(reinterpret_cast<const char*>(&value), sizeof(value));
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
                         size_t species, size_t genotype, size_t drug, int outId,
                         Deploy::Method method)
{
    if (survey == NOT_USED) return;
    assert(measure < measureToStates.size());
    for (size_t idx : measureToStates[measure]) {
        MeasureStore& store = stores[idx];
        if (outId != 0 && store.layout.outMeasure != outId) continue;
        if (method == Deploy::NA && store.layout.deployMask != Deploy::NA) continue;
        if (method != Deploy::NA && (store.layout.deployMask & method) == Deploy::NA) continue;
        const size_t offset = survey * store.layout.size();
        const size_t index = offset + store.layout.index(ageIndex, cohortSet, species, genotype, drug);
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
    MeasureLayout& layout = store.layout;
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

void writeRow(ostream& stream, bool binary, int survey, int column, int measure, double value, bool isDouble)
{
    if (binary) {
        writeBinaryValue<int32_t>(stream, survey);
        writeBinaryValue<int32_t>(stream, column);
        writeBinaryValue<int32_t>(stream, measure);
        writeBinaryValue<double>(stream, value);
        return;
    }
    stream << survey << '\t' << column << '\t' << measure << '\t';
    if (isDouble) {
        stream << value;
    } else {
        assert(std::trunc(value) == value);
        stream << static_cast<long long>(value);
    }
    stream << lineEnd;
}

void writeMeasure(ostream& stream, bool binary, size_t survey, const OutMeasure& om)
{
    assert(om.m < runtime.surveyStore.measureToStates.size());
    const int surveyNum = static_cast<int>(survey + 1);
    const internal::MeasureStore* found = nullptr;
    for (size_t idx : runtime.surveyStore.measureToStates[om.m]) {
        const internal::MeasureStore& store = runtime.surveyStore.stores[idx];
        if (store.layout.outMeasure != om.outId) continue;
        found = &store;
        break;
    }
    if (!found) {
        assert(false && "measure not found in records");
        return;
    }

    const internal::MeasureStore& store = *found;
    const internal::MeasureLayout& layout = store.layout;
    const vector<double>& results = store.reports;
    const size_t surveyStart = survey * layout.size();
    assert(results.size() >= surveyStart + layout.size());
    const int ageGroupAdd = hasDim(om.dims, Dim::Age) ? 1 : 0;
    const size_t nAgeCats = layout.nAges == 1 ? 1 : layout.nAges - 1;

    auto writeCell = [&](size_t ageGroup, size_t cohortSet, size_t species, size_t genotype, size_t drug, int col2) {
        const double value = results[surveyStart + layout.index(ageGroup, cohortSet, species, genotype, drug)];
        writeRow(stream, binary, surveyNum, col2, om.outId, value, om.isDouble);
    };

    if (hasDim(om.dims, Dim::Species)) {
        assert(layout.nAges == 1 && layout.nCohorts == 1 && layout.nDrugs == 1);
        for (size_t species = 0; species < layout.nSpecies; ++species) {
        for (size_t genotype = 0; genotype < layout.nGenotypes; ++genotype) {
            writeCell(0, 0, species, genotype, 0, species + 1 + 1000000 * genotype);
        } }
        return;
    }

    if (hasDim(om.dims, Dim::Drug)) {
        assert(layout.nSpecies == 1 && layout.nGenotypes == 1);
        for (size_t cohortSet = 0; cohortSet < layout.nCohorts; ++cohortSet) {
        for (size_t ageGroup = 0; ageGroup < nAgeCats; ++ageGroup) {
        for (size_t drug = 0; drug < layout.nDrugs; ++drug) {
            const int col2 = ageGroup + ageGroupAdd +
                1000 * cohortSetOutputId(cohortSet) +
                1000000 * (drug + 1);
            writeCell(ageGroup, cohortSet, 0, 0, drug, col2);
        } } }
        return;
    }

    assert(layout.nSpecies == 1 && layout.nDrugs == 1);
    for (size_t cohortSet = 0; cohortSet < layout.nCohorts; ++cohortSet) {
    for (size_t ageGroup = 0; ageGroup < nAgeCats; ++ageGroup) {
    for (size_t genotype = 0; genotype < layout.nGenotypes; ++genotype) {
        const int col2 = ageGroup + ageGroupAdd +
            1000 * cohortSetOutputId(cohortSet) +
            1000000 * genotype;
        writeCell(ageGroup, cohortSet, 0, genotype, 0, col2);
    } } }
}

void write(ostream& stream, bool binary)
{
    if (binary) {
        static_assert(std::endian::native == std::endian::little, "output.bin uses little-endian records");
        static_assert(sizeof(double) == 8 && std::numeric_limits<double>::is_iec559, "output.bin needs IEEE-754 doubles");
        const char magic[8] = {'O', 'M', 'O', 'U', 'T', 'B', '1', '\0'};
        const uint32_t version = 1;
        const uint32_t rowSize = 3 * sizeof(int32_t) + sizeof(double);
        stream.write(magic, sizeof(magic));
        writeBinaryValue(stream, version);
        writeBinaryValue(stream, rowSize);
    }
    for (size_t survey = 0; survey < runtime.nSurveys; ++survey) {
        for (const OutMeasure& om : runtime.reportedMeasures) {
            writeMeasure(stream, binary, survey, om);
        }
    }
    if (runtime.reportIMR >= 0) {
        writeRow(stream, binary, 1, 1, runtime.reportIMR, Clinical::InfantMortality::allCause(), true);
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
            size_t species, size_t genotype, size_t drug, double val, int outId, Deploy::Method method)
{
    runtime.surveyStore.record(val, measure, survey, age, cohort, species, genotype, drug, outId, method);
}

void recordStat(Measure measure, const Host::Human& human, double val, size_t species, size_t genotype, size_t drug, int outId)
{
    record(measure, statSurveyNumber(), human.monitoringAgeGroup, human.getCohortSet(), species, genotype, drug, val, outId);
}

void recordEvent(Measure measure, const Host::Human& human, double val)
{
    record(measure, eventSurveyNumber(), human.monitoringAgeGroup, human.getCohortSet(), 0, 0, 0, val);
}

void recordDeploy(Measure measure, const Host::Human& human, Deploy::Method method, double val)
{
    record(measure, eventSurveyNumber(), human.monitoringAgeGroup, human.getCohortSet(), 0, 0, 0, val, 0, method);
    const Measure treatDeployments = ::OM::mon::measure("nTreatDeployments");
    if (measure != treatDeployments)
        record(treatDeployments, eventSurveyNumber(), human.monitoringAgeGroup, human.getCohortSet(), 0, 0, 0, val, 0, method);
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

    if (util::CommandLine::getOutputFormat() == util::CommandLine::OutputFormat::BIN) {
        ofstream stream(filename, mode);
        write(stream, true);
        return;
    }

    if (util::CommandLine::option(util::CommandLine::COMPRESS_OUTPUT)) {
        filename.append(".gz");
        ogzstream stream(filename.c_str(), mode);
        write(stream, false);
    } else {
        ofstream stream(filename, mode);
        write(stream, false);
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
