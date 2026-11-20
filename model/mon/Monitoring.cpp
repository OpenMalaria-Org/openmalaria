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
#include <type_traits>

namespace OM {
namespace mon {

using internal::Condition;
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

template <typename T>
void writeBinaryValue(ostream& stream, const T& value)
{
    stream.write(reinterpret_cast<const char*>(&value), sizeof(value));
}

} // namespace

namespace internal {

size_t MeasureStore::index(size_t a, size_t c, size_t sp, size_t g, size_t d) const
{
    assert(nAges == 1 || a < nAges);
    assert(nCohorts == 1 || c < nCohorts);
    assert(nSpecies == 1 || sp < nSpecies);
    assert(nGenotypes == 1 || g < nGenotypes);
    assert(nDrugs == 1 || d < nDrugs);
    return (d % nDrugs) + nDrugs *
        ((g % nGenotypes) + nGenotypes *
        ((sp % nSpecies) + nSpecies *
        ((c % nCohorts) + nCohorts *
        (a % nAges))));
}

namespace {

void addStore(const OutMeasure& om, size_t nSpecies, size_t nDrugs, bool forceNoCategories)
{
    assert(om.measure < MeasureCount);
    MeasureStore store;
    store.output = om;
    if (forceNoCategories) store.output.outId = -1;
    if (!forceNoCategories && (om.dims & Dim::Age)) {
        assert(!runtime.ageGroupUpperBound.empty());
    }
    store.nAges = forceNoCategories ? 1 : (om.dims & Dim::Age ? runtime.ageGroupUpperBound.size() : 1);
    store.nCohorts = forceNoCategories ? 1 : (om.dims & Dim::Cohort ? runtime.nCohorts : 1);
    store.nSpecies = forceNoCategories ? 1 : (om.dims & Dim::Species ? nSpecies : 1);
    store.nGenotypes = forceNoCategories ? 1 : (om.dims & Dim::Genotype ? WithinHost::Genotypes::N() : 1);
    store.nDrugs = forceNoCategories ? 1 : (om.dims & Dim::Drug ? nDrugs : 1);
    const size_t reportCount = store.size() * runtime.nSurveys;
    if (om.isDouble) store.reports.emplace<vector<double>>(reportCount, 0.0);
    else store.reports.emplace<vector<int>>(reportCount, 0);
    runtime.storesByMeasure[om.measure].push_back(runtime.stores.size());
    runtime.stores.push_back(std::move(store));
}

} // namespace

void initStores(const vector<OutMeasure>& enabledMeasures, size_t nSpecies, size_t nDrugs)
{
    runtime.stores.clear();
    runtime.storesByMeasure.assign(MeasureCount, {});
    for (const OutMeasure& om : enabledMeasures) addStore(om, nSpecies, nDrugs, false);
}

void ensureConditionStore(const OutMeasure& om)
{
    assert(om.measure < MeasureCount);
    if (runtime.storesByMeasure[om.measure].empty()) addStore(om, 1, 1, true);
}

double surveySum(Measure measure, size_t survey)
{
    assert(measure < runtime.storesByMeasure.size());
    if (runtime.storesByMeasure[measure].empty()) throw SWITCH_DEFAULT_EXCEPTION;
    const MeasureStore& store = runtime.stores[runtime.storesByMeasure[measure].front()];
    const size_t begin = survey * store.size();
    return std::visit([&](const auto& reports) {
        return std::accumulate(reports.begin() + begin, reports.begin() + begin + store.size(), 0.0);
    }, store.reports);
}

RuntimeState runtime;

} // namespace internal

namespace {

void updateConditions()
{
    const size_t survey = statSurveyNumber();
    if (survey == NOT_USED) return;
    for (Condition& cond : runtime.conditions) {
        const double val = internal::surveySum(cond.measure, survey);
        cond.value = (val >= cond.min && val <= cond.max);
    }
}

void writeTextRow(ostream& stream, int survey, int column, int measure, double value, bool isDouble)
{
    stream << survey << '\t' << column << '\t' << measure << '\t';
    if (isDouble) {
        stream << value;
    } else {
        assert(std::trunc(value) == value);
        stream << static_cast<long long>(value);
    }
    stream << '\n';
}

size_t outputAgeGroups(const internal::MeasureStore& store)
{
    return store.nAges == 1 ? 1 : store.nAges - 1;
}

template<class Writer>
void writeMeasure(ostream& stream, Writer writeRow, size_t survey, const internal::MeasureStore& store)
{
    const OutMeasure& om = store.output;
    const size_t surveyStart = survey * store.size();
    const size_t reportCount = std::visit([](const auto& reports) { return reports.size(); }, store.reports);
    assert(reportCount >= surveyStart + store.size());
    const bool bySpecies = om.dims & Dim::Species;
    const bool byDrug = om.dims & Dim::Drug;
    const int ageGroupAdd = om.dims & Dim::Age ? 1 : 0;

    if (bySpecies) assert(store.nAges == 1 && store.nCohorts == 1 && store.nDrugs == 1);
    if (byDrug) assert(store.nSpecies == 1 && store.nGenotypes == 1);
    for (size_t cohortSet = 0; cohortSet < store.nCohorts; ++cohortSet) {
    for (size_t ageGroup = 0; ageGroup < outputAgeGroups(store); ++ageGroup) {
    for (size_t species = 0; species < store.nSpecies; ++species) {
    for (size_t genotype = 0; genotype < store.nGenotypes; ++genotype) {
    for (size_t drug = 0; drug < store.nDrugs; ++drug) {
        const int col2 = bySpecies
            ? species + 1 + 1000000 * genotype
            : ageGroup + ageGroupAdd + 1000 * cohortSetOutputId(cohortSet) + 1000000 * (byDrug ? drug + 1 : genotype);
        const size_t index = surveyStart + store.index(ageGroup, cohortSet, species, genotype, drug);
        const double value = std::visit([index](const auto& reports) {
            return static_cast<double>(reports[index]);
        }, store.reports);
        writeRow(stream, static_cast<int>(survey + 1), col2, om.outId, value, om.isDouble);
    } } } } }
}

enum class ValueType { All, Integer, Double };

template<class Writer>
void writeRows(ostream& stream, Writer writeRow, ValueType type = ValueType::All)
{
    for (size_t survey = 0; survey < runtime.nSurveys; ++survey) {
        for (const internal::MeasureStore& store : runtime.stores) {
            if (store.output.outId < 0 ||
                    (type == ValueType::Integer && store.output.isDouble) ||
                    (type == ValueType::Double && !store.output.isDouble)) continue;
            writeMeasure(stream, writeRow, survey, store);
        }
    }
    if (runtime.reportIMR >= 0 && type != ValueType::Integer) {
        writeRow(stream, 1, 1, runtime.reportIMR, Clinical::InfantMortality::allCause(), true);
    }
}

uint64_t countRows(bool isDouble)
{
    uint64_t count = isDouble && runtime.reportIMR >= 0 ? 1 : 0;
    for (const internal::MeasureStore& store : runtime.stores) {
        if (store.output.outId < 0 || store.output.isDouble != isDouble) continue;
        count += static_cast<uint64_t>(runtime.nSurveys) * outputAgeGroups(store) *
            store.nCohorts * store.nSpecies * store.nGenotypes * store.nDrugs;
    }
    return count;
}

void writeBinary(ostream& stream)
{
    static_assert(std::endian::native == std::endian::little, "output.bin uses little-endian values");
    static_assert(sizeof(double) == 8 && std::numeric_limits<double>::is_iec559, "output.bin needs IEEE-754 doubles");
    static_assert(std::numeric_limits<int>::digits <= std::numeric_limits<int32_t>::digits,
        "integer monitoring values must fit in int32_t");
    const char magic[8] = {'O', 'M', 'O', 'U', 'T', 'B', 'I', 'N'};
    const uint32_t version = 2;
    const uint64_t integerRows = countRows(false);
    const uint64_t doubleRows = countRows(true);
    stream.write(magic, sizeof(magic));
    writeBinaryValue(stream, version);
    writeBinaryValue(stream, integerRows);
    writeBinaryValue(stream, doubleRows);

    auto writeIdColumn = [&](auto select) {
        for (ValueType type : {ValueType::Integer, ValueType::Double}) {
            writeRows(stream, [&](ostream& output, int survey, int column, int measure, double, bool) {
                writeBinaryValue<int32_t>(output, select(survey, column, measure));
            }, type);
        }
    };
    writeIdColumn([](int survey, int, int) { return survey; });
    writeIdColumn([](int, int column, int) { return column; });
    writeIdColumn([](int, int, int measure) { return measure; });
    writeRows(stream, [](ostream& output, int, int, int, double value, bool) {
        assert(std::trunc(value) == value);
        writeBinaryValue<int32_t>(output, static_cast<int32_t>(value));
    }, ValueType::Integer);
    writeRows(stream, [](ostream& output, int, int, int, double value, bool) {
        writeBinaryValue<double>(output, value);
    }, ValueType::Double);
}

void writeText(ostream& stream)
{
    writeRows(stream, writeTextRow);
}

template<class Stream, class Writer>
void writeOutput(const string& filename, ios::openmode mode, Writer writer)
{
    Stream stream(filename.c_str(), mode);
    if (!stream.rdbuf()->is_open()) {
        throw util::base_exception("Unable to open monitoring output file \"" + filename + "\"", util::Error::FileIO);
    }
    writer(stream);
    stream.close();
    if (!stream) {
        throw util::base_exception("Unable to write monitoring output file \"" + filename + "\"", util::Error::FileIO);
    }
}

void updateSurveyNumbers()
{
    runtime.survNumEvent = NOT_USED;
    for (size_t i = runtime.surveyIndex; i < runtime.surveyDates.size(); ++i) {
        runtime.survNumEvent = runtime.surveyDates[i].num;
        if (runtime.survNumEvent != NOT_USED) break;
    }
}

} // namespace
size_t setupCondition(const string& measureName, double minValue, double maxValue, bool initialState)
{
    const OutMeasure om = findOutMeasure(measureName);
    if (om.measure == invalidMeasure) {
        throw util::xml_scenario_error("unrecognised measure: " + measureName);
    }

    if (!isValidCondition(om.measure)) {
        throw util::xml_scenario_error("cannot use measure " + measureName + " as condition of deployment");
    }
    internal::ensureConditionStore(om);

    runtime.conditions.push_back({initialState, om.measure, minValue, maxValue});
    return runtime.conditions.size() - 1;
}

bool checkCondition(size_t conditionKey)
{
    assert(conditionKey < runtime.conditions.size());
    return runtime.conditions[conditionKey].value;
}

static void recordValue(Measure measure, size_t survey, size_t age, uint32_t cohort,
                        size_t species, size_t genotype, size_t drug, double val, int outId = 0)
{
    if (survey == NOT_USED) return;
    assert(measure < runtime.storesByMeasure.size());
    for (size_t idx : runtime.storesByMeasure[measure]) {
        internal::MeasureStore& store = runtime.stores[idx];
        if (outId != 0 && store.output.outId != outId) continue;
        const size_t index = survey * store.size() + store.index(age, cohort, species, genotype, drug);
        std::visit([index, val](auto& reports) {
            using Value = typename std::decay_t<decltype(reports)>::value_type;
            assert(index < reports.size());
            if constexpr (std::is_integral_v<Value>) assert(std::trunc(val) == val);
            reports[index] += static_cast<Value>(val);
        }, store.reports);
    }
}

void recordStat(Measure measure, double val, size_t species, size_t genotype)
{ recordValue(measure, statSurveyNumber(), 0, 0, species, genotype, 0, val); }

void recordStat(Measure measure, const Host::Human& human, double val, size_t species, size_t genotype, size_t drug)
{ recordValue(measure, statSurveyNumber(), human.monitoringAgeGroup, human.getCohortSet(), species, genotype, drug, val); }

void recordEvent(Measure measure, size_t survey, size_t age, uint32_t cohort)
{ recordValue(measure, survey, age, cohort, 0, 0, 0, 1.0); }

void recordEvent(Measure measure, const Host::Human& human, double val, int outId)
{ recordValue(measure, eventSurveyNumber(), human.monitoringAgeGroup, human.getCohortSet(), 0, 0, 0, val, outId); }

void recordDeploy(Measure timedMeasure, Measure ctsMeasure, const Host::Human& human, Deploy::Method method, double val)
{
    assert(method == Deploy::TIMED || method == Deploy::CTS || method == Deploy::TREAT);
    const Measure selected = method == Deploy::TIMED ? timedMeasure :
        method == Deploy::CTS ? ctsMeasure : ::OM::mon::measure("nTreatDeployments");
    recordEvent(selected, human, val);
}

template <typename Stream>
void checkpoint(Stream& stream)
{
    runtime.isInit & stream;
    runtime.surveyIndex & stream;
    runtime.survNumEvent & stream;
    for (internal::MeasureStore& store : runtime.stores) {
        std::visit([&stream](auto& reports) {
            const size_t expectedSize = reports.size();
            size_t storedSize = expectedSize;
            storedSize & stream;
            if (storedSize != expectedSize) {
                throw util::checkpoint_error("Monitoring report size mismatch: checkpoint has " +
                    std::to_string(storedSize) + ", expected " + std::to_string(expectedSize));
            }
            for (auto& report : reports) report & stream;
        }, store.reports);
    }
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
    const bool binary = util::CommandLine::getOutputFormat() == util::CommandLine::OutputFormat::BIN;
    const auto writer = binary ? writeBinary : writeText;
    const ios::openmode mode = ios::out | ios::binary;

    if (util::CommandLine::option(util::CommandLine::COMPRESS_OUTPUT)) {
        writeOutput<ogzstream>(filename, mode, writer);
    } else {
        writeOutput<ofstream>(filename, mode, writer);
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
