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
#include "mon/OutputMeasures.h"
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
#include <type_traits>
#include <variant>

namespace OM {
namespace mon {

NamedMeasureMapT namedOutMeasures;
set<Measure> validCondMeasures;

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
        size_t nAgeCats = nAges == 1 ? 1 : nAges - 1;
        if( om.bySpecies ){
            assert( nAges == 1 && nCohorts == 1 && nDrugs == 1 );
            for( size_t species = 0; species < nSpecies; ++species ){
            for( size_t genotype = 0; genotype < nGenotypes; ++genotype ){
                const int col2 = species + 1 +
                    1000000 * genotype;
                T value = results[surveyStart + index(0, 0, species, genotype, 0)];
                stream << surveyNum << '\t' << col2 << '\t' << om.outId
                    << '\t' << value << lineEnd;
            } }
        }else if( om.byDrug ){
            assert( nSpecies == 1 && nGenotypes == 1 );
            for( size_t cohortSet = 0; cohortSet < nCohorts; ++cohortSet ){
            // Last age category is not reported
            for( size_t ageGroup = 0; ageGroup < nAgeCats; ++ageGroup ){
            for( size_t drug = 0; drug < nDrugs; ++drug ){
                // Yeah, >999 age groups clashes with cohort sets, but unlikely a real issue
                const int col2 = ageGroup + ageGroupAdd +
                    1000 * internal::cohortSetOutputId( cohortSet ) +
                    1000000 * (drug + 1);
                T value = results[surveyStart + index(ageGroup, cohortSet, 0, 0, drug)];
                stream << surveyNum << '\t' << col2 << '\t' << om.outId
                    << '\t' << value << lineEnd;
            } } }
        }else{
            assert( nSpecies == 1 && nDrugs == 1 );
            for( size_t cohortSet = 0; cohortSet < nCohorts; ++cohortSet ){
            // Last age category is not reported
            for( size_t ageGroup = 0; ageGroup < nAgeCats; ++ageGroup ){
            for( size_t genotype = 0; genotype < nGenotypes; ++genotype ){
                // Yeah, >999 age groups clashes with cohort sets, but unlikely a real issue
                const int col2 = ageGroup + ageGroupAdd +
                    1000 * internal::cohortSetOutputId( cohortSet ) +
                    1000000 * genotype;
                T value = results[surveyStart + index(ageGroup, cohortSet, 0, genotype, 0)];
                stream << surveyNum << '\t' << col2 << '\t' << om.outId
                    << '\t' << value << lineEnd;
            } } }
        }
    }
};

struct MeasureState {
    MonIndex layout;
    std::variant<vector<int>, vector<double>> reports;
};

inline size_t stateSurveySize(const MeasureState& state)
{
    return state.layout.size();
}

inline size_t stateSize(const MeasureState& state)
{
    return stateSurveySize(state) * impl::nSurveys;
}

vector<MeasureState> measureStates;
vector<vector<size_t>> measureToStates;

inline bool isDoubleState(const MeasureState& state)
{
    return std::holds_alternative<vector<double>>(state.reports);
}

MeasureState makeState(const OutMeasure& om, size_t nSp, size_t nD, bool forceNoCategories)
{
    MeasureState state;
    state.layout.outMeasure = om.outId;
    state.layout.nAges = forceNoCategories ? 1 : (om.byAge ? AgeGroup::numGroups() : 1);
    state.layout.nCohorts = forceNoCategories ? 1 : (om.byCohort ? impl::nCohorts : 1);
    state.layout.nSpecies = forceNoCategories ? 1 : (om.bySpecies ? nSp : 1);
    state.layout.nGenotypes = forceNoCategories ? 1 : (om.byGenotype ? WithinHost::Genotypes::N() : 1);
    state.layout.nDrugs = forceNoCategories ? 1 : (om.byDrug ? nD : 1);
    state.layout.deployMask = om.method;

    const size_t n = stateSize(state);
    if (om.isDouble) state.reports = vector<double>(n, 0.0);
    else state.reports = vector<int>(n, 0);
    return state;
}

void addState(const OutMeasure& om, size_t nSp, size_t nD, bool forceNoCategories)
{
    assert(om.m < MeasureCount);
    size_t idx = measureStates.size();
    measureStates.push_back(makeState(om, nSp, nD, forceNoCategories));
    measureToStates[om.m].push_back(idx);
}

void initStates(const vector<OutMeasure>& enabledMeasures, size_t nSp, size_t nD)
{
    measureStates.clear();
    measureToStates.assign(MeasureCount, {});
    for (const OutMeasure& om : enabledMeasures)
    {
        if (om.m >= MeasureCount) continue;
        addState(om, nSp, nD, false);
    }
}

void ensureConditionState(const OutMeasure& om)
{
    assert(om.m < MeasureCount);
    for (size_t idx : measureToStates[om.m])
    {
        const MeasureState& state = measureStates[idx];
        if (isDoubleState(state) == om.isDouble && state.layout.deployMask == om.method) return;
    }
    addState(om, 1, 1, true);
}

inline bool acceptsValueType(const MeasureState& state, int) { return std::holds_alternative<vector<int>>(state.reports); }
inline bool acceptsValueType(const MeasureState& state, double) { return std::holds_alternative<vector<double>>(state.reports); }

template <typename T>
inline void addValue(MeasureState& state, size_t index, T val)
{
    auto* values = std::get_if<vector<T>>(&state.reports);
    assert(values != nullptr);
    assert(index < values->size());
    (*values)[index] += val;
}

template <typename T>
void recordValue(T val, Measure measure, size_t survey, size_t ageIndex,
                 uint32_t cohortSet, size_t species, size_t genotype, size_t drug, int outId = 0)
{
    if (survey == NOT_USED) return;
    assert(measure < measureToStates.size());
    for (size_t idx : measureToStates[measure])
    {
        MeasureState& state = measureStates[idx];
        if (!acceptsValueType(state, val)) continue;
        if (state.layout.deployMask != Deploy::NA) continue;
        if (outId != 0 && state.layout.outMeasure != outId) continue;
        size_t i = survey * stateSurveySize(state) + state.layout.index(ageIndex, cohortSet, species, genotype, drug);
        addValue(state, i, val);
    }
}

void recordDeployValue(int val, Measure measure, size_t survey, size_t ageIndex,
                       uint32_t cohortSet, Deploy::Method method)
{
    if (survey == NOT_USED) return;
    assert(method == Deploy::TIMED || method == Deploy::CTS || method == Deploy::TREAT);
    assert(measure < measureToStates.size());
    for (size_t idx : measureToStates[measure])
    {
        MeasureState& state = measureStates[idx];
        if (isDoubleState(state)) continue;
        if ((state.layout.deployMask & method) == Deploy::NA) continue;
        assert(state.layout.nSpecies == 1 && state.layout.nGenotypes == 1);
        size_t i = survey * stateSurveySize(state) + state.layout.index(ageIndex, cohortSet, 0, 0, 0);
        addValue<int>(state, i, val);
    }
}

double sumStateMeasure(Measure measure, bool isDouble, uint8_t method, size_t survey)
{
    assert(survey != NOT_USED);
    assert(measure < measureToStates.size());
    for (size_t idx : measureToStates[measure])
    {
        const MeasureState& state = measureStates[idx];
        if (isDoubleState(state) != isDouble) continue;
        if (state.layout.deployMask != method) continue;
        const size_t off = survey * stateSurveySize(state);
        const size_t end = off + stateSurveySize(state);
        if (isDouble) {
            const auto& values = std::get<vector<double>>(state.reports);
            return std::accumulate(values.begin() + off, values.begin() + end, 0.0);
        } else {
            const auto& values = std::get<vector<int>>(state.reports);
            return std::accumulate(values.begin() + off, values.begin() + end, 0.0);
        }
    }
    throw SWITCH_DEFAULT_EXCEPTION;
}

void writeMeasureState(ostream& stream, size_t survey, const OutMeasure& om)
{
    assert(om.m < measureToStates.size());
    for (size_t idx : measureToStates[om.m])
    {
        const MeasureState& state = measureStates[idx];
        if (state.layout.outMeasure != om.outId) continue;
        std::visit([&](const auto& values) {
            state.layout.write(stream, survey + 1, om, values, survey * stateSurveySize(state));
        }, state.reports);
        return;
    }
    assert(false && "measure not found in records");
}

bool isMeasureUsed(Measure measure, bool isDouble)
{
    assert(measure < measureToStates.size());
    for (size_t idx : measureToStates[measure])
        if (isDoubleState(measureStates[idx]) == isDouble) return true;
    return false;
}

template <typename Stream, typename T>
void checkpointVec(Stream& stream, vector<T>& values, size_t expectedSize)
{
    if constexpr (std::is_base_of_v<std::ostream, std::remove_reference_t<Stream>>) {
        values.size() & stream;
        for (T& y : values) y & stream;
    } else {
        size_t l;
        l & stream;
        if (l != expectedSize) throw util::checkpoint_error("mon::reports: invalid list size");
        values.resize(l);
        for (T& y : values) y & stream;
    }
}

template <typename Stream>
void checkpointStates(Stream& stream)
{
    for (MeasureState& state : measureStates)
    {
        std::visit([&](auto& values) {
            checkpointVec(stream, values, stateSize(state));
        }, state.reports);
    }
}

// Enabled measures:
vector<OutMeasure> reportedMeasures;
int reportIMR = -1; // special output for fitting

struct MeasureByOutId{
    bool operator() (const OutMeasure& i,const OutMeasure& j) {
        return i.outId < j.outId;
    }
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
    recordValue(val, measure, survey, age, cohort, species, genotype, drug, outId);
}

void record(Measure measure, size_t survey, size_t age, uint32_t cohort, size_t species, size_t genotype, size_t drug, double val)
{
    recordValue(val, measure, survey, age, cohort, species, genotype, drug);
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

bool isUsed( Measure measure ){
    return isMeasureUsed(measure, false) || isMeasureUsed(measure, true);
}

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

void updateConditions();        // defined in this file

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
