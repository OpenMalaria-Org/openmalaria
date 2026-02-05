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

#include "mon/info.h"
#include "mon/reporting.h"
#include "mon/management.h"
#define H_OM_mon_cpp
#include "mon/OutputMeasures.h"
#include "Host/WithinHost/Diagnostic.h"
#include "Host/WithinHost/Genotypes.h"
#include "Clinical/ClinicalModel.h"
#include "Host/Human.h"
#include "util/errors.h"
#include "schema/scenario.h"

#include <typeinfo>
#include <iostream>

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
    // Measure which this accepts.
    Measure measure;
    // Measure number used in output file
    int outMeasure;
    // Index of first item in result array
    size_t offset;
    // Number of categories. Must be > 0. If 1, index is set to zero, otherwise
    // indices *should* be less than this.
    // 
    // nAges may include a final, unreported category.
    size_t nAges, nCohorts, nSpecies, nGenotypes, nDrugs;
    // Either Deploy::NA (not tracking deployments) or a binary 'or' of at
    // least one of Deploy::TIMED, Deploy::CTS, Deploy::TREAT.
    uint8_t deployMask;
    
    // Used to calculate next offset. This is max output of `index(...)` + 1.
    inline size_t size() const{
        return nAges * nCohorts * nSpecies * nGenotypes * nDrugs;
    }
    // Get the index in the result array to store this data at
    // (age group, cohort, species, genotype, drug).
    // 
    // First index is `self.offset`, last is `self.offset + self.size() - 1`.
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
        return offset +
            (d % nDrugs) + nDrugs *
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

struct MonIndByMeasure{
    bool operator() (const MonIndex& i,const MonIndex& j) {
        if( i.measure != j.measure ) return i.measure < j.measure;
        else return i.size() < j.size();
    }
} monIndByMeasure;

// Store data of type T which is to be reported
template<typename T>
class Store{
public:
    Store() : surveySize(0) {}
    
private:
    // This lists all enabled outputs, sorted by `measure` (first field, of
    // type `mon::Measure`), ideally with the least-subcategorised first.
    vector<MonIndex> measures;
    
    // This provides a fast way of finding items in `measures` which
    // accepts data of this measure type (MHR_HOSTS, etc.). Is set to the pair
    // `(0, 0)` for unused measures.
    typedef pair<uint16_t, uint16_t> MeasureRange;
    vector<MeasureRange> measure_map;
    
    // Number of indices in `reports` used by a single survey
    size_t surveySize;
    // These are the stored reports (multidimensional; size is `size()` and
    // indices are `survey * surveySize + measures[m].index(...)` for some `m`).
    vector<T> reports;
    
    // get size of reports
    inline size_t size(){ return surveySize * impl::nSurveys; }
    
public:
    // Set up ready to accept reports. The passed list includes all measures
    // used; we ignore those of the wrong type.
    void init( const vector<OutMeasure>& enabledMeasures, size_t nSp, size_t nD ){
        for( const OutMeasure& om : enabledMeasures ){
            // Two types: double and int. Skip if type is wrong.
            if( om.isDouble != (typeid(T) == typeid(double)) ) continue;
            // We don't track weird measures like IMR here:
            if( om.m >= M_NUM ) continue;
            
            MonIndex m;
            m.measure = om.m;
            m.outMeasure = om.outId;
            m.offset = 0;
            m.nAges = om.byAge ? (AgeGroup::numGroups()) : 1;
            m.nCohorts = om.byCohort ? impl::nCohorts : 1;
            m.nSpecies = om.bySpecies ? nSp : 1;
            m.nGenotypes = om.byGenotype ? WithinHost::Genotypes::N() : 1;
            m.nDrugs = om.byDrug ? nD : 1;
            m.deployMask = om.method;
            measures.push_back(m);
        }
        
        sortEnabledMeasures();
        
        // Leave a few spare slots for potential conditions using variables not already reported:
        reports.reserve(size() + 12);
        reports.assign(size(), 0);
    }
    
    // Enable reporting by an additional measure, which does not categorise.
    // (Called after init(); does nothing if this measure is already enabled.)
    // 
    // This *must* be called before any reporting takes place, since it adjusts
    // offsets in `reports`.
    void enableCondition( const OutMeasure& om ){
        assert( om.isDouble ? typeid(T) == typeid(double) : typeid(T) == typeid(int) );
        assert( om.m < M_NUM );
        
        // Skip if we already track this measure:
        for( size_t i = measure_map[om.m].first, end = measure_map[om.m].second;
            i < end; ++i )
        {
            const MonIndex& ind = measures[i];
            assert(ind.measure == om.m);
            if( ind.deployMask == om.method ){
                return;
            }
        }
        
        MonIndex m;
        m.measure = om.m;
        m.outMeasure = om.outId;
        m.offset = 0;
        m.nAges = 1;
        m.nCohorts = 1;
        m.nSpecies = 1;
        m.nGenotypes = 1;
        m.nDrugs = 1;
        m.deployMask = om.method;
        measures.push_back(m);
        
        sortEnabledMeasures();
        reports.resize(size(), 0);
    }
    
    // Sort measures, then fix the offsets and surveySize, then set measure_map
    void sortEnabledMeasures() {
        std::sort( measures.begin(), measures.end(), monIndByMeasure );
        measure_map.assign(M_NUM, make_pair(0, 0));
        surveySize = 0;
        for( size_t i = 0; i < measures.size(); ++i ){
            measures[i].offset = surveySize;
            surveySize += measures[i].size();
            
            Measure m = measures[i].measure;

            assert(m < measure_map.size());
            // Indices are initialised to zero, but zero may also be correct; measures[0] will be valid:
            if( measures.at(measure_map[m].first).measure != m )
                measure_map[m].first = i;
            measure_map[m].second = i + 1;
        }
    }
    
    // Take a reported value and either store it or forget it.
    // If some of ageIndex, cohortSet, species are not applicable, use 0.
    void report( T val, Measure measure, size_t survey, size_t ageIndex,
                 uint32_t cohortSet, size_t species, size_t genotype, size_t drug, int outId = 0)
    {
        if( survey == NOT_USED ) return; // pre-main-sim & unit tests we ignore all reports
        assert(measure < measure_map.size());
        for( size_t i = measure_map[measure].first, end = measure_map[measure].second;
            i < end; ++i )
        {
            assert(i < measures.size());
            const MonIndex& ind = measures[i];
            assert(ind.measure == measure);
            if( ind.deployMask != Deploy::NA ) continue;        // skip measures tracking deployments
            if( outId != 0 && ind.outMeasure != outId) continue;     // skip if supplied outID is different
            size_t index = survey * surveySize +
                    ind.index(ageIndex, cohortSet, species, genotype, drug);
            assert( index < reports.size() );
            reports[index] += val;
        }
    }
    
    // Take a deployment report and potentially store it in one or more places
    // If some ageIndex or cohortSet are not applicable, use 0.
    void deploy( T val, Measure measure, size_t survey, size_t ageIndex,
                 uint32_t cohortSet, Deploy::Method method )
    {
        if( survey == NOT_USED ) return; // pre-main-sim & unit tests we ignore all reports
        assert( method == Deploy::TIMED ||
            method == Deploy::CTS || method == Deploy::TREAT );
        assert(measure < measure_map.size());
        for( size_t i = measure_map[measure].first, end = measure_map[measure].second;
            i < end; ++i )
        {
            assert(i < measures.size());
            const MonIndex& ind = measures[i];
            assert(ind.measure == measure);
            // skip measures not tracking deployments or not tracking this type of deployment
            if( (ind.deployMask & method) == Deploy::NA ) continue;
            assert( ind.nSpecies == 1 && ind.nGenotypes == 1 );     // never used for deployments
            
            size_t index = survey * surveySize +
                    ind.index(ageIndex, cohortSet, 0, 0, 0);
            assert( index < reports.size() );
            reports[index] += val;
        }
    }
    
    /// Get the sum of all reported values for some measure, method and survey.
    /// 
    /// Method may be a bit-or-ed combination of Deploy flags, but must exactly
    /// match some measure being recorded.
    T get_sum( Measure measure, uint8_t method, size_t survey ){
        assert( survey != NOT_USED );
        // We use the first compatible measure
        assert(measure < measure_map.size());
        for( size_t i = measure_map[measure].first, end = measure_map[measure].second;
            i < end; ++i )
        {
            assert(i < measures.size());
            const MonIndex& ind = measures[i];
            assert(ind.measure == measure);
            if( ind.deployMask != method ) continue;    // incompatible deployment mode: skip
            
            const size_t off = survey * surveySize + ind.offset;
            T sum = 0;
            size_t end2 = off + ind.size();
            assert(end2 <= reports.size());
            for( size_t i = off; i < end2; ++i ){
                sum += reports[i];
            }
            return sum;
        }
        // Whoops, didn't find a compatible measure! Program *should* ensure there is always one.
        throw SWITCH_DEFAULT_EXCEPTION;
    }
    
    // Return true if reports by this measure are recorded, false if they are discarded.
    bool isUsed( Measure measure ){
        assert( measure < measure_map.size() );
        return measure_map[measure].second > measure_map[measure].first;
    }
    
    // Write stored values to stream for some output measure, om
    void write( ostream& stream, size_t survey, const OutMeasure& om ){
        assert(om.m < measure_map.size());
        for( size_t i = measure_map[om.m].first, end = measure_map[om.m].second;
            i < end; ++i )
        {
            assert(i < measures.size());
            if( measures[i].outMeasure == om.outId ){
                measures[i].write( stream, survey + 1, om, reports, survey * surveySize );
                return;
            }
        }
        assert(false && "measure not found in records");
    }
    
    // Checkpointing
    void checkpoint( ostream& stream ){
        reports.size() & stream;
        for (T& y : reports) {
            y & stream;
        }
        // reports is the only field which changes after initialisation
    }
    void checkpoint( istream& stream ){
        size_t l;
        l & stream;
        if( l != size() ){
            throw util::checkpoint_error( "mon::reports: invalid list size" );
        }
        reports.resize (l);
        for (T& y : reports) {
            y & stream;
        }
        // reports is the only field which changes after initialisation
    }
};

// Enabled measures:
vector<OutMeasure> reportedMeasures;
// Stores of reported data by two different types:
Store<int> storeI;
Store<double> storeF;
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
    
    set<int> outIds;    // all measure numbers used in output
    for( const scnXml::MonitoringOption& optElt : optsElt.getOption() ){
        if( optElt.getValue() == false ) continue;      // option is disabled
        
        auto it = namedOutMeasures.find( optElt.getName() );
        if( it == namedOutMeasures.end() ){
            throw util::xml_scenario_error("unrecognised survey option: " + string(optElt.getName()));
        }
        OutMeasure om = it->second;     // copy; we may modify below
        
        if( om.m >= M_NUM ){
            if( om.m == M_ALL_CAUSE_IMR ){
                if( om.isDouble && !om.byAge && !om.byCohort && !om.bySpecies ){
                    reportIMR = om.outId;
                }else{
                    throw util::xml_scenario_error( "measure allCauseIMR does not support any categorisation" );
                }
            } else if( om.m == M_OBSOLETE ){
                throw util::xml_scenario_error("obsolete survey option: " + string(optElt.getName()));
            } else TRACED_EXCEPTION_DEFAULT("invalid measure code");
        }
        
        if( om.m == MHF_LOG_DENSITY || om.m == MHF_LOG_DENSITY_GENOTYPE){
            if( WithinHost::diagnostics::monitoringDiagnostic().allowsFalsePositives() ){
                throw util::xml_scenario_error("measure " + string(optElt.getName()) + " may not be used when monitoring diagnostic sensitivity < 1");
            }
        }
        
        // Categorisation can be disabled but not enabled. We check each now.
        if( optElt.getByAge().present() ){
            if( om.byAge ){
                om.byAge = optElt.getByAge().get();     // disable or keep
            }else if( optElt.getByAge().get() ){
                throw util::xml_scenario_error("measure " + string(optElt.getName()) + " does not support categorisation by age group");
            }
        }
        if( optElt.getByCohort().present() ){
            if( om.byCohort ){
                om.byCohort = optElt.getByCohort().get();   // disable or keep
            }else if( optElt.getByCohort().get() ){
                throw util::xml_scenario_error("measure " + string(optElt.getName()) + " does not support categorisation by cohort");
            }
        }
        if( optElt.getBySpecies().present() ){
            if( om.bySpecies ){
                om.bySpecies = optElt.getBySpecies().get(); // disable or keep
            }else if( optElt.getBySpecies().get() ){
                throw util::xml_scenario_error("measure " + string(optElt.getName()) + " does not support categorisation by species");
            }
        }
        if( optElt.getByGenotype().present() ){
            if( om.byGenotype ){
                om.byGenotype = optElt.getByGenotype().get();   // disable or keep
            }else if( optElt.getByGenotype().get() ){
                throw util::xml_scenario_error("measure " + string(optElt.getName()) + " does not support categorisation by genotype");
            }
        }
        if( optElt.getByDrugType().present() ){
            if( om.byDrug ){
                om.byDrug = optElt.getByDrugType().get();   // disable or keep
            }else if( optElt.getByDrugType().get() ){
                throw util::xml_scenario_error("measure " + string(optElt.getName()) + " does not support categorisation by drug type");
            }
        }
        
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
    
    storeI.init( reportedMeasures, nSpecies, nDrugs );
    storeF.init( reportedMeasures, nSpecies, nDrugs );
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
    if( om.isDouble ) storeF.enableCondition(om);
    else storeI.enableCondition(om);
    
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
        double val = cond.isDouble ?
            storeF.get_sum( cond.measure, cond.method, impl::survNumStat ) :
            storeI.get_sum( cond.measure, cond.method, impl::survNumStat );
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
            if( om.m >= M_NUM ){
                // "Special" measures are not reported this way. The only such measure is IMR.
                assert( om.m == M_ALL_CAUSE_IMR && reportIMR >= 0 );
                continue;
            } else if( om.isDouble ) {
                storeF.write( stream, survey, om );
            } else {
                storeI.write( stream, survey, om );
            }
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

// Report functions: each reports to all usable stores (i.e. correct data type
// and where parameters don't have to be fabricated).
// void reportMI( Measure measure, int val ){
//     storeI.report( val, measure, impl::currentSurvey, 0, 0, 0, 0, 0 );
// }
void reportEventMHI( Measure measure, const Host::Human& human, int val ){
    const size_t survey = impl::survNumEvent;
    const size_t ageIndex = human.monitoringAgeGroup.i();
    storeI.report( val, measure, survey, ageIndex, human.getCohortSet(), 0, 0, 0 );
}
void reportStatMHI( Measure measure, const Host::Human& human, int val ){
    const size_t survey = impl::survNumStat;
    const size_t ageIndex = human.monitoringAgeGroup.i();
    storeI.report( val, measure, survey, ageIndex, human.getCohortSet(), 0, 0, 0 );
}
void reportEventMHI_CMDT( Measure measure, const Host::Human& human, int val, int outId ){
    const size_t survey = impl::survNumStat;
    const size_t ageIndex = human.monitoringAgeGroup.i();
    storeI.report( val, measure, survey, ageIndex, human.getCohortSet(), 0, 0, 0, outId );
}
void reportMSACI( Measure measure, size_t survey,
                  AgeGroup ageGroup, uint32_t cohortSet, int val )
{
    storeI.report( val, measure, survey, ageGroup.i(), cohortSet, 0, 0, 0 );
}
void reportStatMHGI( Measure measure, const Host::Human& human, size_t genotype,
                 int val )
{
    const size_t survey = impl::survNumStat;
    const size_t ageIndex = human.monitoringAgeGroup.i();
    storeI.report( val, measure, survey, ageIndex, human.getCohortSet(), 0, genotype, 0 );
}
void reportStatMHPI( Measure measure, const Host::Human& human, size_t drugIndex,
                int val )
{
    const size_t survey = impl::survNumStat;
    const size_t ageIndex = human.monitoringAgeGroup.i();
    storeI.report( val, measure, survey, ageIndex, human.getCohortSet(), 0, 0, drugIndex );
}
// Deployment reporting uses a different function to handle the method
// (mostly to make other types of report faster).
void reportEventMHD( Measure measure, const Host::Human& human,
                Deploy::Method method )
{
    const int val = 1;  // always report 1 deployment
    const size_t survey = impl::survNumEvent;
    size_t ageIndex = human.monitoringAgeGroup.i();
    storeI.deploy( val, measure, survey, ageIndex, human.getCohortSet(), method );
    // This is for nTreatDeployments:
    measure = MHD_ALL_DEPLOYS;
    storeI.deploy( val, measure, survey, ageIndex, human.getCohortSet(), method );
}

void reportStatMF( Measure measure, double val ){
    storeF.report( val, measure, impl::survNumStat, 0, 0, 0, 0, 0 );
}
void reportStatMHF( Measure measure, const Host::Human& human, double val ){
    const size_t survey = impl::survNumStat;
    const size_t ageIndex = human.monitoringAgeGroup.i();
    storeF.report( val, measure, survey, ageIndex, human.getCohortSet(), 0, 0, 0 );
}
void reportStatMACGF( Measure measure, size_t ageIndex, uint32_t cohortSet,
                  size_t genotype, double val )
{
    const size_t survey = impl::survNumStat;
    storeF.report( val, measure, survey, ageIndex, cohortSet, 0, genotype, 0 );
}
void reportStatMACSGF( Measure measure, size_t ageIndex, uint32_t cohortSet,
                  size_t species, size_t genotype, double val )
{
    const size_t survey = impl::survNumStat;
    storeF.report( val, measure, survey, ageIndex, cohortSet, species, genotype, 0 );
}
void reportStatMHPF( Measure measure, const Host::Human& human, size_t drug, double val ){
    const size_t survey = impl::survNumStat;
    const size_t ageIndex = human.monitoringAgeGroup.i();
    storeF.report( val, measure, survey, ageIndex, human.getCohortSet(), 0, 0, drug );
}
void reportStatMHGF( Measure measure, const Host::Human& human, size_t genotype,
                 double val )
{
    reportStatMACGF( measure, human.monitoringAgeGroup.i(), human.getCohortSet(),
                 genotype, val );
}
void reportStatMSF( Measure measure, size_t species, double val ){
    const size_t survey = impl::survNumStat;
    storeF.report( val, measure, survey, 0, 0, species, 0, 0 );
}
void reportStatMSGF( Measure measure, size_t species, size_t genotype, double val ){
    const size_t survey = impl::survNumStat;
    storeF.report( val, measure, survey, 0, 0, species, genotype, 0 );
}

bool isUsedM( Measure measure ){
    return storeI.isUsed(measure) || storeF.isUsed(measure);
}

void checkpoint( ostream& stream ){
    impl::isInit & stream;
    impl::surveyIndex & stream;
    impl::survNumEvent & stream;
    impl::survNumStat & stream;
    impl::nextSurveyDate & stream;
    
    storeI.checkpoint(stream);
    storeF.checkpoint(stream);
}
void checkpoint( istream& stream ){
    impl::isInit & stream;
    impl::surveyIndex & stream;
    impl::survNumEvent & stream;
    impl::survNumStat & stream;
    impl::nextSurveyDate & stream;
    
    storeI.checkpoint(stream);
    storeF.checkpoint(stream);
}

}
}
