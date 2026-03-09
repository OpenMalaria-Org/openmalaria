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

#ifndef H_OM_mon_Monitoring
#define H_OM_mon_Monitoring

#include "Global.h"
#include "mon/OutMeasures.h"
#include <cassert>
#include <cstdint>
#include <iosfwd>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <vector>

namespace scnXml {
class Scenario;
class Monitoring;
}

namespace OM {
namespace Host {
class Human;
}
namespace interventions {
struct ComponentId;
}
namespace mon {

inline constexpr size_t NOT_USED = std::numeric_limits<size_t>::max();
inline constexpr char lineEnd = '\n';

namespace internal {

struct Condition {
    bool value;
    Measure measure;
    uint8_t method;
    double min, max;
};

struct SurveyDate {
    SimTime date = sim::never();
    size_t num = NOT_USED;

    bool isReported() const { return num != NOT_USED; }
};

struct MeasureLayout {
    int outMeasure = 0;
    size_t nAges = 1, nCohorts = 1, nSpecies = 1, nGenotypes = 1, nDrugs = 1;
    uint8_t deployMask = Deploy::NA;

    size_t size() const;
    size_t index(size_t a, size_t c, size_t sp, size_t g, size_t d) const;
    void write(std::ostream& stream, int surveyNum, const OutMeasure& om,
               const std::vector<double>& results, size_t surveyStart) const;
};

struct MeasureStore {
    MeasureLayout layout;
    std::vector<double> reports;
};

struct SurveyStore {
    std::vector<MeasureStore> stores;
    std::vector<std::vector<size_t>> measureToStates;

    void init(const std::vector<OutMeasure>& enabledMeasures, size_t nSpecies, size_t nDrugs);
    void ensureConditionState(const OutMeasure& om);
    void record(double val, Measure measure, size_t survey, size_t ageIndex, uint32_t cohortSet,
                size_t species, size_t genotype, size_t drug, int outId = 0);
    void recordDeploy(int val, Measure measure, size_t survey, size_t ageIndex,
                      uint32_t cohortSet, Deploy::Method method);
    double sum(Measure measure, uint8_t method, size_t survey) const;
    void write(std::ostream& stream, size_t survey, const OutMeasure& om) const;
    bool uses(Measure measure) const;
    void checkpoint(std::ostream& stream);
    void checkpoint(std::istream& stream);

private:
    void add(const OutMeasure& om, size_t nSpecies, size_t nDrugs, bool forceNoCategories);
};

struct RuntimeState {
    bool isInit = false;
    size_t surveyIndex = 0;
    size_t survNumEvent = NOT_USED, survNumStat = NOT_USED;
    SimTime nextSurveyDate = sim::future();
    size_t nSurveys = 0;
    size_t nCohorts = 1;
    NamedMeasureMapT namedOutMeasures;
    std::set<Measure> validCondMeasures;
    std::vector<OutMeasure> reportedMeasures;
    int reportIMR = -1;
    std::vector<Condition> conditions;
    std::vector<SurveyDate> surveyDates;
    std::vector<SimTime> ageGroupUpperBound;
    std::vector<uint32_t> cohortSubPopNumbers;
    std::map<interventions::ComponentId, uint32_t> cohortSubPopIds;
    SurveyStore surveyStore;
};

extern RuntimeState runtime;

} // namespace internal

// ----- info API -----

inline size_t eventSurveyNumber() { return internal::runtime.survNumEvent; }
inline size_t statSurveyNumber() { return internal::runtime.survNumStat; }
inline bool isReported() { return !internal::runtime.isInit || internal::runtime.survNumStat != NOT_USED; }
inline SimTime nextSurveyDate() { return internal::runtime.nextSurveyDate; }
inline size_t numCohortSets() { return internal::runtime.nCohorts; }

size_t setupCondition(const std::string& measureName, double minValue, double maxValue, bool initialState);
bool checkCondition(size_t conditionKey);
uint32_t updateCohortSet(uint32_t old, interventions::ComponentId subPop, bool isMember);

// ----- management API -----

void initAgeGroups(const scnXml::Monitoring& monitoring);
void updateAgeGroup(size_t& index, SimTime age);
void initMainSim();
void concludeSurvey();
void writeSurveyData();

template <typename Stream>
void checkpoint(Stream& stream);
extern template void checkpoint<std::ostream>(std::ostream& stream);
extern template void checkpoint<std::istream>(std::istream& stream);

// ----- direct recording API -----

void record(Measure measure, size_t survey, size_t age, uint32_t cohort, size_t species, size_t genotype, size_t drug, int val, int outId = 0);
void record(Measure measure, size_t survey, size_t age, uint32_t cohort, size_t species, size_t genotype, size_t drug, double val);
void recordStat(Measure measure, const Host::Human& human, int val, size_t species = 0, size_t genotype = 0, size_t drug = 0, int outId = 0);
void recordStat(Measure measure, const Host::Human& human, double val, size_t species = 0, size_t genotype = 0, size_t drug = 0);
void recordEvent(Measure measure, const Host::Human& human, int val);
void recordDeploy(Measure measure, const Host::Human& human, Deploy::Method method, int val = 1);
inline bool isUsed(Measure measure)
{
    return internal::runtime.surveyStore.uses(measure);
}

} // namespace mon
} // namespace OM

#endif
