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
#include <string>
#include <variant>
#include <vector>

namespace OM {
namespace Host {
class Human;
}
namespace interventions {
struct ComponentId;
}
namespace mon {

inline constexpr size_t NOT_USED = std::numeric_limits<size_t>::max();

namespace Deploy {
enum Method { TIMED, CTS, TREAT };
}

namespace internal {

struct Condition {
    bool value;
    Measure measure;
    double min, max;
};

struct SurveyDate {
    SimTime date = sim::never();
    size_t num = NOT_USED;

    bool isReported() const { return num != NOT_USED; }
};

struct MeasureStore {
    OutMeasure output;
    size_t nAges = 1, nCohorts = 1, nSpecies = 1, nGenotypes = 1, nDrugs = 1;
    std::variant<std::vector<int>, std::vector<double>> reports;

    size_t size() const { return nAges * nCohorts * nSpecies * nGenotypes * nDrugs; }
    size_t index(size_t a, size_t c, size_t sp, size_t g, size_t d) const;
};

struct RuntimeState {
    bool isInit = false;
    size_t surveyIndex = 0;
    size_t survNumEvent = NOT_USED;
    size_t nSurveys = 0;
    size_t nCohorts = 1;
    int reportIMR = -1;
    std::vector<Condition> conditions;
    std::vector<SurveyDate> surveyDates;
    std::vector<SimTime> ageGroupUpperBound;
    std::vector<uint32_t> cohortSubPopNumbers;
    std::map<interventions::ComponentId, uint32_t> cohortSubPopIds;
    std::vector<MeasureStore> stores;
    std::vector<std::vector<size_t>> storesByMeasure;
};

extern RuntimeState runtime;
void initStores(const std::vector<OutMeasure>& enabledMeasures, size_t nSpecies, size_t nDrugs);

} // namespace internal

// ----- info API -----

inline size_t eventSurveyNumber() { return internal::runtime.survNumEvent; }
inline size_t statSurveyNumber() {
    const auto& r = internal::runtime;
    return r.isInit && r.surveyIndex < r.surveyDates.size() ? r.surveyDates[r.surveyIndex].num : NOT_USED;
}
inline bool isReported() { return !internal::runtime.isInit || statSurveyNumber() != NOT_USED; }
inline SimTime nextSurveyDate() {
    const auto& r = internal::runtime;
    return r.isInit && r.surveyIndex < r.surveyDates.size() ? r.surveyDates[r.surveyIndex].date : sim::future();
}
size_t setupCondition(const std::string& measureName, double minValue, double maxValue, bool initialState);
bool checkCondition(size_t conditionKey);
uint32_t updateCohortSet(uint32_t old, interventions::ComponentId subPop, bool isMember);

// ----- management API -----

void updateAgeGroup(size_t& index, SimTime age);
void initMainSim();
void concludeSurvey();
void writeSurveyData();

template <typename Stream>
void checkpoint(Stream& stream);
extern template void checkpoint<std::ostream>(std::ostream& stream);
extern template void checkpoint<std::istream>(std::istream& stream);

// ----- direct recording API -----

void recordStat(Measure measure, double val, size_t species = 0, size_t genotype = 0);
void recordStat(Measure measure, const Host::Human& human, double val = 1.0, size_t species = 0, size_t genotype = 0, size_t drug = 0);
void recordEvent(Measure measure, size_t survey, size_t age, uint32_t cohort);
void recordEvent(Measure measure, const Host::Human& human, double val = 1.0, int outId = 0);
void recordDeploy(Measure timedMeasure, Measure ctsMeasure, const Host::Human& human, Deploy::Method method, double val = 1.0);
inline bool isUsed(Measure measure)
{
    assert(measure < MeasureCount);
    return !internal::runtime.storesByMeasure[measure].empty();
}

} // namespace mon
} // namespace OM

#endif
