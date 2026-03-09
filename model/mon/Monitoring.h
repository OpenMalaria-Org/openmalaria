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
#include <cstdint>
#include <iosfwd>
#include <limits>
#include <string>

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

// ----- info API -----

inline constexpr size_t NOT_USED = std::numeric_limits<size_t>::max();
inline constexpr char lineEnd = '\n';

size_t eventSurveyNumber();
size_t statSurveyNumber();
bool isReported();
SimTime nextSurveyDate();
size_t numCohortSets();

size_t setupCondition(const std::string& measureName, double minValue, double maxValue, bool initialState);
bool checkCondition(size_t conditionKey);
uint32_t updateCohortSet(uint32_t old, interventions::ComponentId subPop, bool isMember);

// ----- management API -----

SimTime readSurveyDates(const scnXml::Monitoring& monitoring);
void initAgeGroups(const scnXml::Monitoring& monitoring);
void updateAgeGroup(size_t& index, SimTime age);
size_t numAgeGroups();
void initReporting(const scnXml::Scenario& scenario);
void initCohorts(const scnXml::Monitoring& monitoring);
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
bool isUsed(Measure measure);

} // namespace mon
} // namespace OM

#endif
