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
#include "mon/AgeGroup.h"
#include <cassert>
#include <cstdint>
#include <iosfwd>
#include <limits>
#include <map>
#include <set>
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

// ----- reporting types -----

enum Measure {
    nHost,
    nInfect, nInfect_Imported, nInfect_Introduced, nInfect_Indigenous,
    nPatent, nPatent_Imported, nPatent_Introduced, nPatent_Indigenous,
    totalInfs, totalInfs_Imported, totalInfs_Introduced, totalInfs_Indigenous,
    totalPatentInf, totalPatentInf_Imported, totalPatentInf_Introduced, totalPatentInf_Indigenous,
    nNewInfections, nNewInfections_Imported, nNewInfections_Introduced, nNewInfections_Indigenous,
    nSubPopRemovalFirstEvent,
    nSubPopRemovalTooOld,
    nInfectByGenotype,
    nPatentByGenotype,
    nHostDrugConcNonZero,

    nTreatments1,
    nTreatments2,
    nTreatments3,
    nNMFTreatments,
    nLiverStageTreatments,
    nTreatDiagnostics,

    nUncomp, nUncomp_Imported, nUncomp_Introduced, nUncomp_Indigenous,
    nSevere,
    nSevereWithoutComorbidities,
    nNMFever,

    nDirDeaths,
    nIndDeaths,
    nSeq,
    nHospitalDeaths,
    nHospitalRecovs,
    nHospitalSeqs,
    nNmfDeaths,
    Clinical_FirstDayDeaths,
    Clinical_HospitalFirstDayDeaths,

    vaccinations,
    pev,
    bsv,
    tbv,
    itn,
    irs,
    gvi,
    treat,
    screen,
    recruit,
    nTreatDeployments,

    nCMDTReport,

    nExpectd,
    sumPyrogenThresh,
    sumLogPyrogenThres,
    sumlogDens,
    logDensByGenotype,
    sumAge,
    sumLogDrugConcNonZero,
    expectedDirectDeaths,
    expectedHospitalDeaths,
    expectedIndirectDeaths,
    expectedSequelae,
    expectedSevere,
    expectedSevereWithoutComorbidities,

    nTransmit,
    annAvgK,
    inputEIR,
    simulatedEIR, simulatedEIR_Introduced, simulatedEIR_Indigenous,
    innoculationsPerAgeGroup,
    Vector_Nv0,
    Vector_Nv,
    Vector_Ov,
    Vector_Sv,

    MeasureCount,
    obsoleteMeasure,
    allCauseIMR
};

namespace Deploy {
enum Method {
    NA = 0,
    TIMED = 1 << 0,
    CTS = 1 << 1,
    TREAT = 1 << 2
};
}

// ----- info API -----

namespace impl {
extern size_t nSurveys;
extern size_t nCohorts;
extern bool isInit;
extern size_t survNumEvent, survNumStat;
extern SimTime nextSurveyDate;
}

const size_t NOT_USED = std::numeric_limits<size_t>::max();
const char lineEnd = '\n';

inline size_t eventSurveyNumber() { return impl::survNumEvent; }
inline size_t statSurveyNumber() { return impl::survNumStat; }
inline bool isReported() { return !impl::isInit || impl::survNumStat != NOT_USED; }
inline SimTime nextSurveyDate() { return impl::nextSurveyDate; }
inline size_t numCohortSets() { return impl::nCohorts; }

size_t setupCondition(const std::string& measureName, double minValue, double maxValue, bool initialState);
bool checkCondition(size_t conditionKey);
uint32_t updateCohortSet(uint32_t old, interventions::ComponentId subPop, bool isMember);

// ----- management API -----

SimTime readSurveyDates(const scnXml::Monitoring& monitoring);
void initReporting(const scnXml::Scenario& scenario);
void initCohorts(const scnXml::Monitoring& monitoring);
void initMainSim();
void concludeSurvey();
void writeSurveyData();
template <typename Stream>
void checkpoint(Stream& stream);
extern template void checkpoint<std::ostream>(std::ostream& stream);
extern template void checkpoint<std::istream>(std::istream& stream);

namespace internal {
void write(std::ostream& stream);
uint32_t cohortSetOutputId(uint32_t cohortSet);
}

// ----- direct recording API -----

void record(Measure measure, size_t survey, size_t age, uint32_t cohort, size_t species, size_t genotype, size_t drug, int val, int outId = 0);
void record(Measure measure, size_t survey, size_t age, uint32_t cohort, size_t species, size_t genotype, size_t drug, double val);
void recordStat(Measure measure, const Host::Human& human, int val, size_t species = 0, size_t genotype = 0, size_t drug = 0, int outId = 0);
void recordStat(Measure measure, const Host::Human& human, double val, size_t species = 0, size_t genotype = 0, size_t drug = 0);
void recordEvent(Measure measure, const Host::Human& human, int val);
void recordDeploy(Measure measure, const Host::Human& human, Deploy::Method method, int val = 1);
bool isUsed(Measure measure);

// ----- output-measure registry -----

struct OutMeasure {
    int outId;
    Measure m;
    bool isDouble;
    bool byAge;
    bool byCohort;
    bool bySpecies;
    bool byGenotype;
    bool byDrug;
    uint8_t method;

    OutMeasure()
        : outId(-1), m(MeasureCount), isDouble(false), byAge(false), byCohort(false), bySpecies(false), byGenotype(false), byDrug(false), method(0)
    {
    }

    OutMeasure(int outId, Measure m, bool isDouble, bool byAge, bool byCohort, bool bySpecies, bool byGenotype, bool byDrug, uint8_t method)
        : outId(outId), m(m), isDouble(isDouble), byAge(byAge), byCohort(byCohort), bySpecies(bySpecies), byGenotype(byGenotype), byDrug(byDrug), method(method)
    {
    }

    static OutMeasure value(int outId, Measure m, bool isDouble) { return OutMeasure(outId, m, isDouble, false, false, false, false, false, Deploy::NA); }
    static OutMeasure humanAC(int outId, Measure m, bool isDouble) { return OutMeasure(outId, m, isDouble, true, true, false, false, false, Deploy::NA); }
    static OutMeasure humanACG(int outId, Measure m, bool isDouble) { return OutMeasure(outId, m, isDouble, true, true, false, true, false, Deploy::NA); }
    static OutMeasure humanACP(int outId, Measure m, bool isDouble) { return OutMeasure(outId, m, isDouble, true, true, false, false, true, Deploy::NA); }
    static OutMeasure species(int outId, Measure m, bool byGenotype) { return OutMeasure(outId, m, true, false, false, true, byGenotype, false, Deploy::NA); }
    static OutMeasure humanDeploy(int outId, Measure m, Deploy::Method method)
    {
        assert(method >= 0 && method <= (Deploy::TIMED | Deploy::CTS | Deploy::TREAT));
        return OutMeasure(outId, m, false, true, true, false, false, false, method);
    }
    static OutMeasure obsolete(int outId) { return OutMeasure(outId, obsoleteMeasure, false, false, false, false, false, false, Deploy::NA); }
};

typedef std::map<std::string, OutMeasure> NamedMeasureMapT;
extern NamedMeasureMapT namedOutMeasures;
extern std::set<Measure> validCondMeasures;

void findNamedMeasuresUsing(Measure m, std::ostream& msg);
void defineOutMeasures();

} // namespace mon
} // namespace OM

#endif
