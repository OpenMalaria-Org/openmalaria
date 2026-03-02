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
    MHR_HOSTS,
    MHR_INFECTED_HOSTS, MHR_INFECTED_HOSTS_IMPORTED, MHR_INFECTED_HOSTS_INTRODUCED, MHR_INFECTED_HOSTS_INDIGENOUS,
    MHR_PATENT_HOSTS, MHR_PATENT_HOSTS_IMPORTED, MHR_PATENT_HOSTS_INTRODUCED, MHR_PATENT_HOSTS_INDIGENOUS,
    MHR_INFECTIONS, MHR_INFECTIONS_IMPORTED, MHR_INFECTIONS_INTRODUCED, MHR_INFECTIONS_INDIGENOUS,
    MHR_PATENT_INFECTIONS, MHR_PATENT_INFECTIONS_IMPORTED, MHR_PATENT_INFECTIONS_INTRODUCED, MHR_PATENT_INFECTIONS_INDIGENOUS,
    MHR_NEW_INFECTIONS, MHR_NEW_INFECTIONS_IMPORTED, MHR_NEW_INFECTIONS_INTRODUCED, MHR_NEW_INFECTIONS_INDIGENOUS,
    MHR_SUB_POP_REM_FIRST_EVENT,
    MHR_SUB_POP_REM_TOO_OLD,
    MHR_INFECTED_GENOTYPE,
    MHR_PATENT_GENOTYPE,
    MHR_HOSTS_POS_DRUG_CONC,

    MHT_TREATMENTS_1,
    MHT_TREATMENTS_2,
    MHT_TREATMENTS_3,
    MHT_NMF_TREATMENTS,
    MHT_LS_TREATMENTS,
    MHT_TREAT_DIAGNOSTICS,

    MHE_UNCOMPLICATED_EPISODES, MHE_UNCOMPLICATED_EPISODES_IMPORTED, MHE_UNCOMPLICATED_EPISODES_INTRODUCED, MHE_UNCOMPLICATED_EPISODES_INDIGENOUS,
    MHE_SEVERE_EPISODES,
    MHE_SEVERE_EPISODES_WITHOUT_COMORBIDITIES,
    MHE_NON_MALARIA_FEVERS,

    MHO_DIRECT_DEATHS,
    MHO_INDIRECT_DEATHS,
    MHO_SEQUELAE,
    MHO_HOSPITAL_DEATHS,
    MHO_HOSPITAL_RECOVERIES,
    MHO_HOSPITAL_SEQUELAE,
    MHO_NMF_DEATHS,
    MHO_FIRST_DAY_DEATHS,
    MHO_HOSPITAL_FIRST_DAY_DEATHS,

    MHD_VACCINATIONS,
    MHD_PEV,
    MHD_BSV,
    MHD_TBV,
    MHD_ITN,
    MHD_IRS,
    MHD_GVI,
    MHD_TREAT,
    MHD_SCREEN,
    MHD_RECRUIT,
    MHD_ALL_DEPLOYS,

    MCD_CMDT_REPORT,

    MHF_EXPECTED_INFECTED,
    MHF_PYROGENIC_THRESHOLD,
    MHF_LOG_PYROGENIC_THRESHOLD,
    MHF_LOG_DENSITY,
    MHF_LOG_DENSITY_GENOTYPE,
    MHF_AGE,
    MHF_LOG_DRUG_CONC,
    MHF_EXPECTED_DIRECT_DEATHS,
    MHF_EXPECTED_HOSPITAL_DEATHS,
    MHF_EXPECTED_INDIRECT_DEATHS,
    MHF_EXPECTED_SEQUELAE,
    MHF_EXPECTED_SEVERE,
    MHF_EXPECTED_SEVERE_WITHOUT_COMORBIDITIES,

    MVF_NUM_TRANSMIT,
    MVF_ANN_AVG_K,
    MVF_INPUT_EIR,
    MVF_SIM_EIR, MVF_SIM_EIR_INTRODUCED, MVF_SIM_EIR_INDIGENOUS,
    MVF_INOCS,
    MVF_LAST_NV0,
    MVF_LAST_NV,
    MVF_LAST_OV,
    MVF_LAST_SV,

    M_NUM,
    M_OBSOLETE,
    M_ALL_CAUSE_IMR
};

namespace measure {
constexpr Measure nHost = MHR_HOSTS;
constexpr Measure nInfect = MHR_INFECTED_HOSTS;
constexpr Measure nInfect_Imported = MHR_INFECTED_HOSTS_IMPORTED;
constexpr Measure nInfect_Introduced = MHR_INFECTED_HOSTS_INTRODUCED;
constexpr Measure nInfect_Indigenous = MHR_INFECTED_HOSTS_INDIGENOUS;
constexpr Measure nPatent = MHR_PATENT_HOSTS;
constexpr Measure nPatent_Imported = MHR_PATENT_HOSTS_IMPORTED;
constexpr Measure nPatent_Introduced = MHR_PATENT_HOSTS_INTRODUCED;
constexpr Measure nPatent_Indigenous = MHR_PATENT_HOSTS_INDIGENOUS;
constexpr Measure totalInfs = MHR_INFECTIONS;
constexpr Measure totalInfs_Imported = MHR_INFECTIONS_IMPORTED;
constexpr Measure totalInfs_Introduced = MHR_INFECTIONS_INTRODUCED;
constexpr Measure totalInfs_Indigenous = MHR_INFECTIONS_INDIGENOUS;
constexpr Measure totalPatentInf = MHR_PATENT_INFECTIONS;
constexpr Measure totalPatentInf_Imported = MHR_PATENT_INFECTIONS_IMPORTED;
constexpr Measure totalPatentInf_Introduced = MHR_PATENT_INFECTIONS_INTRODUCED;
constexpr Measure totalPatentInf_Indigenous = MHR_PATENT_INFECTIONS_INDIGENOUS;
constexpr Measure nNewInfections = MHR_NEW_INFECTIONS;
constexpr Measure nNewInfections_Imported = MHR_NEW_INFECTIONS_IMPORTED;
constexpr Measure nNewInfections_Introduced = MHR_NEW_INFECTIONS_INTRODUCED;
constexpr Measure nNewInfections_Indigenous = MHR_NEW_INFECTIONS_INDIGENOUS;
constexpr Measure nSubPopRemovalFirstEvent = MHR_SUB_POP_REM_FIRST_EVENT;
constexpr Measure nSubPopRemovalTooOld = MHR_SUB_POP_REM_TOO_OLD;
constexpr Measure nInfectByGenotype = MHR_INFECTED_GENOTYPE;
constexpr Measure nPatentByGenotype = MHR_PATENT_GENOTYPE;
constexpr Measure nHostDrugConcNonZero = MHR_HOSTS_POS_DRUG_CONC;

constexpr Measure nTreatments1 = MHT_TREATMENTS_1;
constexpr Measure nTreatments2 = MHT_TREATMENTS_2;
constexpr Measure nTreatments3 = MHT_TREATMENTS_3;
constexpr Measure nLiverStageTreatments = MHT_LS_TREATMENTS;
constexpr Measure nTreatDiagnostics = MHT_TREAT_DIAGNOSTICS;

constexpr Measure nUncomp = MHE_UNCOMPLICATED_EPISODES;
constexpr Measure nUncomp_Imported = MHE_UNCOMPLICATED_EPISODES_IMPORTED;
constexpr Measure nUncomp_Introduced = MHE_UNCOMPLICATED_EPISODES_INTRODUCED;
constexpr Measure nUncomp_Indigenous = MHE_UNCOMPLICATED_EPISODES_INDIGENOUS;
constexpr Measure nSevere = MHE_SEVERE_EPISODES;
constexpr Measure nSevereWithoutComorbidities = MHE_SEVERE_EPISODES_WITHOUT_COMORBIDITIES;
constexpr Measure nNMFever = MHE_NON_MALARIA_FEVERS;

constexpr Measure nDirDeaths = MHO_DIRECT_DEATHS;
constexpr Measure nIndDeaths = MHO_INDIRECT_DEATHS;
constexpr Measure nSeq = MHO_SEQUELAE;
constexpr Measure nHospitalDeaths = MHO_HOSPITAL_DEATHS;
constexpr Measure nHospitalRecovs = MHO_HOSPITAL_RECOVERIES;
constexpr Measure nHospitalSeqs = MHO_HOSPITAL_SEQUELAE;
constexpr Measure nNmfDeaths = MHO_NMF_DEATHS;
constexpr Measure Clinical_FirstDayDeaths = MHO_FIRST_DAY_DEATHS;
constexpr Measure Clinical_HospitalFirstDayDeaths = MHO_HOSPITAL_FIRST_DAY_DEATHS;

constexpr Measure vaccinations = MHD_VACCINATIONS;
constexpr Measure pev = MHD_PEV;
constexpr Measure bsv = MHD_BSV;
constexpr Measure tbv = MHD_TBV;
constexpr Measure itn = MHD_ITN;
constexpr Measure irs = MHD_IRS;
constexpr Measure gvi = MHD_GVI;
constexpr Measure treat = MHD_TREAT;
constexpr Measure screen = MHD_SCREEN;
constexpr Measure recruit = MHD_RECRUIT;
constexpr Measure nTreatDeployments = MHD_ALL_DEPLOYS;

constexpr Measure nCMDTReport = MCD_CMDT_REPORT;

constexpr Measure nExpectd = MHF_EXPECTED_INFECTED;
constexpr Measure sumPyrogenThresh = MHF_PYROGENIC_THRESHOLD;
constexpr Measure sumLogPyrogenThres = MHF_LOG_PYROGENIC_THRESHOLD;
constexpr Measure sumlogDens = MHF_LOG_DENSITY;
constexpr Measure logDensByGenotype = MHF_LOG_DENSITY_GENOTYPE;
constexpr Measure sumAge = MHF_AGE;
constexpr Measure sumLogDrugConcNonZero = MHF_LOG_DRUG_CONC;
constexpr Measure expectedDirectDeaths = MHF_EXPECTED_DIRECT_DEATHS;
constexpr Measure expectedHospitalDeaths = MHF_EXPECTED_HOSPITAL_DEATHS;
constexpr Measure expectedIndirectDeaths = MHF_EXPECTED_INDIRECT_DEATHS;
constexpr Measure expectedSequelae = MHF_EXPECTED_SEQUELAE;
constexpr Measure expectedSevere = MHF_EXPECTED_SEVERE;
constexpr Measure expectedSevereWithoutComorbidities = MHF_EXPECTED_SEVERE_WITHOUT_COMORBIDITIES;

constexpr Measure nTransmit = MVF_NUM_TRANSMIT;
constexpr Measure annAvgK = MVF_ANN_AVG_K;
constexpr Measure inputEIR = MVF_INPUT_EIR;
constexpr Measure simulatedEIR = MVF_SIM_EIR;
constexpr Measure simulatedEIR_Introduced = MVF_SIM_EIR_INTRODUCED;
constexpr Measure simulatedEIR_Indigenous = MVF_SIM_EIR_INDIGENOUS;
constexpr Measure innoculationsPerAgeGroup = MVF_INOCS;
constexpr Measure Vector_Nv0 = MVF_LAST_NV0;
constexpr Measure Vector_Nv = MVF_LAST_NV;
constexpr Measure Vector_Ov = MVF_LAST_OV;
constexpr Measure Vector_Sv = MVF_LAST_SV;
}

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
void checkpoint(std::ostream& stream);
void checkpoint(std::istream& stream);

namespace internal {
void write(std::ostream& stream);
uint32_t cohortSetOutputId(uint32_t cohortSet);
}

// ----- direct recording API -----

struct Key {
    size_t survey = NOT_USED;
    size_t age = 0;
    uint32_t cohort = 0;
    size_t species = 0;
    size_t genotype = 0;
    size_t drug = 0;

    inline Key& withAge(size_t v) { age = v; return *this; }
    inline Key& withCohort(uint32_t v) { cohort = v; return *this; }
    inline Key& withSpecies(size_t v) { species = v; return *this; }
    inline Key& withGenotype(size_t v) { genotype = v; return *this; }
    inline Key& withDrug(size_t v) { drug = v; return *this; }
};

Key statKey();
Key eventKey();
Key surveyKey(size_t survey);
Key humanStatKey(const Host::Human& human);
Key humanEventKey(const Host::Human& human);

void record(Measure measure, const Key& key, int val, int outId = 0);
void record(Measure measure, const Key& key, double val);
void recordDeploy(Measure measure, const Key& key, Deploy::Method method, int val = 1);
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
        : outId(-1), m(M_NUM), isDouble(false), byAge(false), byCohort(false), bySpecies(false), byGenotype(false), byDrug(false), method(0)
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
    static OutMeasure obsolete(int outId) { return OutMeasure(outId, M_OBSOLETE, false, false, false, false, false, false, Deploy::NA); }
};

typedef std::map<std::string, OutMeasure> NamedMeasureMapT;
extern NamedMeasureMapT namedOutMeasures;
extern std::set<Measure> validCondMeasures;

void findNamedMeasuresUsing(Measure m, std::ostream& msg);
void defineOutMeasures();

} // namespace mon
} // namespace OM

#endif
