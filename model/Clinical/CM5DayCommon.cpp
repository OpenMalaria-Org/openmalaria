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

#include "Clinical/DecisionTree5Day.h"
#include "Clinical/CMDecisionTree.h"
#include "Host/WithinHost/WHInterface.h"
#include "Host/WithinHost/WHVivax.h"
#include "util/errors.h"
#include "util/ModelOptions.h"
#include "util/random.h"

namespace OM {
namespace Clinical {
using namespace ::OM::util;

// ———  static pararmeters and methods  ———

// These parameters are set by setHealthSystem() and do not need checkpointing.
mon::Measure CM5DayCommon::measures[NumCaseTypes] = {
    mon::MHT_TREATMENTS_1,  // first line official
    mon::MHT_TREATMENTS_2   // second line official
};
double CM5DayCommon::accessUCAny[NumCaseTypes];
double CM5DayCommon::accessUCSelfTreat[NumCaseTypes];
double CM5DayCommon::accessSevere;
double CM5DayCommon::cureRateSevere;
WithinHost::TreatmentId CM5DayCommon::treatmentSevere;

bool cfr_pf_use_hospital = false;

void CM5DayCommon::init(){
    cfr_pf_use_hospital = util::ModelOptions::option( util::CFR_PF_USE_HOSPITAL );
}


// ———  per-human, construction and destruction  ———

CM5DayCommon::CM5DayCommon (double tSF) :
        m_treatmentSeekingFactor (tSF)
{}


// ———  per-human, update  ———

void CM5DayCommon::doClinicalUpdate (Human& human, double ageYears, WithinHost::Pathogenesis::StatePair &pg) {
    Episode::State pgState = static_cast<Episode::State>( pg.state );

    if (pgState & Episode::MALARIA) {
        if (pgState & Episode::COMPLICATED){
            severeMalaria (human, pgState, ageYears, doomed);
        }else if (indirectMortBugfix || !pg.indirectMortality) {
            // NOTE: the "not indirect mortality" bit is a historical accident.
            // Validity is debatable, but there's no point changing now.
            // (This does affect tests.)
            uncomplicatedEvent (human, pgState);
        }
    } else if (pgState & Episode::SICK) { // sick but not from malaria
        uncomplicatedEvent (human, pgState);
    }
    
    if (pg.indirectMortality && doomed == NOT_DOOMED)
        doomed = -sim::oneTS();
    
    if( m_tLastTreatment == sim::ts0() )
        human.removeFirstEvent(interventions::SubPopRemove::ON_FIRST_TREATMENT);
    
    if( pgState & Episode::SICK )
        human.removeFirstEvent(interventions::SubPopRemove::ON_FIRST_BOUT);
}

void CM5DayCommon::severeMalaria (
    Human &human,
    Episode::State pgState,
    double ageYears,
    int& doomed
) {
    // Probability of getting treatment (only part which is case managment):
    double p2 = accessSevere * m_treatmentSeekingFactor;
    // Probability of getting cured after getting treatment:
    double p3 = cureRateSevere;
    // p4 is the hospital case-fatality rate from Tanzania
    double p4 = caseFatalityRate.eval (ageYears);
    // p5a here is the community threshold case-fatality rate
    double p5a = getCommunityCFR (p4);
    // p5b here is the in-hospital treatment-failure case-fatality rate
    // In our code it was originally the community CFR, but in the published model description it should be the hospital CFR
    double p5b = (cfr_pf_use_hospital) ? p4 : p5a;
    // p6 is P(seq) for treated patients
    double p6 = pSequelaeInpatient.eval (ageYears);
    // p7 is P(seq) when parasites aren't cleared
    // note: if p7 were not equal p6, we might have to do something similar to p5b
    double p7 = p6;

    double q[9];
    // Community deaths
    q[0] = (1 - p2) * p5a;
    // Community sequelae
    q[1] = q[0] + (1 - p2) * (1 - p5a) * p7;
    // Community survival
    q[2] = q[1] + (1 - p2) * (1 - p5a) * (1 - p7);
    // In-hospital parasitological failure deaths
    q[3] = q[2] + p2 * (1 - p3) * p5b;
    // In-hospital parasitological failure sequelae
    q[4] = q[3] + p2 * (1 - p3) * (1 - p5b) * p7;
    // In-hospital parasitological failure survivors
    q[5] = q[4] + p2 * (1 - p3) * (1 - p5b) * (1 - p7);
    // In-hospital parasitological success deaths
    q[6] = q[5] + p2 * p3 * p4;
    // In-hospital parasitological success sequelae
    q[7] = q[6] + p2 * p3 * (1 - p4) * p6;
    // In-hospital parasitological success survival
    q[8] = q[7] + p2 * p3 * (1 - p4) * (1 - p6);
    
    // Expectation of death is:
    // double exDeath = (q[6] - q[5]) + (q[3] - q[2]) + q[0];
    // expanded, separating in-hospital and without (q[0]):
    const double exHospitalDeath = p2 * (p3 * p4 + (1 - p3) * p5b);
    const double exDeath = exHospitalDeath + (1 - p2) * p5a;
    mon::reportStatMHF( mon::MHF_EXPECTED_HOSPITAL_DEATHS, human, exHospitalDeath );
    mon::reportStatMHF( mon::MHF_EXPECTED_DIRECT_DEATHS, human, exDeath );
    
    // Expectation of sequelae is:
    // double exSeq = (q[7] - q[6]) + (q[4] - q[3]) + (q[1] - q[0]);
    // expanded and simplified, noting that p7 == p6 :
    const double exSeq = (p2 * (p3 * (1 - p4) + (1 - p3) * (1 - p5b)) + (1 - p2) * (1 - p5a)) * p6;
    mon::reportStatMHF( mon::MHF_EXPECTED_SEQUELAE, human, exSeq );

    double prandom = human.rng.uniform_01();
    
    //NOTE: we do not model diagnostics in this case
    if( prandom >= q[2] ){      // treated in hospital
        m_tLastTreatment = sim::ts0();
        mon::reportEventMHI( mon::MHT_TREATMENTS_3, human, 1 );
        Episode::State stateTreated = Episode::State (pgState | Episode::EVENT_IN_HOSPITAL);
        
        if( prandom >= q[5] ){  // treatment successful at clearing parasites
            // this actually means successful treatment in this case:
            human.withinHostModel->treatment( human, treatmentSevere );
            
            if( prandom < q[6] ){       // death (despite success in clearing parasites)
                latestReport.update (human, Episode::State (stateTreated | Episode::DIRECT_DEATH));
                doomed  = DOOMED_COMPLICATED;
            }else if( prandom < q[7] ){ // Patient recovers, but with sequelae (don't report full recovery)
                latestReport.update (human, Episode::State (stateTreated | Episode::SEQUELAE));
            }else /*prandom < q[8]*/{   // patient recovers fully
                latestReport.update (human, Episode::State (stateTreated | Episode::RECOVERY));
            }
        }else /*prandom < q[5]*/{  // Treated but parasites not cleared (in hospital)
            // No change in parasitological status among in-hospital patients
            
            if( prandom < q[3] ){       // death
                latestReport.update (human, Episode::State (stateTreated | Episode::DIRECT_DEATH));
                doomed  = DOOMED_COMPLICATED;
            }else if( prandom < q[4] ){ // sequelae without parasite clearance
                latestReport.update (human, Episode::State (stateTreated | Episode::SEQUELAE));
            }else /*prandom < q[5]*/{   // full recovery from episode
                latestReport.update (human, pgState);
            }
        }
    }else /*prandom < q[2]*/{
        // Not treated, thus no change in parasitological status
        
        if( prandom < q[0] ){   // outcome 0: death in community
            latestReport.update (human, Episode::State (pgState | Episode::DIRECT_DEATH));
            doomed  = DOOMED_COMPLICATED;
        }else if( prandom < q[1] ){     // outcome 1: sequelae in community
            latestReport.update (human, Episode::State (pgState | Episode::SEQUELAE));
        }else /*prandom < q[2]*/{       // outcome 2: full recovery from episode in community
            latestReport.update (human, pgState);
        }
    }
}


// ———  per-human, intervention & checkpointing  ———

void CM5DayCommon::checkpoint (istream& stream) {
    ClinicalModel::checkpoint (stream);
    m_tLastTreatment & stream;
    m_treatmentSeekingFactor & stream;
}
void CM5DayCommon::checkpoint (ostream& stream) {
    ClinicalModel::checkpoint (stream);
    m_tLastTreatment & stream;
    m_treatmentSeekingFactor & stream;
}

}
}
