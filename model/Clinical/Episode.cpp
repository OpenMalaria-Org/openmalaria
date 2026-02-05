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

#include "Clinical/Episode.h"
#include "Clinical/ClinicalModel.h"
#include "Host/Human.h"
#include "Host/WithinHost/WHInterface.h"

namespace OM {
namespace Clinical {

Episode::~Episode ()
{
    report ();
}

void Episode::flush() {
    report();
    time = sim::never();
}


void Episode::update (const Host::Human& human, Episode::State newState)
{
    bool toReport = false;

    if (healthSystemMemoryFix && (time + ClinicalModel::hsMemory() <= sim::ts0()))
        toReport = true;
    else if (!healthSystemMemoryFix && (time + ClinicalModel::hsMemory() < sim::ts0()))
        toReport = true;

    if(toReport) {
        infectionType = human.withinHostModel->getInfectionOrigin();

        report ();

        time = sim::ts0();
        surveyPeriod = mon::eventSurveyNumber();
        ageGroup = human.monitoringAgeGroup;
        cohortSet = human.getCohortSet();
        state = newState;
    } else {
        state = Episode::State (state | newState);
    }
}

void Episode::report () {
    if (time < sim::zero())        // Nothing to report
        return;
    
    // Reports malarial/non-malarial UC fever dependent on cause, not diagnosis.
    if (state & Episode::MALARIA) {
        // Malarial fevers: report bout
        if (state & Episode::COMPLICATED) {
            mon::reportMSACI( mon::MHE_SEVERE_EPISODES, surveyPeriod, ageGroup, cohortSet, 1 );
            if (state & Episode::SEVERE) {
                mon::reportMSACI( mon::MHE_SEVERE_EPISODES_WITHOUT_COMORBIDITIES, surveyPeriod, ageGroup, cohortSet, 1 );
            }
        } else { // UC or UC2
            mon::reportMSACI( mon::MHE_UNCOMPLICATED_EPISODES, surveyPeriod, ageGroup, cohortSet, 1 );
            if(infectionType == WithinHost::InfectionOrigin::Indigenous)
                mon::reportMSACI( mon::MHE_UNCOMPLICATED_EPISODES_INDIGENOUS, surveyPeriod, ageGroup, cohortSet, 1 );
            else if(infectionType == WithinHost::InfectionOrigin::Introduced)
                mon::reportMSACI( mon::MHE_UNCOMPLICATED_EPISODES_INTRODUCED, surveyPeriod, ageGroup, cohortSet, 1 );
            else
                mon::reportMSACI( mon::MHE_UNCOMPLICATED_EPISODES_IMPORTED, surveyPeriod, ageGroup, cohortSet, 1 );
        }

        // Report outcomes of malarial fevers
        if (state & Episode::EVENT_IN_HOSPITAL) {
            if (state & Episode::DIRECT_DEATH) {
                mon::reportMSACI( mon::MHO_DIRECT_DEATHS, surveyPeriod, ageGroup, cohortSet, 1 );
                mon::reportMSACI( mon::MHO_HOSPITAL_DEATHS, surveyPeriod, ageGroup, cohortSet, 1 );
                if (state & Episode::EVENT_FIRST_DAY){
                    mon::reportMSACI( mon::MHO_FIRST_DAY_DEATHS, surveyPeriod, ageGroup, cohortSet, 1 );
                    mon::reportMSACI( mon::MHO_HOSPITAL_FIRST_DAY_DEATHS, surveyPeriod, ageGroup, cohortSet, 1 );
                }
            }
            else if (state & Episode::SEQUELAE) {
                mon::reportMSACI( mon::MHO_SEQUELAE, surveyPeriod, ageGroup, cohortSet, 1 );
                mon::reportMSACI( mon::MHO_HOSPITAL_SEQUELAE, surveyPeriod, ageGroup, cohortSet, 1 );
            }
            else if (state & Episode::RECOVERY){
                mon::reportMSACI( mon::MHO_HOSPITAL_RECOVERIES, surveyPeriod, ageGroup, cohortSet, 1 );
            }
        } else {
            if (state & Episode::DIRECT_DEATH) {
                mon::reportMSACI( mon::MHO_DIRECT_DEATHS, surveyPeriod, ageGroup, cohortSet, 1 );
                if (state & Episode::EVENT_FIRST_DAY){
                    mon::reportMSACI( mon::MHO_FIRST_DAY_DEATHS, surveyPeriod, ageGroup, cohortSet, 1 );
                }
            }
            else if (state & Episode::SEQUELAE){
                mon::reportMSACI( mon::MHO_SEQUELAE, surveyPeriod, ageGroup, cohortSet, 1 );
            }
            // Don't care about out-of-hospital recoveries
        }
    } else if (state & Episode::SICK) {
        // Report non-malarial fever and outcomes
        mon::reportMSACI( mon::MHE_NON_MALARIA_FEVERS, surveyPeriod, ageGroup, cohortSet, 1 );

        if (state & Episode::DIRECT_DEATH) {
            mon::reportMSACI( mon::MHO_NMF_DEATHS, surveyPeriod, ageGroup, cohortSet, 1 );
        }
    }
}

void Episode::operator& (istream& stream) {
    time & stream;
    if (time > sim::zero()) {
        surveyPeriod & stream;
        ageGroup & stream;
        cohortSet & stream;
        int s;
        s & stream;
        state = Episode::State(s);
    }
}
void Episode::operator& (ostream& stream) {
    time & stream;
    if (time >= sim::zero()) {
        surveyPeriod & stream;
        ageGroup & stream;
        cohortSet & stream;
        state & stream;
    }
}

}
}
