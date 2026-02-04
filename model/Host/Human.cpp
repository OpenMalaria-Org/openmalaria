/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2026 Swiss Tropical and Public Health Institute
 * Copyright (C) 2005-2015 Liverpool School Of Tropical Medicine
 * Copyright (C) 2020-2026 University of Basel
 * Copyright (C) 2026 The Kids Research Institute Australia
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

#include "Host/Human.h"
#include "Host/InfectionIncidenceModel.h"
#include "Clinical/ClinicalModel.h"
#include "Host/WithinHost/WHInterface.h"

#include "Transmission/TransmissionModel.h"
#include "util/ModelOptions.h"
#include "util/vectors.h"
#include "util/StreamValidator.h"
#include "Population.h"
#include "interventions/InterventionManager.h"
#include "mon/reporting.h"
#include "schema/scenario.h"

namespace OM { namespace Host {
    using namespace OM::util;
    using interventions::ComponentId;

struct HumanHet {
    /* Human heterogeneity; affects:
     * comorbidityFactor (stored in PathogenesisModel)
     * treatmentSeekingFactor (stored in CaseManagementModel)
     * availabilityFactor (stored in Transmission::PerHost) */
    double comorbidityFactor;
    double treatmentSeekingFactor;
    double availabilityFactor;

    HumanHet() : comorbidityFactor(1.0), treatmentSeekingFactor(1.0), availabilityFactor(1.0) {}
};

HumanHet hetSample(util::LocalRng& rng)
{
    HumanHet het;
    if(util::ModelOptions::option (util::TRANS_HET)){
        het.availabilityFactor = 0.2;
        if( rng.bernoulli(0.5) ){
            het.availabilityFactor = 1.8;
        }
    }
    if(util::ModelOptions::option (util::COMORB_HET)){
        het.comorbidityFactor = 0.2;
        if( rng.bernoulli(0.5) ){
            het.comorbidityFactor = 1.8;
        }
    }
    if(util::ModelOptions::option (util::TREAT_HET)){
        het.treatmentSeekingFactor = 0.2;
        if( rng.bernoulli(0.5) ){
            het.treatmentSeekingFactor = 1.8;
        }
    }
    if(util::ModelOptions::option (util::TRANS_TREAT_HET)){
        het.treatmentSeekingFactor = 0.2;
        het.availabilityFactor = 1.8;
        if( rng.bernoulli(0.5) ){
            het.treatmentSeekingFactor = 1.8;
            het.availabilityFactor = 0.2;
        }
    }else if(util::ModelOptions::option (util::COMORB_TREAT_HET)){
        if( rng.bernoulli(0.5) ){
            het.comorbidityFactor = 1.8;
            het.treatmentSeekingFactor = 0.2;
        }else{
            het.comorbidityFactor = 0.2;
            het.treatmentSeekingFactor = 1.8;
        }
    }else if(util::ModelOptions::option (util::COMORB_TRANS_HET)){
        het.availabilityFactor = 1.8;
        het.comorbidityFactor = 1.8;
        if( rng.bernoulli(0.5) ){
            het.availabilityFactor = 0.2;
            het.comorbidityFactor = 0.2;
        }
    }else if(util::ModelOptions::option (util::TRIPLE_HET)){
        het.availabilityFactor = 1.8;
        het.comorbidityFactor = 1.8;
        het.treatmentSeekingFactor = 0.2;
        if( rng.bernoulli(0.5) ){
            het.availabilityFactor = 0.2;
            het.comorbidityFactor = 0.2;
            het.treatmentSeekingFactor = 1.8;
        }
    }
    return het;
}

Human::Human(SimTime dateOfBirth) :
    infIncidence(InfectionIncidenceModel::createModel()),
    rng(util::master_RNG),
    dateOfBirth(dateOfBirth)
{
    // Initial humans are created at time 0 and may have dateOfBirth in past. Otherwise dateOfBirth must be now.
    assert( dateOfBirth == sim::nowOrTs1() || (sim::now() == sim::zero() && dateOfBirth < sim::now()) );

    HumanHet het = hetSample(rng);
    withinHostModel = WithinHost::WHInterface::createWithinHostModel( rng, het.comorbidityFactor );
    double iiFactor = infIncidence->getAvailabilityFactor(rng, 1.0);
    perHostTransmission.initialise (rng, het.availabilityFactor * iiFactor);
    clinicalModel = Clinical::ClinicalModel::createClinicalModel(het.treatmentSeekingFactor);
}

void Human::addToCohort(ComponentId id, SimTime duration )
{
    if( duration <= sim::zero() ) return; // nothing to do
    subPopExp[id] = sim::nowOrTs1() + duration;
    cohortSet = mon::updateCohortSet( cohortSet, id, true );
}

void Human::removeFromCohort(interventions::ComponentId id)
{
    subPopExp.erase(id);
}

void Human::removeFirstEvent(interventions::SubPopRemove::RemoveAtCode code )
{
    const vector<ComponentId>& removeAtList = interventions::removeAtIds[code];
    for( auto it = removeAtList.begin(), end = removeAtList.end(); it != end; ++it ){
        auto expIt = subPopExp.find( *it );
        if( expIt != subPopExp.end() ){
            if( expIt->second >= sim::nowOrTs0() ){
                // removeFirstEvent() is used for onFirstBout, onFirstTreatment
                // and onFirstInfection cohort options. Health system memory must
                // be reset for this to work properly; in theory the memory should
                // be independent for each cohort, but this is a usable approximation.
                clinicalModel->flushReports();     // reset HS memory
                
                // report removal due to first infection/bout/treatment
                mon::reportEventMHI( mon::MHR_SUB_POP_REM_FIRST_EVENT, *this, 1 );
            }
            cohortSet = mon::updateCohortSet( cohortSet, expIt->first, false );
            // remove (affects reporting, restrictToSubPop and cumulative deployment):
            subPopExp.erase( expIt );
        }
    }
}

void Human::updateCohortSet()
{
    // check sub-pop expiry
    for( auto expIt = subPopExp.begin(), expEnd = subPopExp.end(); expIt != expEnd; ) {
        if( !(expIt->second >= sim::ts0()) ){       // membership expired
            // don't flush reports
            // report removal due to expiry
            mon::reportEventMHI( mon::MHR_SUB_POP_REM_TOO_OLD, *this, 1 );
            cohortSet = mon::updateCohortSet( cohortSet, expIt->first, false );
            // erase element, but continue iteration
            expIt = subPopExp.erase( expIt );
        }else{
            ++expIt;
        }
    }
}

void Human::checkpoint(istream &stream)
{
    perHostTransmission & stream;
    infIncidence & stream;
    withinHostModel & stream;
    clinicalModel & stream;
    rng.checkpoint(stream);
    dateOfBirth & stream;
    vaccine & stream;
    monitoringAgeGroup & stream;
    cohortSet & stream;
    nextCtsDist & stream;
    subPopExp & stream;
}

void Human::checkpoint(ostream &stream)
{
    perHostTransmission & stream;
    infIncidence & stream;
    withinHostModel & stream;
    clinicalModel & stream;
    rng.checkpoint(stream);
    dateOfBirth & stream;
    vaccine & stream;
    monitoringAgeGroup & stream;
    cohortSet & stream;
    nextCtsDist & stream;
    subPopExp & stream;
}

// -----  Non-static functions: per-time-step update  -----
void summarize(Human &human, bool surveyOnlyNewEp) {
    if( surveyOnlyNewEp && human.clinicalModel->isExistingCase() ){
        // This modifies the denominator to treat the health-system-memory
        // period immediately after a bout as 'not at risk'.
        return;
    }
    
    mon::reportStatMHI( mon::MHR_HOSTS, human, 1 );
    mon::reportStatMHF( mon::MHF_AGE, human, sim::inYears(human.age(sim::now())) );
    bool patent = human.withinHostModel->summarize (human);
    human.infIncidence->summarize (human);
    
    if( patent && mon::isReported() ){
        // this should happen after all other reporting!
        human.removeFirstEvent(interventions::SubPopRemove::ON_FIRST_INFECTION);
    }
}

void update(Human &human, Transmission::TransmissionModel& transmission)
{
    // For integer age checks we use age0 to e.g. get 73 steps comparing less than 1 year old
    SimTime age0 = human.age(sim::ts0());
    if (human.clinicalModel->isDead(age0)) {
        human.kill();
        return;
    }
    
    util::streamValidate( age0 );

    // monitoringAgeGroup is the group for the start of the time step.
    human.monitoringAgeGroup.update( age0 );
    
    human.updateCohortSet();

    // Age at  the end of the update period. In most cases
    // the difference between this and age at the start is not especially
    // important in the model design, but since we parameterised with
    // age1 we should stick with it.
    double age1 = sim::inYears(human.age(sim::ts1()));

    // age1 used only in PerHost::relativeAvailabilityAge(); difference to age0 should be minor
    vector<double> EIR_per_genotype_i, EIR_per_genotype_l;
    transmission.getEIR(human, age0, age1, EIR_per_genotype_i, EIR_per_genotype_l);

    double EIR_i = util::vectors::sum(EIR_per_genotype_i);
    double EIR_l = util::vectors::sum(EIR_per_genotype_l);
    double EIR = EIR_i + EIR_l;

    double x = 0.0;
    
    if(EIR > 0)
        x = EIR_i / EIR;

    double expectedNNewInfs = human.infIncidence->expectedNumNewInfections( human, EIR);
    double expectedNNewInfs_i = x * expectedNNewInfs;
    double expectedNNewInfs_l = (1.0-x) * expectedNNewInfs;

    // One more rng call if imported infection
    int nNewInfs_l = human.infIncidence->numNewInfections( human, expectedNNewInfs_l);
    int nNewInfs_i = human.infIncidence->numNewInfections( human, expectedNNewInfs_i);

    // age1 used when medicating drugs (small effect) and in immunity model (which was parameterised for it)
    human.withinHostModel->update(human, human.rng, nNewInfs_i, nNewInfs_l, EIR_per_genotype_i, EIR_per_genotype_l, age1);
    human.infIncidence->reportNumNewInfections(human, nNewInfs_i, nNewInfs_l);

    // age1 used to get case fatality and sequelae probabilities, determine pathogenesis
    human.clinicalModel->update( human, age1, age0 == sim::zero() );
    human.clinicalModel->updateInfantDeaths( age0 );
}

} }
