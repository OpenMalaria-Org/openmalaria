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

#ifndef Hmod_human
#define Hmod_human
#include "Global.h"
#include "Transmission/PerHost.h"
#include "InfectionIncidenceModel.h"
#include "mon/AgeGroup.h"
#include "interventions/HumanComponents.h"
#include "util/checkpoint_containers.h"
#include <map>

namespace OM {
namespace Clinical {
    class ClinicalModel;
}
namespace WithinHost {
    class WHInterface;
}
namespace Transmission {
    class TransmissionModel;
}

namespace Host {

/** Interface to all sub-models storing data per-human individual.
 *
 * Still contains some data, but most is now contained in sub-models. */
class Human {
public:
    /** Initialise all variables of a human datatype.
    * 
    * @param dateOfBirth date of birth (usually start of next time step) */
    Human(SimTime dateOfBirth);

    /// Allow move construction
    Human(Human&&) = default;
    Human& operator=(Human&&) = default;

    /// Disable copying
    Human(const Human&) = delete;
    Human& operator=(const Human&) = delete;

    /** Get human's age with respect to some time. */
    inline SimTime age( SimTime time ) const;

    /** Return true if human is a member of the sub-population. */
    inline bool isInSubPop( interventions::ComponentId id ) const;

    inline uint32_t getCohortSet() const;

    /** Return date of birth */
    SimTime getDOB() const;

    /** Return true if dead, false otherwise */
    bool isDead() const;

    /** After this call, isDead() will be true */
    void kill();

    /** Add the human to an intervention component's sub-population for the given duration. 
     * A duration of zero implies no effect. */
    void addToCohort(interventions::ComponentId id, SimTime duration);

    /** Remove the human from an intervention component's sub-population */
    void removeFromCohort(interventions::ComponentId id);

    /** Act on "remove from sub-population on first ..." events.
     * This is only for use during a human update. */
    void removeFirstEvent(interventions::SubPopRemove::RemoveAtCode code);

    /** Remove host from expired cohorts */
    void updateCohortSet();

    /** Checkpoint (read) */
    void checkpoint(istream &stream);

    /** Checkpoint (write) */
    void checkpoint(ostream &stream);

    /** Contains per-species vector data (VectorModel only). */
    Transmission::PerHost perHostTransmission;

    /** The WithinHostModel models parasite density and immunity */
    unique_ptr<WithinHost::WHInterface> withinHostModel;

    /** The InfectionIncidenceModel translates per-host EIR into new infections */
    unique_ptr<InfectionIncidenceModel> infIncidence;

    /** The ClinicalModel encapsulates pathogenesis (sickness status),
    * case management (medicating drugs)
    * and clinical outcomes (morbidity, reporting). */
    unique_ptr<Clinical::ClinicalModel> clinicalModel;

    util::LocalRng rng;

    /** Vaccines */
    interventions::PerHumanVaccine vaccine;

    /** Made persistant to save a lookup each time step (significant performance improvement) */
    mon::AgeGroup monitoringAgeGroup;

    /** The next continuous distribution in the series */
    uint32_t nextCtsDist = 0;

private:
    SimTime dateOfBirth = sim::never();        // date of birth; humans are always born at the end of a time step

    /** Cache, updated when human is added to or removed from a sub-population */
    uint32_t cohortSet = 0;

    /** The state of he human. A human cannot be revived. */
    bool dead = false;

    /** This lists sub-populations of which the human is a member together with
    * expiry time.
    * 
    * Definition: a human is in sub-population p if subPopExp.contains(p) and,
    * for t=subPopExp[p], t > sim::now() (at the time of intervention
    * deployment) or t > sim::ts0() (equiv t >= sim::ts1()) during human update.
    * 
    * NOTE: this discrepancy is because intervention deployment effectively
    * happens at the end of a time step and we want a duration of 1 time step to
    * mean 1 intervention deployment (that where the human becomes a member) and
    * 1 human update (the next). */
    //TODO(optimisation): it might be better to instead store for each
    // ComponentId of interest the set of humans who are members
    std::map<interventions::ComponentId,SimTime> subPopExp;
};

void summarize(Human &human, bool surveyOnlyNewEp);

void update(Human &human, Transmission::TransmissionModel& transmission);

inline SimTime Human::age( SimTime time ) const { 
    return time - dateOfBirth; 
}

inline bool Human::isInSubPop( interventions::ComponentId id ) const {
    auto it = subPopExp.find( id );
    if( it == subPopExp.end() ) return false;   // no history of membership
    else return it->second >= sim::nowOrTs0();   // added: has expired?
}

inline uint32_t Human::getCohortSet() const { 
    return cohortSet; 
}

inline SimTime Human::getDOB() const {
    return dateOfBirth;
}

inline bool Human::isDead() const {
    return dead;
}

inline void Human::kill() {
    dead = true;
}

} }
#endif
