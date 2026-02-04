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

#ifndef Hmod_interventions_HumanComponents
#define Hmod_interventions_HumanComponents

#include "Global.h"
#include "interventions/Interfaces.h"
#include "util/sampler.h"
#include "util/DecayFunction.h"

namespace OM {
namespace Host {
    class Human;
}
namespace interventions {
    class VaccineComponent;
    using util::LocalRng;

// ———  sub-population removal  ———

namespace SubPopRemove{
    /** Identifiers for types of sub-population removal option, followed by NUM
     * in last place (as a counter).
     * 
     * ON_FIRST_BOUT: remove humans from some sub-populations at the start of
     * their first clinical events following inclusion in the sub-population.
     * 
     * ON_FIRST_INFECTION: remove patent humans from some sub-population(s)
     * during each survey.
     * 
     * ON_FIRST_TREATMENT: remove humans from some sub-populations when
     * assigning each course of treatment.
     */
    enum RemoveAtCode {
        ON_FIRST_BOUT,
        ON_FIRST_INFECTION,
        ON_FIRST_TREATMENT,
        NUM };
}

/** For each RemoveAtCode (excluding NUM), this is a list of all sub-population
 * identifiers for which the option is enabled. */
extern vector<ComponentId> removeAtIds[SubPopRemove::NUM];

// ———  vaccines  ———

namespace Vaccine{
    enum Types { PEV, BSV, TBV, NumVaccineTypes };
}

/** Per vaccine effect (type), per human details. */
class PerEffectPerHumanVaccine {
public:
    //Note: this constructor is only for checkpointing
    PerEffectPerHumanVaccine();
    
    /// Checkpointing
    template<class S>
    void operator& (S& stream) {
        component & stream;
        numDosesAdministered & stream;
        timeLastDeployment & stream;
        perGenotypeInitialEfficacy & stream;
        hetSample & stream;
    }
    
private:
    PerEffectPerHumanVaccine( LocalRng& rng, ComponentId id, const VaccineComponent& params );
    
    ComponentId component;      // id of component (for finding parameters)
    // Type of vaccine
    /** Number of vaccine doses this individual has received.
     *
     * If an individual misses one EPI (continuous) vaccine dose, it's
     * intentional that they also miss following EPI doses (unless a timed mass
     * vaccination reintroduces them to the EPI schedule). */
    uint32_t numDosesAdministered;
    /// Time of last vaccination with this vaccine type
    SimTime timeLastDeployment = sim::never();
    /// Efficacy at last deployment (undecayed)
    vector<double> perGenotypeInitialEfficacy;
    unique_ptr<util::DecayFunction> hetSample;
    
    friend class PerHumanVaccine;
};

/** Per-human vaccine code. */
class PerHumanVaccine {
public:
    PerHumanVaccine() {}
    
    /** Get one minus the efficacy of the vaccine (1 for no effect, 0 for full effect). */
    double getFactor( Vaccine::Types type, uint32_t genotype = 0) const;
    
    /** Vaccinate unless the passed VaccineLimits specify not to.
     * 
     * @returns true when the vaccine is administered */
    bool possiblyVaccinate( Host::Human& human, ComponentId componentId,
                           interventions::VaccineLimits vaccLimits );

#if 0
    /// Hack for R_0 experiment: make current human the infection source
    inline void specialR_0(){
        // At this point special TBV should have been given to everyone but the PEV to no-one
        effects[R_0_PEV].initialEfficacy = 1.0;
        effects[R_0_TBV].initialEfficacy = 0.0;
    }
#endif
    
    /// Checkpointing
    template<class S>
    void operator& (S& stream) {
        effects & stream;
    }

private:
    /// Details for each deployed vaccine for this human
    typedef std::vector<PerEffectPerHumanVaccine> EffectList;
    EffectList effects;
};

}
}
#endif
