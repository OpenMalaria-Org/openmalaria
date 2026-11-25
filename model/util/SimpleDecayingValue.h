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

#ifndef Hmod_SimpleDecayingValue
#define Hmod_SimpleDecayingValue

#include "Global.h"
#include "util/DecayFunction.h"

namespace OM {
namespace util {

/** Parameters to for a decaying intervention whose effect can be described
 * by a single scalar output.
 * 
 * This class outputs a scalar which defaults to and decays to 0, but
 * immediately after deployment will be the set _initial value_.
 *****************************************************************************/
class SimpleDecayingValue
{
public:
    /** Default construction: always return 0. */
    SimpleDecayingValue() :
            initial(0.0),
            deploy_t(sim::never()) {}
    
    /** Configure from an XML element. */
    inline void set (double initial_value, const scnXml::DecayFunction& elt, const char* name){
        decay = DecayFunction::makeObject(elt, name);
        initial = initial_value;
    }

    /** Configure from an XML element. */
    inline void set (double initial_value){
        initial = initial_value;
    }
    
    /** Trigger a deployment, if decay function was set. */
    inline void deploy (LocalRng& rng, SimTime time){
        if( decay.get() == 0 ) return;  // cannot deploy
        deploy_t = time;
        het = decay->hetSample(rng);
    }
    
    /** Get the value (0 if before any deployment or after complete decay,
     * also 0 if no decay function or initial value was set,
     * otherwise between zero and the inital value). */
    inline double current_value (SimTime time) const{
        if( decay.get() == 0 || het == nullptr) return 0.0;  // decay wasn't set: the always return 0
        else return initial * het->eval( time - deploy_t );
    }
    
    /** Checkpointing: only checkpoint parameters which change after initial
     * set-up. */
    template<class S>
    void operator& (S& stream) {
        if (decay) decay->checkpointSample(het, stream);
        deploy_t & stream;
    }
    
private:
    /** Description of decay of effects on emergence. */
    unique_ptr<util::DecayFunction> decay;
    
    /** Initial value. Is initialised to 0. */
    double initial;
    
    /** Description of larviciding decay. */
    unique_ptr<util::DecayFunction> het;
    
    /** Time of larviciding deployment. */
    SimTime deploy_t = sim::never();
};

/** A time series of multipliers: one factor per time step, relative to a
 * deployment time (used by the changeEmergence intervention).
 * 
 * Unlike SimpleDecayingValue, the output is a factor rather than an effect
 * size: 1 means no change. The factor in effect is the value of the series at
 * the index given by the number of time steps since deployment (the series
 * becomes active on the time step in which the deployment happens). Before
 * the first deployment and once the series is exhausted the factor is 1 (no
 * change). Re-deployment restarts the series. No heterogeneity factor applies
 * and no random numbers are used.
 *****************************************************************************/
class TimestepFactors
{
public:
    /** Default construction: no series, always return 1. */
    TimestepFactors() : deploy_t(sim::never()) {}
    
    /** Set the series: one factor per time step after deployment.
     * epsilon is a global multiplier applied to every value (the factor in
     * effect is epsilon times values[n]; outside the window the factor stays
     * 1). The factors should have been validated (non-empty and finite;
     * negative factors are accepted for stress testing only), as should
     * epsilon (finite). */
    void set (const std::vector<double>& factors, double scale = 1.0) {
        values = factors;
        epsilon = scale;
    }
    
    /** Trigger a deployment: the series restarts from this time. */
    void deploy (SimTime time) {
        deploy_t = time;
    }
    
    /** The factor to apply at the given time: epsilon * values[n], where n is
     * the number of time steps since deployment; 1 before any deployment or
     * once the series is exhausted. */
    double current_factor (SimTime time) const {
        if (deploy_t == sim::never()) return 1.0;
        const int n = sim::inSteps(time - deploy_t);
        if (n < 0 || n >= static_cast<int>(values.size())) return 1.0;
        return epsilon * values[n];
    }
    
    /** Checkpointing: only the deployment time changes after set-up; the
     * series itself is reloaded from the scenario file. */
    template<class S>
    void operator& (S& stream) {
        deploy_t & stream;
    }
    
private:
    std::vector<double> values;
    double epsilon = 1.0; ///< global multiplier for the values
    SimTime deploy_t;
};

} }
#endif
