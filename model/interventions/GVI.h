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

#ifndef OM_INTERVENTIONS_GVI
#define OM_INTERVENTIONS_GVI

#include "Transmission/PerHost.h"
#include "interventions/Interfaces.h"
#include "util/DecayFunction.h"
#include "util/sampler.h"
#include "schema/interventions.h"

namespace OM { namespace interventions {
    using util::DecayFunction;
    using util::LognormalSampler;
    using util::LocalRng;
    using Transmission::PerHostInterventionData;

/** Constant parameters for generic vector intervention model. */
class GVIComponent : public Transmission::HumanVectorInterventionComponent {
public:
    /** Initialise parameters
     * 
     * @param elt Component description from XML
     * @param species_name_map Map of species names to indices.
     */
    GVIComponent( ComponentId, const scnXml::GVIDescription& elt,
               const map<string,size_t>& species_name_map );
    
    void deploy( Host::Human& human, mon::Deploy::Method method, VaccineLimits )const;
    
    virtual Component::Type componentType() const;
    
    virtual void print_details( std::ostream& out )const;
    
    virtual unique_ptr<PerHostInterventionData> makeHumanPart(LocalRng& rng) const;
    virtual unique_ptr<PerHostInterventionData> makeHumanPart( istream& stream, ComponentId ) const;
    
private:
    /** Per mosquito-species parameters for generic vector intervention model. */
    class GVIAnopheles {
    public:
        GVIAnopheles() :
            proportionProtected( numeric_limits<double>::signaling_NaN() ),
            proportionUnprotected( numeric_limits<double>::signaling_NaN() ),
            deterrency( numeric_limits<double>::signaling_NaN() ),
            preprandialKilling( numeric_limits<double>::signaling_NaN() ),
            postprandialKilling( numeric_limits<double>::signaling_NaN() ),
            fecundityReduction( numeric_limits<double>::signaling_NaN() )
        {}
        void init(const scnXml::GVIDescription::AnophelesParamsType& elt, double proportionUse);
        
        /// Return x*proportionProtected + proportionUnprotected
        inline double byProtection(double x) const{
            return x*proportionProtected + proportionUnprotected;
        }
        
        double proportionProtected;
        double proportionUnprotected;
        double deterrency;
        double preprandialKilling;
        double postprandialKilling;
        double fecundityReduction;
        
        friend class HumanGVI;
    };
    
    unique_ptr<DecayFunction> decay;
    vector<GVIAnopheles> species;  // vector specific params
    
    // This is sparse vector: only indexes corresponding to a GVI component are used
    // No memory management
    static vector<GVIComponent*> componentsByIndex;
    
    friend class HumanGVI;
};

/** Low-level (generic) vector intervention model. Has three effects:
 * deterrency, pre-prandial killing and post-prandial killing.
 * 
 * This is the per-host (but not per vector) part.
 */
class HumanGVI : public PerHostInterventionData {
public:
    HumanGVI( LocalRng& rng, const GVIComponent& params );
    HumanGVI( istream& stream, ComponentId );
    
    virtual void redeploy( LocalRng& rng, const Transmission::HumanVectorInterventionComponent& );
    
    /** This is the survival factor of the effect. */
    inline double getEffectSurvival(const GVIComponent& params)const{
        SimTime age = sim::nowOrTs1() - deployTime;  // implies age 1 TS on first use
        return decayHet->eval( age );
    }
    
    virtual void update(Host::Human& human);
    
    /// Get deterrency. See ComponentParams::effect for a more detailed description.
    virtual double relativeAttractiveness(size_t speciesIndex) const;
    /// Get killing effect on mosquitoes before they've eaten.
    /// See ComponentParams::effect for a more detailed description.
    virtual double preprandialSurvivalFactor(size_t speciesIndex) const;
    /// Get killing effect on mosquitoes after they've eaten.
    /// See ComponentParams::effect for a more detailed description.
    virtual double postprandialSurvivalFactor(size_t speciesIndex) const;
    /// Get the mosquito fecundity multiplier (1 for no effect).
    virtual double relFecundity(size_t speciesIndex) const;
    
protected:
    virtual void checkpoint( ostream& stream );
    
private:
    // this parameter is sampled on first deployment, but never resampled for the same human:
    unique_ptr<DecayFunction> decayHet;
};

} }

#endif
