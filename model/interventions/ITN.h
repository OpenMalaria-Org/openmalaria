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

#ifndef OM_INTERVENTIONS_ITN
#define OM_INTERVENTIONS_ITN

#include "util/DecayFunction.h"
#include "Transmission/PerHost.h"
#include "util/sampler.h"
#include "schema/interventions.h"

namespace OM {
namespace interventions {
    using util::DecayFunction;
    using util::NormalSampler;
    using util::LognormalSampler;
    using util::LocalRng;
    using Transmission::PerHostInterventionData;


// —————  utility classes (internal use only)  —————

namespace factors {
    class SurvivalFactor {
    public:
        /// Set parameters.
        /// 
        /// It is checked that parameters lie in a suitible range, giving a
        /// survival factor between 0 and 1.
        /// 
        /// @param elt  Element from XML
        /// @param maxInsecticide
        /// @param eltName  Element name (used in error messages)
        /// @param isDeterrent True when used to describe a deterrency effect, false when used to
        ///     describe a killing effect.
        void init(const scnXml::ITNKillingEffect& elt, double maxInsecticide,
                    const char* eltName, bool attractivenessConstraints);
        
        /// Variant to set paramters for the logit model.
        /// 
        /// @param elt  Element from XML
        /// @param maxHoleIndex Parameter from parent element in XML
        /// @param isDeterrent True when used to describe a deterrency effect, false when used to
        ///     describe a killing effect.
        void initLogit(const scnXml::ITNEffectLogit& elt,
                                                   double maxHoleIndex,
                                                   bool isDeterrent);
        
        /// Initialise the model to always return no effect (factor 1).
        void init1();
        
        /** Part of survival factor, used by new ITN deterrency model. */
        double rel_pAtt( double holeIndex, double insecticideContent )const;
        /** Calculate additional survival factor imposed by nets on pre-/post-
        * prandial killing. Should be bounded to [0,1] and tend to 1 as the
        * net ages. */
        double survivalFactor( double holeIndex, double insecticideContent )const;
        
    private:
        union {
            struct {
                double BF, HF, PF, IF;  // base, hole, insecticide and interaction factors
                double holeScaling, insecticideScaling;
                double invBaseSurvival; // stored for performance only; ≥1
            } a;
            struct {
                double BF, HF, PF, IF;
                double hMax;
                double invBaseEffect;   // 1 over base effect
            } b;
        };
        bool useLogitEqns;
    };

    class RelativeAttractiveness {
        enum DeterrencyModel {
            NO_MODEL, SINGLE_STAGE, TWO_STAGE, TWO_STAGE_LOGIT
        };
    public:
        RelativeAttractiveness() : model(0) {}
        
        /** Set parameters for single-stage deterrency model.
         * 
         * It is checked that input parameters lie in a range such that
         * the relative availability is always in the range (0,1] — that is,
         * the deterrent can never be perfect, but can have zero effect. */
        void initSingleStage(const scnXml::ITNDeterrency& elt, double maxInsecticide);
        
        /** Set parameters for two-stage deterrency model.
         * 
         * It is checked that input parameters lie in a range such that
         * the relative availability is always in the range (0,1] — that is,
         * the deterrent can never be perfect, but can have zero effect. */
        void initTwoStage(const scnXml::TwoStageDeterrency& elt,
                          double maxInsecticide,
                          double holeIndexMax, bool holeIndexMaxPresent = false);
        
        /** Calculate effect. Range of output is any value ≥ 0.
         * 
         * 0 implies a fully effective deterrent, 0.5 a 50% effective
         * deterrent, 1 has no effect, >1 attracts extra mosquitoes. */
        double relativeAttractiveness (double holeIndex, double insecticideContent) const;
        
    private:
        union {
            struct {
                double lHF, lPF, lIF;      // logs of hole, insecticide and interaction factors
                double holeScaling, insecticideScaling;
            } a;
            struct {
                double lPFEntering;      // log of insecticide factor
                double insecticideScalingEntering;
                double _pad;
                SurvivalFactor pAttacking;
            } b;
            struct {
                double entBaseFactor;
                double entInsecticideFactor;
                double pEnt0Inv;
                SurvivalFactor pAttacking;
            } c;
        };
        int model;
    };
}


// —————  main, public classes  —————

class ITNComponent : public Transmission::HumanVectorInterventionComponent {
public:
    ITNComponent( ComponentId id, const scnXml::ITNDescription& elt,
               const map< string, size_t >& species_name_map );
    
    virtual void deploy( Host::Human& human, mon::Deploy::Method method, VaccineLimits )const;
    
    virtual Component::Type componentType() const;
    
    virtual void print_details( std::ostream& out )const;
    
    virtual unique_ptr<PerHostInterventionData> makeHumanPart(LocalRng& rng) const;
    virtual unique_ptr<PerHostInterventionData> makeHumanPart( istream& stream, ComponentId id ) const;
    
private:
    /** Per mosquito-species parameters for extended ITN model. */
    class ITNAnopheles {
    public:
        ITNAnopheles() :
            proportionProtected( numeric_limits<double>::signaling_NaN() ),
            proportionUnprotected( numeric_limits<double>::signaling_NaN() )
        {}
        void init(const scnXml::ITNDescription::AnophelesParamsType& elt,
                  double proportionUse, double maxInsecticide);
        
        /// Get deterrency. See ComponentParams::effect for a more detailed description.
        /// Range: ≥0 where 0=fullly deter, 1=no effect, >1 = attract
        inline double relativeAttractiveness( double holeIndex, double insecticideContent )const{
            return byProtection( relAttractiveness.relativeAttractiveness( holeIndex, insecticideContent ) );
        }
        /// Get killing effect on mosquitoes before feeding.
        /// See ComponentParams::effect for a more detailed description.
        inline double preprandialSurvivalFactor( double holeIndex, double insecticideContent )const{
            return byProtection( preprandialKillingEffect.survivalFactor( holeIndex, insecticideContent ) );
        }
        /// Get killing effect on mosquitoes after they've eaten.
        /// See ComponentParams::effect for a more detailed description.
        inline double postprandialSurvivalFactor( double holeIndex, double insecticideContent )const{
            return byProtection( postprandialKillingEffect.survivalFactor( holeIndex, insecticideContent ) );
        }
        /// Get mosquito fecundity multiplier
        inline double relFecundity( double holeIndex, double insecticideContent )const{
            return byProtection( relFecundityEffect.survivalFactor( holeIndex, insecticideContent ) );
        }
        
        /// Return x*proportionProtected + proportionUnprotected
        inline double byProtection(double x) const{
            return x*proportionProtected + proportionUnprotected;
        }
        
    private:
        double proportionProtected;
        double proportionUnprotected;
        factors::RelativeAttractiveness relAttractiveness;
        factors::SurvivalFactor preprandialKillingEffect;
        factors::SurvivalFactor postprandialKillingEffect;
        factors::SurvivalFactor relFecundityEffect;
        
        friend class HumanITN;
    };

    NormalSampler initialInsecticide;
    unique_ptr<LognormalSampler> holeRate;	// holes per step
    unique_ptr<LognormalSampler> ripRate;	// rips per hole per step
    double maxInsecticide;		// maximum initial insecticide
    double ripFactor;			// factor expressing how significant rips are in comparison to holes
    unique_ptr<DecayFunction> insecticideDecay;
    unique_ptr<DecayFunction> attritionOfNets;
    vector<ITNAnopheles> species; // vector specific params
    
    // This is sparse vector: only indexes corresponding to ITN components are
    // used. No memory management.
    static vector<ITNComponent*> componentsByIndex;
    
    friend class HumanITN;
};

/** Extended ITN model by OB.
 * 
 * Each instance describes a hypothetical net (or no net).
 */
class HumanITN : public PerHostInterventionData {
public:
    HumanITN( LocalRng& rng, const ITNComponent& params );
    HumanITN( istream& stream, ComponentId id );
    
    virtual void redeploy(LocalRng& rng, const Transmission::HumanVectorInterventionComponent& params);
    
    inline double getHoleIndex()const{
        return holeIndex;
    }
    inline double getInsecticideContent(const ITNComponent& params)const{
            SimTime age = sim::nowOrTs1() - deployTime;  // implies age 1 TS on first use
        double effectSurvival = insecticideDecayHet->eval( age );
        return initialInsecticide * effectSurvival;
    }
    
    /// Call once per time step to update holes
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
    // these parameters express the current state of the net:
    SimTime disposalTime = sim::never();	// time at which net will be disposed of (if it's not already been replaced)
    int nHoles;				// total number of holes
    double holeIndex;		// a measure of both the number and size of holes
    double initialInsecticide;	// units: mg/m²
    
    // these parameters are sampled from log-normal per net, but thereafter constant:
    double holeRate;	// rate at which new holes are created (holes/time-step)
    double ripRate;		// rate at which holes are enlarged (rips/hole/time-step)
    unique_ptr<DecayFunction> insecticideDecayHet;

};

} }

#endif
