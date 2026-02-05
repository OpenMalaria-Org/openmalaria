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

#ifndef Hmod_PerHost
#define Hmod_PerHost

#include "Global.h"
#include "interventions/Interfaces.h"
#include "util/AgeGroupInterpolation.h"
#include "util/DecayFunction.h"
#include "util/checkpoint_containers.h"

#include "util/random.h"

namespace OM {
namespace Transmission {

using util::AgeGroupInterpolator;
using util::DecayFunction;
using util::LocalRng;

class HumanVectorInterventionComponent;

/** Stores vector model data applicable between a category of host and a
 * mosquito species: intervention descriptions and model parameters.
 *
 * Parameters are read from XML, and the availability rate is adjusted. */
class PerHostAnophParams {
public:
    PerHostAnophParams (const scnXml::Mosq& mosq);

    PerHostAnophParams(PerHostAnophParams&&) noexcept = default;
    PerHostAnophParams& operator=(PerHostAnophParams&&) noexcept = default;

    PerHostAnophParams(const PerHostAnophParams&) = delete;
    PerHostAnophParams& operator=(const PerHostAnophParams&) = delete;

    static inline void init (const scnXml::Mosq& mosq) {
        params.emplace_back(mosq);
    }
    
    /// Get the number of vector species
    static inline size_t numSpecies() {
        return params.size();
    }
    
    /// Get parameters for the given vector species
    static inline const PerHostAnophParams& get(size_t species) {
        return params[species];
    }
    
    /** entoAvailability is calculated externally, then set after other
     * parameters have been initialised.
     * 
     * This function doesn't need to exist, but helps make this fact obvious.
     * 
     * It should be called exactly once.
      */
    inline static void setEntoAvailabilityFactor(size_t species, double entoAvailability){
        params[species].entoAvailabilityFactor = entoAvailability;
    }

    inline static void calcAvailabilityPercentiles(){
        int nSamples = 100000;
        vector<double> samples(nSamples, 0.0);

        // Using a different seed to make sure old results will not change
        // Does it really matter? 
        util::MasterRng seed(0, 0);
        seed.seed(0, 0);
        
        util::LocalRng rng(seed);

        // Sample
        for(int i=0; i<nSamples; i++)
        {
            double avail = 0.0;
            for(size_t s = 0; s < Transmission::PerHostAnophParams::numSpecies(); ++s)
                avail += Transmission::PerHostAnophParams::get(s).entoAvailabilityFactor * Transmission::PerHostAnophParams::get(s).entoAvailability->sample(rng);
            samples[i] = avail;
        }

        // Sort samples
        sort(samples.begin(), samples.end());

        entoAvailabilityPercentiles.resize(100);

        // Calc percetiles threshold values
        for(int i=0; i<100; i++)
            entoAvailabilityPercentiles[i] = samples[int(i*nSamples/100)];
        entoAvailabilityPercentiles[0] = 0.0;
    }

    inline static double getEntoAvailabilityPercentile(int p){
        if(p >= 0 && p < 100)
            return entoAvailabilityPercentiles[p];
        else if(p == 100)
            return std::numeric_limits<double>::max();
        else throw util::xml_scenario_error("Availability percentiles in intervention specifications can only be in [0 and 100]");
    }

    /** @brief Probabilities of finding a host and surviving a feeding cycle
     * 
     * These parameters describe the mean and heterogeneity of α_i, P_B_i,
     * P_C_i and P_D_i across the human population. */
    //@{
    /** Availability rate (α_i) */
    double entoAvailabilityFactor = 1.0;

    /** Availability heterogeneity factor
     * This multiplies the entoAvailabilityFactor. */
    unique_ptr<util::Sampler> entoAvailability;

    /** Probability of mosquito successfully biting host (P_B_i) */
    util::BetaSampler probMosqBiting;

    /** Probability of mosquito escaping human and finding a resting site without
     * dying, after biting the human (P_C_i). */
    util::BetaSampler probMosqFindRestSite;

    /** Probability of mosquito successfully resting after finding a resting site
     * (P_D_i). */
    util::BetaSampler probMosqSurvivalResting;
    //@}
    
private:
    static vector<PerHostAnophParams> params;
    static vector<double> entoAvailabilityPercentiles;

};

/**
 * A base class for interventions affecting human-vector interaction.
 * 
 * The constructor should initialise the data to represent an intervention
 * deployed at this time (sim::now()).
 * 
 * redeploy() should reset the intervention to a freshly deployed state. If
 * necessary, PerHost::deployComponent can be updated to make it create a new
 * instance instead of calling redeploy.
 */
class PerHostInterventionData {
public:
    virtual ~PerHostInterventionData() {}
    
    /** Deploy an intervention. */
    virtual void redeploy( LocalRng& rng, const HumanVectorInterventionComponent& params ) =0;
    
    /** Per time step update. Used by ITNs to update hole decay. */
    virtual void update(Host::Human& human) =0;
    
    /** Get effect of deterrencies of interventions, as an attractiveness multiplier.
     * 
     * @return a value describing effect on attractiveness. Must not be
     * negative. 0 means mosquitoes are fully deterred, 1 that the intervention
     * has no effect, 2 that the intervention attracts twice as many mosquitoes
     * as would otherwise come. */
    virtual double relativeAttractiveness(size_t species) const =0;
    /** Get the killing effect on mosquitoes before they've eaten as a survival
     * multiplier. */
    virtual double preprandialSurvivalFactor(size_t species) const =0;
    /** Get the killing effect on mosquitoes after they've eaten as a survival
     * multiplier. */
    virtual double postprandialSurvivalFactor(size_t species) const =0;
    /// Get the mosquito fecundity multiplier (1 for no effect).
    virtual double relFecundity(size_t species) const =0;
    
    /// Index of effect describing the intervention
    inline interventions::ComponentId id() const { return m_id; }
    
    /// Return true if this component is deployed (i.e. currently active)
    inline bool isDeployed() const{ return deployTime != sim::never(); }
    
    /// Checkpointing (write only)
    void operator& (ostream& stream) {
        m_id & stream; // must be first; read externally so that the correct
                // makeHumanPart function can be found
        checkpoint( stream );
    }
    
protected:
    /// Set the component id
    explicit PerHostInterventionData( interventions::ComponentId id ) :
            deployTime( sim::now() ),
            m_id(id) {}
    
    /// Checkpointing: write
    virtual void checkpoint( ostream& stream ) =0;
    
    SimTime deployTime = sim::never();        // time of deployment or sim::never()
    interventions::ComponentId m_id;       // component id; don't change
};

/** A base class for human vector intervention parameters. */
class HumanVectorInterventionComponent : public interventions::HumanInterventionComponent {
public:
    virtual ~HumanVectorInterventionComponent() {}
    
    /** Create a new object to store human-specific details of deployment. */
    virtual unique_ptr<PerHostInterventionData> makeHumanPart(LocalRng&) const =0;
    virtual unique_ptr<PerHostInterventionData> makeHumanPart( istream& stream,
            interventions::ComponentId id ) const =0;
protected:
    explicit HumanVectorInterventionComponent(interventions::ComponentId id) :
            HumanInterventionComponent(id) {}
};

/** Contains TransmissionModel parameters which need to be stored per host.
 *
 * Currently many members are public and directly accessed. */
class PerHost
{
public:
    /// @brief Static member functions
    //@{
    /** Static initialisation. */
    static void init (const scnXml::AgeGroupValues& availabilityToMosquitoes);
    //@}
    
    ///@brief Initialisation / checkpionting
    //@{
    PerHost ();
    void initialise (LocalRng& rng, double availabilityFactor);
    //@}
    
    /// Call once per time step. Updates net holes.
    void update(Host::Human& human);
  
    /// Deploy some intervention component
    void deployComponent( LocalRng& rng, const HumanVectorInterventionComponent& params );
    //@}
    
    /** Calculates the adjustment for body size in exposure to mosquitoes,
     * relative to an average adult.
     * 
     * The bites are assumed proportional to average surface area for hosts of
     * the given age. Linear interpolation is used to calculate this from the
     * input array of surface areas. 
     * 
     * @param ageYears Age of host
     * @return the ratio of bites received by the host to the average for an adult 
     *
     * This is the age factor of availiability; mean output should be
     * mean population availability (that is, 1.0/invMeanPopAvail).
     * 
     * Also has a switch to put individuals entirely outside transmission. */
    inline double relativeAvailabilityAge (double ageYears) const {
        return outsideTransmission ? 0.0 : relAvailAge.eval( ageYears );
    }
    
    /** Get the availability of this host to mosquitoes relative to an average
     * adult (including heterogeneity and age effects).
     *
     * Used to drive a simulation from an input EIR.
     * Is relativeAvailabilityHet()*relativeAvailabilityAge(ageYears).
     * 
     * Mean output is less than 1.0 (roughly 1.0/invMeanPopAvail).
     */
    inline double relativeAvailabilityHetAge (double ageYears) const {
        return relativeAvailabilityHet * relativeAvailabilityAge (ageYears);
    }
    
    /** Availability rate of human to mosquitoes (α_i). Equals 
     * entoAvailabilityHetVecItv()*getRelativeAvailability().
     *
     * To be clear, this includes effects from HetVecItc (het, interv, availability
     * rate) as well as age (avail. relative to an adult). It does not divide
     * by the average availability of the population, which was incorrectly done
     * in the past. */
    inline double entoAvailabilityFull (size_t species, double ageYears) const {
        return entoAvailabilityHetVecItv (species) * relativeAvailabilityAge (ageYears);
    }
    //@}

    /** Availability of host to mosquitoes (α_i) excluding age factor
     * 
     * (Includes heterogeneity, intervention, and human-to-vector availability
     * rate factors.)
     * 
     * Assume mean is human-to-vector availability rate factor. */
    double entoAvailabilityHetVecItv (size_t species) const;
    
    ///@brief Get effects of interventions pre/post biting
    //@{
    /** Probability of a mosquito succesfully biting a host (P_B_i). */
    double probMosqBiting (size_t species) const;
    /** Probability of a mosquito succesfully finding a resting
     * place after biting and then resting (P_C_i * P_D_i). */
    double probMosqResting (size_t species) const;
    /** Multiplicative factor for the number of fertile eggs laid by mosquitoes
     * after feeding on this host. Should be 1 normally, less than 1 to reduce
     * fertility, greater than 1 to increase. */
    double relMosqFecundity (size_t species) const;
    //@}
    
    ///@brief Convenience wrappers around several functions
    //@{
    /// entoAvailabilityHetVecItv * probMosqBiting
    inline double availBite (size_t species) const{
        return entoAvailabilityHetVecItv(species) * probMosqBiting(species);
    }
    //@}
    
    ///@brief Miscellaneous
    //@{
    /** Get the age at which individuals are considered adults (i.e. where
     * availability to mosquitoes reaches its maximum). */
    static inline SimTime adultAge() {
        return sim::fromYearsD( relAvailAge.firstGlobalMaximum() );
    }
    
    /** Get whether the user has any active deployments of interventions of
     * the given type, where type is one of those interventions deriving the
     * PerHostInterventionData class (in other cases this will always return
     * false). */
    bool hasActiveInterv( interventions::Component::Type type ) const;
    
    /// Checkpointing
    template<class S>
    void operator& (S& stream) {
        anophEntoAvailabilityRaw & stream;
        anophEntoAvailability & stream;
        anophProbMosqBiting & stream;
        anophProbMosqResting & stream;
        relativeAvailabilityHet & stream;
        outsideTransmission & stream;
        copulaGaussianSamples & stream;
        checkpointIntervs( stream );
    }
    //@}
    
    // Determines whether human is outside transmission
    bool outsideTransmission;

    // Heterogeneity factor in availability; this is already multiplied into the
    // entoAvailability param stored in HostMosquitoInteraction.
    double relativeAvailabilityHet;

    /** Species availability rate of human to mosquitoes */
    vector<double> anophEntoAvailabilityRaw;

    /** Species availability rate of human to mosquitoes, including heterogeneity factor
     * and base rate, but excluding age and intervention factors. */
    vector<double> anophEntoAvailability;
    
    /** Species probability of mosquito successfully biting host (P_B_i) in the absense of
     * interventions. */
    vector<double> anophProbMosqBiting;
    
    /** Species probability of mosquito escaping human and finding a resting site, then
     * resting without dying, after biting the human (P_C_i * P_D_i) in the
     * absense of interventions. */
    vector<double> anophProbMosqResting;

    /** @brief gaussian sample used in the deployment of intervention 
     * components using the gaussian copula transformation. The factors are 
     * component specific and are sampled once per human during the first 
     * deployment attempt. */
    std::map<interventions::ComponentId, double> copulaGaussianSamples; 

private:
    void checkpointIntervs( ostream& stream );
    void checkpointIntervs( istream& stream );

    vector<unique_ptr<PerHostInterventionData>> activeComponents;
    
    static AgeGroupInterpolator relAvailAge;
};

}
}
#endif
