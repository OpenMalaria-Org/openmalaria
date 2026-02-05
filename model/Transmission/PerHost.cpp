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

#include "Transmission/PerHost.h"
#include "Transmission/VectorModel.h"
#include "interventions/InterventionManager.h"
#include "util/errors.h"
#include "util/checkpoint.h"

namespace OM {
namespace Transmission {
using namespace OM::util;

// -----  PerHost static  -----

AgeGroupInterpolator PerHost::relAvailAge;
vector<PerHostAnophParams> PerHostAnophParams::params;
vector<double> PerHostAnophParams::entoAvailabilityPercentiles;

PerHostAnophParams::PerHostAnophParams (const scnXml::Mosq& mosq) {
    if(mosq.getAvailability().getDistr() != "const" && (util::ModelOptions::option(NEGATIVE_BINOMIAL_MASS_ACTION) || util::ModelOptions::option(LOGNORMAL_MASS_ACTION) || util::ModelOptions::option(COMORB_TRANS_HET)))
        throw util::xml_scenario_error( "PerHostAnophParams::PerHostAnophParams(): ModelOptions NEGATIVE_BINOMIAL_MASS_ACTION, LOGNORMAL_MASS_ACTION, COMORB_TRANS_HET are not compatible with ento-availability heterogeneity distributions other than \"const\"");

    entoAvailability = util::createSampler(1.0, mosq.getAvailability());
    probMosqBiting.setParams( mosq.getMosqProbBiting() );
    probMosqFindRestSite.setParams( mosq.getMosqProbFindRestSite() );
    probMosqSurvivalResting.setParams( mosq.getMosqProbResting() );
}

void PerHost::init ( const scnXml::AgeGroupValues& availabilityToMosquitoes ) {
    relAvailAge.set( availabilityToMosquitoes, "availabilityToMosquitoes" );
}

// -----  PerHost non-static -----

PerHost::PerHost () :
        outsideTransmission(false),
        relativeAvailabilityHet(numeric_limits<double>::signaling_NaN())
{
}

void PerHost::initialise (LocalRng& rng, double availabilityFactor) {
    relativeAvailabilityHet = availabilityFactor;
    anophEntoAvailabilityRaw.resize(PerHostAnophParams::numSpecies());
    anophEntoAvailability.resize(PerHostAnophParams::numSpecies());
    anophProbMosqBiting.resize(PerHostAnophParams::numSpecies());
    anophProbMosqResting.resize(PerHostAnophParams::numSpecies());

    for(size_t i = 0; i < PerHostAnophParams::numSpecies(); ++i) {
        const PerHostAnophParams& base = PerHostAnophParams::get(i);
        anophEntoAvailabilityRaw[i] = base.entoAvailability->sample(rng);
        anophEntoAvailability[i] = base.entoAvailabilityFactor * anophEntoAvailabilityRaw[i] * availabilityFactor;
        anophProbMosqBiting[i] = base.probMosqBiting.sample(rng);
        auto pRest1 = base.probMosqFindRestSite.sample(rng);
        auto pRest2 = base.probMosqSurvivalResting.sample(rng);
        anophProbMosqResting[i] = pRest1 * pRest2;
    }
}

void PerHost::update(Host::Human& human){
    for( auto iter = activeComponents.begin(); iter != activeComponents.end(); ++iter ){
        (*iter)->update(human);
    }
}

void PerHost::deployComponent( LocalRng& rng, const HumanVectorInterventionComponent& params ){
    // This adds per-host per-intervention details to the host's data set.
    // This data is never removed since it can contain per-host heterogeneity samples.
    for( auto iter = activeComponents.begin(); iter != activeComponents.end(); ++iter ){
        if( (*iter)->id() == params.id() ){
            // already have a deployment for that description; just update it
            (*iter)->redeploy( rng, params );
            return;
        }
    }
    // no deployment for that description: must make a new one
    activeComponents.push_back( params.makeHumanPart(rng) );
}

// Note: in the case an intervention is not present, we can use the approximation
// of Weibull decay over the time span now - sim::never()
// (easily large enough for conceivable Weibull params that the value is 0.0 when
// rounded to a double. Performance-wise it's perhaps slightly slower than using
// an if() when interventions aren't present.
double PerHost::entoAvailabilityHetVecItv (size_t species) const {
    double alpha_i = anophEntoAvailability[species];
    for( auto iter = activeComponents.begin(); iter != activeComponents.end(); ++iter ){
        alpha_i *= (*iter)->relativeAttractiveness( species );
    }
    return alpha_i;
}
double PerHost::probMosqBiting (size_t species) const {
    double P_B_i = anophProbMosqBiting[species];
    for( auto iter = activeComponents.begin(); iter != activeComponents.end(); ++iter ){
        P_B_i *= (*iter)->preprandialSurvivalFactor( species );
    }
    return P_B_i;
}
double PerHost::probMosqResting (size_t species) const {
    double pRest = anophProbMosqResting[species];
    for( auto iter = activeComponents.begin(); iter != activeComponents.end(); ++iter ){
        pRest *= (*iter)->postprandialSurvivalFactor( species );
    }
    return pRest;
}

double PerHost::relMosqFecundity (size_t species) const {
    double relFecundity = 1.0;
    for( auto iter = activeComponents.begin(); iter != activeComponents.end(); ++iter ){
        relFecundity *= (*iter)->relFecundity( species );
    }
    return relFecundity;
}

bool PerHost::hasActiveInterv(interventions::Component::Type type) const{
    for( auto iter = activeComponents.begin(); iter != activeComponents.end(); ++iter ){
        if( (*iter)->isDeployed() ){
            if( interventions::InterventionManager::getComponent( (*iter)->id() ).componentType() == type )
                return true;
        }
    }
    return false;
}

void PerHost::checkpointIntervs( ostream& stream ){
    activeComponents.size() & stream;
    for( auto iter = activeComponents.begin(); iter != activeComponents.end(); ++iter ){
        *iter & stream;
    }
}
void PerHost::checkpointIntervs( istream& stream ){
    size_t l;
    l & stream;
    validateListSize(l);
    activeComponents.clear();
    for( size_t i = 0; i < l; ++i ){
        interventions::ComponentId id( stream );
        try{
            const interventions::HumanInterventionComponent& gen_params =
                interventions::InterventionManager::getComponent( id );   // may throw
            const HumanVectorInterventionComponent *params =
                dynamic_cast<const HumanVectorInterventionComponent*>( &gen_params );
            if( params == 0 )
                throw util::base_exception( "" );       // see catch block below
            activeComponents.push_back( params->makeHumanPart( stream, id ) );
        }catch( util::base_exception& e ){
            // two causes, both boil down to index being wrong
            throw util::checkpoint_error( "bad value in checkpoint file" );
        }
    }
}

}
}
