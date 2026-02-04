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

#include "Global.h"
#include "interventions/GVI.h"
#include "Host/Human.h"
#include "util/SpeciesIndexChecker.h"
#include "util/errors.h"
#include "util/CommandLine.h"
#include <cmath>

namespace OM { namespace interventions {

vector<GVIComponent*> GVIComponent::componentsByIndex;

GVIComponent::GVIComponent( ComponentId id, const scnXml::GVIDescription& elt,
        const map<string,size_t>& species_name_map ) :
        Transmission::HumanVectorInterventionComponent(id)
{
    decay = DecayFunction::makeObject( elt.getDecay(), "interventions.human.vector.decay" );
    // assume usage modifier is 100% if none is specified
    double propUse;
    if (elt.getUsage().present()) {
        propUse = elt.getUsage().get().getValue();
    }
    else {
        propUse = 1.0;
    }
    if( !( propUse >= 0.0 && propUse <= 1.0 ) ){
        throw util::xml_scenario_error("GVI.description.usage: must be within range [0,1]");
    }
    
    typedef scnXml::GVIDescription::AnophelesParamsSequence AP;
    const AP& ap = elt.getAnophelesParams();
    species.resize( species_name_map.size() );
    util::SpeciesIndexChecker checker( "GVI intervention", species_name_map );
    for( AP::const_iterator it = ap.begin(); it != ap.end(); ++it ) {
        species[checker.getIndex(it->getMosquito())].init (*it, propUse);
    }
    checker.checkNoneMissed();
    
    if( componentsByIndex.size() <= id.id ) componentsByIndex.resize( id.id+1, 0 );
    componentsByIndex[id.id] = this;
}

void GVIComponent::deploy( Host::Human& human, mon::Deploy::Method method, VaccineLimits )const{
    human.perHostTransmission.deployComponent(human.rng, *this);
    mon::reportEventMHD( mon::MHD_GVI, human, method );
}

Component::Type GVIComponent::componentType()const{ return Component::GVI; }

void GVIComponent::print_details( std::ostream& out )const{
    out << id().id << "\tGVI";
}

unique_ptr<PerHostInterventionData> GVIComponent::makeHumanPart(LocalRng& rng) const{
    return unique_ptr<PerHostInterventionData>(new HumanGVI( rng, *this ));
}
unique_ptr<PerHostInterventionData> GVIComponent::makeHumanPart( istream& stream, ComponentId id ) const{
    return unique_ptr<PerHostInterventionData>(new HumanGVI( stream, id ));
}

void GVIComponent::GVIAnopheles::init(const scnXml::GVIDescription::AnophelesParamsType& elt,
                                          double proportionUse)
{
    assert( (std::isnan)(deterrency) ); // double init
    deterrency = elt.getDeterrency().present() ? elt.getDeterrency().get().getValue() : 0.0;
    preprandialKilling = elt.getPreprandialKillingEffect().present() ? elt.getPreprandialKillingEffect().get().getValue() : 0.0;
    postprandialKilling = elt.getPostprandialKillingEffect().present() ? elt.getPostprandialKillingEffect().get().getValue() : 0.0;
    fecundityReduction = elt.getFecundityReduction().present() ? elt.getFecundityReduction().get().getValue() : 0.0;
    // Simpler version of ITN usage/action:
    double propActive = elt.getPropActive();
    if(propActive != 1.0 && propActive != 0.0)
        throw util::xml_scenario_error("GVI.description.propActive: must be either 0 or 1");

    assert( proportionUse >= 0.0 && proportionUse <= 1.0 );
    assert( propActive >= 0.0 && propActive <= 1.0 );
    proportionProtected = proportionUse * propActive;
    proportionUnprotected = 1.0 - proportionProtected;
}


// ———  per-human data  ———
HumanGVI::HumanGVI ( LocalRng& rng, const GVIComponent& params ) :
    PerHostInterventionData( params.id() )
{
    // Varience factor of decay is sampled once per human: human is assumed
    // to account for most variance.
    decayHet = params.decay->hetSample(rng);
}

void HumanGVI::redeploy(LocalRng& rng, const Transmission::HumanVectorInterventionComponent&) {
    deployTime = sim::nowOrTs1();
}

void HumanGVI::update(Host::Human& human){
}

double HumanGVI::relativeAttractiveness(size_t speciesIndex) const{
    const GVIComponent& params = *GVIComponent::componentsByIndex[m_id.id];
    const GVIComponent::GVIAnopheles& anoph = params.species[speciesIndex];
    double effect = (1.0 - anoph.deterrency * getEffectSurvival(params));
    return anoph.byProtection( effect );
}

double HumanGVI::preprandialSurvivalFactor(size_t speciesIndex) const{
    const GVIComponent& params = *GVIComponent::componentsByIndex[m_id.id];
    const GVIComponent::GVIAnopheles& anoph = params.species[speciesIndex];
    double effect = (1.0 - anoph.preprandialKilling * getEffectSurvival(params));
    return anoph.byProtection( effect );
}

double HumanGVI::postprandialSurvivalFactor(size_t speciesIndex) const{
    const GVIComponent& params = *GVIComponent::componentsByIndex[m_id.id];
    const GVIComponent::GVIAnopheles& anoph = params.species[speciesIndex];
    double effect = (1.0 - anoph.postprandialKilling * getEffectSurvival(params));
    return anoph.byProtection( effect );
}
double HumanGVI::relFecundity(size_t speciesIndex) const{
    const GVIComponent& params = *GVIComponent::componentsByIndex[m_id.id];
    const GVIComponent::GVIAnopheles& anoph = params.species[speciesIndex];
    double effect = (1.0 - anoph.fecundityReduction * getEffectSurvival(params));
    return anoph.byProtection( effect );
}

void HumanGVI::checkpoint( ostream& stream ){
    deployTime & stream;
    decayHet & stream;
}
HumanGVI::HumanGVI( istream& stream, ComponentId id ) : PerHostInterventionData( id )
{
    deployTime & stream;
    decayHet & stream;
}

} }
