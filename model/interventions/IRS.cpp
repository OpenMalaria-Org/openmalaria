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
#include "interventions/IRS.h"
#include "Host/Human.h"
#include "util/errors.h"
#include "util/SpeciesIndexChecker.h"
#include "util/CommandLine.h"
#include "R_nmath/qnorm.h"
#include <cmath>

namespace OM { namespace interventions {
vector<IRSComponent*> IRSComponent::componentsByIndex;

IRSComponent::IRSComponent( ComponentId id, const scnXml::IRSDescription& elt,
        const map<string,size_t>& species_name_map ) :
        Transmission::HumanVectorInterventionComponent(id)
{
    initialInsecticide.setParams( elt.getInitialInsecticide() );
    const double maxProp = 0.999;       //NOTE: this could be exposed in XML, but probably doesn't need to be
    maxInsecticide = R::qnorm5(maxProp, initialInsecticide.getMu(), initialInsecticide.getSigma(), true, false);
    insecticideDecay = DecayFunction::makeObject( elt.getInsecticideDecay(),
                                                  "interventions.human.IRS.description.insecticideDecay" );
    // assume usage modifier is 100% if none is specified
    double propUse;
    if (elt.getUsage().present()) {
        propUse = elt.getUsage().get().getValue();
    }
    else {
        propUse = 1.0;
    }
    if( !( propUse >= 0.0 && propUse <= 1.0 ) ){
        throw util::xml_scenario_error("ITN.description.usage: must be within range [0,1]");
    }
    
    typedef scnXml::IRSDescription::AnophelesParamsSequence AP;
    const AP& ap = elt.getAnophelesParams();
    species.resize( species_name_map.size() );
    util::SpeciesIndexChecker checker( "IRS intervention", species_name_map );
    for( AP::const_iterator it = ap.begin(); it != ap.end(); ++it ) {
        species[checker.getIndex(it->getMosquito())].init (*it, propUse, maxInsecticide);
    }
    checker.checkNoneMissed();
    
    if( componentsByIndex.size() <= id.id ) componentsByIndex.resize( id.id+1, 0 );
    componentsByIndex[id.id] = this;
}

void IRSComponent::deploy( Host::Human& human, mon::Deploy::Method method, VaccineLimits )const{
    human.perHostTransmission.deployComponent(human.rng, *this);
    mon::reportEventMHD( mon::MHD_IRS, human, method );
}

Component::Type IRSComponent::componentType()const{ return Component::IRS; }

void IRSComponent::print_details( std::ostream& out )const{
    out << id().id << "\tIRS";
}

unique_ptr<PerHostInterventionData> IRSComponent::makeHumanPart(LocalRng& rng) const{
    return unique_ptr<PerHostInterventionData>(new HumanIRS( rng, *this ));
}
unique_ptr<PerHostInterventionData> IRSComponent::makeHumanPart( istream& stream, ComponentId id ) const{
    return unique_ptr<PerHostInterventionData>(new HumanIRS( stream, id ));
}

void IRSComponent::IRSAnopheles::init(
    const scnXml::IRSDescription::AnophelesParamsType& elt,
    double proportionUse,
    double maxInsecticide)
{
    _relativeAttractiveness.init( elt.getDeterrency() );
    _preprandialKillingEffect.init( elt.getPreprandialKillingEffect(), false, maxInsecticide );
    _postprandialKillingEffect.init( elt.getPostprandialKillingEffect(), true, maxInsecticide );
    if (elt.getFecundityReduction().present()) {
        _fecundityEffect.init( elt.getFecundityReduction().get(), false/*TODO: err msg*/, maxInsecticide );
    } else {
        _fecundityEffect.init1();
    }
    assert( proportionUse >= 0.0 && proportionUse <= 1.0 );
    // Simpler version of ITN usage/action:
    double propActive = elt.getPropActive();
    if(propActive != 1.0 && propActive != 0.0)
        throw util::xml_scenario_error("GVI.description.propActive: must be either 0 or 1");

    assert( propActive >= 0.0 && propActive <= 1.0 );
    proportionProtected = proportionUse * propActive;
    proportionUnprotected = 1.0 - proportionProtected;
}

IRSComponent::IRSAnopheles::RelativeAttractiveness::RelativeAttractiveness() :
    lPF( numeric_limits< double >::signaling_NaN() ),
    insecticideScaling( numeric_limits< double >::signaling_NaN() )
{}
void IRSComponent::IRSAnopheles::RelativeAttractiveness::init(const scnXml::IRSDeterrency& elt){
    double PF = elt.getInsecticideFactor();
    insecticideScaling = elt.getInsecticideScalingFactor();
    if( !( PF > 0.0) ){
        ostringstream msg;
        msg << "IRS.description.anophelesParams.relativeAvailability: expected insecticideFactor to be positive.";
        throw util::xml_scenario_error( msg.str() );
    }
    
    /* We need to ensure the relative availability is non-negative. However,
     * since it's an exponentiated value, it always will be.
     * 
     * If we don't want IRS to be able to increase transmission, the following
     * limits could also be applied. In general, however, there is no reason
     * IRS couldn't make individuals more attractive to mosquitoes.
     * 
     * To ensure relative availability is at most one: relative availability is
     *  exp( log(PF)*p ) = PF^p
     * where PF is the insecticide factor, with p∈[0,1] defined as:
     *  p=1−exp(-insecticideContent*insecticideScalingFactor).
     * We therefore just need PF ≤ 1. */
//     if( !( PF <= 1.0 ) ) {
        // Potentially warn about this... but not necessary since making humans
        // more attractive isn't really an issue.
//     }
    assert( (std::isnan)(lPF) ); // double init
    lPF = log( PF );
}
IRSComponent::IRSAnopheles::SurvivalFactor::SurvivalFactor() :
    BF( numeric_limits< double >::signaling_NaN() ),
    PF( numeric_limits< double >::signaling_NaN() ),
    insecticideScaling( numeric_limits< double >::signaling_NaN() ),
    invBaseSurvival( numeric_limits< double >::signaling_NaN() )
{}
void IRSComponent::IRSAnopheles::SurvivalFactor::init(const scnXml::IRSKillingEffect& elt,
                                                   bool postPrandial, double maxInsecticide){
    BF = elt.getBaseFactor();
    PF = elt.getInsecticideFactor();
    insecticideScaling = elt.getInsecticideScalingFactor();
    invBaseSurvival = 1.0 / (1.0 - BF);
    if( !( BF >= 0.0 && BF < 1.0) ){
        ostringstream msg;
        msg << "IRS.description.anophelesParams." << (postPrandial?"post":"pre")
            << "killingFactor: expected baseFactor to be in range [0,1]";
        throw util::xml_scenario_error( msg.str() );
    }
    if( !(insecticideScaling>=0.0) ){
        ostringstream msg;
        msg << "IRS.description.anophelesParams." << (postPrandial?"post":"pre")
            << "killingFactor: expected scaling factor to be non-negative";
        throw util::xml_scenario_error( msg.str() );
    }
    
    /* We want the calculated survival factor (1−K)/(1−BF) to be in the range
    [0,1] where K is the killing factor: K=BF+PF×p, with p defined as:
    p=1−exp(-insecticideContent×insecticideScalingFactor).
    
    By its nature, insecticideContent ≥ 0. We restrict
        insecticideScalingFactor ≥ 0
    which implies that p lies in the range [0,1]. We also know the base
    survival factor, 1−BF, is in the range [0,1].
    
    To make sure the survival factor is not negative we need (1−K)/(1−BF) ≥ 0.
    Since 1−BF > 0 we need 1−K ≥ 0, which, substituting K, gives us
        BF + PF×p ≤ 1	(1)
    We also want to make sure the survival factor is not greater than one (since
    IRS shouldn't increase mosquito survival), (1−K)/(1−BF) ≤ 1 or equivalently
    1-K ≤ 1-BF or K ≥ BF, which, substituting K, yields
        PF×p ≥ 0		(2)
    Since p ≥ 0, ensuring PF is non-negative is enough for (2).
    
    The maximum value for p depends on the maximum insecticide content; denote
    pmax = max(p). In this extreme case (1) becomes:
        BF + PF×pmax ≤ 1	(3)
    
    As with the ITN model, we impose a maximum value on the initial insecticide
    content, Pmax, such that the probability of sampling a value from our
    parameterise normal distribution greater than Pmax is 0.001. */
    double pmax = 1.0-exp(-maxInsecticide*insecticideScaling);
    if( !( PF >= 0.0 && BF+PF*pmax <= 1.0 ) ){
        ostringstream msg;
        msg << "IRS.description.anophelesParams." << (postPrandial?"post":"pre")
            << "killingFactor: expected insecticideFactor≥0, baseFactor+"
            <<pmax<<"×insecticideFactor≤1";
        throw util::xml_scenario_error( msg.str() );
    }
}
void IRSComponent::IRSAnopheles::SurvivalFactor::init1(){
    BF = 0.0;
    PF = 0.0;
    insecticideScaling = 0.0;
    invBaseSurvival = 1.0;
}
double IRSComponent::IRSAnopheles::RelativeAttractiveness::relativeAttractiveness(
        double insecticideContent
) const {
    double insecticideComponent = 1.0 - exp(-insecticideContent*insecticideScaling);
    double relAvail = exp( lPF*insecticideComponent );
    assert( relAvail>=0.0 );
    return relAvail;
}
double IRSComponent::IRSAnopheles::SurvivalFactor::survivalFactor(
        double insecticideContent
) const {
    double insecticideComponent = 1.0 - exp(-insecticideContent*insecticideScaling);
    double killingEffect = BF + PF*insecticideComponent;
    double survivalFactor = (1.0 - killingEffect) * invBaseSurvival;
    assert( survivalFactor >= 0.0 );
    assert( survivalFactor <= 1.0 );
    return survivalFactor;
}

double IRSComponent::sampleInitialInsecticide(LocalRng& rng) const{
    double sample = initialInsecticide.sample(rng);
    if( sample < 0.0 )
        sample = 0.0;       // avoid negative samples
    if( sample > maxInsecticide )
        sample = maxInsecticide;
    return sample;
}


// ———  per-human data  ———
HumanIRS::HumanIRS( LocalRng& rng, const IRSComponent& params ) :
    PerHostInterventionData( params.id() ),
    initialInsecticide( params.sampleInitialInsecticide(rng) )
{
    // Varience factor of decay is sampled once per human: human is assumed
    // to account for most variance.
    insecticideDecayHet = params.insecticideDecay->hetSample(rng);
}

void HumanIRS::redeploy( LocalRng& rng, const OM::Transmission::HumanVectorInterventionComponent& params ) {
    deployTime = sim::nowOrTs1();
    const IRSComponent* irsParams = dynamic_cast<const IRSComponent*>(&params);
    assert( irsParams != 0 );   // code error if this fails
    initialInsecticide = irsParams->sampleInitialInsecticide(rng);
}

void HumanIRS::update(Host::Human& human){
}

double HumanIRS::relativeAttractiveness(size_t speciesIndex) const{
    const IRSComponent& params = *IRSComponent::componentsByIndex[m_id.id];
    const IRSComponent::IRSAnopheles& anoph = params.species[speciesIndex];
    double effect = anoph.relativeAttractiveness( getInsecticideContent(params) );
    return anoph.byProtection( effect );
}

double HumanIRS::preprandialSurvivalFactor(size_t speciesIndex) const{
    const IRSComponent& params = *IRSComponent::componentsByIndex[m_id.id];
    const IRSComponent::IRSAnopheles& anoph = params.species[speciesIndex];
    double effect = anoph.preprandialSurvivalFactor( getInsecticideContent(params) );
    return anoph.byProtection( effect );
}

double HumanIRS::postprandialSurvivalFactor(size_t speciesIndex) const{
    const IRSComponent& params = *IRSComponent::componentsByIndex[m_id.id];
    const IRSComponent::IRSAnopheles& anoph = params.species[speciesIndex];
    double effect = anoph.postprandialSurvivalFactor( getInsecticideContent(params) );
    return anoph.byProtection( effect );
}
double HumanIRS::relFecundity(size_t speciesIndex) const{
    const IRSComponent& params = *IRSComponent::componentsByIndex[m_id.id];
    const IRSComponent::IRSAnopheles& anoph = params.species[speciesIndex];
    double effect = anoph.fecundityEffect( getInsecticideContent(params) );
    return anoph.byProtection( effect );
}

void HumanIRS::checkpoint( ostream& stream ){
    deployTime & stream;
    initialInsecticide & stream;
    insecticideDecayHet & stream;
}
HumanIRS::HumanIRS( istream& stream, ComponentId id ) : PerHostInterventionData( id )
{
    deployTime & stream;
    initialInsecticide & stream;
    insecticideDecayHet & stream;
}

} }
