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

#include "interventions/ITN.h"
#include "util/random.h"
#include "util/errors.h"
#include "util/SpeciesIndexChecker.h"
#include "Host/Human.h"
#include "R_nmath/qnorm.h"
#include "util/CommandLine.h"
#include <cmath>

namespace OM { namespace interventions {

vector<ITNComponent*> ITNComponent::componentsByIndex;

// —————  utility classes (internal use only)  —————


void factors::SurvivalFactor::init(const scnXml::ITNKillingEffect& elt,
                                                   double maxInsecticide, const char* eltName,
                                                   bool isDeterrent){
    useLogitEqns = false;
    a.BF = elt.getBaseFactor();
    a.HF = elt.getHoleFactor();
    a.PF = elt.getInsecticideFactor();
    a.IF = elt.getInteractionFactor();
    a.holeScaling = elt.getHoleScalingFactor();
    a.insecticideScaling = elt.getInsecticideScalingFactor();
    a.invBaseSurvival = 1.0 / (1.0 - a.BF);
    if( !( a.BF >= 0.0 && a.BF < 1.0) ){
        ostringstream msg;
        msg << eltName << ": expected baseFactor to be in range [0,1]";
        throw util::xml_scenario_error( msg.str() );
    }
    if( !(a.holeScaling>=0.0 && a.insecticideScaling>=0.0) ){
        ostringstream msg;
        msg << eltName << ": expected scaling factors to be non-negative";
        throw util::xml_scenario_error( msg.str() );
    }
    // see below
    double pmax = 1.0-exp(-maxInsecticide * a.insecticideScaling);
    
    if( isDeterrent ){
        // Note: the following argument is a modification of the one below
        // (when !isDeterrent). The original may make more sense.
    /* We want K ≥ 0 where K is the killing factor:
    K=BF+HF×h+PF×p+IF×h×p, with h and p defined as:
    h=exp(-holeIndex×holeScalingFactor),
    p=1−exp(-insecticideContent×insecticideScalingFactor). 
    
    By their nature, holeIndex ≥ 0 and insecticideContent ≥ 0. We restrict:
        holeScalingFactor ≥ 0
        insecticideScalingFactor ≥ 0
    Which implies both h and p lie in the range [0,1]. We also know 0 ≤ BF ≤ 1.
    
    We need K ≥ 0 or:
        BF + HF×h + PF×p + IF×h×p ≥ 0   (1)
    
    Lets derive some limits on HF, PF and IF such that the above inequality (1)
    is satisfied.
    
    A net can theoretically be unholed (holeIndex=0 ⇒ h=1) and have no
    insecticide (thus have p=0). Substituting these values in (1) yields:
        BF + HF ≥ 0     (2)
    
    The maximum value for p depends on the maximum insecticide content; denote
    pmax = max(p). Note that holeIndex has no finite maximum; thus, although
    for any finite value of holeIndex, h > 0, there is no h₀ > 0 s.t. for all
    values of holeIndex h ≥ h₀. For the limiting case of a tattered but
    insecticide-saturated net our parameters are therefore p=pmax, h=0:
        BF + PF×pmax ≥ 0        (3)
    
    Consider a net saturated with insecticide (p=pmax) and without holes (h=1):
        BF + HF + (PF+IF)×pmax ≥ 0      (4)
    
    The opposite extreme (the limiting case of a decayed net with no remaining
    insecticide and a large number of holes) yields only BF ≥ 0.
    
    Some of the above examples of nets may be unlikely, but there is only one
    restriction in our model making any of these cases impossible: some
    insecticide must have been lost by the time any holes occur. We ignore this
    since its effect is likely small, and thus all of the above are required to
    keep the factor in the range [0,1]. Further, the inequalities (2) - (3) are
    sufficient to keep the factor within [0,1] since h and p are non-negative
    and act linearly in (1).
    
    From the definition of p, we always have pmax ≤ 1, so substituting pmax=1
    in (5) - (8) gives us bounds which imply our requirement, however if pmax
    is finite they are stricter than necessary. Since insecticideScalingFactor
    is constant, max(p) coincides with max(insecticideContent) which, since
    insecticide content only decays over time, coincides with the maximum
    initial insecticide content, Pmax. Since the initial insecticide content is
    sampled from a normal distribution in our model it should have no finite
    maximum, thus implying we cannot achieve more relaxed bounds than (5) - (8)
    when pmax=1 (unless the standard deviation of our normal distribution is 1).
    We would however like to impose less strict bounds than these, thus we
    impose a maximum value on the initial insecticide content, Pmax, such that
    the probability of sampling a value from our parameterise normal
    distribution greater than Pmax is 0.001. */
        if( !( a.BF + a.HF >= 0.0
            && a.BF + a.PF * pmax >= 0.0
            && a.BF + a.HF + (a.PF + a.IF) * pmax >= 0.0 ) )
        {
            ostringstream msg;
            msg << eltName << ": bounds not met:";
            if( !(a.BF + a.HF >= 0.0) )
                msg << " baseFactor+holeFactor≥0";
            if( !(a.BF + a.PF * pmax >= 0.0) )
                msg << " baseFactor+"<<pmax<<"×insecticideFactor≥0";
            if( !(a.PF + a.HF + (a.PF + a.IF) * pmax >= 0.0) )
                msg << " baseFactor+holeFactor+"<<pmax<<"×(insecticideFactor+interactionFactor)≥0";
            throw util::xml_scenario_error( msg.str() );
        }
    }else{
    /* We want the calculated survival factor (1−K)/(1−BF) to be in the range
    [0,1] where K is the killing factor: K=BF+HF×h+PF×p+IF×h×p, with h and p
    defined as: h=exp(-holeIndex×holeScalingFactor),
    p=1−exp(-insecticideContent×insecticideScalingFactor). 
    
    By their nature, holeIndex ≥ 0 and insecticideContent ≥ 0. We restrict:
        holeScalingFactor ≥ 0
        insecticideScalingFactor ≥ 0
    Which implies both h and p lie in the range [0,1]. We also know the base
    survival factor, 1−BF, is in the range [0,1].
    
    To make sure the survival factor is not negative we need (1−K)/(1−BF) ≥ 0.
    Since 1−BF > 0 we need 1−K ≥ 0, which, substituting K, gives us
        BF + HF×h + PF×p + IF×h×p ≤ 1	(1)
    We also want to make sure the survival factor is not greater than one (since
    nets shouldn't increase mosquito survival), (1−K)/(1−BF) ≤ 1 or equivalently
    1-K ≤ 1-BF or K ≥ BF, which, substituting K, yields
        HF×h + PF×p + IF×h×p ≥ 0		(2)
    
    Lets derive some limits on HF, PF and IF such that the above inequalities
    (1) and (2) are satisfied.
    
    A net can theoretically be unholed (holeIndex=0 ⇒ h=1) and have no
    insecticide (thus have p=0). Substituting these values in (1) and (2) yields:
        BF + HF ≤ 1	(3)
        HF ≥ 0		(4)
    
    The maximum value for p depends on the maximum insecticide content; denote
    pmax = max(p). Note that holeIndex has no finite maximum; thus, although
    for any finite value of holeIndex, h > 0, there is no h₀ > 0 s.t. for all
    values of holeIndex h ≥ h₀. For the limiting case of a tattered but
    insecticide-saturated net our parameters are therefore p=pmax, h=0:
        BF + PF×pmax ≤ 1	(5)
        PF×pmax ≥ 0			(6)
    (Assuming pmax > 0, (6) is equivalent to PF ≥ 0.)
    
    Consider a net saturated with insecticide (p=pmax) and without holes (h=1):
        BF + HF + (PF+IF)×pmax ≤ 1	(7)
        HF + (PF+IF)×pmax ≥ 0		(8)
    
    The opposite extreme (the limiting case of a decayed net with no remaining
    insecticide and a large number of holes) yields only BF ≤ 1 which we already
    know.
    
    Some of the above examples of nets may be unlikely, but there is only one
    restriction in our model making any of these cases impossible: some
    insecticide must have been lost by the time any holes occur. We ignore this
    since its effect is likely small, and thus all of the above are required to
    keep the survival factor in the range [0,1]. Further, these six inequalities
    (3) - (8) are sufficient to keep the survival factor within [0,1] since h
    and p are non-negative and act linearly in (1) and (2).
    
    From the definition of p, we always have pmax ≤ 1, so substituting pmax=1
    in (5) - (8) gives us bounds which imply our requirement, however if pmax
    is finite they are stricter than necessary. Since insecticideScalingFactor
    is constant, max(p) coincides with max(insecticideContent) which, since
    insecticide content only decays over time, coincides with the maximum
    initial insecticide content, Pmax. Since the initial insecticide content is
    sampled from a normal distribution in our model it should have no finite
    maximum, thus implying we cannot achieve more relaxed bounds than (5) - (8)
    when pmax=1 (unless the standard deviation of our normal distribution is 1).
    We would however like to impose less strict bounds than these, thus we
    impose a maximum value on the initial insecticide content, Pmax, such that
    the probability of sampling a value from our parameterise normal
    distribution greater than Pmax is 0.001. */
        if( !( a.BF + a.HF <= 1.0 && a.HF >= 0.0
            && a.BF + a.PF * pmax <= 1.0 && a.PF >= 0.0
            && a.BF + a.HF + (a.PF + a.IF) * pmax <= 1.0 && a.HF + (a.PF + a.IF) * pmax >= 0.0 ) )
        {
            ostringstream msg;
            msg << eltName << ": bounds not met:";
            if( !(a.BF+a.HF<=1.0) )
                msg << " baseFactor+holeFactor≤1";
            if( !(a.HF>=0.0) )
                msg << " holeFactor≥0";
            if( !(a.BF+a.PF*pmax<=1.0) )
                msg << " baseFactor+"<<pmax<<"×insecticideFactor≤1";
            if( !(a.PF>=0.0) )
                msg << " insecticideFactor≥0";      // if this fails, we know pmax>0 (since it is in any case non-negative) — well, or an NaN
            if( !(a.PF+a.HF+(a.PF+a.IF)*pmax<=1.0) )
                msg << " baseFactor+holeFactor+"<<pmax<<"×(insecticideFactor+interactionFactor)≤1";
            if( !(a.HF+(a.PF+a.IF)*pmax>=0.0) )
                msg << " holeFactor+"<<pmax<<"×(insecticideFactor+interactionFactor)≥0";
            throw util::xml_scenario_error( msg.str() );
        }
    }
}

void factors::SurvivalFactor::initLogit(const scnXml::ITNEffectLogit& elt,
                                                   double holeIndexMax,
                                                   bool isDeterrent)
{
    useLogitEqns = true;
    b.BF = elt.getBaseFactor();
    b.HF = elt.getHoleFactor();
    b.PF = elt.getInsecticideFactor();
    b.IF = elt.getInteractionFactor();
    b.hMax = log(holeIndexMax + 1.0);
    const double x = exp(b.BF + b.HF * b.hMax);
    if (isDeterrent) {
        // pAtt0 =x/(x+1); this is 1/pAtt0:
        b.invBaseEffect = (x + 1.0) / x;
    } else {
        // K0 = x/(x+1), so 1/(1-K0) = x+1:
        b.invBaseEffect = x + 1.0;
        
        // We expect K >= K0 and deduce these "advisory limits":
        if (b.PF < 0.0 || b.PF + b.IF * b.hMax < 0.0) {
            cerr << "ITN.description.anophelesParams.*KillingEffectLogit: expected"
                << "\n(insecticide) P >= 0, found P = " << b.PF
                << "\n(interaction) P+I*log(holeIndexMax+1) >= 0, found " << b.PF + b.IF * b.hMax
                << endl;
        }
    }
}

void factors::SurvivalFactor::init1() {
    useLogitEqns = false;
    a.invBaseSurvival = 1.0;
    a.BF = 0.0;
    a.HF = 0.0;
    a.PF = 0.0;
    a.IF = 0.0;
    a.insecticideScaling = 0.0;
    a.holeScaling = 0.0;
}

double factors::SurvivalFactor::rel_pAtt( double holeIndex, double insecticideContent )const {
    if (!useLogitEqns) {
        const double holeComponent = exp(-holeIndex * a.holeScaling);
        const double insecticideComponent = 1.0 - exp(-insecticideContent * a.insecticideScaling);
        const double pAtt = a.BF
                + a.HF * holeComponent
                + a.PF * insecticideComponent
                + a.IF * holeComponent * insecticideComponent;
        assert( pAtt >= 0.0 );
        return pAtt / a.BF;
    } else {
        const double holeComponent = min(log(holeIndex + 1.0), b.hMax);
        const double insecticideComponent = log(insecticideContent + 1.0);
        const double x = exp(b.BF
                + b.HF * holeComponent
                + b.PF * insecticideComponent
                + b.IF * holeComponent * insecticideComponent);
        const double pAtt = x / (x + 1.0);
        
        // By construction, this is non-negative:
        return pAtt * b.invBaseEffect;
    }
}

double factors::SurvivalFactor::survivalFactor( double holeIndex, double insecticideContent )const {
    if (!useLogitEqns) {
        double holeComponent = exp(-holeIndex * a.holeScaling);
        double insecticideComponent = 1.0 - exp(-insecticideContent * a.insecticideScaling);
        double killingEffect = a.BF
                + a.HF * holeComponent
                + a.PF * insecticideComponent
                + a.IF * holeComponent * insecticideComponent;
        
        double survivalFactor = (1.0 - killingEffect) * a.invBaseSurvival;
        assert( killingEffect <= 1.0 );
        // survivalFactor might be out of bounds due to precision error, see #49
        if (survivalFactor < 0.0)
            return 0.0;
        else if (survivalFactor > 1.0)
            return 1.0;
        return survivalFactor;
    } else {
        const double holeComponent = min(log(holeIndex + 1.0), b.hMax);
        const double insecticideComponent = log(insecticideContent + 1.0);
        const double x = exp(b.BF
                + b.HF * holeComponent
                + b.PF * insecticideComponent
                + b.IF * holeComponent * insecticideComponent);
        // K = x/(x+1), so (1-K) = 1/(x+1):
        const double surv = 1.0 / (x + 1.0);
        
        return surv * b.invBaseEffect;
    }
}


void factors::RelativeAttractiveness::initSingleStage(
        const scnXml::ITNDeterrency& elt, double maxInsecticide)
{
    model = SINGLE_STAGE;
    double HF = elt.getHoleFactor();
    double PF = elt.getInsecticideFactor();
    double IF = elt.getInteractionFactor();
    a.holeScaling = elt.getHoleScalingFactor();
    a.insecticideScaling = elt.getInsecticideScalingFactor();
    if( !(a.holeScaling>=0.0 && a.insecticideScaling>=0.0) ){
        throw util::xml_scenario_error("ITN.description.anophelesParams.deterrency: expected scaling factors to be non-negative");
    }
    
    /* We need to ensure the relative availability is non-negative. However,
     * since it's an exponentiated value, it always will be.
     * 
     * If don't want nets to be able to increase transmission, the following
     * limits could also be applied. In general, however, there is no reason
     * nets couldn't make individuals more attractive to mosquitoes.
     * 
     * To ensure relative availability is at most one: relative availability is
     *  exp( log(HF)*h + log(PF)*p + log(IF)*h*p )
     * where HF, PF and IF are the hole, insecticide and interaction factors
     * respectively, with h and p defined as:
     *  h=exp(-holeIndex*holeScalingFactor),
     *  p=1−exp(-insecticideContent*insecticideScalingFactor).
     * We therefore need to ensure that:
     *  log(HF)*h + log(PF)*p + log(IF)*h*p ≤ 0
     * 
     * As with the argument below concerning limits of the killing effect
     * parameters, h and p will always be in the range [0,1] and p ≤ pmax.
     * We can then derive some bounds for HF and PF:
     *  log(HF) ≤ 0
     *  log(PF)×pmax ≤ 0
     *  log(HF) + (log(PF)+log(IF))×pmax = log(HF×(PF×IF)^pmax) ≤ 0
     * or equivalently (assuming pmax>0):
     *  HF ∈ (0,1]
     *  PF ∈ (0,1]
     *  HF×(PF×IF)^pmax ∈ (0,1]
     *
     * Weaker limits would not be sufficient, as with the argument for the
     * limits of killing effect arguments below. */
    // Print out a warning if nets may increase transmission.
    double pmax = 1.0 - exp(-maxInsecticide * a.insecticideScaling);
    if( !( HF > 0.0 && PF > 0.0 && IF > 0.0 &&
            HF <= 1.0 && PF <= 1.0 && HF*pow(PF*IF,pmax) <= 1.0 ) )
    {
        cerr << "Note: since the following bounds are not met, the ITN could make humans more\n";
        cerr << "attractive to mosquitoes than they would be without a net.\n";
        cerr << "ITN.description.anophelesParams.deterrency: bounds not met:\n";
        if( !(HF>0.0) )
            cerr << "  holeFactor>0\n";
        if( !(PF>0.0) )
            cerr << "  insecticideFactor>0\n";
        if( !(IF>0.0) )
            cerr << "  interactionFactor>0\n";
        if( !(HF<=1.0) )
            cerr << "  holeFactor≤1\n";
        if( !(PF<=1.0) )
            cerr << "  insecticideFactor≤1\n";
        if( !(HF*pow(PF*IF,pmax)<=1.0) )
            cerr << "  holeFactor×(insecticideFactor×interactionFactor)^"<<pmax<<"≤1\n";
        cerr.flush();
    }
    a.lHF = log( HF );
    a.lPF = log( PF );
    a.lIF = log( IF );
}

void factors::RelativeAttractiveness::initTwoStage (
        const scnXml::TwoStageDeterrency& elt, double maxInsecticide,
        double holeIndexMax, bool holeIndexMaxPresent)
{
    b.lPFEntering = numeric_limits<double>::quiet_NaN();
    if (elt.getEntering().present()) {
        model = TWO_STAGE;
        const double PF = elt.getEntering().get().getInsecticideFactor();
        b.insecticideScalingEntering = elt.getEntering().get().getInsecticideScalingFactor();
        if( !( PF > 0.0) ){
            // we take the log of PF, so it must be positive
            ostringstream msg;
            msg << "ITN.description.anophelesParams.twoStageDeterrency.entering: insecticideFactor must be positive since we take its logarithm.";
            throw util::xml_scenario_error( msg.str() );
        }
        
        /* We need to ensure the relative availability is non-negative. However,
        * since it's an exponentiated value, it always will be.
        * 
        * If we don't want ITNs to be able to increase transmission, the following
        * limits could also be applied. In general, however, there is no reason
        * ITNs couldn't make individuals more attractive to mosquitoes.
        * 
        * To ensure relative availability is at most one: relative availability is
        *  exp( log(PF)*p ) = PF^p
        * where PF is the insecticide factor, with p∈[0,1] defined as:
        *  p=1−exp(-insecticideContent*insecticideScalingFactor).
        * We therefore just need PF ≤ 1. */
        // Print out a warning if ITNs may increase transmission.
        if( !( PF <= 1.0 ) ) {
            cerr << "Note: since the following bounds are not met, the IRS could make humans more\n";
            cerr << "attractive to mosquitoes than they would be without IRS.\n";
            cerr << "IRS.description.anophelesParams.deterrency: bounds not met:\n";
            cerr << "  0<insecticideFactor≤1\n";
            cerr.flush();
        }
        b.lPFEntering = log( PF );
    } else {
        assert( elt.getEnteringLogit().present() );
        model = TWO_STAGE_LOGIT;
        c.entBaseFactor = elt.getEnteringLogit().get().getBaseFactor();
        c.entInsecticideFactor = elt.getEnteringLogit().get().getInsecticideFactor();
        if (c.entInsecticideFactor > 0.0) {
            cerr << "ITN.description.anophelesParams.twoStageDeterrency.enteringLogit: \n"
                << "insecticideFactor should be negative to for nets to reduce chance of entering."
                << endl;
        }
        // pre-calculate for efficieny:
        c.pEnt0Inv = (exp(c.entBaseFactor) + 1.0) / exp(c.entBaseFactor);
    }
    
    // Note: b.pAttacking and c.pAttacking overlap, so we can use either:
    if (elt.getAttacking().present()) {
        b.pAttacking.init(elt.getAttacking().get(),
                          maxInsecticide,
                          "ITN.description.anophelesParams.twoStageDeterrency.attacking",
                          true);
    } else {
        assert (elt.getAttackingLogit().present());
        if (!holeIndexMax) {
            throw util::xml_scenario_error("ITN.description.anophelesParams: holeIndexMax required when using logit attacking deterrency");
        }
        b.pAttacking.initLogit(elt.getAttackingLogit().get(), holeIndexMax, true);
    }
}

double factors::RelativeAttractiveness::relativeAttractiveness (
        double holeIndex, double insecticideContent) const
{
    if (model == SINGLE_STAGE) {
        double holeComponent = exp(-holeIndex * a.holeScaling);
        double insecticideComponent = 1.0 - exp(-insecticideContent * a.insecticideScaling);
        double relAvail = exp(a.lHF * holeComponent
                + a.lPF * insecticideComponent
                + a.lIF * holeComponent * insecticideComponent );
        assert( relAvail>=0.0 );
        return relAvail;
    }
    
    double factor = numeric_limits<double>::quiet_NaN();
    if (model == TWO_STAGE) {
        // This is essentially a combination of the relative attractiveness as used
        // by IRS and a killing factor.
        
        // Note that an alternative, simpler, model could have been used, but was
        // not for consistency with other models. Alternative (here we don't take
        // the logarithm of PF):
        // pEnt = 1 - PFEntering × insecticideComponent
        
        const double insecticideComponent = 1.0 - exp(-insecticideContent * b.insecticideScalingEntering);
        const double pEnt = exp(b.lPFEntering * insecticideComponent);
        // In this model, effect with 0 insectice, pEnt0 = exp(0) = 1, hence we don't need to
        // divide by a denominator like in the logit model:
        factor = pEnt;
    } else {
        assert (model == TWO_STAGE_LOGIT);
        const double p = log(insecticideContent + 1.0);
        // We directly take the exponential (this is exp(logit.pEnt0):
        const double q = exp(c.entBaseFactor + c.entInsecticideFactor * p);
        const double pEnt = q / (q + 1.0);
        factor = pEnt * c.pEnt0Inv;
    }
    assert( factor >= 0.0 );
    
    // Note: b.pAttacking and c.pAttacking overlap, so we can use either:
    const double rel_pAtt = b.pAttacking.rel_pAtt( holeIndex, insecticideContent );
    return factor * rel_pAtt;
}


// —————  main, public classes  —————

ITNComponent::ITNComponent( ComponentId id, const scnXml::ITNDescription& elt,
        const map<string, size_t>& species_name_map ) :
        Transmission::HumanVectorInterventionComponent(id),
        ripFactor( numeric_limits<double>::signaling_NaN() )
{
    initialInsecticide.setParams( elt.getInitialInsecticide() );
    const double maxProp = 0.999;       //NOTE: this could be exposed in XML, but probably doesn't need to be
    maxInsecticide = R::qnorm5(maxProp, initialInsecticide.getMu(), initialInsecticide.getSigma(), true, false);
    holeRate = util::createSampler(elt.getHoleRate());    // per year
    holeRate->scaleMean( sim::yearsPerStep() );  // convert to per step
    ripRate = util::createSampler(elt.getRipRate());
    ripRate->scaleMean( sim::yearsPerStep() );
    ripFactor = elt.getRipFactor().getValue();
    insecticideDecay = DecayFunction::makeObject( elt.getInsecticideDecay(), "ITNDescription.insecticideDecay" );
    attritionOfNets = DecayFunction::makeObject( elt.getAttritionOfNets(), "ITNDescription.attritionOfNets" );
    // assume usage modifier is 100% if none is specified
    double propUse;
    if (elt.getUsage().present()) {
        propUse = elt.getUsage().get().getValue();
    }
    else {
        propUse = 1.0;
    }
    if( !( propUse >= 0.0 && propUse <= 1.0 ) ){
        throw util::xml_scenario_error("ITN.description.proportionUse: must be within range [0,1]");
    }
    
    typedef scnXml::ITNDescription::AnophelesParamsSequence AP;
    const AP& ap = elt.getAnophelesParams();
    species.resize(species_name_map.size());
    util::SpeciesIndexChecker checker( "ITN", species_name_map );
    for( AP::const_iterator it = ap.begin(); it != ap.end(); ++it ){
        species[checker.getIndex(it->getMosquito())].init (*it, propUse, maxInsecticide);
    }
    checker.checkNoneMissed();
    
    if( componentsByIndex.size() <= id.id ) componentsByIndex.resize( id.id+1, 0 );
    componentsByIndex[id.id] = this;
}

void ITNComponent::deploy( Host::Human& human, mon::Deploy::Method method, VaccineLimits )const{
    human.perHostTransmission.deployComponent( human.rng, *this );
    mon::reportEventMHD( mon::MHD_ITN, human, method );
}

Component::Type ITNComponent::componentType() const{
    return Component::ITN;
}

void ITNComponent::print_details( std::ostream& out )const{
    out << id().id << "\tITN";
}

unique_ptr<PerHostInterventionData> ITNComponent::makeHumanPart(LocalRng& rng) const{
    return unique_ptr<PerHostInterventionData>(new HumanITN( rng, *this ));
}
unique_ptr<PerHostInterventionData> ITNComponent::makeHumanPart( istream& stream, ComponentId id ) const{
    return unique_ptr<PerHostInterventionData>(new HumanITN( stream, id ));
}

void ITNComponent::ITNAnopheles::init(
    const scnXml::ITNDescription::AnophelesParamsType& elt,
    double proportionUse,
    double maxInsecticide)
{   bool holeIndexMaxPresent = false;

    double holeIndexMax = 0.0;
    if (elt.getHoleIndexMax().present()) {
        holeIndexMax = elt.getHoleIndexMax().get().getValue();
        holeIndexMaxPresent = true;
    }
    
    if (elt.getDeterrency().present()) {
        relAttractiveness.initSingleStage(elt.getDeterrency().get(), maxInsecticide);
    } else {
        assert (elt.getTwoStageDeterrency().present());
        relAttractiveness.initTwoStage(elt.getTwoStageDeterrency().get(), maxInsecticide, holeIndexMax, holeIndexMaxPresent);
    }
    if (elt.getPreprandialKillingEffect().present()) {
        preprandialKillingEffect.init(elt.getPreprandialKillingEffect().get(),
                                      maxInsecticide,
                                      "ITN.description.anophelesParams.preprandialKillingFactor",
                                      false);
    } else {
        assert (elt.getPreprandialKillingEffectLogit().present());
        if (!holeIndexMax) {
            throw util::xml_scenario_error("ITN.description.anophelesParams: holeIndexMax required when using logit killing effect");
        }
        preprandialKillingEffect.initLogit(elt.getPreprandialKillingEffectLogit().get(),
                                           holeIndexMax,
                                           false);
    }
    if (elt.getPostprandialKillingEffect().present()) {
        postprandialKillingEffect.init(elt.getPostprandialKillingEffect().get(),
                                       maxInsecticide,
                                       "ITN.description.anophelesParams.postprandialKillingFactor",
                                       false);
    } else {
        assert (elt.getPostprandialKillingEffectLogit().present());
        if (!holeIndexMax) {
            throw util::xml_scenario_error("ITN.description.anophelesParams: holeIndexMax required when using logit killing effect");
        }
        postprandialKillingEffect.initLogit(elt.getPostprandialKillingEffectLogit().get(),
                                            holeIndexMax,
                                            false);
    }
    if (elt.getFecundityReduction().present()) {
        relFecundityEffect.init(elt.getFecundityReduction().get(),
                                maxInsecticide,
                                "ITN.description.anophelesParams.fecundityReduction",
                                false);
    } else if (elt.getFecundityReductionLogit().present()) {
        if (!holeIndexMax) {
            throw util::xml_scenario_error("ITN.description.anophelesParams: holeIndexMax required when using logit fecundity effect");
        }
        relFecundityEffect.initLogit(elt.getFecundityReductionLogit().get(),
                                     holeIndexMax,
                                     false);
    } else {
        relFecundityEffect.init1();
    }
    // Nets only affect people while they're using the net. NOTE: we may want
    // to revise this at some point (heterogeneity, seasonal usage patterns).
    double propActive = elt.getPropActive();
    if(propActive != 1.0 && propActive != 0.0)
        throw util::xml_scenario_error("GVI.description.propActive: must be either 0 or 1");
    assert( proportionUse >= 0.0 && proportionUse <= 1.0 );
    assert( propActive >= 0.0 && propActive <= 1.0 );
    proportionProtected = proportionUse * propActive;
    proportionUnprotected = 1.0 - proportionProtected;
}

HumanITN::HumanITN( LocalRng& rng, const ITNComponent& params ) :
        PerHostInterventionData( params.id() ),
        nHoles( 0 ),
        holeIndex( 0.0 )
{
    // Net rips and insecticide loss are assumed to co-vary dependent on
    // handling of net. They are sampled once per human: human handling is
    // presumed to be the largest cause of variance.
    util::NormalSample x = util::NormalSample::generate(rng);
    holeRate = params.holeRate->sample(x);
    ripRate = params.ripRate->sample(x);
    insecticideDecayHet = params.insecticideDecay->hetSample(x);

    // Sample per-deployment variables as in redeploy:
    disposalTime = sim::now() + params.attritionOfNets->sampleAgeOfDecay(rng);
    // this is sampled independently: initial insecticide content doesn't depend on handling
    initialInsecticide = params.initialInsecticide.sample(rng);
    if( initialInsecticide < 0.0 )
        initialInsecticide = 0.0;       // avoid negative samples
    if( initialInsecticide > params.maxInsecticide )
        initialInsecticide = params.maxInsecticide;
}

void HumanITN::redeploy(LocalRng& rng, const OM::Transmission::HumanVectorInterventionComponent& params0) {
    const ITNComponent& params = *dynamic_cast<const ITNComponent*>(&params0);
    
    deployTime = sim::nowOrTs1();
    disposalTime = sim::nowOrTs1() + params.attritionOfNets->sampleAgeOfDecay(rng);
    nHoles = 0;
    holeIndex = 0.0;
    // this is sampled independently: initial insecticide content doesn't depend on handling
    initialInsecticide = params.initialInsecticide.sample(rng);
    if( initialInsecticide < 0.0 )
        initialInsecticide = 0.0;	// avoid negative samples
    if( initialInsecticide > params.maxInsecticide )
        initialInsecticide = params.maxInsecticide;
}

void HumanITN::update(Host::Human& human){
    const ITNComponent& params = *ITNComponent::componentsByIndex[m_id.id];
    if( deployTime != sim::never() ){
        // First use is at age 0 relative to ts0()
        if( sim::ts0() >= disposalTime ){
            deployTime = sim::never();
            human.removeFromCohort(id());
            return;
        }
        
        int newHoles = 0;
        if( holeRate > 0)
            newHoles = human.rng.poisson( holeRate );

        nHoles += newHoles;
        holeIndex += newHoles;

        double lambda = nHoles * ripRate;
        if (lambda > 0)
            holeIndex += params.ripFactor * human.rng.poisson( lambda );
    }
}

double HumanITN::relativeAttractiveness(size_t speciesIndex) const{
    if( deployTime == sim::never() ) return 1.0;
    const ITNComponent& params = *ITNComponent::componentsByIndex[m_id.id];
    const ITNComponent::ITNAnopheles& anoph = params.species[speciesIndex];
    return anoph.relativeAttractiveness( holeIndex, getInsecticideContent(params) );
}

double HumanITN::preprandialSurvivalFactor(size_t speciesIndex) const{
    if( deployTime == sim::never() ) return 1.0;
    const ITNComponent& params = *ITNComponent::componentsByIndex[m_id.id];
    const ITNComponent::ITNAnopheles& anoph = params.species[speciesIndex];
    return anoph.preprandialSurvivalFactor( holeIndex, getInsecticideContent(params) );
}

double HumanITN::postprandialSurvivalFactor(size_t speciesIndex) const{
    if( deployTime == sim::never() ) return 1.0;
    const ITNComponent& params = *ITNComponent::componentsByIndex[m_id.id];
    const ITNComponent::ITNAnopheles& anoph = params.species[speciesIndex];
    return anoph.postprandialSurvivalFactor( holeIndex, getInsecticideContent(params) );
}
double HumanITN::relFecundity(size_t speciesIndex) const{
    if( deployTime == sim::never() ) return 1.0;
    const ITNComponent& params = *ITNComponent::componentsByIndex[m_id.id];
    const ITNComponent::ITNAnopheles& anoph = params.species[speciesIndex];
    return anoph.relFecundity( holeIndex, getInsecticideContent(params) );
}

void HumanITN::checkpoint( ostream& stream ){
    deployTime & stream;
    disposalTime & stream;
    nHoles & stream;
    holeIndex & stream;
    initialInsecticide & stream;
    holeRate & stream;
    ripRate & stream;
    insecticideDecayHet & stream;
}
HumanITN::HumanITN( istream& stream, ComponentId id ) : PerHostInterventionData( id )
{
    deployTime & stream;
    disposalTime & stream;
    nHoles & stream;
    holeIndex & stream;
    initialInsecticide & stream;
    holeRate & stream;
    ripRate & stream;
    insecticideDecayHet & stream;
}

} }
