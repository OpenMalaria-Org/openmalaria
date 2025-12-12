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

#include <limits>
#include <memory>

#include "Global.h"
#include "util/sampler.h"
#include "util/UnitParse.h"
#include "util/DecayFunction.h"

namespace OM {
namespace util {
    
inline double readLToDays( const scnXml::DecayFunction& elt ){
    if( !elt.getL().present() ){
        throw xml_scenario_error( "decay function: attribute L required" );
    }
    return UnitParse::durationToDays(elt.getL().get(), UnitParse::YEARS);
}

class ConstantDecayFunction : public DecayFunction {
public:
    ConstantDecayFunction( const scnXml::DecayFunction& elt ) :
        DecayFunction(elt.getIncreasing(), elt.getInitialEfficacy(), elt.getCV()),
        hetFactor(0.0)
    {}
    
    double compute(double effectiveAge) const {
        // Note: we now require all decay functions to return 0 when time > 0
        // and the DecayFunction is default-constructed. So const *after deployment*.
        if( effectiveAge * hetFactor == numeric_limits<double>::infinity() )
            return 0.0;
        return 1.0;
    }
    SimTime sampleAgeOfDecay (LocalRng& rng) const {
        return sim::future();        // decay occurs "in the future" (don't use sim::never() because that is interpreted as being in the past)
    }

    unique_ptr<DecayFunction> hetSample(double hetFactor) const {
        unique_ptr<ConstantDecayFunction> copy = std::make_unique<ConstantDecayFunction>(*this);
        copy->hetFactor = hetFactor;
        return std::move(copy);
    }

private:
    double hetFactor;
};

class StepDecayFunction : public DecayFunction {
public:
    StepDecayFunction( const scnXml::DecayFunction& elt ) :
        DecayFunction(elt.getIncreasing(), elt.getInitialEfficacy(), elt.getCV()),
        invL( 1.0 / readLToDays(elt) ),
        hetFactor(0.0)
    {}
    
    double compute(double effectiveAge) const{
        if( effectiveAge * invL * hetFactor < 1.0 ){
            return 1.0;
        }else{
            return 0.0;
        }
    }
    
    SimTime sampleAgeOfDecay (LocalRng& rng) const{
        return sim::roundToTSFromDays( 1.0 / invL );
    }

    unique_ptr<DecayFunction> hetSample(double hetFactor) const {
        unique_ptr<StepDecayFunction> copy = make_unique<StepDecayFunction>(*this);
        copy->hetFactor = hetFactor;
        return std::move(copy);
    }
    
private:
    double invL, hetFactor;
};

class LinearDecayFunction : public DecayFunction {
public:
    LinearDecayFunction( const scnXml::DecayFunction& elt ) :
        DecayFunction(elt.getIncreasing(), elt.getInitialEfficacy(), elt.getCV()),
        invL( 1.0 / readLToDays(elt) ),
        hetFactor(0.0)
    {}
    
    double compute(double effectiveAge) const{
        if( effectiveAge * invL * hetFactor < 1.0 )
            return 1.0 - effectiveAge * invL * hetFactor;
        else
            return 0.0;
    }
    
    SimTime sampleAgeOfDecay (LocalRng& rng) const{
        // Note: rounds to nearest. Object may decay instantly or at time L.
        return sim::roundToTSFromDays(rng.uniform_01() / invL);
    }

    unique_ptr<DecayFunction> hetSample(double hetFactor) const {
        unique_ptr<LinearDecayFunction> copy = make_unique<LinearDecayFunction>(*this);
        copy->hetFactor = hetFactor;
        return std::move(copy);
    }
    
private:
    double invL, hetFactor;
};

class ExponentialDecayFunction : public DecayFunction {
public:
    ExponentialDecayFunction( const scnXml::DecayFunction& elt ) :
        DecayFunction(elt.getIncreasing(), elt.getInitialEfficacy(), elt.getCV()),
        invL( log(2.0) / readLToDays(elt) ),
        hetFactor(0.0)
    {}
    
    double compute(double effectiveAge) const{
        return exp( -effectiveAge * invL * hetFactor);
    }
    
    SimTime sampleAgeOfDecay (LocalRng& rng) const{
        return sim::roundToTSFromDays( -log(rng.uniform_01()) / invL );
    }

    unique_ptr<DecayFunction> hetSample(double hetFactor) const {
        unique_ptr<ExponentialDecayFunction> copy = make_unique<ExponentialDecayFunction>(*this);
        copy->hetFactor = hetFactor;
        return std::move(copy);
    }
    
private:
    double invL, hetFactor;
};

class WeibullDecayFunction : public DecayFunction {
public:
    WeibullDecayFunction( const scnXml::DecayFunction& elt ) :
        DecayFunction(elt.getIncreasing(), elt.getInitialEfficacy(), elt.getCV()),
        constOverLambda( pow(log(2.0),1.0/elt.getK()) / readLToDays(elt) ),
        k( elt.getK() ),
        hetFactor(0.0)
    {}
    
    double compute(double effectiveAge) const{
        double p = -pow(effectiveAge * constOverLambda * hetFactor, k);
        if(p < -700.0)
            return 0.0;
        else
            return exp(p);
    }
    
    SimTime sampleAgeOfDecay (LocalRng& rng) const{
        return sim::roundToTSFromDays( pow( -log(rng.uniform_01()), 1.0/k ) / constOverLambda );
    }

    unique_ptr<DecayFunction> hetSample(double hetFactor) const {
        unique_ptr<WeibullDecayFunction> copy = make_unique<WeibullDecayFunction>(*this);
        copy->hetFactor = hetFactor;
        return std::move(copy);
    }
    
private:
    double constOverLambda, k, hetFactor;
};

class HillDecayFunction : public DecayFunction {
public:
    HillDecayFunction( const scnXml::DecayFunction& elt ) :
        DecayFunction(elt.getIncreasing(), elt.getInitialEfficacy(), elt.getCV()),
        invL( 1.0 / readLToDays(elt) ),
        k( elt.getK() ),
        hetFactor(0.0)
    {}
    
    double compute(double effectiveAge) const{
        return 1.0 / (1.0 + pow(effectiveAge * invL * hetFactor, k));
    }
    
    SimTime sampleAgeOfDecay (LocalRng& rng) const{
        return sim::roundToTSFromDays( pow( 1.0 / rng.uniform_01() - 1.0, 1.0/k ) / invL );
    }

    unique_ptr<DecayFunction> hetSample(double hetFactor) const {
        unique_ptr<HillDecayFunction> copy = make_unique<HillDecayFunction>(*this);
        copy->hetFactor = hetFactor;
        return std::move(copy);
    }
    
private:
    double invL, k, hetFactor;
};

class SmoothCompactDecayFunction : public DecayFunction {
public:
    SmoothCompactDecayFunction( const scnXml::DecayFunction& elt ) :
        DecayFunction(elt.getIncreasing(), elt.getInitialEfficacy(), elt.getCV()),
        invL( 1.0 / readLToDays(elt) ),
        k( elt.getK() ),
        hetFactor(0.0)
    {}

    double compute(double effectiveAge) const{
        if( effectiveAge * invL * hetFactor < 1.0 )
            return exp( k - k / (1.0 - pow(effectiveAge * invL * hetFactor, 2.0)) );
        else
            return 0.0;
    }
    
    SimTime sampleAgeOfDecay (LocalRng& rng) const{
        return sim::roundToTSFromDays( sqrt( 1.0 - k / (k - log( rng.uniform_01() )) ) / invL );
    }

    unique_ptr<DecayFunction> hetSample(double hetFactor) const {
        unique_ptr<SmoothCompactDecayFunction> copy = make_unique<SmoothCompactDecayFunction>(*this);
        copy->hetFactor = hetFactor;
        return std::move(copy);
    }
    
private:
    double invL, k, hetFactor;
};

template<class T>
class OperatorDecayFunction : public DecayFunction {
public:
    OperatorDecayFunction( const scnXml::DecayFunction& elt ) : 
        DecayFunction(elt.getIncreasing(), elt.getInitialEfficacy(), elt.getCV()) {
        const scnXml::DecayFunction::DecaySequence &decaySequence = elt.getDecay();
        if(decaySequence.size() != 2)
            throw xml_scenario_error("Operator decay function expects two decay functions, " + to_string(decaySequence.size()) +"  were given.");

        f1 = DecayFunction::makeObject(decaySequence[0], "Operator::f1");
        f2 = DecayFunction::makeObject(decaySequence[1], "Operator::f2");
    }

    OperatorDecayFunction(const OperatorDecayFunction &copy, unique_ptr<DecayFunction> f1, unique_ptr<DecayFunction> f2) : 
        DecayFunction(copy),
        f1(std::move(f1)), f2(std::move(f2)) {}

    double compute(double effectiveAge) const {
        return max(min(op(f1->eval(effectiveAge), f2->eval(effectiveAge)), 1.0), 0.0);
    }
    
    SimTime sampleAgeOfDecay (LocalRng& rng) const {
        return sim::roundToTSFromDays( max(f1->sampleAgeOfDecay(rng), f2->sampleAgeOfDecay(rng)) );
    }

    unique_ptr<DecayFunction> hetSample(double hetFactor) const {
        unique_ptr<DecayFunction> f1hetSample = f1->hetSample(hetFactor);
        unique_ptr<DecayFunction> f2hetSample = f2->hetSample(hetFactor);
        unique_ptr<OperatorDecayFunction> copy = make_unique<OperatorDecayFunction>(*this, std::move(f1hetSample), std::move(f2hetSample));
        return std::move(copy);
    }

private:
    unique_ptr<DecayFunction> f1, f2;
    T op;
};

class EmaxFunction : public DecayFunction {
public:
    EmaxFunction( const scnXml::DecayFunction& elt ) : 
        DecayFunction(elt.getIncreasing(), elt.getInitialEfficacy(), elt.getCV()) {
        const scnXml::DecayFunction::DecaySequence &decaySequence = elt.getDecay();
        if(decaySequence.size() != 1)
            throw util::xml_scenario_error("The Emax function expects exactly one user function, but " + to_string(decaySequence.size()) +"  were given.");
        if(!elt.getEmax().present())
            throw util::xml_scenario_error("Emax function: the Emax parameter is not declared");
        if(!elt.getIC50().present())
            throw util::xml_scenario_error("Emax function: the IC50 parameter is not declared");
        if(!elt.getSlope().present())
            throw util::xml_scenario_error("Emax function: the Slope parameter is not declared");
        if(!elt.getInitialConcentration().present())
            throw util::xml_scenario_error("Emax function: the InitialConcentration parameter is not declared");
        
        f = makeObject(decaySequence[0], "Emax::f");

        Emax = elt.getEmax().get();
        IC50 = elt.getIC50().get();
        slope = elt.getSlope().get();
        initialConcentration = elt.getInitialConcentration().get();

        IC50 = pow(IC50, slope); // compute once
    }

    EmaxFunction(const EmaxFunction &copy, unique_ptr<DecayFunction> f) : 
        DecayFunction(copy),
        f(move(f)),
        Emax(copy.Emax), 
        IC50(copy.IC50),
        slope(copy.slope),
        initialConcentration(copy.initialConcentration) {}

    double compute(double effectiveAge) const {
        double fcomp = pow(initialConcentration * f->eval(effectiveAge), slope);
        // cout << "t: " << effectiveAge << endl;
        // cout << "AB(t): " << f->eval(effectiveAge) << endl;
        // cout << "rho: " << initialConcentration << endl;
        // cout << "slope: " << slope << endl;
        // cout << "IC50^slope: " << IC50 << endl;
        // cout << "rho * AB(t)^slope: " << fcomp << endl;
        // cout << "Emax(t) = " << Emax * fcomp / (fcomp + IC50) << endl;
        // cout << "Emax(t) = " << Emax * fcomp / (fcomp + IC50) << "\t" << max(min(Emax * fcomp / (fcomp + IC50), 1.0), 0.0) << endl;
        return max(min(Emax * fcomp / (fcomp + IC50), 1.0), 0.0);
    }
    
    SimTime sampleAgeOfDecay (LocalRng& rng) const {
        return sim::roundToTSFromDays( f->sampleAgeOfDecay(rng) );
    }

    unique_ptr<DecayFunction> hetSample(double hetFactor) const {
        unique_ptr<DecayFunction> fhetSample = f->hetSample(hetFactor);
        unique_ptr<EmaxFunction> copy = make_unique<EmaxFunction>(*this, move(fhetSample));
        return move(copy);
    }

private:
    unique_ptr<DecayFunction> f;
    double Emax, IC50, slope, initialConcentration;
};

// -----  interface / static functions  -----
unique_ptr<DecayFunction> DecayFunction::makeObject(
    const scnXml::DecayFunction& elt, const char* eltName
){
    // Type mostly equivalent to a std::string:
    const scnXml::Function& func = elt.getFunction();
    if( func == "constant" )
        return make_unique<ConstantDecayFunction>( elt );
    else if( func == "step" )
        return make_unique<StepDecayFunction>( elt );
    else if( func == "linear" )
        return make_unique<LinearDecayFunction>( elt );
    else if( func == "exponential" )
        return make_unique<ExponentialDecayFunction>( elt );
    else if( func == "weibull" )
        return make_unique<WeibullDecayFunction>( elt );
    else if( func == "hill" )
        return make_unique<HillDecayFunction>( elt );
    else if( func == "smooth-compact" )
        return make_unique<SmoothCompactDecayFunction>( elt );
    else if( func == "plus" )
        return make_unique<OperatorDecayFunction<std::plus<double>>>( elt );
    else if( func == "minus" )
        return make_unique<OperatorDecayFunction<std::minus<double>>>( elt );
    else if( func == "divides" )
        return make_unique<OperatorDecayFunction<std::divides<double>>>( elt );
    else if( func == "multiplies" )
        return make_unique<OperatorDecayFunction<std::multiplies<double>>>( elt );
    else if( func == "emax" )
        return make_unique<EmaxFunction>( elt );
    else
        throw xml_scenario_error("decay function " + string(eltName) + " of type " + string(func) + " unrecognized");
}

} }
