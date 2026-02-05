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

#ifndef Hmod_DecayFunction
#define Hmod_DecayFunction

#include <limits>
#include <memory>

#include "Global.h"
#include "util/sampler.h"
#include "util/UnitParse.h"

namespace OM {
namespace util {

/** An interface for a few types of decay function (some of which may also be
 * suitible survival functions).
 *
 * Heterogeneity is implemented by passing a DecayFunction object to the eval
 * function; this should be sampled by the same DecayFunction.
 *****************************************************************************/
class DecayFunction
{
public:
    DecayFunction(bool increasing = false, double initialEfficacy = 1.0, double CV = 0.0) : 
        increasing(increasing), initialEfficacy(initialEfficacy), het(std::move(LognormalSampler::fromMeanCV(1.0, CV))) {}

    DecayFunction(const DecayFunction& other) : 
        increasing(other.increasing),
        initialEfficacy(other.initialEfficacy),
        het(other.het ? other.het->clone() : nullptr) {}

    virtual ~DecayFunction() {}

    /** Return a new decay function, constructed from an XML element.
     * 
     * @param elt XML element specifying which function to use and parameters
     * @param eltName Name of XML element (for reasonable error reporting)
     */
    static unique_ptr<DecayFunction> makeObject(
        const scnXml::DecayFunction& elt, const char* eltName
    );

    /** Sample a DecayFunction value (should be stored per individual).
     * 
     * Note that a DecayFunction is needed to call eval() even if heterogeneity
     * is not wanted. If sigma = 0 then the random number stream will not be
     * touched. */
    unique_ptr<DecayFunction> hetSample(LocalRng& rng) const {
        return hetSample(het->sample(rng));
    }
    
    /** Generate a DecayFunction value from an existing sample. */
    std::unique_ptr<DecayFunction> hetSample(NormalSample sample) const {
        return hetSample(het->sample(sample));
    }
    
    /** Generate a DecayFunction value from an existing sample. */
    virtual unique_ptr<DecayFunction> hetSample(double hetFactor) const =0;

    /** Say you have a population of objects which each have two states:
     * decayed and not decayed. If you want to use a DecayFunction to model
     * the proportion of objects which have decayed, you need to work out per
     * object the age of decay. This function does that.
     * 
     * This is only valid where mu and sigma parameters are zero.
     * 
     * @returns Age at which an object should decay. */
    virtual SimTime sampleAgeOfDecay (LocalRng& rng) const =0;
    
    inline double eval(double ageDays) const {
        if(increasing)
            return 1.0 - compute( ageDays ) * initialEfficacy;
        else
            return compute( ageDays ) * initialEfficacy;
    }

    virtual double compute(double ageDays) const = 0;

    /// Checkpointing
    template<class S>
    void operator& (S& stream) {
        // timeFactorHet & stream;
    }

private:
    bool increasing;
    double initialEfficacy;
    std::unique_ptr<util::Sampler> het;
};

} }
#endif
