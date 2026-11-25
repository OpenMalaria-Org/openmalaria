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

#ifndef Hmod_SamplerSuite
#define Hmod_SamplerSuite

#include <cxxtest/TestSuite.h>
#include "util/sampler.h"

class SamplerSuite : public CxxTest::TestSuite {
public:
    void testLognormalVarianceMatchesCV() {
        for (double mean : {0.5, 1.0, 3.0}) {
            for (double variance : {0.0, 0.25, 1.0, 4.0}) {
                auto fromVariance = OM::util::LognormalSampler::fromMeanVariance(mean, variance);
                auto fromCV = OM::util::LognormalSampler::fromMeanCV(mean, std::sqrt(variance) / mean);
                TS_ASSERT_DELTA(fromVariance->mean(), mean, 1e-12);
                for (double value : {0.1, 0.5, 1.0, 3.0, 10.0}) {
                    TS_ASSERT_DELTA(fromVariance->cdf(value), fromCV->cdf(value), 1e-12);
                }
            }
        }
    }

    void testLognormalVarianceFromXml() {
        scnXml::SampledValueLN description(1.0);
        description.setDistr("lognormal");
        description.setVariance(0.25);
        auto sampler = OM::util::createSampler(description);

        // A lognormal with arithmetic mean 1 and variance 0.25 has this median.
        const double median = 1.0 / std::sqrt(1.25);
        TS_ASSERT_DELTA(sampler->cdf(median), 0.5, 1e-12);
        TS_ASSERT_DELTA(sampler->mean(), 1.0, 1e-12);
    }
};

#endif
