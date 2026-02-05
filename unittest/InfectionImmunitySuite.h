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

#ifndef Hmod_InfectionImmunitySuite
#define Hmod_InfectionImmunitySuite

#include <cxxtest/TestSuite.h>
#include "Host/WithinHost/DescriptiveWithinHost.h"
#include <limits>

using namespace OM::WithinHost;

class InfectionImmunitySuite : public CxxTest::TestSuite
{
public:
    void setUp () {
        UnittestUtil::initTime(5);
        
	// Note: these values were pulled from one source and shouldn't be taken as authoritative
	double invCumulativeYstar = 1.0 / 68564384.7102;
	double invCumulativeHstar = 1.0 / 71.676733;
	double alpha_m = 1.0 - exp(- 2.411434);
	double decayM = 2.717773;
        WHFalciparum::setParams(invCumulativeYstar, invCumulativeHstar, alpha_m, decayM);
        
        // We need a concrete class deriving from WHFalciparum; this will do
        LocalRng rng(0, 721347520444481703);
	wh = new DescriptiveWithinHostModel{ rng, numeric_limits<double>::quiet_NaN() };
    }
    void tearDown () {
	delete wh;
    }
    
    void testImmunity () {
	// Base case: virtually no immunity due to mother's immunity or past infections
        wh->m_cumulative_h = 0.0;
        wh->m_cumulative_Y = 0.0;
        double cumExposureJ = 0.0;
	TS_ASSERT_APPROX (wh->immunitySurvivalFactor (100., cumExposureJ), 1.0);
        
	// maternal immunity
	TS_ASSERT_APPROX (wh->immunitySurvivalFactor (0.1, cumExposureJ), 0.30631938551881299);
        
	// past infections, no density
        wh->m_cumulative_h = 100.0;
	TS_ASSERT_APPROX (wh->immunitySurvivalFactor (100., cumExposureJ), 0.41995608739475931);
        
	// previous with cum. density of 1e8 (but none from current infection)
        wh->m_cumulative_Y = 1e8;
	TS_ASSERT_APPROX (wh->immunitySurvivalFactor (100., cumExposureJ), 0.17081918453312689);
    }
    
private:
    WHFalciparum* wh;
};

#endif
