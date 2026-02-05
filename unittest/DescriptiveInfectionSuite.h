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

#ifndef Hmod_DescriptiveInfectionSuite
#define Hmod_DescriptiveInfectionSuite

#include <cxxtest/TestSuite.h>
#include "Host/WithinHost/Infection/DescriptiveInfection.h"
#include "util/random.h"
#include "UnittestUtil.h"
#include <limits>

using namespace OM::WithinHost;

class DescriptiveInfectionSuite : public CxxTest::TestSuite
{
public:
    void setUp () {
        UnittestUtil::initTime( 5 );
	UnittestUtil::Infection_init ();
	UnittestUtil::DescriptiveInfection_init ();
	DescriptiveInfection::initParameters();	//FIXME: this gets data from InputData
	gsl::setUp (83);	// seed is unimportant, but must be fixed
	infection = new DescriptiveInfection ();
    }
    void tearDown () {
	delete infection;
	gsl::tearDown ();
    }
    
    // Should test DescriptiveInfection::determineDensities (and maybe others), but due to large number of inputs & 2 outputs isn't so easy.
    
private:
    DescriptiveInfection* infection;
};

#endif
