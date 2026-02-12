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

#ifndef Hmod_DecayFunctionSuite
#define Hmod_DecayFunctionSuite

#include <cxxtest/TestSuite.h>
#include "UnittestUtil.h"
#include "ExtraAsserts.h"
#include "util/DecayFunction.h"

using ::OM::util::DecayFunction;

class DecayFunctionSuite : public CxxTest::TestSuite
{
public:
    DecayFunctionSuite () :
        m_rng(0, 0), dfElt( "" )
    {
        dfElt.setL( "10y" );
        dfElt.setK( 1.6 );
    }
    
    void setUp() {
        m_rng.seed(0, 721347520444481703);
        UnittestUtil::initTime(5);
    }
    
    void testBad() {
        dfElt.setFunction( "unknown" );
        string errMsg = "decay function DecayFunctionSuite of type unknown unrecognized";
        TS_ASSERT_THROWS_EQUALS( df = DecayFunction::makeObject( dfElt, "DecayFunctionSuite" ),
                          const std::runtime_error &e, e.what(), errMsg );
    }
    
    void testConstant () {
        dfElt.setFunction( "constant" );
        df = DecayFunction::makeObject( dfElt, "DecayFunctionSuite" );
        // Second test: with an appropriately sampled helper value, we should get the results we want
        unique_ptr<DecayFunction> dHet = df->hetSample(m_rng);
        TS_ASSERT_APPROX( dHet->eval( sim::zero() ) , 1.0 );
        TS_ASSERT_APPROX( dHet->eval( sim::fromYearsI(10 ) ), 1.0 );
        // Third test: time of decay (age plus now) should always be in the future
        TS_ASSERT_EQUALS( df->sampleAgeOfDecay(m_rng), sim::future() );
    }
    
    void testStep () {
        dfElt.setFunction( "step" );
        df = DecayFunction::makeObject( dfElt, "DecayFunctionSuite" );
        unique_ptr<DecayFunction> dHet = df->hetSample(m_rng);
        TS_ASSERT_APPROX( dHet->eval( sim::zero() ), 1.0 );
        TS_ASSERT_APPROX( dHet->eval( sim::fromYearsI(6 ) ), 1.0 );
        TS_ASSERT_APPROX( dHet->eval( sim::fromYearsI(20 ) ), 0.0 );
        TS_ASSERT_EQUALS( df->sampleAgeOfDecay(m_rng), sim::fromYearsI(10) );
    }
    
    void testLinear () {
        dfElt.setFunction( "linear" );
        df = DecayFunction::makeObject( dfElt, "DecayFunctionSuite" );
        unique_ptr<DecayFunction> dHet = df->hetSample(m_rng);
        TS_ASSERT_APPROX( dHet->eval( sim::zero() ), 1.0 );
        TS_ASSERT_APPROX( dHet->eval( sim::fromYearsI(6 ) ), 0.4 );
        TS_ASSERT_APPROX( dHet->eval( sim::fromYearsI(20 ) ), 0.0 );
    }
    
    void testExponential () {
        dfElt.setFunction( "exponential" );
        df = DecayFunction::makeObject( dfElt, "DecayFunctionSuite" );
        unique_ptr<DecayFunction> dHet = df->hetSample(m_rng);
        TS_ASSERT_APPROX( dHet->eval( sim::zero() ), 1.0 );
        TS_ASSERT_APPROX( dHet->eval( sim::fromYearsI(6 ) ), 0.65975394736842108 );
        TS_ASSERT_APPROX( dHet->eval( sim::fromYearsI(20 ) ), 0.25 );
    }
    
    void testWeibull () {
        dfElt.setFunction( "weibull" );
        df = DecayFunction::makeObject( dfElt, "DecayFunctionSuite" );
        unique_ptr<DecayFunction> dHet = df->hetSample(m_rng);
        TS_ASSERT_APPROX( dHet->eval( sim::zero() ), 1.0 );
        TS_ASSERT_APPROX( dHet->eval( sim::fromYearsI(6 ) ), 0.73631084210526321 );
        TS_ASSERT_APPROX( dHet->eval( sim::fromYearsI(20 ) ), 0.122306 );
    }
    
    void testHill () {
        dfElt.setFunction( "hill" );
        df = DecayFunction::makeObject( dfElt, "DecayFunctionSuite" );
        unique_ptr<DecayFunction> dHet = df->hetSample(m_rng);
        TS_ASSERT_APPROX( dHet->eval( sim::zero() ), 1.0 );
        TS_ASSERT_APPROX( dHet->eval( sim::fromYearsI(6 ) ), 0.6936673684210527 );
        TS_ASSERT_APPROX( dHet->eval( sim::fromYearsI(20 ) ), 0.24805074736842106 );
    }
    
    void testSmoothCompact () {
        dfElt.setFunction( "smooth-compact" );
        df = DecayFunction::makeObject( dfElt, "DecayFunctionSuite" );
        unique_ptr<DecayFunction> dHet = df->hetSample(m_rng);
        TS_ASSERT_APPROX( dHet->eval( sim::zero() ), 1.0 );
        TS_ASSERT_APPROX( dHet->eval( sim::fromYearsI(6 ) ), 0.40656965789473687 );
        TS_ASSERT_APPROX( dHet->eval( sim::fromYearsI(20 ) ), 0.0 );
    }
    
private:
    LocalRng m_rng;
    scnXml::DecayFunction dfElt;
    unique_ptr<DecayFunction> df;
};

#endif
