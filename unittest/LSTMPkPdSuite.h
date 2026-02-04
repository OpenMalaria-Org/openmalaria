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

// Unittest for the LSTM drug model

#ifndef Hmod_LSTMPkPdSuite
#define Hmod_LSTMPkPdSuite

#include <cxxtest/TestSuite.h>
#include "PkPd/LSTMModel.h"
#include "Host/WithinHost/Infection/DummyInfection.h"
#include "UnittestUtil.h"
#include "ExtraAsserts.h"
#include <limits>

using namespace OM;
using namespace OM::PkPd;

/// There is probably little value in this unit-test now that PkPdComplianceSuite exists
class LSTMPkPdSuite : public CxxTest::TestSuite
{
public:
    LSTMPkPdSuite() :
            m_rng(0, 0), proxy(0)
    {
        // This is what a previously-used weight distribution gave us,
        // and is good enough for testing purposes:
        massAt21 = 55.4993;
    }
    
    void setUp () {
        m_rng.seed(0, 721347520444481703);
        UnittestUtil::initTime(1);
	UnittestUtil::PkPdSuiteSetup();
	proxy = new LSTMModel ();
        inf = createDummyInfection(m_rng, 0, InfectionOrigin::Indigenous);
        MQ_index = LSTMDrugType::findDrug( "MQ" );
    }
    void tearDown () {
        delete inf;
	delete proxy;
        LSTMDrugType::clear();
    }
    
    void testNone () {
	TS_ASSERT_EQUALS (proxy->getDrugFactor (m_rng, inf, massAt21), 1.0);
    }
    
    void testOral () {
	UnittestUtil::medicate( m_rng, *proxy, MQ_index, 3000, 0 );
	TS_ASSERT_APPROX (proxy->getDrugFactor (m_rng, inf, massAt21), 0.03174563638523168);
    }
    
    void testOralHalves () {	// the point being: check it can handle two doses at the same time-point correctly
	UnittestUtil::medicate( m_rng, *proxy, MQ_index, 1500, 0 );
	UnittestUtil::medicate( m_rng, *proxy, MQ_index, 1500, 0 );
	TS_ASSERT_APPROX (proxy->getDrugFactor (m_rng, inf, massAt21), 0.03174563638523168);
    }
    
    void testOralSplit () {
	UnittestUtil::medicate( m_rng, *proxy, MQ_index, 3000, 0 );
	UnittestUtil::medicate( m_rng, *proxy, MQ_index, 0, 0.5 );	// insert a second dose half way through the day: forces drug calculation to be split into half-days but shouldn't affect result
	TS_ASSERT_APPROX (proxy->getDrugFactor (m_rng, inf, massAt21), 0.03174563639140275);
    }
    
    void testOralDecayed () {
	UnittestUtil::medicate( m_rng, *proxy, MQ_index, 3000, 0 );
	proxy->decayDrugs (massAt21);
	TS_ASSERT_APPROX (proxy->getDrugFactor (m_rng, inf, massAt21), 0.03174563639501896);
    }
    
    void testOral2Doses () {
	UnittestUtil::medicate( m_rng, *proxy, MQ_index, 3000, 0 );
	proxy->decayDrugs (massAt21);
	UnittestUtil::medicate( m_rng, *proxy, MQ_index, 3000, 0 );
	TS_ASSERT_APPROX (proxy->getDrugFactor (m_rng, inf, massAt21), 0.03174563637686205);
    }
    
private:
    LocalRng m_rng;
    LSTMModel *proxy;
    CommonInfection *inf;
    double massAt21;
    size_t MQ_index;
};

#endif
