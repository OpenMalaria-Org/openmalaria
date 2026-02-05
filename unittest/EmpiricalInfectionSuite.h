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

#ifndef Hmod_EmpiricalInfectionSuite
#define Hmod_EmpiricalInfectionSuite

#include <cxxtest/TestSuite.h>
#include "UnittestUtil.h"
#include "ExtraAsserts.h"
#include "Host/WithinHost/Infection/EmpiricalInfection.h"
#include "Host/WithinHost/CommonWithinHost.h"
#include "util/random.h"
#include <limits>

using namespace OM::WithinHost;

const bool dump_emprical = false;

class EmpiricalInfectionSuite : public CxxTest::TestSuite
{
public:
    EmpiricalInfectionSuite() : m_rng(0, 0) {}
    
    void setUp () {
        m_rng.seed(3978236241, 721347520444481703);
        UnittestUtil::initTime(1);
        UnittestUtil::Infection_init_latentP_and_NaN ();
        EmpiricalInfection::init();
        // pkpdID (1st value) isn't important since we're not using drug model here:
        infection = CommonWithinHost::createInfection( m_rng, 0xFFFFFFFF, InfectionOrigin::Indigenous);
        for( SimTime d = sim::ts1(), end = sim::ts1() + sim::fromDays(15); d < end; d = d + sim::oneDay() ){
            // blood stage starts 15 days after creation
            UnittestUtil::incrTime( sim::oneDay() );
            infection->update( m_rng, 1.0, d, numeric_limits<double>::quiet_NaN() );
        }
    }
    void tearDown () {
        delete infection;
    }

    void testNewInf () {
        TS_ASSERT_APPROX (infection->getDensity(), 0.00000000000000000);
    }

    // Parasite growth is stochastic, so there's not a lot we can test, except for reproducability
    void testUpdatedInf () {
        UnittestUtil::incrTime( sim::oneTS() );
        infection->update (m_rng, 1.0, sim::ts1(), numeric_limits<double>::quiet_NaN());
        if (dump_emprical) cout << setprecision(8) << infection->getDensity() << endl;
        TS_ASSERT_APPROX (infection->getDensity(), 2.3427793);
    }
    void testUpdated2Inf () {
        UnittestUtil::incrTime( sim::oneTS() );
        infection->update (m_rng, 1.0, sim::ts1(), numeric_limits<double>::quiet_NaN());
        UnittestUtil::incrTime( sim::oneTS() );
        infection->update (m_rng, 1.0, sim::ts1(), numeric_limits<double>::quiet_NaN());
        if (dump_emprical) cout << setprecision(8) << infection->getDensity() << endl;
        TS_ASSERT_APPROX (infection->getDensity(), 0.42220637);
    }
    void testUpdated3Inf () {
        UnittestUtil::incrTime( sim::oneTS() );
        infection->update (m_rng, 1.0, sim::ts1(), numeric_limits<double>::quiet_NaN());
        UnittestUtil::incrTime( sim::oneTS() );
        infection->update (m_rng, 1.0, sim::ts1(), numeric_limits<double>::quiet_NaN());
        UnittestUtil::incrTime( sim::oneTS() );
        infection->update (m_rng, 1.0, sim::ts1(), numeric_limits<double>::quiet_NaN());
        if (dump_emprical) cout << setprecision(8) << infection->getDensity() << endl;
        TS_ASSERT_APPROX (infection->getDensity(), 1.2290052);
    }
    void testUpdated4Inf () {
        UnittestUtil::incrTime( sim::oneTS() );
        infection->update (m_rng, 1.0, sim::ts1(), numeric_limits<double>::quiet_NaN());
        UnittestUtil::incrTime( sim::oneTS() );
        infection->update (m_rng, 1.0, sim::ts1(), numeric_limits<double>::quiet_NaN());
        UnittestUtil::incrTime( sim::oneTS() );
        infection->update (m_rng, 1.0, sim::ts1(), numeric_limits<double>::quiet_NaN());
        UnittestUtil::incrTime( sim::oneTS() );
        infection->update (m_rng, 1.0, sim::ts1(), numeric_limits<double>::quiet_NaN());
        if (dump_emprical) cout << setprecision(8) << infection->getDensity() << endl;
        TS_ASSERT_APPROX (infection->getDensity(), 0.4272347);
    }
    void testUpdatedInf1 () {
        UnittestUtil::incrTime( sim::oneTS() );
        infection->update (m_rng, 1.0, sim::ts1(), numeric_limits<double>::quiet_NaN());
        if (dump_emprical) cout << setprecision(8) << infection->getDensity() << endl;
        TS_ASSERT_APPROX (infection->getDensity(), 2.3427793);
    }

    void testUpdatedReducedInf () {
        UnittestUtil::incrTime( sim::oneTS() );
        infection->update (m_rng, 1.0, sim::ts1(), numeric_limits<double>::quiet_NaN());
        UnittestUtil::incrTime( sim::oneTS() );
        infection->update (m_rng, 0.1, sim::ts1(), numeric_limits<double>::quiet_NaN());
        // This is, as expected, 1/10th of that in testUpdated2Inf
        if (dump_emprical) cout << setprecision(8) << infection->getDensity() << endl;
        TS_ASSERT_APPROX (infection->getDensity(), 0.042220637);
    }
    void testUpdatedReducedInf2 () {
        UnittestUtil::incrTime( sim::oneTS() );
        infection->update (m_rng, 0.1, sim::ts1(), numeric_limits<double>::quiet_NaN());
        UnittestUtil::incrTime( sim::oneTS() );
        infection->update (m_rng, 1.0, sim::ts1(), numeric_limits<double>::quiet_NaN());
        // This is completely different due to stochasitic effects
        if (dump_emprical) cout << setprecision(8) << infection->getDensity() << endl;
        TS_ASSERT_APPROX (infection->getDensity(), 0.16787395);
    }

private:
    CommonInfection* infection;
    LocalRng m_rng;
};

#endif
