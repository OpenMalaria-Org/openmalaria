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

#ifndef Hmod_XoshiroSuite
#define Hmod_XoshiroSuite

#include <cxxtest/TestSuite.h>
#include "util/xoshiro.h"

class XoshiroSuite : public CxxTest::TestSuite
{
public:
    XoshiroSuite () {
    }
    
    void setUp () {
    }
    void tearDown () {
    }
    
    void testXoshiro () {
        // These values were produced with the reference implementation:
        // http://xoshiro.di.unimi.it/xoshiro256plus.c
        uint64_t vector[] = {
            5ull, 211106232532999ull, 211106635186183ull, 9223759065350669058ull,
            9250833439874351877ull, 13862484359527728515ull, 2346507365006083650ull,
            1168864526675804870ull, 34095955243042024ull, 3466914240207415127ull,
        };
        
        Xoshiro256P rng(1, 2, 3, 4);
        
        for (int n = 0; n < 10; n++) {
            uint64_t x = rng();
            TS_ASSERT_EQUALS(x, vector[n]);
        }
    }
};

#endif
