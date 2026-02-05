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

#ifndef Hmod_UtilVectorsSuite
#define Hmod_UtilVectorsSuite

#include <cxxtest/TestSuite.h>
#include "ExtraAsserts.h"

#include "util/vectors.h"

using namespace OM::util;
using OM::sim;

class UtilVectorsSuite : public CxxTest::TestSuite
{
public:
    void testInversion() {
        // Some test data. Could be any sequence.
        const double data[] = { 10, 8, 6, 6, 4, 3, 2, 1, 1, 1, 1, 1 };
        vector<double> input( data, data+6 );
//         cout << '\n' << input << endl;
        vector<double> freqDomain( input.size() * 2 - 1 );
        vectors::logDFT( input, freqDomain );
        std::vector<double> result( input.size() );
        vectors::expIDFT( result, freqDomain, 0 );
//         cout << freqDomain << endl;
//         cout << result << endl;
        ETS_ASSERT_EQUALS( input.size(), result.size() );
        for( size_t i=0; i<result.size(); ++i )
            TS_ASSERT_APPROX( input[i], result[i] );
    }
};

#endif
