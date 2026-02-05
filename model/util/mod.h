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

#ifndef Hmod_util_mod
#define Hmod_util_mod

#include <cassert>

namespace OM {
namespace util {

// We want a real modular arithmatic operator (one supporting
// "a + n mod n = a mod n" even where a is negative).
// General form as follows (don't use templates because implicit casts don't
// work with them).
inline int mod (int a, int b) {
    assert (b>0);
    // http://stackoverflow.com/questions/12089514/real-modulo-operator-in-c-c
    const int result = a % b;
    return result >= 0 ? result : result + b;
}
// Optimisation valid when we know a is non-negative.
// Tested; with NamawalaArabiensis scenario this yields a 0.1% improvement in
// release mode. In RelWithRelease (-O2 -g) mode it's 0.03% slower on this
// scenario and 0.015% slower on the Penny scenario (less usage).
// Not highly significant, but the optimisation is implemented so it can stay.
inline int mod_nn (int a, int b) {
    assert (a>=0 && b>0);
    return a % b;
}

}
}
#endif
