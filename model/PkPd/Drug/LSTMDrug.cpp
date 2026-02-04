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

#include "PkPd/Drug/LSTMDrug.h"
#include "util/errors.h"
#include "util/StreamValidator.h"
#include "util/vectors.h"

#include <assert.h>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <sstream>
#include <algorithm>

using namespace std;

namespace OM {
namespace PkPd {

LSTMDrug::LSTMDrug(double Vd): vol_dist(Vd) {}
LSTMDrug::~LSTMDrug() {}

bool comp(pair<double,double> lhs, pair<double,double> rhs){
    return lhs.first < rhs.first;
}

// Create our list of doses. Optimise for the case where we only have 1 or less
// per day, but be able to handle more.
// If two doses coincide, they can be combined but doing so is not essential.

void LSTMDrug::medicate(double time, double qty){
    // Insert in the right position to maintain sorting:
    auto elt = make_pair (time, qty);
    auto pos = lower_bound(doses.begin(), doses.end(), elt, comp);
    doses.insert(pos, std::move(elt));
    assert(is_sorted(doses.begin(), doses.end(), comp));
}

}
}
