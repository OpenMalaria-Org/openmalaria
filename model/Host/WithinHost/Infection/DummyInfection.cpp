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

#include "Host/WithinHost/Infection/DummyInfection.h"
#include "Host/WithinHost/CommonWithinHost.h"
#include "util/ModelOptions.h"

#include <algorithm>
#include <sstream>
#include <string.h>

namespace OM { namespace WithinHost {
    
CommonInfection* createDummyInfection (LocalRng& rng, uint32_t protID, int origin) {
    return new DummyInfection (rng, protID, origin);
}
CommonInfection* checkpointedDummyInfection (istream& stream) {
    return new DummyInfection (stream);
}

void DummyInfection::init () {
    CommonWithinHost::createInfection = &createDummyInfection;
    CommonWithinHost::checkpointedInfection = &checkpointedDummyInfection;
}

DummyInfection::DummyInfection (LocalRng& rng, uint32_t protID, int origin) :
    CommonInfection(protID, origin)
{
    m_density=16;	// increased by DH to avoid zeros in initialKappa
}

bool DummyInfection::updateDensity( LocalRng& rng, double survivalFactor, SimTime, double ){
    const double GROWTH_RATE = 8.0;
    const double PARASITE_THRESHOLD = 1;
    
    m_density = (mod_nn(int(m_density*GROWTH_RATE), 20000)) * survivalFactor;
    m_cumulativeExposureJ += m_density;
    
    if (m_density < PARASITE_THRESHOLD) {
        return true;
    } else {
        return false;
    }
}


DummyInfection::DummyInfection (istream& stream) :
    CommonInfection (stream)
{}

} }
