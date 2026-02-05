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

#ifndef OM_UTIL_SPECIES_INDEX_CHECKER
#define OM_UTIL_SPECIES_INDEX_CHECKER

#include "Global.h"

#include <set>

namespace OM { namespace util {
using std::set;

/** Checks keys exist and that all are used. Fails with an xml_scenario_error
 * and a suitable error message. */
class SpeciesIndexChecker{
public:
    SpeciesIndexChecker( const string& intervName, const map<string,size_t>& indices ) :
        _intervName( intervName ), _indices( indices ){}
    
    /** Return the index in speciesIndex of mosquito, throwing if not found. */
    size_t getIndex( string species );
    
    /** Throw if some species was missed. */
    void checkNoneMissed() const;
    
private:
    const string& _intervName;
    const map<string,size_t>& _indices;
    set<string> found;
};

} }
#endif
