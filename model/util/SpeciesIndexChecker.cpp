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

#include "util/SpeciesIndexChecker.h"
#include "util/errors.h"
#include <sstream>

namespace OM { namespace util {
using std::ostringstream;

size_t SpeciesIndexChecker::getIndex( string species ){
    if( found.count( species ) > 0 ){
        ostringstream msg;
        msg << "Intervention \"" << _intervName
            << "\" has multiple descriptions for vector species \""
            << species << "\"";
        throw util::xml_scenario_error( msg.str() );
    }
    auto sIndex = _indices.find( species );
    if( sIndex == _indices.end() ){
        ostringstream msg;
        msg << "Intervention \"" << _intervName
            << "\" has a description for vector species \"" << species
            << "\", but this species is not mentioned in the entomology section";
        throw util::xml_scenario_error( msg.str() );
    }
    found.insert( species );
    return sIndex->second;
}

void SpeciesIndexChecker::checkNoneMissed() const{
    for( auto it = _indices.begin(); it != _indices.end(); ++it ){
        if( found.count( it->first ) == 0 ){
            ostringstream msg;
            msg << "Intervention \"" << _intervName
                << "\" has a no description for vector species \""
                << it->first << "\"";
            throw util::xml_scenario_error( msg.str() );
        }
    }
}

 } }
 
