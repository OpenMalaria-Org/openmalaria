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

#include "Host/WithinHost/Treatments.h"
#include "util/errors.h"
#include "util/UnitParse.h"
#include "schema/healthSystem.h"

namespace OM {
namespace WithinHost {

// ———  static  ———

vector<Treatments> Treatments::treatments;

TreatmentId Treatments::addTreatment( const scnXml::TreatmentOption& desc ){
    Treatments treatment( desc );
    TreatmentId id( treatments.size() );
    treatments.push_back( std::move(treatment) );
    return id;
}


// ———   non-static  ———

Treatments::Stages stageFromString( const std::string& str ){
    if( str == "liver" ) return Treatments::LIVER;
    if( str == "blood" ) return Treatments::BLOOD;
    if( str == "both" ) return Treatments::BOTH;
    throw util::format_error( std::string("stage must be liver, blood or both, not ").append(str) );
}
Treatments::Treatments( const scnXml::TreatmentOption& elt ) :
    TriggeredDeployments(elt), timeLiver(sim::zero()), timeBlood(sim::zero())
{
    for( auto it = elt.getClearInfections().begin(), end = elt.getClearInfections().end();
        it != end; ++it )
    {
        try{
            //NOTE: if changing XSD, this should not be called "timesteps" or have a default unit
            SimTime len = UnitParse::readShortDuration( it->getTimesteps(), UnitParse::STEPS );
            if( len < -sim::oneTS() || len == sim::zero() ){ 
                throw util::format_error( "timesteps must be ≥ 1 or have special value -1" );
            }
            Stages stage = stageFromString( it->getStage() );
            if( stage & LIVER ){
                if( timeLiver != sim::zero() )   // existing treatment configuration
                    throw util::format_error( "multiple specification of liver stage effect" );
                timeLiver = len;
            }
            if( stage & BLOOD ){
                if( timeBlood != sim::zero() )   // existing treatment configuration
                    throw util::format_error( "multiple specification of blood stage effect" );
                timeBlood = len;
            }
        }catch( const util::format_error& e ){
            throw util::xml_scenario_error( string(".../clearInfections/timesteps: ").append(e.message()) );
        }
    }
}

}
}
