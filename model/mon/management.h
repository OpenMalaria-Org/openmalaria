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

#ifndef H_OM_mon_management
#define H_OM_mon_management

#include <fstream>

namespace scnXml{
    class Scenario;
    class Monitoring;
}

/** This header manages monitoring: it reads configuration and writes output.
 *
 * It does not store reported data (directly) and does not handle reports. */
namespace OM {
    class Parameters;
namespace mon {

/// Read survey times from XML. Return the date of the final survey.
SimTime readSurveyDates( const scnXml::Monitoring& monitoring );

/// Call before start of simulation to set up outputs. Call readSurveyDates first.
void initReporting( const scnXml::Scenario& scenario );

/// Call after initialising interventions
void initCohorts( const scnXml::Monitoring& monitoring );

/// Call just before the start of the intervention period
void initMainSim();

/// Call after all data for some survey number has been provided
void concludeSurvey();

/// Write survey data to output.txt (or configured file)
void writeSurveyData();

// Checkpointing
void checkpoint( std::ostream& stream );
void checkpoint( std::istream& stream );

// Functions for internal use (within mon package)
namespace internal{
    // Write results to stream
    void write( std::ostream& stream );
    
    /** Get the output cohort set numeric identifier given the internal one
     * (as returned by Survey::updateCohortSet()). */
    uint32_t cohortSetOutputId( uint32_t cohortSet );
}

}
}
#endif
