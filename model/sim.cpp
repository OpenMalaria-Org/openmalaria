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

#include "Global.h"
#include "util/errors.h"
#include "util/ModelOptions.h"
#include "util/UnitParse.h"
#include "schema/scenario.h"
#include "mon/management.h"

namespace OM {

void sim::init( const scnXml::Scenario& scenario, util::ModelNameProvider mnp ){
    initInterval(scenario, mnp);
    SimData::steps_per_year = sim::inSteps(sim::oneYear());
    SimData::years_per_step = 1.0 / SimData::steps_per_year;
    sim::s_max_human_age = sim::fromYearsD( scenario.getDemography().getMaximumAgeYrs() );
    
    sim::s_start = sim::origin();
    auto mon = scenario.getMonitoring();
    if( mon.getStartDate().present() ){
        try{
            // on failure, this throws or returns sim::never()
            sim::s_start = UnitParse::parseDate( mon.getStartDate().get() );
            if( sim::s_start == sim::never() ){
                throw util::format_error( "invalid format (expected YYYY-MM-DD)" );
            }
        }catch( const util::format_error& e ){
            throw util::xml_scenario_error( string("monitoring/startDate: ").append(e.message()) );
        }
    }
    
    sim::s_interv = sim::never(); // large negative number
    sim::s_end = mon::readSurveyDates( mon );
}

void sim::initInterval(const scnXml::Scenario& scenario, util::ModelNameProvider mnp)
{
    // If scenario states to use a named model, set the relevant time step duration for
    // that named model here.  It may possibly get overwritten later in this method if
    // the scenario also states an explicit interval value - a behaviour that exists so
    // that users have the ability to override the time step duration if they need to.
    const util::ModelNames namedModelToUse = mnp.GetModelName();
    const bool useNamedModel = namedModelToUse != util::ModelNames::none;
    if (useNamedModel)
    {
        if (namedModelToUse == util::ModelNames::base)
        {
            // This is where we define the length of the default time step, in days,
            // for the base model.
            SimData::interval = 5;
        }
        else
        {
            throw util::xml_scenario_error(
                "No pre-set interval value is available for the specified model name.");
        }
    }

    // Apply the interval duration explicitly written in the scenario, if one exists.
    const bool parameterBlockExistsInScenario = scenario.getModel().getParameters().present();
    if (parameterBlockExistsInScenario)
    {
        // Schema dictates that interval attribute is required if parameter element exists.
        // So no need to check for existsence of interval attribute in this case.
        SimData::interval = scenario.getModel().getParameters().get().getInterval();
    }
}

}
