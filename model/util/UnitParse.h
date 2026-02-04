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

#ifndef Hmod_OM_util_timeConversions
#define Hmod_OM_util_timeConversions

#include "Global.h"

namespace OM {

/** Encapsulation of code to parse times with units from strings.
 * 
 * The goal of this is (a) to make entering values in the scenario documents
 * easier and (b) to avoid confusion.
 * 
 * Use is being implemented by (a) changing existing values in the schema to
 * use this without breaking backwards compatibility (by allowing the unit to
 * be omitted), and (b) by requiring new values be specified with unit. */
namespace UnitParse {
    enum DefaultUnit {
        NONE,   // error if not specified
        DAYS,   // assume days
        YEARS,  // assume years
        STEPS,  // assume steps
    };
    
    static const int monthLen[] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
    
    static const int monthStart[] = { // offsets of start of month from start of year
        0 /*Jan*/, 31, 59, 90,
        120 /*May*/, 151, 181, 212,
        243 /*Sept*/, 273, 304, 334,
        365 /*next year: stop condition in formatDate*/
    };

    /** Parse a short duration from a string found in the input document.
     * 
     * Supports units of days (5d), and steps (2t) with integer values,
     * rounding inputs to the nearest time step.
     * 
     * Call sim::init() first. */
    SimTime readShortDuration( const std::string& str, DefaultUnit defUnit );
    
    /** Like readShortDuration(), but also allow inputs to be in fractional
     * years.
     * 
     * Call sim::init() first. */
    SimTime readDuration( const std::string& str, DefaultUnit defUnit );
    
    /** Like readDuration(), but without rounding. Output is in days and may
     * not be an integer.
     * 
     * Call sim::init() first. */
    double durationToDays( const std::string& str, DefaultUnit defUnit );
    
    /** Returns sim::never() when it doesn't recognise a date. Throws when it does
     * but encounters definite format errors. */
    SimTime parseDate( const std::string& str );

    /** Read a date or relative time specifier found in the XML; dates are
     * rebased relative to a starting date so that they work the same as other
     * ways of specifying intervention-period time from the point of view of
     * code using this function.
     * 
     * Supports dates (e.g. 2015-10-08) as well as times relative to the start
     * of the intervention period (as readDuration will parse). Returns a time
     * to be compared against sim::intervDate(). Again, the result is rounded to
     * the nearest time step. */
    SimTime readDate( const std::string& str, DefaultUnit defUnit );
}

}

#endif
