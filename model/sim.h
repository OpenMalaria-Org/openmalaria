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

#ifndef Hmod_OM_SimTime
#define Hmod_OM_SimTime

#ifndef Hmod_Global
#error "Please include Global.h not SimTime.h directly."
// because checkpoint.h needs this, and we need checkpoint.h
#endif

#include "util/mod.h"
#include "util/checkpoint.h"
#include "util/ModelNameProvider.h"
#include <memory>

class UnittestUtil;

namespace scnXml {
    class Scenario;
}

namespace OM {

inline int floorToInt( double x ){
	return static_cast<int>(std::floor(x));
}

/// Internal use only
// Required to avoid use-before-definition
class SimData {
    static inline int interval;        // days per time step
    static inline size_t steps_per_year;
    static inline double years_per_step;
    friend class sim;
};


/******************************************************************************
 * Since version 43.1, SimTime is an integer instead of a class. Note: a 
 * class can be used again simply by replacing the following typedef 
 * definition by the old SimTime class.
 * 
 * SimTime is now always a number of days.
 * 
 * Time steps, days and dates are derived from this. The values and units of
 * internal variables are an implementation detail (i.e. code outside of
 * sim should not need to know).
 * 
 * The simulation always starts at time zero. "Intervention time" is the 
 * beginning of the scenario and monitoring period.
 * 
 * Granularity: 1 day.
 *****************************************************************************/
typedef int SimTime;

/** Encapsulates static variables: sim time. */
class sim {
public:
    ///@brief Simulation constants
    //@{
    /// Number of days in a year; defined as 365 (leap years are not simulated).
    static const int DAYS_IN_YEAR = 365;
    
    ///@brief Conversions to other types/units
    //@{
    /// Convert to years
    static inline double inYears(SimTime d) { return d * (1.0 / DAYS_IN_YEAR); }
    
    /// Convert to time steps (rounding down)
    static inline int inSteps(SimTime d) { return d / SimData::interval; }
    //@}
    
    /** Return this time in time steps modulo some positive integer. */
    static inline int moduloSteps(SimTime d, int denominator){
        return util::mod_nn(d / SimData::interval, denominator);
    }
    
    /** Return this time in time steps modulo some positive integer. */
    static inline int moduloYearSteps(SimTime d){
        return util::mod_nn(d / SimData::interval, SimData::steps_per_year);
    }

    /** Duration zero and the time at the start of the simulation. */
    static inline SimTime zero(){ return SimTime(0); }

    static inline SimTime origin(){ return SimTime(0); }

        /** Special value representing a time point always in the past, such that
     * never() + x < zero() and x - never() will not to overflow for all valid
     * simulation times x (including any value now() may take as well as
     * never() and future()). */
    static inline SimTime never(){ return SimTime(-0x3FFFFFFF); }
    
    /** Special value representing a time point always in the future, such that
     * now() < future() and now() + future() does not overflow. */
    static inline SimTime future(){ return SimTime(0x3FFFFFFF); }
    
    /** One day. */
    static inline SimTime oneDay(){ return SimTime(1); }
    
    /** One year. See DAYS_IN_YEAR. */
    static inline SimTime oneYear(){ return SimTime(DAYS_IN_YEAR); }
    
    /** One time step (currently either one or five days). */
    static inline SimTime oneTS(){ return SimTime(SimData::interval); }
    //@}
    
    ///@brief Parameterised constructors
    //@{
    /** Convert. */
    static inline SimTime fromTS(int ts){ return SimTime(int(oneTS()) * ts); }
    
    /** Duration in days. Should be fast (currently no conversion required). */
    static inline SimTime fromDays(int days){ return SimTime(days); }
    
    /** Convert from a whole number of years. */
    static inline SimTime fromYearsI(int years){
        return SimTime(DAYS_IN_YEAR * years);
    }
    
    /** Convert from years to nearest time step. */
    static inline SimTime fromYearsN(double years){
        return roundToTSFromDays(DAYS_IN_YEAR * years);
    }
    
    /** Convert from years, rounding down to the next time step. */
    static inline SimTime fromYearsD(double years){
        return fromTS( floorToInt(SimData::steps_per_year * years) );
    }
    
    /** Round to the nearest time-step, where input is in days. */
    static inline SimTime roundToTSFromDays(double days){
        return fromTS( floorToInt( days / SimData::interval + 0.5 ));
    }
    //@}

    /** The number of time steps in one year. */
    static inline size_t stepsPerYear(){ return SimData::steps_per_year; }
    
    /** A cached value: one year divided by one time step. */
    static inline double yearsPerStep(){ return SimData::years_per_step; }
    
    /// Maximum possible age of a human.
    static inline SimTime maxHumanAge(){ return s_max_human_age; }
    
    /// The starting date of the simulation
    static inline SimTime startDate() { return s_start; }
    
    /// The ending date of the simulation
    static inline SimTime endDate() { return s_end; }
    //@}
    
    ///@brief Access simulation time variables
    //@{
    /** Time at the beginning of a time step update.
     *
     * This is what is mostly used during an update. It is never negative and
     * increases throughout the simulation. */
    static inline SimTime ts0(){
        assert(in_update);      // should only be used during updates
        return s_t0;
    }
    /** Time at the end of a time step update.
     * 
     * During an update, ts0() + oneTS() = ts1(). Neither this nor ts0 should
     * be used outside of updates. */
    static inline SimTime ts1(){
        assert(in_update);      // should only be used during updates
        return s_t1;
    }
    /**
     * Time steps are mid-day to mid-day, and this is the time at mid-day (i.e.
     * this equals ts1 from the last step and ts0 from the next one).
     * 
     * This is for monitoring and intervention deployment which happens between
     * updates. Cannot be used during human or vector update.
     */
    static inline SimTime now(){
        assert(!in_update);     // only for use outside of step updates
        return s_t0;   // which is equal to s_t1 outside of updates, but that's a detail
    }
    /** During updates, this is ts0; between, this is now. */
    static inline SimTime nowOrTs0(){ return s_t0; }
    /** During updates, this is ts1; between, this is now. */
    static inline SimTime nowOrTs1(){ return s_t1; }
    /** During updates, this is ts0; between, it is now - 1. */
    static inline SimTime latestTs0(){ return s_t1 - sim::oneTS(); }
    //@}
    
    ///@brief Access intervention-time variables
    //@{
    /// Time relative to the start of the intervention period.
    /// 
    /// This equals (intervDate() - startDate()), but happens to be the most
    /// common way that intervention-period dates are used.
    static inline SimTime intervTime() {
        return s_interv;
    }
    
    /// The current date.
    /// 
    /// Only valid during the intervention phase, since the duration required
    /// for warm-up is not known in advance. (In prior phases, this function
    /// returns a large negative value.)
    /// 
    /// Intervention deployment times are relative to this date.
    static inline SimTime intervDate(){ return s_start + s_interv; }
    //@}
    
// private:
    // Initial set-up: called by Simulator
    static void init( const scnXml::Scenario& scenario, util::ModelNameProvider mnp);

    // Sets up time step duration.
    // Assumes that the scenario either:
    // (1) states to use a named model which has some value for interval, or
    // (2) includes a <parameters> element, which sets an interval value.
    static void initInterval(const scnXml::Scenario& scenario, util::ModelNameProvider mnp);
    
    // Start of update: called by Simulator
    static inline void start_update(){
        s_t1 = s_t1 + sim::oneTS();
#ifndef NDEBUG
        in_update = true;
#endif
    }
    // End of update: called by Simulator
    static inline void end_update(){
#ifndef NDEBUG
        in_update = false;
#endif
        s_t0 = s_t1;
        s_interv = s_interv + sim::oneTS();
    }
    
    // Scenario constants
    static inline SimTime s_start = sim::never();
    static inline SimTime s_end = sim::never();
    static inline SimTime s_max_human_age = sim::never();
    static inline SimTime s_t0 = sim::never();
    static inline SimTime s_t1 = sim::never();
    static inline SimTime s_interv = sim::never();

#ifndef NDEBUG
    static inline bool in_update = false;       // only true during human/population/transmission update
#endif
};

}
#endif
