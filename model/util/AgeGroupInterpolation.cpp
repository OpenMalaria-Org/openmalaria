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

#include <fstream>

#include "util/AgeGroupInterpolation.h"
#include "util/CommandLine.h"
#include "util/errors.h"
#include "schema/healthSystem.h"

/* Required by AgeGroupSplineInterpolation
#include <gsl_interp.h>
#include <gsl_spline.h>
*/

namespace OM { namespace util {

// ———  AgeGroupInterpolation derivatives  ———
    
    void AgeGroupInterpolation::outputSamples( const string name ){
        if( !util::CommandLine::option(util::CommandLine::SAMPLE_INTERPOLATIONS) ){
            return;
        }
        ofstream fstream( (name+".csv").c_str() );
        
        double max = sim::inYears(sim::maxHumanAge());
        for( double age = 0.0; age < max; age += 0.1 ){
            fstream << age << "," << this->eval( age ) << endl;
        }
    }
    
    /** Dummy object, for initialization. Potentially makes code safer.
     ********************************************/
    class AgeGroupDummy : public AgeGroupInterpolation
    {
    public:
        static AgeGroupDummy singleton;
        
        virtual double eval( double ageYears ) const {
            throw TRACED_EXCEPTION_DEFAULT( "AgeGroupDummy" );
        }
        virtual void scale( double factor ) {
            throw TRACED_EXCEPTION_DEFAULT( "AgeGroupDummy" );
        }
        virtual double firstGlobalMaximum() {
            throw TRACED_EXCEPTION_DEFAULT( "AgeGroupDummy" );
        }
    };
    AgeGroupDummy AgeGroupDummy::singleton;
    
    
    /** This class gives direct access to input age-group
    * data (discontinuous).
    ********************************************/
    class AgeGroupPiecewiseConstant : public AgeGroupInterpolation
    {
    public:
        AgeGroupPiecewiseConstant ( const scnXml::AgeGroupValues& ageGroups,
                                    const char* eltName )
        {
            if( ageGroups.getGroup().size() == 0 ){
                throw util::xml_scenario_error(string(eltName) + ": at least one age group required" );
            }

            auto pos = dataGroups.begin();
            assert(pos == dataGroups.end());

            double greatestLbound = -1.0;
            for( const scnXml::Group& group : ageGroups.getGroup() ){
                double lbound = group.getLowerbound();
                if( lbound >= greatestLbound ){
                    greatestLbound = lbound;
                    pos = dataGroups.insert( pos, make_pair( lbound, group.getValue() ) );
                } else
                    throw util::xml_scenario_error(string(eltName) + ": lower bound " + to_string(lbound) + " not greater than previous " + to_string(greatestLbound));
            }

            // first lower-bound must be 0
            if( dataGroups.begin()->first != 0.0 )
                throw util::xml_scenario_error(string(eltName) + ": first lower-bound must be 0");

            outputSamples( eltName );
        }

        virtual double eval( double ageYears ) const {
            assert ( ageYears >= 0.0 );
            auto it = dataGroups.upper_bound( ageYears );
            --it;   // previous point when ordered by age (it->first <= ageYears < (it+1)->first)
            return it->second;
        }
        
        virtual void scale( double factor ){
            for( auto it = dataGroups.begin(); it != dataGroups.end(); ++it ){
                it->second *= factor;
            }
        }

        virtual double firstGlobalMaximum() {
            double age = 0.0, max = numeric_limits<double>::min();
            for( auto it = dataGroups.begin(); it != dataGroups.end(); ++it ){
                if( it->second > max ){
                    max = it->second;
                    age = it->first;
                }
            }
            return age;
        }

    protected:
        // All data groups as (lower-age-bound, value) pairs.
        map<double,double> dataGroups;
    };

    
    /** Filter class to convert data groups into points in middle of groups plus
     * stabilization points at ends.
     ********************************************/
    class AgeGroupInterpolationPoints : public AgeGroupInterpolation
    {
    public:
        AgeGroupInterpolationPoints( const scnXml::AgeGroupValues& ageGroups,
                                     const char* eltName )
        {
            // Our read iterator
            auto it = ageGroups.getGroup().begin();

            if( it == ageGroups.getGroup().end() )
                throw util::xml_scenario_error(string(eltName) + ": at least one age group required");

            // first lower-bound must be 0
            if( it->getLowerbound() != 0.0 )
                throw util::xml_scenario_error(string(eltName) + ": first lower-bound must be 0");

            // Our insert iterator (used for performance)
            auto pos = dataPoints.begin();
            assert(pos == dataPoints.end());

            double greatestLbound = 0.0;    // also last lbound: age groups must be ordered lowest to highest
            double lastValue =it->getValue();       // first value is repeated for const start
            for( ; it!=ageGroups.getGroup().end(); ++it ){
                double lbound = it->getLowerbound();
                if( lbound >= greatestLbound ){
                    double insBound = 0.5*(greatestLbound + lbound);
                    pos = dataPoints.insert( pos, make_pair( insBound, lastValue ) );
                    greatestLbound = lbound;
                    lastValue = it->getValue();     // always insert previous value
                } else
                    throw util::xml_scenario_error(string(eltName) + ": lower bound " + to_string(lbound) + " less than previous " + to_string(greatestLbound));
            }

            // add a point in middle of last age group (taking upper bound as max-age-years:
            dataPoints[ 0.5*(greatestLbound + sim::inYears(sim::maxHumanAge())) ] = lastValue;
        }

        virtual void scale( double factor ){
            for( auto it = dataPoints.begin(); it != dataPoints.end(); ++it ){
                it->second *= factor;
            }
        }
        
    protected:
        // Points to interpolate between in the middle of input age groups. Extra
        // points at zero and infinity are added with equal value to first and last
        // points respectively.
        map<double,double> dataPoints;
    };

    
    /** This class gives piecewise linear inpolation on top of input age-group
     * data (continuous but with discontinuous derivative).
     ********************************************/
    class AgeGroupPiecewiseLinear : public AgeGroupInterpolationPoints
    {
    public:
        AgeGroupPiecewiseLinear (
            const scnXml::AgeGroupValues& ageGroups, const char* eltName
        ) : AgeGroupInterpolationPoints(ageGroups, eltName)
        {
            // Add first point at zero
            dataPoints[ 0.0 ] = dataPoints.begin()->second;
            // add a point at infinity to catch remaining values
            dataPoints[ numeric_limits<double>::infinity() ] = dataPoints.rbegin()->second;
            
            outputSamples( eltName );
        }
        
        virtual double eval( double ageYears ) const{
            assert ( ageYears >= 0.0 );
            //TODO: optimise by putting (f1-f0)/(a1-a0) into the array and checking usage
            // Get the first point with age greater than ageYears:
            auto it = dataPoints.upper_bound( ageYears );
            assert( it != dataPoints.end() );
            double a1 = it->first;  // a1 > ageYears
            double f1 = it->second;
            --it;   // previous point when ordered by age
            double a0 = it->first;  // a0 ≤ ageYears
            double f0 = it->second;
            return (ageYears - a0) / (a1 - a0) * (f1 - f0) + f0;
        }

        virtual double firstGlobalMaximum() {
            double age = 0.0, max = numeric_limits<double>::min();
            for( auto it = dataPoints.begin(); it != dataPoints.end(); ++it ){
                if( it->second > max ){
                    max = it->second;
                    age = it->first;
                }
            }
            return age;
        }
    };
    
    
    /* * Spline interpolation using GSL.
     * 
     * GSL has a few interpolation algorithms: linear, polynomial,
     * cspline, akima. They fit a single curve to the entire interval which
     * looks nothing like the data set in some cases (e.g. CFR).
     ******************************************** /
    class AgeGroupSplineInterpolation : public AgeGroupInterpolationPoints
    {
    public:
        AgeGroupSplineInterpolation (
            const scnXml::AgeGroupValues& ageGroups, const char* eltName
        ) : AgeGroupInterpolationPoints(ageGroups, eltName)
        {
            vector<double> x, y;
            x.reserve( dataPoints.size() );
            y.reserve( dataPoints.size() );
            for( auto it=dataPoints.begin(); it!=dataPoints.end(); ++it ){
                x.push_back( it->first );
                y.push_back( it->second );
            }
            
            acc = gsl_interp_accel_alloc ();
            spline = gsl_spline_alloc (gsl_interp_cspline, 10);
            gsl_spline_init (spline, x.data(), y.data(), 10);
            
            outputSamples( eltName );
        }
        virtual ~AgeGroupSplineInterpolation() {
            gsl_spline_free (spline);
            gsl_interp_accel_free (acc);
        }
        
        virtual double eval( double ageYears ) const{
            return gsl_spline_eval( spline, ageYears, acc );
        }
        
    protected:
        gsl_interp_accel *acc;
        gsl_spline *spline;
    };
    */
    
    
    // -----  AgeGroupInterpolator  -----
    
    AgeGroupInterpolator::AgeGroupInterpolator() :
        obj(&AgeGroupDummy::singleton) {}
    
    void AgeGroupInterpolator::set(
        const scnXml::AgeGroupValues& ageGroups, const char* eltName
    ){
        reset();
        // Type mostly equivalent to a std::string:
        const scnXml::AgeGroupValues::InterpolationOptional& interp = ageGroups.getInterpolation();
        if( !interp.present() || interp.get() == "linear" ){
            obj = new AgeGroupPiecewiseLinear( ageGroups, eltName );
        }else if( interp.get() == "none" ){
            obj = new AgeGroupPiecewiseConstant( ageGroups, eltName );
        }else
            throw util::xml_scenario_error(string("age group interpolation ") + interp.get() + " not implemented" );
    }
    void AgeGroupInterpolator::reset(){
        assert( obj != nullptr );  // should not do that
        if( obj != &AgeGroupDummy::singleton ){
            delete obj;
            obj = &AgeGroupDummy::singleton;
        }
    }
    bool AgeGroupInterpolator::isSet()    {
        return obj != &AgeGroupDummy::singleton;
    }
} }
