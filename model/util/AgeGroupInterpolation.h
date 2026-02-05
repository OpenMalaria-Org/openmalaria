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

#ifndef Hmod_AgeDistribution
#define Hmod_AgeDistribution

#include "Global.h"

namespace scnXml {
    class AgeGroupValues;
}
namespace OM { namespace util {

/** This struct is not for external use. Use AgeGroupInterpolator instead.
 * 
 * It is here to allow the trivial functions below to be inlined. */
struct AgeGroupInterpolation{
    virtual ~AgeGroupInterpolation() {}
    
    virtual double eval (double ageYears) const =0;
    virtual void scale( double factor ) =0;
    virtual double firstGlobalMaximum() =0;
    
protected:
    /** Sample interpolator between 0 and max age, outputting to a csv file
     * called name.csv. */
    void outputSamples( const string name );
};

/** A class representing deterministic interpolation of data collected
 * according to age groups. Derived classes implement the actual interpolation.
 * 
 * Current version does not store an age index, thus an order log(n) lookup
 * must occur each time a value is looked up.
 ********************************************/
struct AgeGroupInterpolator
{
    /** Create. set() must be called before further use to avoid exceptions. */
    AgeGroupInterpolator();
    ~AgeGroupInterpolator(){ reset(); }
    
    /** Set age-group interpolation data from an XML fragment.
     * 
     * @param ageGroups XML element of per-age-group values
     * @param eltName Name of XML element (for reasonable error reporting)
     */
    void set( const scnXml::AgeGroupValues& ageGroups, const char* eltName );
    
    /** Set back to initial state. */
    void reset();
    
    /// Return true if set() was ever called.
    bool isSet();
    
    /** Return a value interpolated for age ageYears. */
    inline double eval( double ageYears )const{
        return obj->eval( ageYears );
    }
    
    /** Scale function by factor. */
    inline void scale( double factor ){
        obj->scale( factor );
    }

    /** Find the youngest age which is the global maximum (i.e. the age at
     * which individuals are considered adults where all adults have the same
     * same value). */
    inline double firstGlobalMaximum() const{
        return obj->firstGlobalMaximum();
    }
    
private:
    AgeGroupInterpolation *obj;
};

} }
#endif
