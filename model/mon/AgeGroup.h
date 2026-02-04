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

#ifndef H_OM_mon_AgeGroup
#define H_OM_mon_AgeGroup

#include "Global.h"
#include <vector>

namespace scnXml{ class Monitoring; }
namespace OM {
namespace mon {
/**
 * Included for type-saftey: don't allow implicit double->int conversions.
 *
 * Incindentally, the constructor can be used implicitly for implicit
 * conversion doing the right thing.
 * 
 * Don't use _this_ class for other index/age-group types. */
class AgeGroup {
  public:
    AgeGroup () : index(0) {}
    
    /** Update age-group. Assumes age only increases (per instance).
     *
     * If called regularly, should be O(1); worst case is O(_upperbound.size()). */
    void update (SimTime age);
    
    /// Checkpointing
    template<class S>
    void operator& (S& stream) {
        index & stream;
    }
    
    /** Get the represented index. */
    inline size_t i () const{
        return index;
    }
    
    /// Read age group bounds from XML data
    static void init (const scnXml::Monitoring& monitoring);
    
    /// Get the total number of age categories (inc. one for indivs. not in any
    /// category given in XML).
    static inline size_t numGroups () {
        assert( upperBound.size() > 0 );      // otherwise not yet initialised
        return upperBound.size();
    }
    
private:
    size_t index;
    
    //BEGIN Static parameters only set by init()
    /** Upper boundaries of age groups.
     * 
     * Converted from years (input) to days, rounding down to the next time
     * step.
     *
     * These are age-groups given in XML plus one with no upper limit for
     * individuals outside other bounds. */
    static std::vector<SimTime> upperBound;
    //END
};

} }
#endif
