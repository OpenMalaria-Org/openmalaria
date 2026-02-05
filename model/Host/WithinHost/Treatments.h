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

#ifndef Hmod_WithinHost_Treatments
#define Hmod_WithinHost_Treatments

#include "Host/WithinHost/WHInterface.h"
#include "interventions/Interfaces.h"

using namespace std;

namespace OM {
namespace WithinHost {

/**
 * Objects of this class describe effects of a treatment, after selecting a
 * compliance/adherence/... group.
 * 
 * The static part provides a way to find objects.
 * 
 * For use within WithinHost only.
 */
class Treatments : public interventions::TriggeredDeployments {
public:
    /// @brief Static methods
    //@{
    /** Configure a new treatment option, and return the code used to select
     * that option later. */
    static TreatmentId addTreatment( const scnXml::TreatmentOption& desc );
    
    /** Return the corresponding treatment description. */
    static inline const Treatments& select( TreatmentId treatId ){
        assert( treatId.id < treatments.size() );
        return treatments[treatId.id];
    }
    //@}
    
    /// @brief Types used by the non-static methods
    //@{
    /// Stages effected
    enum Stages{
        NONE /*i.e. no effect*/,
        LIVER = 1,
        BLOOD = 2,
        BOTH = LIVER | BLOOD
    };
    //@}
    
    /// @brief Non-static methods
    //@{
    /** Get the liver stage action.
     *
     * 0 implies no action, <0 implies retrospective action, and n>0
     * implies treatment for next n time steps. */
    inline SimTime liverEffect()const{ return timeLiver; }
    /** Get the blood stage action.
     *
     * 0 implies no action, <0 implies retrospective action, and n>0
     * implies treatment for next n time steps. */
    inline SimTime bloodEffect()const{ return timeBlood; }
    //@}
    
private:
    // static:
    static vector<Treatments> treatments;
    
    // non-static:
    Treatments( const scnXml::TreatmentOption& elt );
    
    SimTime timeLiver, timeBlood;
    
    friend class ::UnittestUtil;
};

}
}
#endif
