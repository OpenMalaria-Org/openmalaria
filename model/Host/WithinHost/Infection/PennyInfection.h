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

/*This model refers to the paper:
 *PENNY et al  (2011). The potential effects of blood stage vaccines on the within-host dynamics of Plasmodium Falciparum */

#ifndef Hmod_PENNYINFECTION_H
#define Hmod_PENNYINFECTION_H

#include "Host/WithinHost/Infection/CommonInfection.h"

class PennyInfectionSuite;

namespace OM { namespace WithinHost {

class PennyInfection : public CommonInfection {
public:
    /// Static initialization (happens once)
    static void init();
    
    /// Constructor
    PennyInfection(LocalRng& rng, uint32_t protID, int origin);
    /// Resume from a checkpoint
    PennyInfection (istream& stream);
    /// Destructor
    virtual ~PennyInfection () {};
    
    virtual bool updateDensity( LocalRng& rng, double survivalFactor, SimTime bsAge, double );
    
    /** Get the density of sequestered parasites. */
    inline double seqDensity(int ageDays){
        size_t todayV = mod_nn(ageDays, delta_V);
        return seqDensities[todayV];
    }

protected:
    virtual void checkpoint (ostream& stream);

private:
    // function to obtain the summation component of variant specific immunity
    double getVariantSpecificSummation(LocalRng& rng, int ageDays);
    // function to obtain the summation component of clonal immunity
    double getClonalSummation(int ageDays);
    
    // delta_C := delay to clonal antibody response (days) (value 7.2038 round to 7)  
    static const int delta_C=7;
    // delta_V := delay to variant specific antibody response in R_V^x (days) (value 6.3572 round to 6)  
    static const int delta_V=6;
    
    static bool immune_threshold_gamma;
    static bool update_density_gamma;
    
    // Circulating densities, 1 to delta_C days ago.
    // index (age mod delta_C) corresponds to delta_C days ago.
    double cirDensities[delta_C];
    // as above, but length delta_V
    double seqDensities[delta_V];
    
    // threshold_N is critical threshold for innate immunity (for sigmoidal immune function)
    double threshold_N;
    // threshold_V is critical threshold for variant specific immunity (for sigmoidal immune function)
    double threshold_V;
    // threshold_C is critical threshold for clonal immunity (for sigmoidal immune function)
    double threshold_C;
    
    // tracked summation of densities with decay for variant specific immunity (sigmoidal function)
    double variantSpecificSummation;
    // tracked summation of densities with decay for clonal immunity (sigmoidal function)
    double clonalSummation;
    
    // allow unittest to access private vars
    friend class ::PennyInfectionSuite;
};

}}
#endif
