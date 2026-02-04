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

#ifndef Hmod_MOLINEAUXINFECTION_H
#define Hmod_MOLINEAUXINFECTION_H

#include "Host/WithinHost/Infection/CommonInfection.h"

class MolineauxInfectionSuite;

namespace OM { namespace WithinHost {

/**
 * Implementation of a P. falciparum model by Molineaux et. al., published in:
 * 
 * L. MOLINEAUX, H. H. DIEBNER, M. EICHNER, W. E. COLLINS, G. M. JEFFERY and
 * K. DIETZ, 2001: Plasmodium falciparum parasitaemia described by a new
 * mathematical model. Parasitology, 122, pp 379-391
 * doi:10.1017/S0031182001007533
 */
class MolineauxInfection : public CommonInfection {
public:
    ///@brief Static class members
    //@{
    static void init(const OM::Parameters& parameters);
    
    // Constants:
    // v: number of variants per clone (one infection = one new clone); 50 in paper
    //TODO(optimisation): v should have a significant effect on performance, but reducing by,
    // say, half, may not have a big effect on the model. Evaluate.
    static const size_t v = 50;
    // taus: used for the variantTranscending and variantSpecific array, 4 Molineaux time steps = 8 days
    static const size_t taus = 4;
    //@}
    
    // Initialise. Samples several parameters.
    MolineauxInfection(LocalRng& rng, uint32_t protID, int origin);
    // Load from a checkpoint:
    MolineauxInfection (istream& stream);
    virtual ~MolineauxInfection () {};
    
    virtual bool updateDensity( LocalRng& rng, double survivalFactor, SimTime bsAge, double body_mass );
    
protected:
    virtual void checkpoint (ostream& stream);
    
private:
    double getVariantSpecificSummation(int i, double P_current);
    
    // Note: we also have inherited parameters:
    // m_startDate is used to give the age here
    // m_density is equivalent to Pc in paper
    // m_cumulativeExposureJ is cumulative parasite density (used by external immunity function)
    
    float mi[v];        // base multiplication factor per two-day cycle of variant i
    float Sm_summation; // sum in eqn 7
    // index: we use ((bsAge/2) mod 4) for τ = t - δ_v respectively τ = t
    float lagged_Pc[taus];       // Pc(τ) for τ ∈ {t - δ_v, ..., t - 2}
    /* Pc_star, Pm_star: two host-specific critical densities.
     * These two values depend on the first local maximum or the difference
     * between the last positive day and the first positive day. */
    float Pc_star, Pm_star;
    
    struct Variant {
        // Initialise all variables to 0
        Variant ();
        
        /// Checkpointing
        void operator& (ostream& stream);
        void operator& (istream& stream);
        
        float Pi1, Pi2;   // Pi(t+1), Pi(t+2): variant's i density (PRBC/μl blood)
        float Si_summation;     // sum in eqn 6
        // index: we use ((bsAge/2) mod 4) for τ = t - δ_v respectively τ = t
        float lagged_Pi[taus];   // Pi(τ) for τ ∈ {t - δ_v, ..., t - 2}
    };
    // variant-specific data; variants[i-1] corresponds to variant i in the paper
    vector<Variant> variants;
    
    // allow unittest to access private vars
    friend class ::MolineauxInfectionSuite;
};

}}
#endif
