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

#ifndef Hmod_WithinHost_Descriptive
#define Hmod_WithinHost_Descriptive

#include "Host/WithinHost/WHFalciparum.h"
#include "Host/WithinHost/Infection/DescriptiveInfection.h"

using namespace std;

namespace OM { namespace WithinHost {

/** Within-host model class for the original (descriptive) infection model.
 *
 * Note: this implementation has a few bugs with (hopefully) small effect
 * conditionally fixed (see MAX_DENS_CORRECTION and INNATE_MAX_DENS).
 * This allows reproduction of old results and is the main reason it cannot be
 * integrated with the CommonWithinHost class. */
class DescriptiveWithinHostModel : public WHFalciparum {
public:
    /// Must run after monitoring is set up.
    static void initDescriptive();
    
    /// Create a new WHM
    DescriptiveWithinHostModel( LocalRng& rng, double comorbidityFactor );
    virtual ~DescriptiveWithinHostModel();
    
    virtual void importInfection(LocalRng& rng);
    
    /// load an infection from a checkpoint
    virtual void loadInfection(istream& stream);
    virtual void clearImmunity();
    
    virtual void update(Host::Human &human, LocalRng& rng, int &nNewInfs_i, int &nNewInfs_l, 
        vector<double>& genotype_weights_i, vector<double>& genotype_weights_l, double ageInYears);
    
    virtual InfectionOrigin getInfectionOrigin() const;

    virtual bool summarize( Host::Human& human ) const;
    
protected:
    virtual void clearInfections( Treatments::Stages stage );
    
    virtual void checkpoint (istream& stream);
    virtual void checkpoint (ostream& stream);
    
    // Doesn't do anything in this model:
    virtual void treatPkPd(size_t schedule, size_t dosages, double age, double delay_d);
    
    /** The list of all infections this human has.
     * 
     * Since infection models and within host models are very much intertwined,
     * the idea is that each WithinHostModel has its own list of infections. */
     std::list<DescriptiveInfection*> infections;

     bool opt_vaccine_genotype = false;
};

} }
#endif
