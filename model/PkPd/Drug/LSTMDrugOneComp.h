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

#ifndef Hmod_LSTMDrugOnComp
#define Hmod_LSTMDrugOnComp

#include "Global.h"
#include "PkPd/Drug/LSTMDrug.h"
#include "PkPd/Drug/LSTMDrugType.h"

using namespace std;

namespace OM {

namespace PkPd {
struct DoseParams;
}

namespace PkPd {

/** A class holding pkpd drug use info.
 * 
 * One compartment model
 *
 * Each human has an instance for each type of drug present in their blood. */
class LSTMDrugOneComp : public LSTMDrug {
public:
    /** Create a new instance. */
    LSTMDrugOneComp (const LSTMDrugType&, LocalRng& rng);
    
    virtual size_t getIndex() const;
    virtual double getConcentration(size_t index) const;
    
    virtual double calculateDrugFactor(LocalRng& rng, WithinHost::CommonInfection *inf, double body_mass) const;
    virtual void updateConcentration (double body_mass);
    
protected:
    virtual void checkpoint (istream& stream);
    virtual void checkpoint (ostream& stream);
    
    /// Always links a drug instance to its drug-type data
    //TODO: does it still make sense to link this?
    const LSTMDrugType& typeData;
    
    /// Concentration in blood; units: mg / l.
    double concentration;
    
    /// Sampled elimination rate constant
    double neg_elim_sample;
};

}
}
#endif
