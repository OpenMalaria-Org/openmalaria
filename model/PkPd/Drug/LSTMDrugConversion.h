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

#ifndef Hmod_LSTMDrugConversion
#define Hmod_LSTMDrugConversion

#include "Global.h"
#include "PkPd/Drug/LSTMDrug.h"
#include "PkPd/Drug/LSTMDrugType.h"

using namespace std;

namespace OM {

namespace PkPd {
struct DoseParams;
}

namespace PkPd {

struct Params_convFactor;

/** A class holding pkpd drug use info.
 * 
 * Conversion model: two 1-compartment models with conversion.
 *
 * Each human has an instance for each type of drug present in their blood. */
class LSTMDrugConversion : public LSTMDrug {
public:
    /** Create a new instance. */
    LSTMDrugConversion (const LSTMDrugType& parent, const LSTMDrugType& metabolite, LocalRng& rng);
    
    virtual size_t getIndex() const;
    virtual double getConcentration(size_t index) const;
    
    virtual double calculateDrugFactor(LocalRng& rng, WithinHost::CommonInfection *inf, double body_mass) const;
    virtual void updateConcentration (double body_mass);
    double getMetaboliteConcentration() const;
    double getParentConcentration() const;
    
protected:
    virtual void checkpoint (istream& stream);
    virtual void checkpoint (ostream& stream);
    
    double calculateFactor(const Params_convFactor& p, double duration) const;
    
    //TODO: do we need to link these?
    const LSTMDrugType &parentType, &metaboliteType;

    /// Amount of parent drug in gut, parent drug in circulation and metabolite
    // in circulation; units: mg. In the paper, these are labelled A, B, C.
    //TODO: move metabolite data to a separate class object? Or not?
    double qtyG, qtyP, qtyM;
    
    /// @brief Sampled elimination and conversion rate constants
    //@{
    /// Absorbtion rate constant (-x)
    double nka;
    /// Elimination rate constant for parent (-y * body_mass ^ m_exponent)
    double nkP_sample;
    /// Conversion rate constant (-z * body_mass ^ m_exponent)
    double nconv_sample;
    /// Elimination rate constant for metabolite (-k * body_mass ^ m_exponent)
    double nkM_sample;
    /// Volume of distribution of metabolite (parent's Vd uses vol_dist)
    double vol_dist_metabolite;
    /// Last body mass (kg)
    double last_bm;
    //@}
    
private:
    void setConversionParameters(Params_convFactor& p, double body_mass) const;
    void setKillingParameters(LocalRng& rng, Params_convFactor& p, WithinHost::CommonInfection *inf) const;
};

}
}
#endif
