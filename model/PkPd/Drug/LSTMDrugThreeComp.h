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

#ifndef Hmod_LSTMDrugThreeComp
#define Hmod_LSTMDrugThreeComp

#include "Global.h"
#include "PkPd/Drug/LSTMDrug.h"
#include "PkPd/Drug/LSTMDrugType.h"

using namespace std;

namespace OM {

namespace PkPd {
struct DoseParams;
}
namespace PkPd {

struct Params_fC;

/** A class holding PK and PD drug info, per human, per drug type.
 * 
 * Three compartment model:
 * 
 * "Mathematical Expressions of the Pharmacokinetic and Pharmacodynamic Models
 * implemented in the Monolix software", Julie Bertrand and France Mentré, Sept
 * 2008. From page 37, in particular equation (1.72) p44. */
class LSTMDrugThreeComp : public LSTMDrug {
public:
    /** Create a new instance. */
    LSTMDrugThreeComp (const LSTMDrugType&, LocalRng& rng);
    
    virtual size_t getIndex() const;
    virtual double getConcentration(size_t index) const;
    
    virtual double calculateDrugFactor(LocalRng& rng, WithinHost::CommonInfection *inf, double body_mass) const;
    virtual void updateConcentration (double body_mass);
    
protected:
    virtual void checkpoint (istream& stream);
    virtual void checkpoint (ostream& stream);
    
    inline double conc() const{
        return concA + concB + concC - concABC;
    }
    
    // Update constants. bm: body mass (kg)
    void updateCached(double bm) const;
    
    /// Always links a drug instance to its drug-type data
    //TODO: does it still make sense to link this?
    const LSTMDrugType& typeData;
    
    // Concentrations in the blood and other compartments
    // See "Permutation of the three-compartment equation",
    // Diggory Hardy, Swiss TPH, February 3, 2015
    double concA, concB, concC, concABC;
    
    // Sampled constants
    double elim_sample;
    double k12, k21, k13, k31;
    double nka;    // -k_a
    
    // Computed parameters, constant except for dependence on body mass
    // (this is essentially a cache updated by updateCached())
    mutable double last_bm;     // body mass at time of last calculation
    mutable double na, nb, ng;      // -α, -β, -γ, -k_a
    mutable double A, B, C;     // A, B, C
    
private:
    double calculateFactor(const Params_fC& p, double duration) const;
    
    friend double func_fC( double t, void* pp );        // function used in calculateFactor
};

}
}
#endif
