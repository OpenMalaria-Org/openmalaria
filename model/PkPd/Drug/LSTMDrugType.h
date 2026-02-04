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

#ifndef Hmod_LSTMDrugType
#define Hmod_LSTMDrugType

#include "Global.h"
#include "util/sampler.h"

#include <string>
#include <deque>
#include <vector>
#include <set>
#include <cassert>
#include <memory>
#include <map>

using namespace std;

namespace scnXml{
    class Drugs;
    class PKPDDrug;
    class Phenotype;
}
namespace OM {
namespace WithinHost {
    class CommonInfection;
}
namespace PkPd {
using util::LognormalSampler;
using util::NormalSample;
using util::LocalRng;

class LSTMDrugType;
class LSTMDrug;

/** Drug PD parameters (specified per phenotype), as well as applicable
 * functions to calculate drug factors and concentrations.
 * 
 * Parasite genetics is specified elsewhere; we simply map genotype to
 * phenotype as well as describing phenotypes here. */
class LSTMDrugPD {
public:
    /// Construct
    LSTMDrugPD( const scnXml::Phenotype& phenotype );
    
    LSTMDrugPD(LSTMDrugPD&&) = default;
    LSTMDrugPD& operator=(LSTMDrugPD&&) = default;
    
    /** Calculate a survival factor induced by a drug already in the blood.
     * It is expected that no drug doses are taken over the period for which
     * this function calculates a drug factor.
     * 
     * @param Kn IC50^slope, sampled per infection
     * @param neg_elim_rate -k (sampled)
     * @param C0 Concentration of drug in blood at start of period. Will be
     *  updated to correct concentration at end of period. Units: mg/l
     * @param duration Timespan over which the factor is being calculated. Units: days.
     * @return survival factor (unitless)
     */
    double calcFactor( double Kn, double neg_elim_rate, double* C0, double duration ) const;
    
    inline double slope() const{ return n; }
    double IC50_pow_slope(LocalRng& rng, size_t index, WithinHost::CommonInfection *inf) const;
    inline double IC50_pow_slope(NormalSample normal) const {
        return pow(IC50->sample(normal), n);
    }
    inline double max_killing_rate() const{ return V; }
    
private:
    /// Slope of the dose response curve (no unit)
    double n;   // slope
    /// Concentration with 50% of the maximal parasite killing
    std::unique_ptr<LognormalSampler> IC50;  // IC50
    /// Maximal drug killing rate per day (units: 1/days)
    double V;   // max_killing_rate;
};

    
/** Drug PK & PD parameters (one instance per drug type).
 * 
 * Static data contains a list of all available drug types.
 * 
 * None of this data is checkpointed, because it is loaded by init() from XML
 * data. (Although if it cannot be reproduced by reloading it should be
 * checkpointed.) */
class LSTMDrugType {
public:
    ///@brief Static functions
    //@{
    /** Initialise the drug model. Called at start of simulation. */
    static void init (const scnXml::Drugs& data);
    /** Clear previous data. Only needed for testing. */
    static void clear();
    
    /** Get the number of drug types. */
    static size_t numDrugTypes();
    
    /** Find a DrugType by its abbreviation, and returns its index.
     *
     * Throws if the drug isn't found, so you can rely on it returning a valid
     * index if it doesn't throw. */
    static size_t findDrug(string abbreviation);
    
    /** Get a reference to drug type data for some index. This can't be used
     * to change the data. */
    static LSTMDrugType& get(size_t index);
    
    /** Get a list of all drug types which are (possibly) being used. */
    static const vector<size_t>& getDrugsInUse();
    
    /** Create a per-human drug module for a given drug index. */
    static unique_ptr<LSTMDrug> createInstance( LocalRng& rng, size_t index );
    //@}
    
    
    ///@brief Non static (per instance) functions
    //@{
    /** Create a new DrugType.
     *
     * @param index     Index of drug in internal list
     * @param drugData Scenario data for this drug (PK params, PD params per phenotype)
     */
    LSTMDrugType (size_t index, const scnXml::PKPDDrug& drugData);
    ~LSTMDrugType ();
    
    LSTMDrugType(LSTMDrugType&&) = default;
    LSTMDrugType& operator=(LSTMDrugType&&) = default;
    
    inline size_t getIndex() const {
        return index;
    }
    inline double getNegligibleConcentration() const{
        return negligible_concentration;
    }
    inline double neg_m_exponent() const{ return neg_m_exp; }
    inline double molecular_weight_ratio() const{ return mwr; }
    
    /// Generate IC50 sample of metabolite
    inline NormalSample IC50_correlated_sample(NormalSample parent, LocalRng& rng) const {
        return NormalSample::generate_correlated(parent, ic50_log_corr, ic50_corr_factor, rng);
    }

    inline double sample_Vd(LocalRng& rng) const{
        return vol_dist->sample(rng);
    }
    inline double sample_elim_rate(LocalRng& rng) const{
        return elimination_rate->sample(rng);
    }
    inline double sample_conv_rate(LocalRng& rng) const{
        return conversion_rate->sample(rng);
    }
    inline double sample_k12(LocalRng& rng) const{ return k12->sample(rng); }
    inline double sample_k21(LocalRng& rng) const{ return k21->sample(rng); }
    inline double sample_k13(LocalRng& rng) const{ return k13->sample(rng); }
    inline double sample_k31(LocalRng& rng) const{ return k31->sample(rng); }
    inline double sample_ka(LocalRng& rng) const{ return absorption_rate->sample(rng); }
    
    /** Return reference to correct drug-phenotype data. */
    const LSTMDrugPD& getPD( uint32_t genotype ) const;
    //@}
  
private:
    // non-copyable (due to allocation of members of drugPhenotype)
    LSTMDrugType( const LSTMDrugType& ) = delete;
    LSTMDrugType& operator= ( const LSTMDrugType& ) = delete;
    
    /** The drug index in drugTypes and a unique identifier for the drug type.
     * Stored here since LTSMDrug stores a pointer to this struct object, not the index. */
    size_t index;
    
    /// Index of metabolite
    size_t metabolite;
    
    /** A mapping from genotype codes to phenotypes for this drug. */
    vector<uint32_t> genotype_mapping;
    
    // TODO: at the moment we're storing all PK & PD data regardless of which
    // model is used. Evaluate whether this is sensible or not.
    
    /* PD parameters (may vary based on infection genotype). */
    vector<LSTMDrugPD> PD;
    
    /*PK parameters required - varies with humans age and severity of disease*/
    /** Concentration, below which drug is deemed not to have an effect and is
     * removed for performance reasons. (mg/l) */
    double negligible_concentration;

    // Used to calculate elimination rate
    double neg_m_exp;
    double mwr;      // set for parent, not metabolite
    double ic50_log_corr;
    double ic50_corr_factor;
    
    /// Volume of distribution (l/kg)
    std::unique_ptr<LognormalSampler> vol_dist;

    /// Absorbtion rate
    std::unique_ptr<LognormalSampler> absorption_rate;

    /** Terminal elimination rate constant (k). Equals ln(2)/half_life.
     * Units: (1 / days) */
    std::unique_ptr<LognormalSampler> elimination_rate;

    /// Convertion rate
    std::unique_ptr<LognormalSampler> conversion_rate;

    /// Parameters for absorption rates in two- and three-compartment models.
    std::unique_ptr<LognormalSampler> k12, k21, k13, k31;
    
    // Allow LSTMDrug to access private members
    friend class LSTMDrugPD;
    friend inline double drugEffect (const LSTMDrugType& drugType, double& concentration, double duration, double weight_kg, double dose_mg);
};

} }
#endif
