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

#ifndef Hmod_EmpiricalInfection
#define Hmod_EmpiricalInfection

#include "Host/WithinHost/Infection/CommonInfection.h"

namespace OM { namespace WithinHost {
    
class EmpiricalInfection : public CommonInfection {
public:
  ///@brief Static methods
  //@{
  /// Static (shared) data initialisation
  static void init();
  
  /// only for parameterisation?
  static void overrideInflationFactors(double inflationMean, double inflationVariance, double extinctionLevel, double overallMultiplier);
  //@}
  
  /// @brief Construction and destruction
  //@{
  /// For checkpointing (don't use for anything else)
  EmpiricalInfection(istream& stream);
  /// Per instance initialisation; create new inf.
  EmpiricalInfection(LocalRng& rng, uint32_t protID, int origin, double growthRateMultiplier);
  /** Destructor
   * 
   * Note: this destructor does nothing in order to allow shallow copying to
   * the population list. */
  virtual ~EmpiricalInfection() {}
  //@}
  
  /// Set patent growth rate multiplier.
  /// This was used for independant parameterization.
  void setPatentGrowthRateMultiplier(double multiplier);
  
    virtual bool updateDensity( LocalRng& rng, double survivalFactor, SimTime bsAge, double );
  
protected:
    virtual void checkpoint (ostream& stream);
    
private:
  double getInflatedDensity(LocalRng& rng, double nonInflatedDensity);
  double sigma_noise(int ageDays);
  double samplePatentValue(LocalRng& rng, double mu, double sigma, double lowerBound);
  double sampleSubPatentValue(LocalRng& rng, double mu, double sigma, double upperBound);
  
  double _laggedLogDensities[3];
  double _patentGrowthRateMultiplier;
  
  ///@brief Static variables
  ///Set by initParameters and some reset by overrideInflationFactors
  //@{
  static const int _maximumDurationInDays=418; 
  static double _maximumPermittedAmplificationPerCycle;
  static double _subPatentLimit;
  static double _alpha1;
  static double _alpha2;	
  static double _alpha3;
  static double _mu1;	
  static double _mu2;	
  static double _mu3;
  static double _sigma0_res;	
  static double _sigmat_res;
  static double _mu_beta1[_maximumDurationInDays];
  static double _sigma_beta1[_maximumDurationInDays];
  static double _mu_beta2[_maximumDurationInDays];
  static double _sigma_beta2[_maximumDurationInDays];
  static double _mu_beta3[_maximumDurationInDays];
  static double _sigma_beta3[_maximumDurationInDays];
  static double _inflationMean;
  static double _inflationVariance;
  static double _extinctionLevel;
  static double _overallMultiplier;
  //@}
};

} }
#endif
