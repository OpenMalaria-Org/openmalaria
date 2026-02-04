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

#ifndef Hmod_ClinicalEventSchduler
#define Hmod_ClinicalEventSchduler

#include "Global.h"
#include "Clinical/ClinicalModel.h"
#include "Clinical/ESCaseManagement.h"
#include "util/AgeGroupInterpolation.h"

#include <list>

namespace OM {
namespace Clinical {

using util::AgeGroupInterpolator;

/** Tracks clinical status (sickness), triggers case management for new events,
 * medicates treatment, determines patient recovery, death and sequelae.
 *
 * Note: there are several variables that only need to be used during a
 * bout. It's possible that memory usage could be reduced by storing them
 * externally in a temporary object during episodes (but unlikely worth doing).
 */
class ClinicalEventScheduler : public ClinicalModel
{
public:
    static void init (const OM::Parameters& parameters, const scnXml::Clinical& clinical);
    static void setParameters (const scnXml::HSEventScheduler& esData);

    ClinicalEventScheduler (double tSF);
    
    virtual bool isExistingCase();

protected:
    virtual void doClinicalUpdate (Human& human, double ageYears, WithinHost::Pathogenesis::StatePair &pg);

    virtual void checkpoint (istream& stream);
    virtual void checkpoint (ostream& stream);

private:
    /// Maximum number of time steps (including first of case) an individual
    /// will remember they are sick before resetting.
    static SimTime maxUCSeekingMemory;
    /// Length of an uncomplicated case
    static SimTime uncomplicatedCaseDuration;
    /// Length of a complicated case
    static SimTime complicatedCaseDuration;
    /// Time-span for which individual is at risk of death in complicated case
    /// minus length of complicated case (must be <= 0)
    static SimTime extraDaysAtRisk;
    /// First value is probability immediate treatment, second is first +
    /// probability 1-day delay to treatment seeking, etc. Last value must be 1.
    /// Index units are days.
    static vector<double> cumDailyPrImmUCTS;

    /// Parameter of S(t) for t > 0
    static double neg_v;
    /// Parameter
    static double alpha;
    
    /** Base log odds of treatment of non-malarial fevers in absense of a
     * malaria diagnostic and irrespective of whether treatment is needed.
     * 
     * In our model, this is logit(P₀), not β₀. */
    static double logOddsAbBase;
    /** Added to log odds treatment when a malaria diagnostic indicates no
     * parasites. Symbol in model: β₁. */
    static double logOddsAbNegTest;
    /** Added to log odds treatment when a malaria diagnostic indicates
     * parasites. Symbol in model: β₂. */
    static double logOddsAbPosTest;
    /** Added to log odds treatment when NMF is categorized as an illness
     * potentially leading to death (Pathogenesis::NEED_ANTIBIOTIC).
     * Symbol in model: β₃. */
    static double logOddsAbNeed;
    /** Added to log odds treatment when given by an informal provider.
     * Symbol in model: β₄. */
    static double logOddsAbInformal;
    /** One minus the efficacy of antibiotic/NMF treatment (i.e. a multiplier
     * for fatality-rate given that the case is treated). */
    static double oneMinusEfficacyAb;
    /** Case fatality rate of non-malaria fevers requiring treatment given that
     * the case is not treated. */
    static AgeGroupInterpolator severeNmfMortality;
    
    /// Probability that an NMF needs antibiotic treatment and could lead to death.
    static AgeGroupInterpolator NMF_need_antibiotic;
    /** Probability that a malarial fever classified as uncomplicated requires
     * antibiotic treatment (and death could occur from non-malarial causes).
     * Unclear whether using this at the same time as comorbidity parameters
     * makes any sense. */
    static AgeGroupInterpolator MF_need_antibiotic;
    
    // Note on memory usage: Pathogenesis::State is and enum (an int), so we
    // have a vtable followed by 3 ints, a double and a list. Alignment probably
    // wastes some space.
    /// Current state of sickness
    Episode::State pgState;

    /** Set to when a bout should start. If sim::ts0() equals this, a bout
     * is started (UC & severe behaviour different).
     *
     * Note: medications are not delayed by this. */
    SimTime caseStartTime = sim::never();

    /** The individual recovers when sim::ts0() >= timeOfRecovery,
     * assuming they didn't die. */
    SimTime timeOfRecovery = sim::never();

    /// Time at which last treatment was recieved (for second-case considerations).
    SimTime timeLastTreatment = sim::never();

    /// Total parasite density at last time step (used during a bout).
    double previousDensity;
};

}
}
#endif
