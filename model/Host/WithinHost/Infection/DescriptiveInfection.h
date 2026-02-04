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

#ifndef Hmod_DescriptiveInfection
#define Hmod_DescriptiveInfection

#include "Host/WithinHost/Infection/Infection.h"
#include "util/random.h"

namespace OM { namespace WithinHost {

using util::LocalRng;

// Max duration of sample data for an infection in intervals. Simulated
// infections may last longer; in this case the last sample data is repeated
// until the infection terminates.
const int32_t maxDurationTS = 83;
const int32_t numDurations = 84;    // Above plus one (include a category for 0)

//The maximum parasite density we allow per DescriptiveInfection. Higher values are set to maxDens.
// NOTE: should this be enforced before or after attenuation (BSV, immunity, IPT)? Was enforced in both places.
const double maxDens = 2000000;


/** A model of Plasmodium Falciparum infections. Described in AJTMH 75(2)
 * pp19-31.
 * 
 * This model was designed primarily for usage with a 5-day time-step, but is
 * mostly applicable to 1-4 day time-steps too. In such cases the indexes used
 * to access meanLogParasiteCount (or the data contained) would need adjusting.
 * 
 * Note that this class models only a single infection; for the associated
 * handling of multiple infections see the DescriptiveWithinHostModel class.
 */
class DescriptiveInfection : public Infection {
public:
    ///@name Static init/cleanup
    //@{
    /**
     * Loads some constants: parameters used by the empirical models.
     * 
     * An old comment said the following, but due to code changes may not be completely accurate now:
     * Init constants common to all Phase A (AJTMH 75(2)) infections (i.e.
     * those models documented in the referenced publication).
     */
    static void init(const OM::Parameters& parameters);
    //@}
    
    ///@name Constructors
    //@{
    /// Default constructor
    DescriptiveInfection (LocalRng& rng, uint32_t genotype, int origin);
    /// Checkpoint loading constructor
    DescriptiveInfection (istream& stream);
    //@}
    
    /** Returns true when age reaches the pre-determined duration (i.e. when
    * this infection terminates). */
    bool expired () {
        return sim::ts0() > m_startDate + m_duration;
    }
    
    /** Determines parasite density of an individual infection (5-day time step
     * update)
     *
     * @param ageInYears Age (of human)
     * @param cumulativeh Cumulative number of infections
     * @param cumulativeY Previous exposure (cumulative parasite density)
     * @param timeStepMaxDensity (In-out param) Used to return the maximum
     *  parasite density over a 5-day interval.
     * @param innateImmSurvFact Density multiplier for innate immunity.
     * @param bsvFactor Density multiplier for Blood-Stage Vaccine effect.
     */
    void determineDensities(
            LocalRng& rng,
            double cumulativeh,
            double &timeStepMaxDensity,
            double immSurvFact,
            double innateImmSurvFact,
            double bsvFactor);
    
    /** Decide on an infection duration and return it.
     * 
     * Parameters for this model are hard-coded.
     * 
     * Description: determines infection duration by sampling from the log
     * normal distribution using parameters for 53 patients from Georgia.
     * Mean log duration of an infection values from AJTM p.9 eq.5.
     */
    static SimTime infectionDuration(LocalRng& rng);
    
    /// Includes the effect of attenuated infections by SP concentrations, when using IPT
    virtual void IPTattenuateAsexualDensity () {}
    
protected:
    virtual void checkpoint (ostream& stream);
    
    // Arbitrary predetermined maximum duration of the infection
    SimTime m_duration = sim::never(); 
    
    bool notPrintedMDWarning;
    
private:
    /// @brief Static parameters set by init()
    //@{
    /* A triangular matrix: meanLogParasiteCount[i][j] is the
     * Mean Log Parasite Count for age i (in time steps) of an infection which
     * lasts j days. Indices with i>j are unused. */
    static double meanLogParasiteCount[numDurations][numDurations];
    
    /// Sigma0^2 from AJTM p.9 eq. 13
    static double sigma0sq;
    /// XNuStar in AJTM p.9 eq. 13
    static double xNuStar;
    //@}
};

} }
#endif
