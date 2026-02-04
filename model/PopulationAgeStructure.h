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

#ifndef Hmod_PopulationAgeStructure
#define Hmod_PopulationAgeStructure

#include <vector>
#include "Global.h"

namespace scnXml{ class Demography; }
namespace OM
{
    /** Encapsulates code just setting up the age structure (i.e.  cumAgeProp). */
    class AgeStructure
    {
    public:
        /** Set up cumAgeProp from XML data. */
        static void init( const scnXml::Demography& demography );
        
        /** Return maximum individual lifetime in intervals that AgeStructure can handle. */
        static inline size_t getMaxTStepsPerLife(){
            return cumAgeProp.size();
        }
        
        /** Return the expected population size of individuals aged ageTSteps or
        * older, based on a total population size of targetPop. */
        static int targetCumPop( size_t ageTSteps, int targetPop );
    
    private:
        /*! Estimates demography parameters to define a smooth curve for the target
        population age-distribution (age in years) */
        static void estimateRemovalRates( const scnXml::Demography& demography );
	
        /** For input values for alpha1 and mu1, the fit to field data (residualSS)
        * is calculated and returned function called iteratively by
        * estimateRemovalRates. */
        static double setDemoParameters (double param1, double param2);

        /** Takes the best-fitting demography parameters estimated by
        * estimateRemovalRates and calculates the age structure (cumAgeProp). */
	static void calcCumAgeProp ();
        
	//BEGIN static parameters only used by estimateRemovalRates(), setDemoParameters() and calcCumAgeProp()
	/** The bounds for each age group and percentage of population in this age
        * group for the field data demography age groups.
        *
        * ageGroupBounds[i] is the lower-bound for group i, ageGroupBounds[i+1] is
        * the group's upper bound. ageGroupPercent[i] is the percentage of the
        * population in age group i.
        *
        * Set by estimateRemovalRates() and used internally (by
        * setDemoParameters()). */
        //@{
        static vector<double> ageGroupBounds;
        static vector<double> ageGroupPercent;
        //@}
        /** Parameters defining smooth curve of target age-distribution.
        *
        * Set by estimateRemovalRates() (via setDemoParameters()) and used by
        * setupPyramid(). */
        //@{
        static double mu0;
        static double mu1;
        static double alpha0;
        static double alpha1;
        // rho is growth rate of human population.
        constexpr static double rho = 0.0;
        //@}
	//END
	
	//BEGIN static parameters set by init() or calcCumAgeProp()
        /** Target cumulative percentage of population by age, from oldest age to youngest.
	*
	* cumAgeProp[size()+1-i] gives the proportion of people aged i time steps or older.
	*/
	static std::vector<double> cumAgeProp;
	//END
    };
}

#endif
