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

#ifndef Hmod_ImmediateOutcomes
#define Hmod_ImmediateOutcomes

#include "Clinical/CM5DayCommon.h"

namespace OM {
namespace Clinical {

/**
 * This models case management at a 5-day time step with all-or-nothing
 * treatment.
 * 
 * Uncomplicated cases: access, otherwise known as "seeking any type of
 * treatment", is determined by a fixed-function decision, which may be
 * modified by a treatment-seeking factor. Treatment decisions are also fixed
 * function.
 * 
 * Severe cases: all decisions and outcomes are calculated via a fixed-function
 * probability tree, using the same logic for handling severe cases as has long
 * been used.
 * 
 * NOTE: this provides a strict subset of the functionality of DecisionTree5Day.
 */
class ImmediateOutcomes : public CM5DayCommon
{
public:
    /** Load health system data from initial data or an intervention's data (both from XML).
     * (Re)loads all data affected by this healthSystem element. */
    static void setHealthSystem (const scnXml::HSImmediateOutcomes& hsioData);
    
    ImmediateOutcomes (double tSF) : CM5DayCommon(tSF) {}
    
protected:
    /** Called when a non-severe/complicated malaria sickness occurs. */
    virtual void uncomplicatedEvent(Human& human, Episode::State pgState);
    
private:
    static double cureRateUCOfficial[CM5DayCommon::NumCaseTypes];
    static double cureRateUCSelfTreat[CM5DayCommon::NumCaseTypes];
    static WithinHost::TreatmentId treatmentUC[CM5DayCommon::NumCaseTypes];
    static bool useDiagnosticUC;
};

}
}
#endif
