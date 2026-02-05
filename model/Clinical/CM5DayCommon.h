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

#ifndef Hmod_CM5DayCommon
#define Hmod_CM5DayCommon

#include "Host/WithinHost/Pathogenesis/State.h"
#include "Clinical/ClinicalModel.h"
#include "Host/Human.h"
#include "Host/WithinHost/WHInterface.h"
#include "mon/reporting.h"

namespace OM {
namespace Clinical {
using Host::Human;

/**
 * Common parts of 5-day case management models.
 */
class CM5DayCommon : public ClinicalModel
{
public:
    static void init();
    
    /** Construct an instance (per-human).
     * 
     * @param tSF Treatment-seeking factor. Normally 1 but allows simple heterogeneity.
     */
    CM5DayCommon (double tSF);
    
    virtual bool isExistingCase () {
        // If treated in the recent past:
        return sim::now() > m_tLastTreatment && sim::now() <= m_tLastTreatment + healthSystemMemory;
    }
    
protected:
    enum CaseType { FirstLine, SecondLine, NumCaseTypes };
    static mon::Measure measures[NumCaseTypes];
    static double accessUCAny[NumCaseTypes];
    static double accessUCSelfTreat[NumCaseTypes];
    static double accessSevere;
    static double cureRateSevere;
    static WithinHost::TreatmentId treatmentSevere;
    
    virtual void doClinicalUpdate (Human& human, double ageYears, WithinHost::Pathogenesis::StatePair &pg);

    virtual void checkpoint (istream& stream);
    virtual void checkpoint (ostream& stream);

    /** Called when a non-severe/complicated malaria sickness occurs. */
    virtual void uncomplicatedEvent(Human& human, Episode::State pgState) =0;

    //! treatment seeking for heterogeneity
    double m_treatmentSeekingFactor;
    
private:
    /** Called when a severe/complicated (with co-infection) malaria sickness occurs.
     *
     * Note: sets doomed = 4 if patient dies. */
    void severeMalaria(Human& human, Episode::State pgState, double ageYears, int& doomed);
};

}
}
#endif
