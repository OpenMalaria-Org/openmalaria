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

#include "Clinical/ESCaseManagement.h"
#include "Clinical/CMDecisionTree.h"
#include "Clinical/EventScheduler.h"
#include "util/errors.h"
#include "util/ModelOptions.h"

#include <set>
#include <sstream>

namespace OM { namespace Clinical {
    using namespace OM::util;

// -----  ESCaseManagement  -----

const CMDecisionTree *escm_uncomplicated = 0,
    *escm_complicated = 0;

void ESCaseManagement::setHealthSystem(const scnXml::HSEventScheduler& esData){
    escm_uncomplicated = &CMDecisionTree::create( esData.getUncomplicated(), true );
    escm_complicated = &CMDecisionTree::create( esData.getComplicated(), false );
    
    // Calling our parent class like this is messy. Changing this would require
    // moving change-of-health-system handling into ClinicalModel.
    ClinicalEventScheduler::setParameters( esData );
}

CMDTOut ESCaseManagement::execute (
        const CMHostData& hostData
) {
    assert (hostData.pgState & Episode::SICK);
    //TODO: should we remove any existing prescriptions?
    //TODO: Note that these trees do both "access" and "case management" decisions.
    //TODO: medicateQueue.clear();
    
    const CMDecisionTree& tree = (hostData.pgState & Episode::COMPLICATED) ?
        *escm_complicated : *escm_uncomplicated;
    return tree.exec( hostData );
}

} }
