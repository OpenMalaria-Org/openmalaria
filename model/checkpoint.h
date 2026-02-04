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

#ifndef Hmod_GeneralCheckpoint_H
#define Hmod_GeneralCheckpoint_H

#include "Global.h"
#include "Population.h"
#include "Transmission/transmission.h"

namespace OM
{
    /** @brief checkpointing functions
    *
    * readCheckpoint/writeCheckpoint prepare to read/write the file,
    * and read/write read and write the actual data. */
    void writeCheckpoint(const bool startedFromCheckpoint, const string &checkpointFileName, SimTime &endTime, SimTime &estEndTime, Population &population, Transmission::TransmissionModel &transmission);

    void readCheckpoint(const string &checkpointFileName, SimTime &endTime, SimTime &estEndTime, Population &population, Transmission::TransmissionModel &transmission);
}

#endif