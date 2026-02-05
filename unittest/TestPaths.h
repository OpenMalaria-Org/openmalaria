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

#ifndef Hmod_TestPaths
#define Hmod_TestPaths
// If this file is included unconfigured, this will fail.
// When configured, @CMAKE_CONFIGURED@ is replaced by 1.
#if @CMAKE_CONFIGURED@ != 1
assert (false);
#endif

const char* UnittestSourceDir = "@CMAKE_CURRENT_SOURCE_DIR@/";	// must end with '/'
const char* UnittestScenario = "@CMAKE_CURRENT_BINARY_DIR@/configured/scenario.xml";
#endif
