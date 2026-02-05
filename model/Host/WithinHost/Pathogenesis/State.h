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

#ifndef MO_PATHOGENISIS_STATE
#define MO_PATHOGENISIS_STATE

namespace OM { namespace WithinHost { namespace Pathogenesis {
  /** Types of sickness; used by case management.
   *
   * Most values are flags which can be combined in any form. A few
   * combinations set follow. */
  enum State {
    /* Values here are written in hexadecimal: http://en.wikipedia.org/wiki/Hexadecimal
     * Many are designed to be "flags", so the value corresponds to a single bit:
     * http://en.wikipedia.org/wiki/Flag_byte
     * Max: 0x80000000
     * (note & | ^ are C++'s binary AND, OR and XOR operators). */
    NONE		= 0,		///< Not sick
    
    // Flags for current state/worst state to report:
    SICK		= 0x1,		///< Sick (may or may not be from malaria)
    MALARIA		= 0x2,		///< Malaria sickness
    SEVERE		= 0x8,		///< Severe malaria case
    COINFECTION		= 0x4,		///< Malaria with a coinfection
    COMPLICATED		= 0x20,		///< Flag used to indicate SEVERE and/or COINFECTION
    
    // Flags for outcome reporting:
    STATE_NMF           = SICK, ///< Alias: SICK
    //STATE_NMF_NEED_AB    = SICK | NEED_ANTIBIOTIC,       ///< Combination: SICK, NEED_ANTIBIOTIC
    STATE_MALARIA	= SICK | MALARIA,	///< Combination: SICK, MALARIA
    STATE_SEVERE	= STATE_MALARIA | COMPLICATED | SEVERE,	///< Combination: SICK, MALARIA, COMPLICATED, SEVERE
    STATE_COINFECTION	= STATE_MALARIA | COMPLICATED | COINFECTION,	///< Combination: SICK, MALARIA, COMPLICATED, COINFECTION
  };
  struct StatePair{
      StatePair() : state(NONE), indirectMortality(false) {}
      State state;
      bool indirectMortality;   // implies death a little later on which may not be seen as a malaria death
  };

} } }
#endif
