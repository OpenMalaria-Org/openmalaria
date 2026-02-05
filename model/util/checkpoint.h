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

#ifndef OM_util_checkpoint
#define OM_util_checkpoint

#ifndef Hmod_Global
#error "Please include Global.h not checkpoint.h directly."
// otherwise "using ..." declaration in Global.h won't work
#endif

#include <iostream>
#include <vector>
#include <list>
#include <map>
#include <set>
using namespace std;

namespace OM {
namespace interventions{
    struct ComponentId;
}
namespace util {
namespace checkpoint {
    
    const long DEFAULT_MAX_LENGTH = 2000;
    
    ///@brief Utility functions
    //@{
    /** Perform important checks on checkpoint file format.
     *
     * Call at beginning of read/write. */
    void header (ostream& stream);
    void header (istream& stream);	///< ditto
    
  /** For checking the number of elements in a list is sensible when loading
   * checkpoints.
   *
   * While loading checkpoints, if memory is allocated in a loop using a
   * variable loaded from the checkpoint file, it should be checked that length
   * is sensible. This function is a convenient way to do that (although the
   * limit isn't large enough to account for population size).
   * 
   * The reason for this check is that if the checkpoint is read wrongly,
   * variables often get completely wrong values which can cause memory
   * allocation to grind the computer to a halt. In this case the values are
   * usually bad enough that a lenient check will still catch them. */
    void validateListSize (long length, long max = DEFAULT_MAX_LENGTH);
    //@}
    
    ///@brief Operator& for simple data-types
    //@{
    void operator& (bool x, ostream& stream);
    void operator& (bool& x, istream& stream);
    
    void operator& (signed char x, ostream& stream);
    void operator& (signed char& x, istream& stream);
    
    void operator& (short x, ostream& stream);
    void operator& (short& x, istream& stream);
    
    void operator& (int x, ostream& stream);
    void operator& (int& x, istream& stream);
    
    void operator& (long x, ostream& stream);
    void operator& (long& x, istream& stream);
    
    void operator& (long long x, ostream& stream);
    void operator& (long long& x, istream& stream);
    
    void operator& (unsigned char x, ostream& stream);
    void operator& (unsigned char& x, istream& stream);
    
    void operator& (unsigned short x, ostream& stream);
    void operator& (unsigned short& x, istream& stream);
    
    void operator& (unsigned int x, ostream& stream);
    void operator& (unsigned int& x, istream& stream);
    
    void operator& (unsigned long x, ostream& stream);
    void operator& (unsigned long& x, istream& stream);
    
    void operator& (unsigned long long x, ostream& stream);
    void operator& (unsigned long long& x, istream& stream);
    
    void operator& (float x, ostream& stream);
    void operator& (float& x, istream& stream);
    
    void operator& (double x, ostream& stream);
    void operator& (double& x, istream& stream);
    
    void operator& (long double x, ostream& stream);
    void operator& (long double& x, istream& stream);
    //@}
    
    /** @brief Operator& for pointers
     *
     * These could be implemented for pointers but require reading objects in the opposite order to
     * writing them.
     * 
     * Matching up pointers with addresses of objects not stored in references could be more
     * challanging though.
     * 
     * Can't really deal with referenced objects properly. */
    //@{
    //@}
    
    /** @brief Operator& for std::string */
    //@{
    void operator& (string x, ostream& stream);
    void operator& (string& x, istream& stream);
    //@}
    
} } }	// end of namespaces
#endif
