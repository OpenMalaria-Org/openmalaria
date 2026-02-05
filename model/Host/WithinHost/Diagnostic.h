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

#ifndef HS_DIAGNOSTIC
#define HS_DIAGNOSTIC

#include "Global.h"
#include "util/random.h"

class UnittestUtil;
namespace scnXml {
    class Diagnostics;
    class Diagnostic;
}

namespace OM {
    class Parameters;
namespace WithinHost {

using util::LocalRng;

class Diagnostic {
public:
    /** Use the test.
     * 
     * @param dens Current parasite density in parasites per µL
     * @param densHRP2 Equivalent density for diagnostics dependent on HRP2
     * @returns True if outcome is positive. */
    bool isPositive( LocalRng& rng, double dens, double densHRP2 ) const;
    
    inline bool operator!=( const Diagnostic& that )const{
        return specificity != that.specificity ||
            dens_lim != that.dens_lim;
    }
    
    /// True if false positives are possible
    bool allowsFalsePositives() const;
    
private:
    /** Construct from XML parameters. */
    Diagnostic( const Parameters& parameters, const scnXml::Diagnostic& elt );
    /** Construct as deterministic. */
    explicit Diagnostic( double minDens );
    
    // switch: either not-a-number indicating a deterministic test, or specificity
    double specificity;
    // depending on model, this is either the minimum detectible density
    // or the density at which test has half a chance of a positive outcome
    double dens_lim;
    
    // if true, the diagnostic is dependent on the HRP2 antigen
    bool uses_hrp2;
    
    friend class diagnostics;
};

/** Static members to do with Diagnostic: library of parameterised Diagnostic
 * objects. */
class diagnostics {
public:
    /** Initialise from input data. */
    static void init( const Parameters& parameters, const scnXml::Scenario& scenario );
    
    /** Look up a diagnostic by name and get a reference to it.
     * 
     * @throws util::xml_scenario_error on failure
     */
    static const Diagnostic& get( const std::string& name );
    
    /** Make a new diagnostic with deterministic density and return a reference. */
    static const Diagnostic& make_deterministic( double minDens );

    /// Static access functions

    static inline const Diagnostic& monitoringDiagnostic(){
        assert( monitoring_diagnostic != 0 );
        return *monitoring_diagnostic;
    }

    
private:
    static void clear();        // for unit tests
    static const Diagnostic* monitoring_diagnostic;
    friend class ::UnittestUtil;
};

} }
#endif
