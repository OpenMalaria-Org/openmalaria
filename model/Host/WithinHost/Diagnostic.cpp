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

#include "Host/WithinHost/Diagnostic.h"
#include "Parameters.h"
#include "util/random.h"
#include "util/errors.h"
#include "util/ModelOptions.h"
#include "util/CommandLine.h"
#include "schema/scenario.h"
#include <limits>

namespace OM { namespace WithinHost {

// ———  Diagnostic (non-static)  ———

Diagnostic::Diagnostic( const Parameters& parameters, const scnXml::Diagnostic& elt ){
    dens_lim = numeric_limits<double>::quiet_NaN();
    specificity = numeric_limits<double>::quiet_NaN();
    
    if( util::ModelOptions::option( util::VIVAX_SIMPLE_MODEL ) ){
        return; // both parameters NaN (set above)
    }
    
    if( elt.getDeterministic().present() ){
        dens_lim = elt.getDeterministic().get().getMinDensity();
    }else if( elt.getStochastic().present() ){
        dens_lim = elt.getStochastic().get().getDens_50();
        if( dens_lim == 0.0 ){
            // The equation used for stochastic diagnostics breaks down when
            // dens=dens_lim=0 and for other cases the deterministic model is
            // the same when dens_lim=0.
        }else{
            specificity = elt.getStochastic().get().getSpecificity();
            if( specificity < 0.0 || specificity > 1.0 ){
                throw util::xml_scenario_error(
                    string("diagnostics/diagnostic(").append(elt.getName())
                    .append("): specificity must be in range [0,1]") );
            }
        }
    }else{
        // This should be impossible since according to schema one of these
        // elements must be present.
        throw SWITCH_DEFAULT_EXCEPTION;
    }
    if( dens_lim < 0.0 ){
        throw util::xml_scenario_error(
            string("diagnostics/diagnostic(").append(elt.getName())
            .append("): must have density ≥ 0") );
    }
    
    // We use a bias factor to adjust the "units" used to specify the density
    // of this diagnostic, since estimates from Garki and the standard
    // non-Garki sources are not equivalent to those from the Malariatherapy
    // data (which is used internally).
    if( !elt.getUnits().present() ){
        if( util::ModelOptions::option(util::GARKI_DENSITY_BIAS) ){
            // User must be explicit in this case, because presumably the Garki
            // bias is to be used for some diagnostics but likely not all
            // (e.g. neonatal mortality).
            throw util::xml_scenario_error( "diagnostics/diagnostic(*)/units: must specify this attribute when GARKI_DENSITY_BIAS is set" );
        }
        // otherwise we assume "Other"
        dens_lim *= parameters[Parameter::DENSITY_BIAS_NON_GARKI];
    }else if( elt.getUnits().get() == "Other" ){
        dens_lim *= parameters[Parameter::DENSITY_BIAS_NON_GARKI];
    }else if( elt.getUnits().get() == "Garki" ){
        dens_lim *= parameters[Parameter::DENSITY_BIAS_GARKI];
    }else{
        assert( elt.getUnits().get() == "Malariatherapy" );
        // in this case we don't need to use a bias factor
    }
    
    uses_hrp2 = elt.getMechanism() == "HRP2";
}

Diagnostic::Diagnostic(double minDens){
    specificity = numeric_limits<double>::quiet_NaN();
    dens_lim = minDens;
    uses_hrp2 = false;
}

bool Diagnostic::isPositive( LocalRng& rng, double dens, double densHRP2 ) const {
    if( uses_hrp2 ){
        assert( densHRP2 == densHRP2 ); // monitoring diagnostic passes NaN; use of HRP2 is not supported
        dens = densHRP2;
    }
    if( (std::isnan)(specificity) ){
        // use deterministic test
        return dens >= dens_lim;
    }else{
        // dens_lim is dens_50 in this case
        double pPositive = 1.0 + specificity * (dens / (dens + dens_lim) - 1.0);
//         double pPositive = (dens + dens_lim - dens_lim * specificity) / (dens + dens_lim);       // equivalent
        return rng.bernoulli( pPositive );
    }
}

bool Diagnostic::allowsFalsePositives() const{
    if( (std::isnan)(specificity) ) return dens_lim <= 0.0;
    return specificity < 1.0;
}


// ———  diagnostics (static)  ———

map<string, Diagnostic> diagnostic_set;
const Diagnostic* diagnostics::monitoring_diagnostic = 0;

void diagnostics::clear(){
    diagnostic_set.clear();
}

void diagnostics::init( const Parameters& parameters, const scnXml::Scenario& scenario ){
    diagnostic_set.clear(); // compatibility with unit tests
    
    if(scenario.getDiagnostics().present()){
        for( const scnXml::Diagnostic& diagnostic : scenario.getDiagnostics().get().getDiagnostic() ){
            string name = diagnostic.getName();     // conversion fails without this extra line
            bool inserted = diagnostic_set.insert( make_pair(std::move(name), Diagnostic(parameters, diagnostic)) ).second;
            if( !inserted ){
                throw util::xml_scenario_error( string("diagnostic with this name already set: ").append(diagnostic.getName()) );
            }
        }
    }

    const scnXml::Surveys& surveys = scenario.getMonitoring().getSurveys();
    if( util::ModelOptions::option( util::VIVAX_SIMPLE_MODEL ) ){
        // So far the implemented Vivax code does not produce parasite
        // densities, thus this diagnostic model cannot be used.
        diagnostics::monitoring_diagnostic = &diagnostics::make_deterministic( numeric_limits<double>::quiet_NaN() );
    }else if( surveys.getDetectionLimit().present() ){
        if( surveys.getDiagnostic().present() ){
            throw util::xml_scenario_error( "monitoring/surveys: do not "
                "specify both detectionLimit and diagnostic" );
        }
        if( util::CommandLine::option( util::CommandLine::DEPRECATION_WARNINGS ) ){
            std::cerr << "Deprecation warning: monitoring/surveys: "
                "specification of \"diagnostic\" is suggested over \"detectionLimit\"" << std::endl;
        }

        // This controls whether the detection limit is specified relative to
        // the Garki or other methods.
        double densitybias = numeric_limits<double>::quiet_NaN();
        if (util::ModelOptions::option (util::GARKI_DENSITY_BIAS)) {
            densitybias = parameters[Parameter::DENSITY_BIAS_GARKI];
        } else {
            if( scenario.getAnalysisNo().present() ){
                int analysisNo = scenario.getAnalysisNo().get();
                if ((analysisNo >= 22) && (analysisNo <= 30)) {
                    cerr << "Warning: these analysis numbers used to mean "
                        "use Garki density bias. If you do want to use this, "
                        "specify the option GARKI_DENSITY_BIAS; if not, nothing's wrong." << endl;
                }
            }
            densitybias = parameters[Parameter::DENSITY_BIAS_NON_GARKI];
        }
        double detectionLimit = surveys.getDetectionLimit().get() * densitybias;
        diagnostics::monitoring_diagnostic = &diagnostics::make_deterministic( detectionLimit );
    }else{
        if( !surveys.getDiagnostic().present() ){
            throw util::xml_scenario_error( "monitoring/surveys: require "
                "either detectionLimit or diagnostic" );
        }
        if( util::ModelOptions::option(util::GARKI_DENSITY_BIAS) ){
            throw util::xml_scenario_error( "Use of GARKI_DENSITY_BIAS is not "
                "appropriate when monitoring/surveys/diagnostic is used." );
        }
        diagnostics::monitoring_diagnostic = &diagnostics::get( surveys.getDiagnostic().get() );
    }
    
    if( diagnostics::monitoring_diagnostic->uses_hrp2 ){
        throw util::xml_scenario_error( "the diagnostic used for "
            "monitoring may not use HRP2 as its mechanism" );
    }
}

const Diagnostic& diagnostics::get( const string& name ){
    auto it = diagnostic_set.find(name);
    if( it == diagnostic_set.end() ){
        throw util::xml_scenario_error( string("diagnostic not found: ").append(name) );
    }
    return it->second;
}

const Diagnostic& diagnostics::make_deterministic(double minDens){
    string name = "";   // anything not matching an existing name is fine
    assert( diagnostic_set.count(name) == 0 );  // we shouldn't need to call make_deterministic more than once, so I think this will do
    auto ret = diagnostic_set.insert( make_pair(std::move(name), Diagnostic( minDens )) );
    assert(ret.second);
    return ret.first->second;
}

} }
