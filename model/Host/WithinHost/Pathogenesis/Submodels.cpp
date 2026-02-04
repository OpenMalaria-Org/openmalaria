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

#include "Host/WithinHost/Pathogenesis/Submodels.h"
#include "Parameters.h"
#include "mon/reporting.h"

#include <cmath>
using namespace std;

namespace OM {
namespace WithinHost {
namespace Pathogenesis {

// ———  Müller presentation model  ———

// Müller model parameters
double rateMultiplier_31;
double densityExponent_32;

void MuellerPathogenesis::init( const Parameters& parameters ){
    rateMultiplier_31 = parameters[Parameter::MUELLER_RATE_MULTIPLIER]
        * sim::yearsPerStep();
    densityExponent_32 = parameters[Parameter::MUELLER_DENSITY_EXPONENT];
}

double MuellerPathogenesis::getPEpisode(double, double totalDensity) {
    double incidenceDensity = rateMultiplier_31
            * (pow(totalDensity, densityExponent_32));
    return 1.0 - exp(-incidenceDensity);
}


// ———  Pyrogenic threshold model  ———

//Pyrogenic threshold at birth (Y*0)
double initPyroThres;
// Ystar2: critical value in determining increase in pyrogenic threshold
double Ystar2_13;
double Ystar1_26;

//Number of categories in the numerical approximation used in updatePyrogenThres()
const size_t n = 11;
// Derived parameters without good names:
double a = numeric_limits<double>::signaling_NaN(),
    b = numeric_limits<double>::signaling_NaN();

void PyrogenPathogenesis::init( const Parameters& parameters ){
    initPyroThres = parameters[Parameter::Y_STAR_0];
    
    double delt = 1.0 / n;
    double smuY = -log(0.5) /
        (sim::stepsPerYear() * parameters[Parameter::Y_STAR_HALF_LIFE]);
    b = -smuY * delt;
    
    Ystar2_13 = parameters[Parameter::Y_STAR_SQ];
    
    //alpha: factor determining increase in pyrogenic threshold
    double alpha14 = parameters[Parameter::ALPHA];
    a = alpha14 * sim::oneTS() * delt;  
    
    //Ystar1: critical value of parasite density in determing increase in pyrog t
    Ystar1_26 = parameters[Parameter::Y_STAR_1];
}

PyrogenPathogenesis::PyrogenPathogenesis(double cF) :
     PathogenesisModel (cF), _pyrogenThres (initPyroThres)
{}


double PyrogenPathogenesis::getPEpisode(double timeStepMaxDensity, double totalDensity) {
    updatePyrogenThres(totalDensity);
    return timeStepMaxDensity / (timeStepMaxDensity + _pyrogenThres);
}

void PyrogenPathogenesis::summarize (const Host::Human& human) {
    mon::reportStatMHF( mon::MHF_PYROGENIC_THRESHOLD, human, _pyrogenThres );
    mon::reportStatMHF( mon::MHF_LOG_PYROGENIC_THRESHOLD, human, log(_pyrogenThres+1.0) );
}

void PyrogenPathogenesis::updatePyrogenThres(double totalDensity){
    // Note: this calculation is slow (something like 5% of runtime)
    
    //Numerical approximation to equation 2, AJTMH p.57
    for( size_t i = 1; i <= n; ++i ){
        _pyrogenThres += totalDensity * a /
            ( (Ystar1_26 + totalDensity) * (Ystar2_13 + _pyrogenThres) )
            + b * _pyrogenThres;
    }
}


void PyrogenPathogenesis::checkpoint (istream& stream) {
    PathogenesisModel::checkpoint (stream);
    _pyrogenThres & stream;
}
void PyrogenPathogenesis::checkpoint (ostream& stream) {
    PathogenesisModel::checkpoint (stream);
    _pyrogenThres & stream;
}


// ———  Predetermined episodes model  ———

double PredetPathogenesis::getPEpisode(double timeStepMaxDensity, double totalDensity) {
    updatePyrogenThres(totalDensity);
    return (timeStepMaxDensity > _pyrogenThres) ? 1 : 0;
}

} } }
