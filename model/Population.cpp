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

#include "Population.h"
#include "mon/Continuous.h"

#include "Host/Human.h"
#include "Host/NeonatalMortality.h"
#include "Host/WithinHost/WHInterface.h"
#include "Host/WithinHost/Genotypes.h"
#include "Host/WithinHost/Diagnostic.h"
#include "Clinical/ClinicalModel.h"
#include "Transmission/TransmissionModel.h"

#include "util/errors.h"
#include "util/random.h"
#include "util/ModelOptions.h"
#include "util/StreamValidator.h"
#include <schema/scenario.h>

#include <cmath>

namespace OM
{
    using namespace std::placeholders;
    using namespace OM::util;
    using Transmission::TransmissionModel;

/** Variables for continuous reporting */
const vector<double> ctsDemogAgeGroups = { 1.0, 5.0, 10.0, 15.0, 25.0 };

// -----  non-static methods: creation/destruction, checkpointing  -----

Population::Population(size_t size) : size (size) {}

void Population::createInitialHumans()
{
    /* We create a whole population here, regardless of whether humans can
    survive until start of vector init (vector model needs a whole population
    structure in any case). However, we don't update humans known not to survive
    until vector init, which saves computation and memory (no infections). */
    
    int cumulativePop = 0;
    for(size_t iage_prev = AgeStructure::getMaxTStepsPerLife(), iage = iage_prev - 1;
         iage_prev > 0; iage_prev = iage, iage -= 1 )
    {
        int targetPop = AgeStructure::targetCumPop( iage, size );
        while (cumulativePop < targetPop) {
            SimTime dob = sim::zero() - sim::fromTS(iage);
            util::streamValidate( dob );
            humans.push_back( Host::Human (dob) );
            ++cumulativePop;
        }
    }
    
    // Vector setup dependant on human population structure (we *want* to
    // include all humans, whether they'll survive to vector init phase or not).
    assert( sim::now() == sim::zero() );      // assumed below
}

void Population::update()
{
    //NOTE: other parts of code are not set up to handle changing population size. Also
    // size is assumed to be the _actual and exact_ population size by other code.
    int cumPop = 0;

    for (auto it = humans.begin(); it != humans.end();) {
        bool isDead = it->isDead();

        // if (Actual number of people so far > target population size for this age)
        // "outmigrate" some to maintain population shape
        //NOTE: better to use age(sim::ts0())? Possibly, but the difference will not be very significant.
        // Also see targetPop = ... comment above
        bool outmigrate = cumPop >= AgeStructure::targetCumPop(sim::inSteps(it->age(sim::ts1())), size);
        
        if( isDead || outmigrate ){
            it = humans.erase (it);
            continue;
        }
        ++cumPop;
        ++it;
    } // end of per-human updates

    // increase population size to targetPop
    recentBirths += (size - cumPop);
    while (cumPop < (int)size) {
        // humans born at end of this time step = beginning of next, hence ts1
        humans.push_back( Host::Human (sim::ts1()) );
        ++cumPop;
    }
}

void Population::checkpoint(istream& stream)
{
    size & stream;
    recentBirths & stream;

    for(size_t i = 0; i < size && !stream.eof(); ++i) {
        humans.push_back( Host::Human (sim::zero()) );
        humans.back().checkpoint(stream);
    }
    if (humans.size() != size)
        throw util::checkpoint_error("Population: out of data (read " + to_string(humans.size()) + " humans)");
}

void Population::checkpoint(ostream& stream)
{
    size & stream;
    recentBirths & stream;
    for(Host::Human &human : humans)
        human.checkpoint(stream);
}

void registerContinousPopulationCallbacks()
{
    ostringstream ctsDemogTitle;
    for( double ubound : ctsDemogAgeGroups )
        ctsDemogTitle << "\thost % ≤ " << ubound;

    mon::Continuous.registerCallback( "hosts", "\thosts", &ctsHosts );
    mon::Continuous.registerCallback( "host demography", ctsDemogTitle.str(), &ctsHostDemography);
    mon::Continuous.registerCallback( "recent births", "\trecent births", &ctsRecentBirths);
    mon::Continuous.registerCallback( "patent hosts", "\tpatent hosts", &ctsPatentHosts);
    mon::Continuous.registerCallback( "immunity h", "\timmunity h", &ctsImmunityh);
    mon::Continuous.registerCallback( "immunity Y", "\timmunity Y", &ctsImmunityY);
    mon::Continuous.registerCallback( "median immunity Y", "\tmedian immunity Y", &ctsMedianImmunityY);
    mon::Continuous.registerCallback( "human age availability", "\thuman age availability", &ctsMeanAgeAvailEffect);
    mon::Continuous.registerCallback( "ITN coverage", "\tITN coverage", &ctsITNCoverage);
    mon::Continuous.registerCallback( "IRS coverage", "\tIRS coverage", &ctsIRSCoverage);
    mon::Continuous.registerCallback( "GVI coverage", "\tGVI coverage", &ctsGVICoverage);
}

void ctsHosts (Population &population, ostream& stream){
    // this option is intended for debugging human initialization; normally this should equal size.
    stream << '\t' << population.getSize();
}

void ctsHostDemography (Population &population, ostream& stream){
    auto iter = population.humans.crbegin();
    int cumCount = 0;
    for( double ubound : ctsDemogAgeGroups ){
        while( iter != population.humans.crend() && sim::inYears(iter->age(sim::now())) < ubound ){
            ++cumCount;
            ++iter;
        }
        stream << '\t' << cumCount;
    }
}

void ctsRecentBirths (Population &population, ostream& stream){
    stream << '\t' << population.getRecentBirths();
    population.resetRecentBirths();
}

void ctsPatentHosts (Population &population, ostream& stream){
    int patent = 0;
    for(Host::Human &human : population.humans) {
        auto diag = WithinHost::diagnostics::monitoringDiagnostic();
        if( human.withinHostModel->diagnosticResult(human.rng, diag) )
            ++patent;
    }
    stream << '\t' << patent;
}

void ctsImmunityh (Population &population, ostream& stream){
    double x = 0.0;
    for(const Host::Human &human : population.humans) {
        x += human.withinHostModel->getCumulative_h();
    }
    x /= population.getSize();
    stream << '\t' << x;
}

void ctsImmunityY (Population &population, ostream& stream){
    double x = 0.0;
    for(const Host::Human &human : population.humans) {
        x += human.withinHostModel->getCumulative_Y();
    }
    x /= population.getSize();
    stream << '\t' << x;
}

void ctsMedianImmunityY (Population &population, ostream& stream){
    vector<double> list;
    list.reserve( population.getSize() );
    for(const Host::Human &human : population.humans) {
        list.push_back( human.withinHostModel->getCumulative_Y() );
    }
    sort( list.begin(), list.end() );
    double x;
    if( mod_nn(population.getSize(), 2) == 0 ){
        size_t i = population.getSize() / 2;
        x = (list[i-1]+list[i])/2.0;
    }else{
        x = list[population.getSize() / 2];
    }
    stream << '\t' << x;
}

void ctsMeanAgeAvailEffect (Population &population, ostream& stream){
    int nHumans = 0;
    double avail = 0.0;
    for(Host::Human &human : population.humans) {
        if( !human.perHostTransmission.outsideTransmission ){
            ++nHumans;
            avail += human.perHostTransmission.relativeAvailabilityAge(sim::inYears(human.age(sim::now())));
        }
    }
    stream << '\t' << avail/nHumans;
}

void ctsITNCoverage (Population &population, ostream& stream){
    int nActive = 0;
    for(const Host::Human &human : population.humans) {
        nActive += human.perHostTransmission.hasActiveInterv( interventions::Component::ITN );
    }
    double coverage = static_cast<double>(nActive) / population.getSize();
    stream << '\t' << coverage;
}

void ctsIRSCoverage (Population &population, ostream& stream){
    int nActive = 0;
    for(const Host::Human &human : population.humans) {
        nActive += human.perHostTransmission.hasActiveInterv( interventions::Component::IRS );
    }
    double coverage = static_cast<double>(nActive) / population.getSize();
    stream << '\t' << coverage;
}

void ctsGVICoverage (Population &population, ostream& stream){
    int nActive = 0;
    for(const Host::Human &human : population.humans) {
        nActive += human.perHostTransmission.hasActiveInterv( interventions::Component::GVI );
    }
    double coverage = static_cast<double>(nActive) / population.getSize();
    stream << '\t' << coverage;
}

}

