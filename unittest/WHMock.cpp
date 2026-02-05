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

#include "WHMock.h"
#include "Host/WithinHost/Diagnostic.h"
#include "util/errors.h"

using namespace std;

namespace OM {
namespace UnitTest {

WHMock::WHMock() :
    totalDensity(numeric_limits<double>::quiet_NaN()),
    nTreatments(0),
    lastTimeLiver(sim::never()), lastTimeBlood(sim::never())
{}
WHMock::~WHMock() {}

double WHMock::probTransmissionToMosquito(vector<double> &, vector<double> &) const{
    throw util::unimplemented_exception( "not needed in unit test" );
}

bool WHMock::summarize(Host::Human& human)const{
    throw util::unimplemented_exception( "not needed in unit test" );
}

void WHMock::importInfection(LocalRng& rng){
    throw util::unimplemented_exception( "not needed in unit test" );
}

void WHMock::optionalPqTreatment(Host::Human& human){
    throw util::unimplemented_exception( "not needed in unit test" );
}

bool WHMock::treatSimple( Host::Human& human, SimTime timeLiver, SimTime timeBlood ){
    nTreatments += 1;
    lastTimeLiver = timeLiver;
    lastTimeBlood = timeBlood;
    return timeBlood != sim::zero();
}

void WHMock::treatPkPd(size_t schedule, size_t dosages, double age, double delay_d){
    nTreatments += 1;
    pkpd.prescribe( schedule, dosages, age, numeric_limits<double>::quiet_NaN(), delay_d );
}

void WHMock::update(Host::Human &human, LocalRng& rng, int &nNewInfs_i, int &nNewInfs_l, vector<double>& genotype_weights_i, vector<double>& genotype_weights_l, double ageInYears) {
    throw util::unimplemented_exception( "not needed in unit test" );
}

inline double WHMock::getTotalDensity() const{
    return totalDensity;
}

bool WHMock::diagnosticResult( LocalRng& rng, const Diagnostic& diagnostic ) const{
    return diagnostic.isPositive( rng, totalDensity, numeric_limits<double>::quiet_NaN() );
}

void WHMock::treatment( Host::Human& human, TreatmentId treatId ){
    nTreatments += 1;
}

Pathogenesis::StatePair WHMock::determineMorbidity( Host::Human& human, double ageYears, bool isDoomed ){
    throw util::unimplemented_exception( "not needed in unit test" );
}

void WHMock::clearImmunity(){
    throw util::unimplemented_exception( "not needed in unit test" );
}

double WHMock::getCumulative_h() const{
    throw util::unimplemented_exception( "not needed in unit test" );
}

double WHMock::getCumulative_Y() const{
    throw util::unimplemented_exception( "not needed in unit test" );
}

InfectionOrigin WHMock::getInfectionOrigin() const{
    return InfectionOrigin::Indigenous;
}

void WHMock::checkpoint (istream& stream){
    throw util::unimplemented_exception( "not needed in unit test" );
}
void WHMock::checkpoint (ostream& stream){
    throw util::unimplemented_exception( "not needed in unit test" );
}

}
}
