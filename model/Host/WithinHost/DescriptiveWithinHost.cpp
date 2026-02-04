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

#include "Global.h"
#include "Host/Human.h"
#include "Host/WithinHost/DescriptiveWithinHost.h"
#include "Host/WithinHost/Diagnostic.h"
#include "Host/WithinHost/Genotypes.h"
#include "Host/WithinHost/Pathogenesis/PathogenesisModel.h"
#include "util/ModelOptions.h"
#include "util/StreamValidator.h"
#include "util/errors.h"

#include <cassert>

using namespace std;

namespace OM {
namespace WithinHost {

extern bool bugfix_max_dens;    // DescriptiveInfection.cpp
bool reportPatentInfected = false;
// -----  Initialization  -----

void DescriptiveWithinHostModel::initDescriptive(){
    reportPatentInfected = (
                                mon::isUsedM(mon::MHR_PATENT_INFECTIONS) 
                                || mon::isUsedM(mon::MHR_PATENT_INFECTIONS_IMPORTED) 
                                || mon::isUsedM(mon::MHR_PATENT_INFECTIONS_INTRODUCED) 
                                || mon::isUsedM(mon::MHR_PATENT_INFECTIONS_INDIGENOUS) 
                            );
}

DescriptiveWithinHostModel::DescriptiveWithinHostModel( LocalRng& rng, double comorbidityFactor ) :
        WHFalciparum( rng, comorbidityFactor )
{
    assert( sim::oneTS() == sim::fromDays(5) );
    opt_vaccine_genotype = util::ModelOptions::option (util::VACCINE_GENOTYPE);
}

DescriptiveWithinHostModel::~DescriptiveWithinHostModel() {
    for( auto inf = infections.begin(); inf != infections.end(); ++inf ){
        delete *inf;
    }
    infections.clear();
}


// -----  Simple infection adders/removers  -----

void DescriptiveWithinHostModel::loadInfection(istream& stream) {
    infections.push_back(new DescriptiveInfection(stream));
}

void DescriptiveWithinHostModel::clearInfections( Treatments::Stages stage ){
    for(auto inf = infections.begin(); inf != infections.end();) {
        if( stage == Treatments::BOTH ||
            (stage == Treatments::LIVER && !(*inf)->bloodStage()) ||
            (stage == Treatments::BLOOD && (*inf)->bloodStage())
        ){
            delete *inf;
            inf = infections.erase( inf );
        }else{
            ++inf;
        }
    }
    numInfs = infections.size();
}

// -----  Interventions  -----

void DescriptiveWithinHostModel::clearImmunity() {
    for(auto inf = infections.begin(); inf != infections.end(); ++inf) {
        (*inf)->clearImmunity();
    }
    m_cumulative_h = 0.0;
    m_cumulative_Y_lag = 0.0;
}

void DescriptiveWithinHostModel::importInfection(LocalRng& rng){
    if( numInfs < MAX_INFECTIONS ){
        m_cumulative_h += 1;
        numInfs += 1;
        // This is a hook, used by interventions. The newly imported infections
        // should use initial frequencies to select genotypes.
        vector<double> weights( 0 );        // zero length: signal to use initial frequencies
        uint32_t genotype = Genotypes::sampleGenotype(rng, weights);
        infections.push_back(new DescriptiveInfection(rng, genotype, WithinHost::InfectionOrigin::Imported));
    }
    assert( numInfs == static_cast<int>(infections.size()) );
}


// -----  Density calculations  -----

void DescriptiveWithinHostModel::update(Host::Human &human, LocalRng& rng, int &nNewInfs_i, int &nNewInfs_l, 
        vector<double>& genotype_weights_i, vector<double>& genotype_weights_l, double ageInYears)
{
    // Note: adding infections at the beginning of the update instead of the end
    // shouldn't be significant since before latentp delay nothing is updated.
    int nNewInfsToBeCreated_i = nNewInfs_i;
    int nNewInfsToBeCreated_l = nNewInfs_l;

    nNewInfs_l = min(nNewInfs_l,MAX_INFECTIONS-numInfs);
    nNewInfs_i = min(nNewInfs_i,MAX_INFECTIONS-numInfs-nNewInfs_l);

    numInfs += nNewInfs_i;
    assert( numInfs>=0 && numInfs<=MAX_INFECTIONS );
    for( int i=0; i<nNewInfs_i; ++i ) {
        uint32_t genotype = Genotypes::sampleGenotype(rng, genotype_weights_i);

        // If opt_vaccine_genotype is true the infection is discarded with probability 1-vaccineFactor
        if( opt_vaccine_genotype )
        {
            double vaccineFactor = human.vaccine.getFactor( interventions::Vaccine::PEV, genotype );
            if(vaccineFactor == 1.0 || human.rng.bernoulli(vaccineFactor))
                infections.push_back(new DescriptiveInfection (rng, genotype, InfectionOrigin::Introduced));
        }
        else if (opt_vaccine_genotype == false)
            infections.push_back(new DescriptiveInfection (rng, genotype, InfectionOrigin::Introduced));
    }
    assert( numInfs == static_cast<int>(infections.size()) );

    numInfs += nNewInfs_l;
    assert( numInfs>=0 && numInfs<=MAX_INFECTIONS );
    for( int i=0; i<nNewInfs_l; ++i ) {
        uint32_t genotype = Genotypes::sampleGenotype(rng, genotype_weights_l);

        // If opt_vaccine_genotype is true the infection is discarded with probability 1-vaccineFactor
        if( opt_vaccine_genotype )
        {
            double vaccineFactor = human.vaccine.getFactor( interventions::Vaccine::PEV, genotype );
            if(vaccineFactor == 1.0 || human.rng.bernoulli(vaccineFactor))
                infections.push_back(new DescriptiveInfection (rng, genotype, InfectionOrigin::Indigenous));
        }
        else if (opt_vaccine_genotype == false)
            infections.push_back(new DescriptiveInfection (rng, genotype, InfectionOrigin::Indigenous));
    }
    assert( numInfs == static_cast<int>(infections.size()) );

    updateImmuneStatus ();

    totalDensity = 0.0;
    hrp2Density = 0.0;
    timeStepMaxDensity = 0.0;

    bool treatmentLiver = treatExpiryLiver > sim::ts0();
    bool treatmentBlood = treatExpiryBlood > sim::ts0();

    for(auto inf = infections.begin(); inf != infections.end();) {
        //NOTE: it would be nice to combine this code with that in
        // CommonWithinHost.cpp, but a few changes would be needed:
        // INNATE_MAX_DENS and MAX_DENS_CORRECTION would need to be required
        // (couldn't support old parameterisations using buggy versions of code
        // any more).
        // SP drug action and the PK/PD model would need to be abstracted
        // behind a common interface.
        if ( (*inf)->expired() /* infection has self-terminated */ ||
            ((*inf)->bloodStage() ? treatmentBlood : treatmentLiver) )
        {
            delete *inf;
            inf=infections.erase(inf);
            numInfs--;
            continue;
        }
        
        // Should be: infStepMaxDens = 0.0, but has some history.
        // See MAX_DENS_CORRECTION in DescriptiveInfection.cpp.
        double infStepMaxDens = timeStepMaxDensity;
        double immSurvFact = immunitySurvivalFactor(ageInYears, (*inf)->cumulativeExposureJ());
        double bsvFactor = human.vaccine.getFactor(interventions::Vaccine::BSV, opt_vaccine_genotype? (*inf)->genotype() : 0);

        (*inf)->determineDensities(rng, m_cumulative_h, infStepMaxDens, immSurvFact, _innateImmSurvFact, bsvFactor);

        if (bugfix_max_dens)
            infStepMaxDens = std::max(infStepMaxDens, timeStepMaxDensity);
        timeStepMaxDensity = infStepMaxDens;

        double density = (*inf)->getDensity();
        totalDensity += density;
        if( !(*inf)->isHrp2Deficient() ){
            hrp2Density += density;
        }

        ++inf;
    }
    
    // As in AJTMH p22, cumulative_h (X_h + 1) doesn't include infections added
    // this time-step and cumulative_Y only includes past densities.
    m_cumulative_h += nNewInfs_i + nNewInfs_l;
    m_cumulative_Y += sim::oneTS() * totalDensity;
    
    util::streamValidate( totalDensity );
    util::streamValidate( hrp2Density );
    assert( (std::isfinite)(totalDensity) );        // inf probably wouldn't be a problem but NaN would be
    
    // Cache total density for infectiousness calculations
    int y_lag_i = sim::moduloSteps(sim::ts1(), y_lag_len);

    for( size_t g = 0; g < Genotypes::N(); ++g )
    {
        m_y_lag_i[y_lag_i * Genotypes::N() + g] = 0.0;
        m_y_lag_l[y_lag_i * Genotypes::N() + g] = 0.0;
    }

    for( auto inf = infections.begin(); inf != infections.end(); ++inf )
    {
        if((*inf)->origin() == InfectionOrigin::Imported)
            m_y_lag_i[y_lag_i * Genotypes::N() + (*inf)->genotype()] += (*inf)->getDensity();
        else
            m_y_lag_l[y_lag_i * Genotypes::N() + (*inf)->genotype()] += (*inf)->getDensity();
    }

    // This is a bug, we keep it this way to be consistent with old simulations
    if(opt_vaccine_genotype == false)
    {
        nNewInfs_i = nNewInfsToBeCreated_i;
        nNewInfs_l = nNewInfsToBeCreated_l;
    }
}

InfectionOrigin DescriptiveWithinHostModel::getInfectionOrigin()const
{
    return get_infection_origin(infections);
}

// -----  Summarize  -----

bool DescriptiveWithinHostModel::summarize( Host::Human& human )const{
    pathogenesisModel->summarize( human );
    
    InfectionOrigin infectionType = get_infection_origin(infections);

    // If the number of infections is 0 and parasite density is positive we default to Indigenous
    if( infections.size() > 0 ){
        mon::reportStatMHI( mon::MHR_INFECTED_HOSTS, human, 1 );
        if(infectionType == InfectionOrigin::Indigenous)
            mon::reportStatMHI( mon::MHR_INFECTED_HOSTS_INDIGENOUS, human, 1 );
        else if(infectionType == InfectionOrigin::Introduced)
            mon::reportStatMHI( mon::MHR_INFECTED_HOSTS_INTRODUCED, human, 1 );
        else
            reportStatMHI( mon::MHR_INFECTED_HOSTS_IMPORTED, human, 1 );

        int nImported = 0, nIntroduced = 0, nIndigenous = 0;
        for( auto inf = infections.begin(); inf != infections.end(); ++inf )
        {
            if((*inf)->origin() == InfectionOrigin::Indigenous) nIndigenous++;
            else if((*inf)->origin() == InfectionOrigin::Introduced) nIntroduced++;
            else nImported++;
        }

        // (patent) infections are reported by genotype, even though we don't have
        // genotype in this model
        mon::reportStatMHGI( mon::MHR_INFECTIONS, human, 0, infections.size() );
        mon::reportStatMHGI( mon::MHR_INFECTIONS_IMPORTED, human, 0, nImported );
        mon::reportStatMHGI( mon::MHR_INFECTIONS_INTRODUCED, human, 0, nIntroduced );
        mon::reportStatMHGI( mon::MHR_INFECTIONS_INDIGENOUS, human, 0, nIndigenous );

        if( reportPatentInfected ){
            for(auto inf = infections.begin(); inf != infections.end(); ++inf)
            {
                if( diagnostics::monitoringDiagnostic().isPositive( human.rng, (*inf)->getDensity(), std::numeric_limits<double>::quiet_NaN() ) )
                {
                    mon::reportStatMHGI( mon::MHR_PATENT_INFECTIONS, human, 0, 1 );
                    if((*inf)->origin() == InfectionOrigin::Indigenous)
                        mon::reportStatMHGI( mon::MHR_PATENT_INFECTIONS_INDIGENOUS, human, 0, 1 );
                    else if((*inf)->origin() == InfectionOrigin::Introduced)
                        mon::reportStatMHGI( mon::MHR_PATENT_INFECTIONS_INTRODUCED, human, 0, 1 );
                    else
                        mon::reportStatMHGI( mon::MHR_PATENT_INFECTIONS_IMPORTED, human, 0, 1 );
                }
            }
        }
        if( reportInfectionsByGenotype ){
            // accumulate total density by genotype
            map<uint32_t, double> dens_by_gtype;
            for(auto inf = infections.begin(); inf != infections.end(); ++inf)
                dens_by_gtype[(*inf)->genotype()] += (*inf)->getDensity();
            
            for( auto gtype: dens_by_gtype ){
                // we had at least one infection of this genotype
                mon::reportStatMHGI( mon::MHR_INFECTED_GENOTYPE, human, gtype.first, 1 );
                if( diagnostics::monitoringDiagnostic().isPositive(human.rng, gtype.second, std::numeric_limits<double>::quiet_NaN()) ){
                    mon::reportStatMHGI( mon::MHR_PATENT_GENOTYPE, human, gtype.first, 1 );
                    mon::reportStatMHGF( mon::MHF_LOG_DENSITY_GENOTYPE, human, gtype.first, log(gtype.second) );
                }
            }
        }
    }
    
    // Some treatments (simpleTreat with steps=-1) clear infections immediately
    // (and are applied after update()), thus infections.size() may be 0 while
    // totalDensity > 0. Here we report the last calculated density.
    if( diagnostics::monitoringDiagnostic().isPositive(human.rng, totalDensity, std::numeric_limits<double>::quiet_NaN()) ){
        mon::reportStatMHI( mon::MHR_PATENT_HOSTS, human, 1 );
        if(infectionType == InfectionOrigin::Imported)
            mon::reportStatMHI( mon::MHR_PATENT_HOSTS_IMPORTED, human, 1 );
        else if(infectionType == InfectionOrigin::Introduced)
            mon::reportStatMHI( mon::MHR_PATENT_HOSTS_INTRODUCED, human, 1 );
        else if(infectionType == InfectionOrigin::Indigenous)
            mon::reportStatMHI( mon::MHR_PATENT_HOSTS_INDIGENOUS, human, 1 );

        if(totalDensity > 1e-10)
            mon::reportStatMHF( mon::MHF_LOG_DENSITY, human, log(totalDensity) );
        return true;    // patent
    }
    return false;       // not patent
}


// -----  Data checkpointing  -----

void DescriptiveWithinHostModel::checkpoint (istream& stream) {
    WHFalciparum::checkpoint (stream);
    for(int i=0; i<numInfs; ++i) {
        loadInfection(stream);  // create infections using a virtual function call
    }
    assert( numInfs == static_cast<int>(infections.size()) );
}
void DescriptiveWithinHostModel::checkpoint (ostream& stream) {
    WHFalciparum::checkpoint (stream);
    for(DescriptiveInfection* inf : infections)
        *inf & stream;
}

char const*const not_impl = "feature not available with the \"descriptive\" within-host model";
void DescriptiveWithinHostModel::treatPkPd(size_t schedule, size_t dosages, double age, double delay_d){
    throw TRACED_EXCEPTION( not_impl, util::Error::WHFeatures ); }

}
}
