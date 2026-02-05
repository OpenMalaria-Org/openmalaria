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

#include "Host/Human.h"
#include "Host/WithinHost/CommonWithinHost.h"
#include "Host/WithinHost/Diagnostic.h"
#include "Host/WithinHost/Genotypes.h"
#include "Host/WithinHost/Pathogenesis/PathogenesisModel.h"
#include "util/errors.h"
#include "util/ModelOptions.h"
#include "util/AgeGroupInterpolation.h"
#include "util/random.h"
#include "util/StreamValidator.h"
#include "schema/scenario.h"

using namespace std;

namespace OM {
namespace WithinHost {

CommonInfection* (* CommonWithinHost::createInfection) (LocalRng& rng, uint32_t protID, int origin);
CommonInfection* (* CommonWithinHost::checkpointedInfection) (istream& stream);

double hetMassMultStdDev = std::numeric_limits<double>::signaling_NaN();
double minHetMassMult = std::numeric_limits<double>::signaling_NaN();
util::AgeGroupInterpolator massByAge;

bool reportInfectedOrPatentInfected = false;

// -----  Initialization  -----

void CommonWithinHost::init( const scnXml::Scenario& scenario ){
    const scnXml::Human human = scenario.getModel().getHuman();
    if( !human.getWeight().present() ){
        // Technically this is needed by the PK/PD and Molineaux models
        throw util::xml_scenario_error( "model->human->weight element required by certain models" );
    }
    massByAge.set( human.getWeight().get(), "weight" );
    hetMassMultStdDev = human.getWeight().get().getMultStdDev();
    // hetWeightMult must be large enough that birth weight is at least 0.5 kg:
    minHetMassMult = 0.5 / massByAge.eval( 0.0 );
    
    reportInfectedOrPatentInfected = mon::isUsedM(mon::MHR_INFECTIONS) ||
        mon::isUsedM(mon::MHR_PATENT_INFECTIONS);
    
    PkPd::LSTMModel::init( scenario );
}

CommonWithinHost::CommonWithinHost( LocalRng& rng, double comorbidityFactor ) :
        WHFalciparum( rng, comorbidityFactor )
{
    assert( sim::oneTS() == sim::fromDays(1) || sim::oneTS() == sim::fromDays(5) );
    
    opt_vaccine_genotype = util::ModelOptions::option (util::VACCINE_GENOTYPE);

    // Sample a weight heterogeneity factor
#ifndef NDEBUG
    int counter = 0;
#endif
    do {
        hetMassMultiplier = rng.gauss( 1.0, hetMassMultStdDev );
#ifndef NDEBUG
        assert( counter < 100 );        // too many resamples: resamples should rarely be needed...
        ++counter;
#endif
    } while( hetMassMultiplier < minHetMassMult );
}

CommonWithinHost::~CommonWithinHost() {
    for( auto inf = infections.begin(); inf != infections.end(); ++inf ){
        delete *inf;
    }
    infections.clear();
}

// -----  Simple infection adders/removers  -----

void CommonWithinHost::clearInfections( Treatments::Stages stage ){
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

// -----  interventions -----

void CommonWithinHost::treatPkPd(size_t schedule, size_t dosage, double age, double delay_d){
    double mass = massByAge.eval( age ) * hetMassMultiplier;
    pkpdModel.prescribe( schedule, dosage, age, mass, delay_d );
}
void CommonWithinHost::clearImmunity() {
    for(auto inf = infections.begin(); inf != infections.end(); ++inf) {
        (*inf)->clearImmunity();
    }
    m_cumulative_h = 0.0;
    m_cumulative_Y_lag = 0.0;
}
void CommonWithinHost::importInfection(LocalRng& rng){
    if( numInfs < MAX_INFECTIONS ){
        m_cumulative_h += 1;
        numInfs += 1;
        // This is a hook, used by interventions. The newly imported infections
        // should use initial frequencies to select genotypes.
        vector<double> weights( 0 );        // zero length: signal to use initial frequencies
        uint32_t genotype = Genotypes::sampleGenotype(rng, weights);
        infections.push_back(createInfection(rng, genotype, WithinHost::InfectionOrigin::Imported));
    }
    assert( numInfs == static_cast<int>(infections.size()) );
}


// -----  Density calculations  -----

void CommonWithinHost::update(Host::Human &human, LocalRng& rng, int &nNewInfs_i, int &nNewInfs_l, 
        vector<double>& genotype_weights_i, vector<double>& genotype_weights_l, double ageInYears)
{
    // Note: adding infections at the beginning of the update instead of the end
    // shouldn't be significant since before latentp delay nothing is updated.
    nNewInfs_l = min(nNewInfs_l,MAX_INFECTIONS-numInfs);
    nNewInfs_i = min(nNewInfs_i,MAX_INFECTIONS-numInfs-nNewInfs_l);

    int nNewInfsIgnored = nNewInfs_i + nNewInfs_l - (MAX_INFECTIONS-numInfs);

    assert( numInfs>=0 && numInfs<=MAX_INFECTIONS );

    int nNewInfsDiscarded = 0;
    for( int i=0; i<nNewInfs_i; ++i ) {
        uint32_t genotype = Genotypes::sampleGenotype(rng, genotype_weights_i);

        // If opt_vaccine_genotype is true the infection is discarded with probability 1-vaccineFactor
        if( opt_vaccine_genotype )
        {
            double vaccineFactor = human.vaccine.getFactor( interventions::Vaccine::PEV, genotype );
            if(vaccineFactor == 1.0 || human.rng.bernoulli(vaccineFactor))
                infections.push_back(createInfection (rng, genotype, InfectionOrigin::Introduced));
            else
                nNewInfsDiscarded++;
        }
        else if (opt_vaccine_genotype == false)
            infections.push_back(createInfection (rng, genotype, InfectionOrigin::Introduced));
    }
    // Update nNewInfs, this is the number that will be reported in Human
    nNewInfs_i -= nNewInfsDiscarded;
    numInfs += nNewInfs_i;

    nNewInfsDiscarded = 0;
    for( int i=0; i<nNewInfs_l; ++i ) {
        uint32_t genotype = Genotypes::sampleGenotype(rng, genotype_weights_l);

        // If opt_vaccine_genotype is true the infection is discarded with probability 1-vaccineFactor
        if( opt_vaccine_genotype )
        {
            double vaccineFactor = human.vaccine.getFactor( interventions::Vaccine::PEV, genotype );
            if(vaccineFactor == 1.0 || human.rng.bernoulli(vaccineFactor))
                infections.push_back(createInfection (rng, genotype, InfectionOrigin::Indigenous));
            else
                nNewInfsDiscarded++;
        }
        else if (opt_vaccine_genotype == false)
            infections.push_back(createInfection (rng, genotype, InfectionOrigin::Indigenous));
    }
    // Update nNewInfs, this is the number that will be reported in Human
    nNewInfs_l -= nNewInfsDiscarded;
    numInfs += nNewInfs_l;

    assert( numInfs == static_cast<int>(infections.size()) );
    
    updateImmuneStatus ();

    totalDensity = 0.0;
    hrp2Density = 0.0;
    timeStepMaxDensity = 0.0;
    
    bool treatmentLiver = treatExpiryLiver > sim::ts0();
    bool treatmentBlood = treatExpiryBlood > sim::ts0();
    
    double body_mass = massByAge.eval( ageInYears ) * hetMassMultiplier;
    
    for( SimTime now = sim::ts0(), end = sim::ts0() + sim::oneTS(); now < end; now = now + sim::oneDay() ){
        // every day, medicate drugs, update each infection, then decay drugs
        pkpdModel.medicate(rng);
        
        double sumLogDens = 0.0;
        
        for(auto inf = infections.begin(); inf != infections.end();) {
            // Note: this is only one treatment model; there is also the PK/PD model
            bool expires = ((*inf)->bloodStage() ? treatmentBlood : treatmentLiver);
            
            if( !expires ){     /* no expiry due to simple treatment model; do update */
                const double drugFactor = pkpdModel.getDrugFactor(rng, *inf, body_mass);
                const double immFactor = immunitySurvivalFactor(ageInYears, (*inf)->cumulativeExposureJ());
                const double bsvFactor = human.vaccine.getFactor(interventions::Vaccine::BSV, opt_vaccine_genotype? (*inf)->genotype() : 0);
                const double survivalFactor = bsvFactor * _innateImmSurvFact * immFactor * drugFactor;
                // update, may result in termination of infection:
                expires = (*inf)->update(rng, survivalFactor, now, body_mass);
            }
            
            if( expires ){
                delete *inf;
                inf = infections.erase(inf);        // inf points to next infection now so don't increment with ++inf
                --numInfs;
            } else {
                double density = (*inf)->getDensity();
                totalDensity += density;
                if( !(*inf)->isHrp2Deficient() ){
                    hrp2Density += density;
                }
                timeStepMaxDensity = max(timeStepMaxDensity, density);
                if( density > 0 ){
                    // Base 10 logarithms are usually used; +1 because it avoids negatives in output while having very little affect on high densities
                    sumLogDens += log10(1.0 + density);
                }
                ++inf;
            }
        }
        pkpdModel.decayDrugs (body_mass);
    }
    
    // As in AJTMH p22, cumulative_h (X_h + 1) doesn't include infections added
    // this time-step and cumulative_Y only includes past densities, thus we
    // increment these after the update.
    m_cumulative_h += nNewInfs_i + nNewInfs_l;
    m_cumulative_Y += totalDensity;
    
    util::streamValidate(totalDensity);
    util::streamValidate(hrp2Density);
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
    if(nNewInfsIgnored > 0)
        nNewInfs_l += nNewInfsIgnored;
}

void CommonWithinHost::addProphylacticEffects(const vector<double>& pClearanceByTime) {
    // this should actually be easy; it just isn't needed yet
    throw util::unimplemented_exception( "prophylactic effects on 1-day time step" );
}

InfectionOrigin CommonWithinHost::getInfectionOrigin()const
{
    return get_infection_origin(infections);
}

// -----  Summarize  -----

// Used in summarizeInfs.
vector<CommonInfection*> sortedInfs;
struct InfGenotypeSorter {
    bool operator() (CommonInfection* i, CommonInfection* j){
        return i->genotype() < j->genotype();
    }
} infGenotypeSorter;

bool CommonWithinHost::summarize( Host::Human& human )const{
    pathogenesisModel->summarize( human );
    pkpdModel.summarize( human );
    
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

        if( reportInfectedOrPatentInfected ){
            for(auto inf = infections.begin(); inf != infections.end(); ++inf)
            {
                uint32_t genotype = (*inf)->genotype();
                mon::reportStatMHGI( mon::MHR_INFECTIONS, human, genotype, 1 );
                if((*inf)->origin() == InfectionOrigin::Indigenous)
                    mon::reportStatMHGI( mon::MHR_INFECTIONS_INDIGENOUS, human, genotype, 1 );
                else if((*inf)->origin() == InfectionOrigin::Introduced)
                    mon::reportStatMHGI( mon::MHR_INFECTIONS_INTRODUCED, human, genotype, 1 );
                else
                    mon::reportStatMHGI( mon::MHR_INFECTIONS_IMPORTED, human, genotype, 1 );

                if( diagnostics::monitoringDiagnostic().isPositive( human.rng, (*inf)->getDensity(), std::numeric_limits<double>::quiet_NaN() ) ){
                    mon::reportStatMHGI( mon::MHR_PATENT_INFECTIONS, human, genotype, 1 );
                    if((*inf)->origin() == InfectionOrigin::Indigenous)
                        mon::reportStatMHGI( mon::MHR_PATENT_INFECTIONS_INDIGENOUS, human, genotype, 1 );
                    else if((*inf)->origin() == InfectionOrigin::Introduced)
                        mon::reportStatMHGI( mon::MHR_PATENT_INFECTIONS_INTRODUCED, human, genotype, 1 );
                    else
                        mon::reportStatMHGI( mon::MHR_PATENT_INFECTIONS_IMPORTED, human, genotype, 1 );
                }
            }
        }
        
        if( reportInfectionsByGenotype ){
            // Instead of storing nInfs and total density by genotype we sort
            // infections by genotype and report each in sequence.
            // We don't sort in place since that would affect random number sampling
            // order when updating, and the monitoring system should not in my
            // opinion affect outputs (since it would make testing harder).
            sortedInfs.assign( infections.begin(), infections.end() );
            sort( sortedInfs.begin(), sortedInfs.end(), infGenotypeSorter );
            auto inf = sortedInfs.begin();
            while( inf != sortedInfs.end() ){
                uint32_t genotype = (*inf)->genotype();
                double dens = 0.0;
                do{     // at start: genotype is that of the current infection, dens is 0
                    dens += (*inf)->getDensity();
                    ++inf;
                }while( inf != sortedInfs.end() && (*inf)->genotype() == genotype );
                // we had at least one infection of this genotype
                mon::reportStatMHGI( mon::MHR_INFECTED_GENOTYPE, human, genotype, 1 );
                if( diagnostics::monitoringDiagnostic().isPositive(human.rng, dens, std::numeric_limits<double>::quiet_NaN()) ){
                    mon::reportStatMHGI( mon::MHR_PATENT_GENOTYPE, human, genotype, 1 );
                    mon::reportStatMHGF( mon::MHF_LOG_DENSITY_GENOTYPE, human, genotype, log(dens) );
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

        if(totalDensity == 0.0)
            mon::reportStatMHF( mon::MHF_LOG_DENSITY, human, 0.0);
        else
            mon::reportStatMHF( mon::MHF_LOG_DENSITY, human, log(totalDensity) );
        return true;    // patent
    }
    return false;       // not patent
}


void CommonWithinHost::checkpoint (istream& stream) {
    WHFalciparum::checkpoint (stream);
    hetMassMultiplier & stream;
    pkpdModel & stream;
    for(int i = 0; i < numInfs; ++i) {
        infections.push_back (checkpointedInfection (stream));
    }
    assert( numInfs == static_cast<int>(infections.size()) );
}

void CommonWithinHost::checkpoint (ostream& stream) {
    WHFalciparum::checkpoint (stream);
    hetMassMultiplier & stream;
    pkpdModel & stream;
    for(auto inf = infections.begin(); inf != infections.end(); ++inf) {
        (**inf) & stream;
    }
}
}
}
