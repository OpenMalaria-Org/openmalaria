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

#include "Host/WithinHost/Infection/PennyInfection.h"
#include "Host/WithinHost/CommonWithinHost.h"
#include "util/random.h"
#include "util/errors.h"
#include "util/CommandLine.h"
#include "util/ModelOptions.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>

namespace OM {
namespace WithinHost {

using namespace OM::util;

/** @brief The static variables (double)

Parameter names match names in Penny et al (2011)

Innate immunity related
 * beta_N := minimum value for R_N^x (value 0.5198)
 * psi_N := minimum value for R_N^y  (value 0.0946)
 * kappa_N := slope of innate immune responses R_N^x and R_N^y (value 2.9506)
 * sigma_epsilon := standard deviation associated with random effects (logNormal) (value 1.4217) 

Clonal immunity related
 * beta_C := minimum value for R_C^x  (value 0.1872)
 * psi_C := minimum value for R_C^y  (value 0.2224)
 * kappa_C := slope of clonal immune responses (value 1.9535) 
 * delta_C := delay to clonal antibody response (days) (value 7.2038 round to 7)  
 * rho_C := decay of clonal antibodies (value 0.1292) 
 * exp_negRho_C := exp of negative rho_C (value exp(-0.1292))

Variant-specific immunity related
 * beta_V := minimum value for R_V^x (value 0.0427)  
 * kappa_V := slope of variant specific immune responses R_V^x (value 4.1529)  
 * delta_V := delay to variant specific antibody response in R_V^x (days) (value 6.3572 round to 6)  
 * rho_V := decay of variant specific antibodies in R_V^x (value 2.5482)  
 * exp_negRho_V := exp of negative rho_V (value exp(-2.5482))
 * lambda_V := arrival times for new dominant variants (days) (value 4.2119) 

Parameters to assign infection dependent parameters
 * mu_Y := mean of log-Normal to assign Y(0) (log(per microL)) (value 3.9700)
 * sigma_Y := standard-deviation of log-Normal to assign Y(0) (value 1.3436)
 * mu_X := mean of log-Normal to assign X(0) (log(per microL)) (value 1.9969)
 * sigma_X := standard-deviation of log-Normal to assign X(0) (value 0.7424)
 * mu_TN := mean of log-Normal to assign T_N (log(per microL)) (value 7.5872)
 * sigma_TN  := standard-deviation of log-Normal to assign T_N (value 2.8977)
 * mu_TC := mean of log-Normal to assign T_C (log(per microL)) (value 5.5573)  
 * sigma_TC  := standard-deviation of log-Normal to assign T_C (value 0.4068) 
 * mu_TV := mean of log-Normal to assign T_V (log(per microL)) (value 6.12898)  
 * sigma_TV  := standard-deviation of log-Normal to assign T_V (value 1.3768) 
 
 Other infection parameters
 * m_rep := replication per cycle (value 16)
 * Omega := critical density when infection ends (per microL) (value 2.5times 10^-4) 
 
 */
//@{
/* Innate immunity related */
const double beta_N=0.5198;
const double psi_N=0.0946;
const double kappa_N=2.9506;
const double sigma_epsilon=1.4217; 

/* Clonal immunity related */
const double beta_C=0.1872;
const double psi_C=0.2224;
const double kappa_C=1.9535; 
const double rho_C=0.1292; 
const double exp_negRho_C=exp(-rho_C);

/* Variant-specific immunity related */
const double beta_V=0.0427;  
const double kappa_V=4.1529;  
const double rho_V=2.5482; 
const double exp_negRho_V=exp(-rho_V);
const double lambda_V=4.2119; 
const double prob_lambda_V = 1.0 / lambda_V;

/* distribution model paramters */
bool PennyInfection::immune_threshold_gamma = false;
bool PennyInfection::update_density_gamma = false;

/* Parameters to assign infection dependent parameters */
const double mu_Y=3.9700;
const double sigma_Y=1.3436;
const double a_Y = 8.7305;
const double b_Y = 0.4547;

const double mu_X=1.9969;
const double sigma_X=0.7424;
const double a_X=7.2350;
const double b_X=0.2760;

const double mu_TN=7.5872;
const double sigma_TN=2.8977;
const double a_TN = 6.8558;
const double b_TN = 1.1067;

const double mu_TC=5.5573;  
const double sigma_TC=0.4068; 
const double a_TC = 186.6233;
const double b_TC = 0.0297;

const double mu_TV=6.12898;  
const double sigma_TV=1.3768;
const double a_TV = 19.8167;
const double b_TV = 0.3093;

/* Other infection parameters */
const double m_rep=16.0; 
const double Omega=0.00025;
//TODO ending of infection is density per microL of ADULTS, we need to relate weight to BV

//@}

CommonInfection* createPennyInfection (LocalRng& rng, uint32_t protID, int origin) {
    return new PennyInfection (rng, protID, origin);
}

CommonInfection* checkpointedPennyInfection (istream& stream) {
    return new PennyInfection (stream);
}

void PennyInfection::init() {
    CommonWithinHost::createInfection = &createPennyInfection;
    CommonWithinHost::checkpointedInfection = &checkpointedPennyInfection;

    if(util::ModelOptions::option (util::IMMUNE_THRESHOLD_GAMMA)) {
        immune_threshold_gamma = true;
    } else {
        immune_threshold_gamma = false;
    }
    
    if(util::ModelOptions::option (util::UPDATE_DENSITY_GAMMA)) {
        update_density_gamma = true;
    } else {
        update_density_gamma = false;
    }
}

PennyInfection::PennyInfection(LocalRng& rng, uint32_t protID, int origin):
        CommonInfection(protID, origin),
        variantSpecificSummation(0),
        clonalSummation(0)
{
    // assign infection dependent immune thresholds
    if( immune_threshold_gamma )/* using gamma distribution */{
        do {
            threshold_N = exp(rng.gamma(a_TN,b_TN));
            threshold_C = exp(rng.gamma(a_TC,b_TC));
            threshold_V = exp(rng.gamma(a_TV,b_TV));
        } while(threshold_N <= threshold_C || threshold_N <= threshold_V);
    }else /* using lognormal distribution */{
        do {
            threshold_N = exp(rng.gauss(mu_TN,sigma_TN));
            threshold_C = exp(rng.gauss(mu_TC,sigma_TC));
            threshold_V = exp(rng.gauss(mu_TV,sigma_TV));
        } while(threshold_N <= threshold_C || threshold_N <= threshold_V);
    }
    
    for(int i=0; i<delta_C; ++i){
        cirDensities[i] = 0.0;
    }
    for(int i=0; i<delta_V; ++i){
        seqDensities[i] = 0.0;
    }
}


bool PennyInfection::updateDensity( LocalRng& rng, double survivalFactor, SimTime bsAge, double ){
    int ageDays = bsAge;       // lazy
    if( bsAge == sim::zero() ){
        // assign initial densities (Y circulating, X sequestered)
        size_t today = mod_nn(ageDays, delta_C);
        
        if(update_density_gamma) {
            cirDensities[today] = exp(rng.gamma(a_Y,b_Y));
        } else {
            cirDensities[today] = exp(rng.gauss(mu_Y,sigma_Y));
        }
        
        m_density = cirDensities[today];
        today = mod_nn(ageDays, delta_V);
        
        if(update_density_gamma){
            seqDensities[today] = exp(rng.gamma(a_X,b_X));
        } else {
            seqDensities[today] = exp(rng.gauss(mu_X,sigma_X));
        }
    }
    else /*not first day*/
    {
        // save yesterday's density (since getVariantSpecificSummation may reset it to zero)
        size_t yesterdayV = mod_nn(ageDays - 1, delta_V);
        double seqDensityYesterday = seqDensities[yesterdayV];
        
        // The immune responses are represented by the variables
        // R_Nx, RNy (probability that a parasite escapes control by clonal immune response)
        // R_Cx, RCy (probability that a parasite escapes control by clonal immune response)
        // R_Vx,    (probability that a parasite escapes control by clonal immune response)
        
        // innate immunity  
        size_t yesterdayC = mod_nn(ageDays - 1, delta_C);
        double base_N = cirDensities[yesterdayC]/threshold_N;
        double base_Npow = pow(base_N,kappa_N);
        double R_Nx = (1.0-beta_N) / (1.0 + base_Npow) + beta_N;
        double R_Ny = (1.0-psi_N) / (1.0 + base_Npow) + psi_N;
        
        // clonal immunity
        double base_C = getClonalSummation(ageDays)/threshold_C;
        double base_Cpow = pow(base_C,kappa_C);
        double R_Cx = (1.0-beta_C) / (1.0 + base_Cpow) + beta_C;
        double R_Cy = (1.0-psi_C) / (1.0 + base_Cpow) + psi_C;
        
        // variant specific immunity
        double base_V = getVariantSpecificSummation(rng, ageDays)/threshold_V;
        double R_Vx = (1.0-beta_V) / (1.0 + pow(base_V,kappa_V)) + beta_V;
        
        // cirDensity: circulating density of circulating at t
        // seqDensity: sequestered density of circulating at t
        // new cirDensity  = (seqDensity(t-1)) * replication * the probability that
        //  seq parasites escape control by immune response.
        // new seqDensity  = (cirDensity(t-1)) * the probability that
        //  circ parasites escape control by immune response.
        double cirDensity_new = seqDensityYesterday * m_rep * R_Vx * R_Cx * R_Nx;
        double seqDensity_new = cirDensities[yesterdayC] * R_Cy * R_Ny;
        
        // end infection if density less than Omega (per microL)
        // add random biologicaleffect to circulating
        if (cirDensity_new < Omega) {
            cirDensity_new = 0.0;
        } else {
            
            if( update_density_gamma ) {
                double a_cirDens = pow(log(cirDensity_new),2)/pow(sigma_epsilon,2);
                double b_cirDens = pow(sigma_epsilon,2)/log(cirDensity_new);
                cirDensity_new = exp(rng.gamma(a_cirDens,b_cirDens) ) * survivalFactor;
            } else {
                cirDensity_new = exp(rng.gauss(log(cirDensity_new),sigma_epsilon)) * survivalFactor;
            }
            // please don't simplify this, we want more chance at ending infection
            if (cirDensity_new < Omega) {
                cirDensity_new = 0.0;
            }
        }
        seqDensity_new *= survivalFactor;
        if (seqDensity_new < Omega) {
            if (cirDensity_new == 0.0){
                // infection is extinct
                return true;
            }
            seqDensity_new = 0.0;
        }
        
        size_t todayC = mod_nn(ageDays, delta_C);
        cirDensities[todayC] = cirDensity_new;
        m_density = cirDensities[todayC];
        
        size_t todayV = mod_nn(ageDays, delta_V);
        seqDensities[todayV] = seqDensity_new;
    }
    
    // used for immunity across infections
    m_cumulativeExposureJ += m_density;
    
    // if we haven't already exited this funciton, the infection is not extinct (so return false)
    return false;
}

double PennyInfection::getVariantSpecificSummation(LocalRng& rng, int ageDays) {
    // check if new dominant variant has arrived
    bool newVarDominant = rng.bernoulli(prob_lambda_V);
    //draw from bernouli distrb, prob 1/lambda_V
    if (newVarDominant) {
        variantSpecificSummation = 0;
        for(int i=0; i<delta_V; ++i){
            seqDensities[i] = 0.0;
        }
    }
    //The effective exposure is computed by adding in the delta_V-day lagged parasite density 
    //and decaying the previous value for the effective exposure with decay parameter rho_V
    size_t index = mod_nn(ageDays, delta_V);
    variantSpecificSummation = (variantSpecificSummation * exp_negRho_V) + seqDensities[index];
    
    return variantSpecificSummation;
}

double PennyInfection::getClonalSummation(int ageDays) {
    //The effective exposure is computed by adding in the delta_C-day lagged parasite density 
    //and decaying the previous value for the effective exposure with decay parameter rho_C
    size_t index = mod_nn(ageDays, delta_C);
    clonalSummation = (clonalSummation * exp_negRho_C) + cirDensities[index];
    
    return clonalSummation;
}


PennyInfection::PennyInfection (istream& stream) :
        CommonInfection(stream)
{
    for(int i=0; i<delta_C; ++i){
        cirDensities[i] & stream;
    }
    for(int i=0; i<delta_V; ++i){
        seqDensities[i] & stream;
    }
    threshold_N & stream;
    threshold_V & stream;
    threshold_C & stream;
    variantSpecificSummation & stream;
    clonalSummation & stream;
}

void PennyInfection::checkpoint (ostream& stream) {
    CommonInfection::checkpoint(stream);
    for(int i=0; i<delta_C; ++i){
        cirDensities[i] & stream;
    }
    for(int i=0; i<delta_V; ++i){
        seqDensities[i] & stream;
    }
    threshold_N & stream;
    threshold_V & stream;
    threshold_C & stream;
    variantSpecificSummation & stream;
    clonalSummation & stream;
}

}
}
