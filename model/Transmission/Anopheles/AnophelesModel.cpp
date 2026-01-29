/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2025 Swiss Tropical and Public Health Institute
 * Copyright (C) 2005-2015 Liverpool School Of Tropical Medicine
 * Copyright (C) 2020-2025 University of Basel
 * Copyright (C) 2025 The Kids Research Institute Australia
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
#include "Transmission/TransmissionModel.h"
#include "Transmission/Anopheles/AnophelesModel.h"
#include "Transmission/PerHost.h"
#include "Population.h"
#include "Host/WithinHost/Genotypes.h"
#include "util/vectors.h"
#include "util/errors.h"
#include "util/ModelOptions.h"
#include "util/StreamValidator.h"

#include <cmath>

#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

struct quadratic_params
{
    double ah; // The availability rate of hosts h. Time-1
    double an; // The availability rate of hosts n. Time-1
    double mu_v_interv; // Mosquito mortality rate while host-seeking. Time-1
    double theta_d; // Proportion of the night that the mosquito spends host-seeking. Time.
    double PAt_wanted; // PAt is the daily probability of a mosquito feeding on an ATSB
};

double quadratic_f (double x, void *params)
{
    struct quadratic_params *p = (struct quadratic_params *) params;

    double ah = p->ah;
    double an = p->an;
    double mu_v_interv = p->mu_v_interv;
    double theta_d = p->theta_d;
    double PAt_wanted = p->PAt_wanted;
    
    double tmp = (x + ah + an + mu_v_interv);

    return (1.0 - exp(-tmp * theta_d)) * x / tmp - PAt_wanted;
}

double quadratic_df (double x, void *params)
{
    struct quadratic_params *p = (struct quadratic_params *) params;

    double ah = p->ah;
    double an = p->an;
    double mu_v_interv = p->mu_v_interv;
    double theta_d = p->theta_d;

    double tmp = (x + ah + an + mu_v_interv);

    double a = (1.0 - exp(-tmp * theta_d)) * x / (tmp * tmp);
    double b = (1.0 - exp(-tmp * theta_d)) / tmp;
    double c = exp(-tmp * theta_d) * x * theta_d / tmp;

    return -a + b + c;
}

void quadratic_fdf (double x, void *params, double *y, double *dy)
{
    struct quadratic_params *p = (struct quadratic_params *) params;

    double ah = p->ah;
    double an = p->an;
    double mu_v_interv = p->mu_v_interv;
    double theta_d = p->theta_d;
    double PAt_wanted = p->PAt_wanted;

    double tmp = (x + ah + an + mu_v_interv);

    double a = (1.0 - exp(-tmp)) * x / (tmp * tmp);
    double b = (1.0 - exp(-tmp)) / tmp;
    double c = exp(-tmp * theta_d) * x * theta_d / tmp;

    *y = (1.0 - exp(-tmp * theta_d)) * x / tmp - PAt_wanted;
    *dy = -a + b + c;
}

#define USE_DERIVATIVE

#ifdef USE_DERIVATIVE
    using gsl_root_solver_type = gsl_root_fdfsolver_type;
    using gsl_root_solver = gsl_root_fdfsolver;
    auto gsl_root_solver_alloc = gsl_root_fdfsolver_alloc;
    auto gsl_root_solver_free = gsl_root_fdfsolver_free;
    auto gsl_root_solver_name = gsl_root_fdfsolver_name;
    auto gsl_root_solver_iterate = gsl_root_fdfsolver_iterate;
    auto gsl_root_solver_root = gsl_root_fdfsolver_root;
    // auto gsl_root_solver_method = gsl_root_fdfsolver_newton;
    auto gsl_root_solver_method = gsl_root_fdfsolver_secant;
    // auto gsl_root_solver_method = gsl_root_fdfsolver_steffenson;
#else
    using gsl_root_solver_type = gsl_root_fsolver_type;
    using gsl_root_solver = gsl_root_fsolver;
    auto gsl_root_solver_alloc = gsl_root_fsolver_alloc;
    auto gsl_root_solver_free = gsl_root_fsolver_free;
    auto gsl_root_solver_name = gsl_root_fsolver_name;
    auto gsl_root_solver_iterate = gsl_root_fsolver_iterate;
    auto gsl_root_solver_root = gsl_root_fsolver_root;
    // auto gsl_root_solver_method = gsl_root_fsolver_brent;
    auto gsl_root_solver_method = gsl_root_fsolver_bisection;
    // auto gsl_root_solver_method = gsl_root_fsolver_falsepos;
#endif


namespace OM
{
namespace Transmission
{
namespace Anopheles
{
using namespace OM::util;
using WithinHost::Genotypes;

// -----  Initialisation of model, done before human warmup  ------

void AnophelesModel::initialise(size_t species, MosquitoParams mosqParams)
{
    mosq = mosqParams;

    N_v_length = (mosq.EIPDuration + mosq.restDuration);

    // -----  allocate memory  -----
    // Set up fArray and ftauArray. Each step, all elements not set here are
    // calculated, even if they aren't directly used in the end;
    // however all calculated values are used in calculating the next value.
    fArray.resize((mosq.EIPDuration - mosq.restDuration + sim::oneDay()));
    fArray[sim::zero()] = 1.0;
    ftauArray.resize(mosq.EIPDuration);
    for (int i = 0; i < mosq.restDuration; i += 1)
    {
        ftauArray[i] = 0.0;
    }
    ftauArray[mosq.restDuration] = 1.0;
    uninfected_v.resize(N_v_length);
    uninfected_v[sim::zero()] = numeric_limits<double>::quiet_NaN(); // index not used
}

void AnophelesModel::initAvailability(size_t species, const vector<NhhParams> &nhhs, int populationSize)
{
    // A: Host seeking
    // Proportion of host-seeking parous mosquitoes (those who have laid eggs)
    // which laid eggs that day:
    const double A0 = mosq.laidEggsSameDayProportion;
    // Probability that the mosquito survives the feeding cycle.
    // Note: Pf = M, the parous rate (prop mosqs which have laid eggs):
    const double Pf = mosq.survivalFeedingCycleProbability;
    const double humanBloodIndex = mosq.humanBloodIndex; // χ (chi)
    // Cycle probabilities, when biting a human:
    // B: Host encountered
    const double P_B1 = mosq.probBiting;
    // C: Fed
    const double P_C1 = mosq.probFindRestSite;
    // D: Resting
    const double P_D1 = mosq.probResting;
    // E: Laying eggs (ovipositing)
    const double P_E1 = mosq.probOvipositing;

    // -----  Calculate P_A, P_A1, P_A_n  -----
    // P_A is prob that a mosq is still host seeking after 1 day. It is also the
    // proportion of parous mosquitoes who have waited at least 1 day since
    // laying, thus 1 - P_A = A0.
    const double initP_A = 1.0 - A0;
    // This is multiplied by the probability of encountering some type of host
    // on a given night is the total availability of this type of host (N_i * α_i).
    // Note: P_Ai = α_i * N_i / availFactor
    const double availFactor = -log(initP_A) / (mosq.seekingDuration * (1.0 - initP_A));
    // Probability that a mosquito encounters a human on a given night
    double P_A1 = numeric_limits<double>::quiet_NaN();
    // Probability that a mosquito encounters any non human host on a given night
    // (confusingly labelled P_Ah in paper)
    double P_Ah = numeric_limits<double>::quiet_NaN();

    if (nhhs.empty())
    {
        // Number of non-human hosts: χ=1
        P_A1 = A0 * Pf / (P_B1 * P_C1 * P_D1 * P_E1);
        P_Ah = 0.0;
    }
    else
    {
        // Have non-human hosts: χ<1

        // let v = χ * P_D_1 * P_E_1; note that this is the average for humans
        const double v = humanBloodIndex * P_D1 * P_E1;
        const double chi1 = 1.0 - humanBloodIndex; // 1 - χ

        double sum_xi = 0.0;  // sum rel. avail across NNHs; should be 1
        double sum_u = 0.0;   // sum u across NNHs where u = xi * P_B * P_C
        double sum_uvw = 0.0; // sum u*(v+w) across NNHs where w = (1-χ)*P_D*P_E

        for (auto &nhh : nhhs)
        {
            // availability population of hosts of this type relative to other non-human hosts:
            const double xi_i = nhh.mosqRelativeEntoAvailability;
            // cycle probabilities, when biting this type of host:
            const double P_B_i = nhh.mosqProbBiting;
            const double P_C_i = nhh.mosqProbFindRestSite;
            const double P_D_i = nhh.mosqProbResting;
            sum_xi += xi_i;
            const double u_i = xi_i * P_B_i * P_C_i;
            sum_u += u_i;
            const double w_i = chi1 * P_D_i * P_E1; // note: we now assume P_E_i = P_E1
            sum_uvw += u_i * (v + w_i);
        }

        if (!(sum_xi > 0.9999 && sum_xi < 1.0001))
        {
            throw xml_scenario_error("The sum of the relative entomological availability (ξ_i) "
                                     "across non-human hosts must be 1!");
        }

        // equations (14), (15) of paper:
        P_A1 = (A0 * Pf * humanBloodIndex * sum_u) / (P_B1 * P_C1 * sum_uvw);
        P_Ah = (A0 * Pf * chi1) / sum_uvw;
    }

    // -----  Calculate availability rate of hosts (α_i) and non-human population data  -----
    PerHostAnophParams::setEntoAvailabilityFactor(species, (P_A1 / populationSize) * availFactor);

    initAvail = (P_A1 / populationSize) * availFactor;

    nhh_avail = 0.0;
    nhh_sigma_df = 0.0;
    nhh_sigma_dff = 0.0;
    for (auto &nhh : nhhs)
    {
        // availability population of hosts of this type relative to other non-human hosts:
        const double xi_i = nhh.mosqRelativeEntoAvailability;
        // cycle probabilities, when biting this type of host:
        const double P_B_i = nhh.mosqProbBiting;
        const double P_C_i = nhh.mosqProbFindRestSite;
        const double P_D_i = nhh.mosqProbResting;
        const double rel_fecundity = nhh.hostFecundityFactor;
        const double P_Ahi = P_Ah * xi_i;           // probability of encountering this type of NNH on a given night
        const double avail_i = P_Ahi * availFactor; // N_i * α_i

        nhh_avail += avail_i;                              // N * α (nhh_avail is not used anymore)
        const double df = avail_i * P_B_i * P_C_i * P_D_i; // term in P_df series
        nhh_sigma_df += df;
        nhh_sigma_dff += df * rel_fecundity;
        // Note: we would do the same for P_dif except that it's multiplied by
        // infectiousness of host to mosquito which is zero.

        string name = nhh.name;

        Nhh nhhi;
        nhhi.avail_i = avail_i;
        nhhi.P_B_I = P_B_i;
        nhhi.P_C_I = P_C_i;
        nhhi.P_D_I = P_D_i;
        nhhi.rel_fecundity = rel_fecundity;
        nhhi.expiry = sim::future();

        nhhInstances[name] = nhhi;
    }

    // ———  set mosqSeekingDeathRate  ———
    // Note: sum_k( P_Ak ) = P_A1 + P_Ah
    const double mu1 = (1.0 - (initP_A + P_A1 + P_Ah)) / (1.0 - initP_A);
    const double mu2 = -log(initP_A) / mosq.seekingDuration;
    mosq.seekingDeathRate = mu1 * mu2;

    if (mosq.seekingDeathRate <= 0)
        throw xml_scenario_error("Anopheles::init mosq.seekingDeathRate is <= 0. \
The probability of a mosquito surviving the entire feeding cycle \
should be less than the product of the probabilities of a \
mosquito surviving biting, finding a resting site, resting and ovipositing.");
}

/** Called by initialise function to init variables directly related to EIR
 *
 * @param anoph Data from XML
 * @param initialisationEIR In/out parameter: TransmissionModel::initialisationEIR
 * @param EIPDuration parameter from MosqTransmission (used for an estimation)
 */
void AnophelesModel::initEIR(const vector<double> &initEIR365, vector<double> FSCoefficInit, double EIRRotateAngleInit, double propInfectious, double propInfected)
{
    FSCoeffic = FSCoefficInit;
    EIRRotateAngle = EIRRotateAngleInit;
    speciesEIR = initEIR365;

    partialInitEIR.resize(sim::stepsPerYear(), 0.0);

    // Add to the TransmissionModel's EIR, used for the initalization phase.
    // Note: sum stays the same, units changes to per-time-step.
    for (SimTime i = sim::zero(); i < sim::oneYear(); i = i + sim::oneDay())
    {
        partialInitEIR[mod_nn(sim::inSteps(i), sim::stepsPerYear())] += initEIR365[i];
    }

    if (util::CommandLine::option(util::CommandLine::PRINT_ANNUAL_EIR))
        cout << "Annual EIR for " << mosq.name << ": " << vectors::sum(speciesEIR) << endl;

    // Set other data used for mosqEmergeRate calculation:
    initNvFromSv = 1.0 / propInfectious;
    initOvFromSv = initNvFromSv * propInfected;
}

// -----  Initialisation of model which is done after creating initial humans  -----

void AnophelesModel::init2(int nHumans, double meanPopAvail, double sum_avail, double sigma_f, double sigma_df, double sigma_dff)
{
    // -----  Calculate P_A, P_Ai, P_df based on pop age structure  -----

    // ν_A: rate at which mosquitoes find hosts or die (i.e. leave host-seeking state)
    double leaveRate = sum_avail + nhh_avail + mosq.seekingDeathRate;
    sigma_df += nhh_sigma_df;
    sigma_dff += nhh_sigma_dff;

    // Probability of a mosquito not finding a host this day:
    double tsP_A = exp(-leaveRate * mosq.seekingDuration);
    double availDivisor = (1.0 - tsP_A) / leaveRate; // α_d

    // Input per-species EIR is the mean EIR experienced by a human adult.
    // We use sumPFindBite below to get required S_v.
    double sumPFindBite = sigma_f * availDivisor;
    double tsP_df = sigma_df * availDivisor * mosq.probMosqSurvivalOvipositing;
    double tsP_dff = sigma_dff * availDivisor * mosq.probMosqSurvivalOvipositing;

    double tsP_Amu = (1 - tsP_A) * mosq.seekingDeathRate / (mosq.seekingDeathRate + sum_avail + nhh_avail);
    double tsP_A1 = (1 - tsP_A) * sum_avail / (mosq.seekingDeathRate + sum_avail + nhh_avail);
    double tsP_Ah = (1 - tsP_A) * nhh_avail / (mosq.seekingDeathRate + sum_avail + nhh_avail);
    // for( auto it = nhhInstances.begin(); it != nhhInstances.end(); ++it){
    //     initialP_Ah += (1-initialP_A) * it->second.avail_i / (mosqSeekingDeathRate + sum_avail + nhh_avail);
    // }

    // -----  Calculate required S_v based on desired EIR  -----
    initNv0FromSv = initNvFromSv * (1.0 - tsP_A - tsP_df);
    initSvFromEIR = nHumans * meanPopAvail / sumPFindBite;

    // We estimate forcedS_v from EIR
    forcedS_v = speciesEIR;
    vectors::scale(forcedS_v, initSvFromEIR);

    N_v.resize(N_v_length, numeric_limits<double>::quiet_NaN());
    O_v_i.resize(N_v_length * Genotypes::N(), 0.0);
    O_v_l.resize(N_v_length * Genotypes::N(), numeric_limits<double>::quiet_NaN());
    S_v_i.resize(N_v_length * Genotypes::N(), 0.0);
    S_v_l.resize(N_v_length * Genotypes::N(), numeric_limits<double>::quiet_NaN()); 
    P_A.resize(N_v_length, numeric_limits<double>::quiet_NaN());
    P_df.resize(N_v_length, numeric_limits<double>::quiet_NaN());
    P_dif_i.resize(N_v_length * Genotypes::N(), 0.0); 
    P_dif_l.resize(N_v_length * Genotypes::N(), 0.0);
    P_dff.resize(N_v_length, numeric_limits<double>::quiet_NaN());
    P_Amu.resize(N_v_length, numeric_limits<double>::quiet_NaN());
    P_A1.resize(N_v_length, numeric_limits<double>::quiet_NaN());
    P_Ah.resize(N_v_length, numeric_limits<double>::quiet_NaN());

    // Initialize per-day variables; S_v, N_v and O_v are only estimated
    assert(N_v_length <= (int)forcedS_v.size());
    for (int t = 0; t < N_v_length; t++)
    {
        P_A[t] = tsP_A;
        P_Amu[t] = tsP_Amu;
        P_A1[t] = tsP_A1;
        P_Ah[t] = tsP_Ah;
        P_df[t] = tsP_df;
        P_dff[t] = tsP_dff;
        N_v[t] = forcedS_v[t] * initNvFromSv;
        for (size_t genotype = 0; genotype < Genotypes::N(); ++genotype)
        {
            S_v_l[t * Genotypes::N() + genotype] = forcedS_v[t] * Genotypes::initialFreq(genotype);
            O_v_l[t * Genotypes::N() + genotype] = S_v_l[t * Genotypes::N() + genotype] * initOvFromSv;
        }
    }

    // Crude estimate of mosqEmergeRate: (1 - P_A(t) - P_df(t)) / (T * ρ_S) * S_T(t)
    mosqEmergeRate = forcedS_v;
    vectors::scale(mosqEmergeRate, initNv0FromSv);

    // All set up to drive simulation from forcedS_v
}

bool AnophelesModel::initIterate() { return true; }
//@}

void AnophelesModel::initVectorInterv(const scnXml::VectorSpeciesIntervention &elt, size_t instance)
{
    if (emergenceReduction.size() <= instance) emergenceReduction.resize(instance + 1);

    if (elt.getEmergenceReduction().present())
    {
        const scnXml::EmergenceReduction &elt2 = elt.getEmergenceReduction().get();
        if (elt2.getInitial() > 1.0) throw util::xml_scenario_error("emergenceReduction intervention: initial effect must be ≤ 1");
        emergenceReduction[instance].set(elt2.getInitial(), elt2.getDecay(), "emergenceReduction");
    }
    if (probAdditionalDeathSugarFeedingIntervs.size() <= instance) probAdditionalDeathSugarFeedingIntervs.resize(instance + 1);
    if (seekingDeathRateIntervs.size() <= instance) seekingDeathRateIntervs.resize(instance + 1);
    if (probDeathOvipositingIntervs.size() <= instance) probDeathOvipositingIntervs.resize(instance + 1);

    if (elt.getProbAdditionalDeathSugarFeeding().present())
    {
        const scnXml::ProbAdditionalDeathSugarFeeding &elt2 = elt.getProbAdditionalDeathSugarFeeding().get();
        if (elt2.getInitial() < -1.0) throw util::xml_scenario_error("probAdditionalDeathSugarFeeding intervention: initial effect must be ≥ -1");
        probAdditionalDeathSugarFeedingIntervs[instance].set(elt2.getInitial(), elt2.getDecay(), "seekingDeathRateIncrease");
    }
    if (elt.getSeekingDeathRateIncrease().present())
    {
        const scnXml::SeekingDeathRateIncrease &elt2 = elt.getSeekingDeathRateIncrease().get();
        if (elt2.getInitial() < -1.0) throw util::xml_scenario_error("seekingDeathRateIncrease intervention: initial effect must be ≥ -1");
        seekingDeathRateIntervs[instance].set(elt2.getInitial(), elt2.getDecay(), "seekingDeathRateIncrease");
    }
    if (elt.getProbDeathOvipositing().present())
    {
        const scnXml::ProbDeathOvipositing &elt2 = elt.getProbDeathOvipositing().get();
        if (elt2.getInitial() < 0.0 || elt2.getInitial() > 1.0)
            throw util::xml_scenario_error("probDeathOvipositing intrevention: initial effect must be in range [0,1]");
        probDeathOvipositingIntervs[instance].set(elt2.getInitial(), elt2.getDecay(), "probDeathOvipositing");
    }
}
void AnophelesModel::initVectorTrap(const scnXml::Description1 &desc, size_t instance)
{
    assert(trapParams.size() == instance); // if triggered, this is a code error not XML
    TrapParams params;
    params.relAvail = desc.getRelativeAvailability().getValue();
    params.availDecay = DecayFunction::makeObject(desc.getDecayOfAvailability(), "decayOfAvailability");
    trapParams.push_back(std::move(params));
}

void AnophelesModel::deployVectorPopInterv(LocalRng &rng, size_t instance)
{
    assert(instance < emergenceReduction.size());
    emergenceReduction[instance].deploy(rng, sim::now());
    // do same as in above function (of EmergenceModel)
    assert(instance < seekingDeathRateIntervs.size() && instance < probDeathOvipositingIntervs.size());
    probAdditionalDeathSugarFeedingIntervs[instance].deploy(rng, sim::now());
    seekingDeathRateIntervs[instance].deploy(rng, sim::now());
    probDeathOvipositingIntervs[instance].deploy(rng, sim::now());
}
void AnophelesModel::deployVectorTrap(LocalRng &rng, size_t species, size_t instance, double popSize, SimTime lifespan)
{
    assert(instance < trapParams.size());
    TrapData data;
    data.instance = instance;
    const PerHostAnophParams &anophParams = PerHostAnophParams::get(species);
    const double adultAvail = anophParams.entoAvailability->mean() * anophParams.entoAvailabilityFactor;
    data.initialAvail = popSize * adultAvail * trapParams[instance].relAvail;
    data.availHet = trapParams[instance].availDecay->hetSample(rng);
    data.deployTime = sim::now();
    data.expiry = sim::now() + lifespan;
    baitedTraps.push_back(std::move(data));
}

// Every sim::oneTS() days:
void AnophelesModel::advancePeriod(double sum_avail, double sigma_df, vector<double> &sigma_dif_i, vector<double> &sigma_dif_l, double sigma_dff, bool isDynamic)
{
    /* Largely equations correspond to Nakul Chitnis's model in
      "A mathematic model for the dynamics of malaria in
      mosquitoes feeding on a heterogeneous host population" [MMDM]
    section 2, 3.5-3.6, plus extensions to a non-autonomous case from
      "Nonautonomous Difference Equations for Malaria Dynamics
                   in a Mosquito Population" [NDEMD]

    We calculate EIR over a time step (one or five days) as:
      sum_{for t over days in step} σ_i[t] * s_v[t]
      = sum_... (N_v[t] * P_Ai[t] * P_B_i[t])/(T*N_i[t]) * S_v[t]/N_v[t]
      = sum_... P_Ai[t] * P_B_i[t] * S_v[t]
    (since T == 1 and N_i[t] == 1 for all t).

      P_Ai[t] = (1 - P_A[t]) α_i[t] / sum_{h in hosts} α_h[t]
    (letting N_h[t] == 1 for all h,t). The only part of this varying per-host is
      α_i[t] = host.entoAvailability (index, human.getAgeInYears())
      Let availDivisor[t] = (1 - P_A[t]) / (sum_{h in hosts} α_h[t] + μ_vA).

    Note that although the model allows α_i and P_B_i to vary per-day, they only
    vary per time step of the main simulation. Hence:
      EIR = (sum_{t=...} S_v[t] * availDivisor[t]) * α_i * P_B_i

    Since S_v[t] * availDivisor[t] does not vary per individual, we calculate this
    per time step of the main simulation as partialEIR:
      partialEIR = (sum_{t=...} S_v[t] * availDivisor[t])

    Hence calculateEIR() only needs to do the following:
      EIR = partialEIR * α_i * P_B_i
    */

    // -----  Calculate P_A, P_Ai, P_df, P_dif based on human pop  -----

    // ν_A: rate at which mosquitoes find hosts or die (i.e. leave host-seeking state
    double leaveRate = mosq.seekingDeathRate;
    for (const util::SimpleDecayingValue &increase : seekingDeathRateIntervs)
    {
        leaveRate *= 1.0 + increase.current_value(sim::ts0());
    }
    double mu_v_interv = leaveRate;
    leaveRate += sum_avail;

    // NON-HUMAN HOSTS INTERVENTIONS
    // Check if some nhh must be removed
    for (auto it = nhhInstances.begin(); it != nhhInstances.end();)
    {
        if (sim::ts0() >= it->second.expiry)
        {
            it = nhhInstances.erase(it);
            continue;
        }
        it++;
    }

    double modified_nhh_avail = 0.0;
    double modified_nhh_sigma_df = 0.0;
    double modified_nhh_sigma_dff = 0.0;

    map<string, Nhh> currentNhh = nhhInstances;

    for (auto it = reduceNhhAvailability.begin(); it != reduceNhhAvailability.end(); ++it)
    {
        for (const auto &decay : it->second)
        {
            if (currentNhh.count(it->first) != 0) // Check that the non-human hosts still exist
                currentNhh[it->first].avail_i *= 1.0 - decay.current_value(sim::ts0());
        }
    }

    for (auto it = reduceP_B_I.begin(); it != reduceP_B_I.end(); ++it)
    {
        for (const auto &decay : it->second)
        {
            if (currentNhh.count(it->first) != 0) // Check that the non-human hosts still exist
                currentNhh[it->first].P_B_I *= 1.0 - decay.current_value(sim::ts0());
        }
    }

    for (auto it = reduceP_C_I.begin(); it != reduceP_C_I.end(); ++it)
    {
        for (const auto &decay : it->second)
        {
            if (currentNhh.count(it->first) != 0) // Check that the non-human hosts still exist
                currentNhh[it->first].P_C_I *= 1.0 - decay.current_value(sim::ts0());
        }
    }

    for (auto it = reduceP_D_I.begin(); it != reduceP_D_I.end(); ++it)
    {
        for (const auto &decay : it->second)
        {
            if (currentNhh.count(it->first) != 0) // Check that the non-human hosts still exist
                currentNhh[it->first].P_D_I *= 1.0 - decay.current_value(sim::ts0());
        }
    }

    for (auto it = reduceFecundity.begin(); it != reduceFecundity.end(); ++it)
    {
        for (const auto &decay : it->second)
        {
            if (currentNhh.count(it->first) != 0) // Check that the non-human hosts still exist
                currentNhh[it->first].rel_fecundity *= 1.0 - decay.current_value(sim::ts0());
        }
    }

    for (auto it = currentNhh.begin(); it != currentNhh.end(); ++it)
    {
        modified_nhh_avail += it->second.avail_i;
        const double df = it->second.avail_i * it->second.P_B_I * it->second.P_C_I * it->second.P_D_I; // term in P_df series
        modified_nhh_sigma_df += df;
        modified_nhh_sigma_dff += df * it->second.rel_fecundity;
    }

    leaveRate += modified_nhh_avail;
    sigma_df += modified_nhh_sigma_df;
    sigma_dff += modified_nhh_sigma_dff;
    // NON-HUMAN HOSTS INTERVENTIONS

    for (auto it = baitedTraps.begin(); it != baitedTraps.end();)
    {
        if (sim::ts0() > it->expiry)
        {
            it = baitedTraps.erase(it);
            continue;
        }
        SimTime age = sim::ts0() - it->deployTime;
        double decayCoeff = it->availHet->eval(age);
        leaveRate += it->initialAvail * decayCoeff;
        // sigma_df doesn't change: mosquitoes do not survive traps
        it++;
    }

    // Calculate alpha_t for ATSB interventions with fixed target PA
    // =============================================================
    double pDeathSeeking = 0.0;
    for (const util::SimpleDecayingValue &pDeath : probAdditionalDeathSugarFeedingIntervs)
        pDeathSeeking += pDeath.current_value(sim::ts0());

    if(pDeathSeeking > 1.0)
         throw xml_scenario_error("VectorPop: the cumulative probability of death wile seeking is greater than 1 during the simulation");

    if(pDeathSeeking > 0.0)
    {
        int status;
        int iter = 0, max_iter = 20;

        double x0, a_t = 0.5; // guess

        struct quadratic_params params = { sum_avail, modified_nhh_avail, mu_v_interv, mosq.seekingDuration, pDeathSeeking};

        const gsl_root_solver_type *T = gsl_root_solver_method;
        gsl_root_solver *s = gsl_root_solver_alloc(T);

        #ifdef USE_DERIVATIVE
            gsl_function_fdf FDF;
            FDF.f = &quadratic_f;
            FDF.df = &quadratic_df;
            FDF.fdf = &quadratic_fdf;
            FDF.params = &params;

            gsl_root_fdfsolver_set (s, &FDF, a_t);
        #else
            gsl_function F;
            F.function = &quadratic_f;
            F.params = &params;

            gsl_root_fsolver_set (s, &F, 0.0, 1.0);
        #endif

        do
        {
            iter++;
            status = gsl_root_solver_iterate (s);
            x0 = a_t;
            a_t = gsl_root_solver_root (s);

            status = gsl_root_test_delta (a_t, x0, 0, 1e-4);
        }
        while (status == GSL_CONTINUE && iter < max_iter);

        gsl_root_solver_free (s);

        leaveRate += a_t;
    }
    // =============================================================

    // Probability of a mosquito not finding a host this day:
    double tsP_A = exp(-leaveRate * mosq.seekingDuration);
    double availDivisor = (1.0 - tsP_A) / leaveRate; // α_d

    // alphaE (α_E) is α_d * P_E, where P_E may be adjusted by interventions
    double alphaE = availDivisor * mosq.probMosqSurvivalOvipositing;
    for (const util::SimpleDecayingValue &pDeath : probDeathOvipositingIntervs)
    {
        alphaE *= 1.0 - pDeath.current_value(sim::ts0());
    }
    double tsP_df = sigma_df * alphaE;
    double tsP_dff = sigma_dff * alphaE;

    // from now, sigma_dif becomes P_dif (but we can't simply rename):
    vectors::scale(sigma_dif_i, alphaE);
    vectors::scale(sigma_dif_l, alphaE);

    // Summed per day:
    partialEIR_i.assign(WithinHost::Genotypes::N(), 0.0);
    partialEIR_l.assign(WithinHost::Genotypes::N(), 0.0);

    resetTSStats();

    // Computing for output only
    double tsP_Amu = (1 - tsP_A) * mosq.seekingDeathRate / (mosq.seekingDeathRate + sum_avail + modified_nhh_avail);
    double tsP_A1 = (1 - tsP_A) * sum_avail / (mosq.seekingDeathRate + sum_avail + modified_nhh_avail);
    double tsP_Ah = (1 - tsP_A) * modified_nhh_avail / (mosq.seekingDeathRate + sum_avail + modified_nhh_avail);

    // The code within the for loop needs to run per-day, wheras the main
    // simulation uses one or five day time steps.
    const SimTime nextTS = sim::ts0() + sim::oneTS();
    for (SimTime d0 = sim::ts0(); d0 < nextTS; d0 = d0 + sim::oneDay())
    {
        update(d0, tsP_A, tsP_Amu, tsP_A1, tsP_Ah, tsP_df, sigma_dif_i, sigma_dif_l, tsP_dff, isDynamic, partialEIR_i, partialEIR_l, availDivisor);
    }
}

void AnophelesModel::update(SimTime d0, double tsP_A, double tsP_Amu, double tsP_A1, double tsP_Ah, double tsP_df,
    const vector<double> &tsP_dif_i, const vector<double> &tsP_dif_l, double tsP_dff, bool isDynamic, 
    vector<double> &partialEIR_i, vector<double> &partialEIR_l, double EIR_factor)
{
    double interventionSurvival = 1.0;
    for (size_t i = 0; i < emergenceReduction.size(); ++i)
        interventionSurvival *= 1.0 - emergenceReduction[i].current_value(sim::ts0());

    int d1 = d0 + 1; //sim::oneDay(); // end of step

    // We add N_v_length so that we can use mod_nn() instead of mod().
    int d1Mod = d1 + N_v_length;
    assert(d1Mod >= N_v_length);
    // Indecies for end time, start time, and mosqRestDuration days before end time:
    int t1 = util::mod_nn(d1, N_v_length);
    int t0 = util::mod_nn(d0, N_v_length);
    int ttau = util::mod_nn(d1Mod - mosq.restDuration, N_v_length);

    // These only need to be calculated once per time step, but should be
    // present in each of the previous N_v_length - 1 positions of arrays.
    P_A[t1] = tsP_A;
    P_Amu[t1] = tsP_Amu;
    P_A1[t1] = tsP_A1;
    P_Ah[t1] = tsP_Ah;
    P_df[t1] = tsP_df;
    P_dff[t1] = tsP_dff;

    const size_t n = Genotypes::N();

    for (size_t i = 0; i < n; ++i)
    {
        P_dif_i[t1 * n + i] = tsP_dif_i[i];
        P_dif_l[t1 * n + i] = tsP_dif_l[i];
    }

    // BEGIN cache calculation: fArray, ftauArray, uninfected_v
    // Set up array with n in 1..θ_s−τ for f(d1Mod-n) (NDEMD eq. 1.6)
    for (int n = 1; n <= mosq.restDuration; n ++)
    {
        const int tn = util::mod_nn(d1Mod - n, N_v_length);
        fArray[n] = fArray[n - 1] * P_A[tn];
    }
    fArray[mosq.restDuration] += P_df[ttau];

    const int fAEnd = mosq.EIPDuration - mosq.restDuration;
    for (int n = mosq.restDuration + 1; n <= fAEnd; n++)
    {
        const int tn = util::mod_nn(d1Mod - n, N_v_length);
        fArray[n] = P_df[tn] * fArray[n - mosq.restDuration] + P_A[tn] * fArray[n - 1];
    }

    // Set up array with n in 1..θ_s−1 for f_τ(d1Mod-n) (NDEMD eq. 1.7)
    const int fProdEnd = mosq.restDuration * 2;
    for (int n = mosq.restDuration + 1; n <= fProdEnd; n++)
    {
        int tn = util::mod_nn(d1Mod - n, N_v_length);
        ftauArray[n] = ftauArray[n - 1] * P_A[tn];
    }
    ftauArray[fProdEnd] += P_df[util::mod_nn(d1Mod - fProdEnd, N_v_length)];

    for (int n = fProdEnd + 1; n < mosq.EIPDuration; n++)
    {
        int tn = util::mod_nn(d1Mod - n, N_v_length);
        ftauArray[n] = P_df[tn] * ftauArray[n - mosq.restDuration] + P_A[tn] * ftauArray[n - 1];
    }

    for (int d = 1; d < N_v_length; d++)
    {
        int t = util::mod_nn(d1Mod - d, N_v_length);
        double sum = N_v[t];
        for (size_t i = 0; i < n; ++i)
            sum -= (O_v_i[t * n + i] + O_v_l[t * n + i]); // .at(t, i);
        uninfected_v[d] = sum;
    }
    // END cache calculation: fArray, ftauArray, uninfected_v

    double total_S_v = 0.0;
    for (size_t g = 0; g < n; ++g)
    {
        // Num infected seeking mosquitoes is the new ones (those who were
        // uninfected tau days ago, started a feeding cycle then, survived and
        // got infected) + those who didn't find a host yesterday + those who
        // found a host tau days ago and survived a feeding cycle.
        // O_v.at(t1, g) = P_dif.at(ttau, g) * uninfected_v[mosq.restDuration] + P_A[t0] * O_v.at(t0, g) +
        //                        P_df[ttau] * O_v.at(ttau, g);
        O_v_i[t1 * n + g] = P_dif_i[ttau * n + g] * uninfected_v[mosq.restDuration] + P_A[t0] * O_v_i[t0 * n + g] + P_df[ttau] * O_v_i[ttau * n + g];
        O_v_l[t1 * n + g] = P_dif_l[ttau * n + g] * uninfected_v[mosq.restDuration] + P_A[t0] * O_v_l[t0 * n + g] + P_df[ttau] * O_v_l[ttau * n + g];            
        
        // BEGIN S_v
        double sum_i = 0.0, sum_l = 0.0;
        const int ts = d1Mod - mosq.EIPDuration;
        for (int l = 1; l < mosq.restDuration; l++)
        {
            const int tsl = util::mod_nn(ts - l, N_v_length); // index d1Mod - theta_s - l
            sum_i += P_dif_i[tsl * n + g] * P_df[ttau] * (uninfected_v[mosq.EIPDuration + l]) * ftauArray[mosq.EIPDuration + l - mosq.restDuration];
            sum_l += P_dif_l[tsl * n + g] * P_df[ttau] * (uninfected_v[mosq.EIPDuration + l]) * ftauArray[mosq.EIPDuration + l - mosq.restDuration];
        }

        const int tsm = util::mod_nn(ts, N_v_length); // index d1Mod - theta_s
        S_v_i[t1 * n + g] = P_dif_i[tsm * n + g] * fArray[mosq.EIPDuration - mosq.restDuration] * (uninfected_v[mosq.EIPDuration]) + sum_i + P_A[t0] * S_v_i[t0 * n + g] + P_df[ttau] * S_v_i[ttau * n + g];
        S_v_l[t1 * n + g] = P_dif_l[tsm * n + g] * fArray[mosq.EIPDuration - mosq.restDuration] * (uninfected_v[mosq.EIPDuration]) + sum_l + P_A[t0] * S_v_l[t0 * n + g] + P_df[ttau] * S_v_l[ttau * n + g];

        if (isDynamic)
        {
            // We cut-off transmission when no more than X mosquitos are infected to
            // allow true elimination in simulations. Unfortunately, it may cause problems with
            // trying to simulate extremely low transmission, such as an R_0 case.
            if (S_v_i[t1 * n + g] + S_v_l[t1 * n + g] <= mosq.minInfectedThreshold)
            {
                S_v_l[t1 * n + g] = 0.0;
                S_v_i[t1 * n + g] = 0.0; // Removing this will break unit tests
            }
        }
        
        partialEIR_i[g] += S_v_i[t1 * n + g] * EIR_factor;
        partialEIR_l[g] += S_v_l[t1 * n + g] * EIR_factor;
        total_S_v += S_v_i[t1 * n + g] + S_v_l[t1 * n + g];
        // END S_v
    }

    // We use time at end of step (i.e. start + 1) in index:
    int d5Year = util::mod_nn(d1, sim::fromYearsI(5));
    quinquennialS_v[d5Year] = total_S_v;

    const double nOvipositing = P_dff[ttau] * N_v[ttau]; // number ovipositing on this step
    const double newAdults = getEmergenceRate(d0, mosqEmergeRate, nOvipositing) * interventionSurvival;
    util::streamValidate(newAdults);

    // num seeking mosquitos is: new adults + those which didn't find a host
    // yesterday + those who found a host tau days ago and survived cycle:
    N_v[t1] = newAdults + P_A[t0] * N_v[t0] + nOvipositing;

    timeStep_N_v0 += newAdults;
}

// -----  Summary and intervention functions  -----
void AnophelesModel::changeEIRIntervention(const scnXml::NonVector &nonVectorData)
{
    // assert(nonVectorData.getEIRDaily().present()); // XML loading code should enforce this
    const scnXml::NonVector::EIRDailySequence &daily = nonVectorData.getEIRDaily();

    // const scnXml::DailyValues &seasM = seasonality.getDailyValues().get();
    // const scnXml::DailyValues::ValueSequence &daily = seasM.getValue();

    // The minimum EIR allowed in the array. The product of the average EIR and a constant.
    double minEIR = 0.01 * averageEIR(nonVectorData);

    if (daily.size() < static_cast<size_t>(sim::oneYear()))
        throw util::xml_scenario_error("entomology.anopheles.seasonality.dailyValues insufficient daily data for a year");

    vector<int> nDays(sim::inSteps(sim::fromDays(daily.size() - 1)) + 1, 0);

    interventionEIR.assign(nDays.size(), 0.0);

    size_t required_days = static_cast<size_t>((sim::endDate() - sim::startDate()) + 1);

    if (daily.size() < required_days)
        throw util::xml_scenario_error("Insufficient intervention phase EIR values provided");

    for (SimTime mpcday = sim::zero(), endDay = sim::fromDays(daily.size()); mpcday < endDay; mpcday = mpcday + sim::oneDay())
    {
        double EIRdaily = std::max(static_cast<double>(daily[mpcday]), minEIR);

        // istep is the time period to which the day is assigned.
        size_t istep = sim::inSteps(mpcday);
        nDays[istep]++;
        interventionEIR[istep] += EIRdaily;
    }
    // divide by number of records assigned to each interval (usually one per day)
    for (size_t i = 0; i < interventionEIR.size(); ++i)
    {
        interventionEIR[i] *= sim::oneTS() / static_cast<double>(nDays[i]);
    }
}

void AnophelesModel::uninfectVectors()
{
    for(size_t i=0; i<N_v_length * Genotypes::N(); i++)
    {
        O_v_i[i] = 0.0;
        O_v_l[i] = 0.0;
        S_v_i[i] = 0.0;
        S_v_l[i] = 0.0;
        P_dif_i[i] = 0.0;
        P_dif_l[i] = 0.0;
    }
}

double sum1(const std::vector<double> &arr, int end, int N_v_length)
{
    double val = 0.0;
    // Last time step ended at sim::now(). Values are stored per day, and for
    // the last time step values at sim::now() and four previos were set.
    for (int d1 = end - sim::oneTS(); d1 < end; d1++)
    {
        val += arr[util::mod_nn(d1, N_v_length)];
    }
    return val / sim::oneTS();
}

double sum2(const std::vector<double> &arr, int end, int N_v_length)
{
    double val = 0.0;
    // Last time step ended at sim::now(). Values are stored per day, and for
    // the last time step values at sim::now() and four previos were set.
    for (int d1 = end - sim::oneTS(); d1 < end; d1++)
    {
        int i1 = util::mod_nn(d1, N_v_length);
        for (size_t g = 0; g < Genotypes::N(); ++g)
        {
            val += arr[i1 * Genotypes::N() + g]; //.at(i1, g);
        }
    }
    return val / sim::oneTS();
}

double sum3(const std::vector<double> &arr, size_t g, int end, int N_v_length)
{
    double val = 0.0;
    // Last time step ended at sim::now(). Values are stored per day, and for
    // the last time step values at sim::now() and four previos were set.
    for (int d1 = end - sim::oneTS(); d1 < end; d1++)
    {
        val += arr[util::mod_nn(d1, N_v_length) * Genotypes::N() + g];// .at(mod_nn(d1, N_v_length), g);
    }
    return val / sim::oneTS();
}

double AnophelesModel::getLastVecStat(VecStat vs) const
{
    // Last time step ended at sim::now(). Values are stored per day, and for
    // the last time step values at sim::now() and four previos were set.
    // One plus last, plus (0 mod N_v_length) to avoid negatives:
    int end = sim::now() + 1 + N_v_length;
    switch (vs)
    {
        case PA: return sum1(P_A, end, N_v_length);
        case PDF: return sum1(P_df, end, N_v_length);
        case PDIF: return sum2(P_dif_i, end, N_v_length) + sum2(P_dif_l, end, N_v_length);
        case NV: return sum1(N_v, end, N_v_length);
        case OV: return sum2(O_v_i, end, N_v_length) + sum2(O_v_l, end, N_v_length);
        case SV: return sum2(S_v_i, end, N_v_length) + sum2(S_v_l, end, N_v_length);
        case PAmu: return sum1(P_Amu, end, N_v_length);
        case PA1: return sum1(P_A1, end, N_v_length);
        case PAh: return sum1(P_Ah, end, N_v_length);
        default: throw SWITCH_DEFAULT_EXCEPTION;
    }
}

void AnophelesModel::summarize(size_t species) const
{
    // Last time step ended at sim::now(). Values are stored per day, and for
    // the last time step values at sim::now() and four previos were set.
    // One plus last, plus (0 mod N_v_length) to avoid negatives:
    int end = sim::now() + 1 + N_v_length;
    mon::reportStatMSF(mon::MVF_LAST_NV0, species, getLastN_v0());
    mon::reportStatMSF(mon::MVF_LAST_NV, species, sum1(N_v, end, N_v_length));
    for (size_t g = 0; g < Genotypes::N(); ++g)
    {
        mon::reportStatMSGF(mon::MVF_LAST_OV, species, g, sum3(O_v_i, g, end, N_v_length) + sum3(O_v_l, g, end, N_v_length));
        mon::reportStatMSGF(mon::MVF_LAST_SV, species, g, sum3(S_v_i, g, end, N_v_length) + sum3(S_v_l, g, end, N_v_length));
    }
}

} // namespace Anopheles
} // namespace Transmission
} // namespace OM
