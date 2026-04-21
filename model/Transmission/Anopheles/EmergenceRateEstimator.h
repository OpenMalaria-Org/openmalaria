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

#ifndef Hmod_EmergenceRateEstimator
#define Hmod_EmergenceRateEstimator

#include "Global.h"
#include "Transmission/PerHost.h"
#include "util/SimpleDecayingValue.h"
#include "util/vectors.h"
#include "util/CommandLine.h"
#include "util/errors.h"

#include "CalcNvO.h"

#include <vector>
#include <limits>

namespace OM {
namespace Transmission {
namespace Anopheles {
    using std::numeric_limits;
    using namespace OM::util;

inline int argmax(const std::vector<double> &vec)
{
    int imax = 0;
    double max = 0.0;
    for(size_t i=0; i<vec.size(); i++)
    {
        double v = vec[i];
        if(v >= max)
        {
            max = v;
            imax = i;
        }
    }
    return imax;
}

inline double findAngle(const double EIRRotageAngle, const vector<double> & FSCoeffic, const std::vector<double> &sim)
{
    std::vector<double> temp(sim.size(), 0.0);

    double delta = 2.0 * M_PI / 365.0;

    double min = std::numeric_limits<double>::infinity();
    double minAngle = 0.0;
    for(double angle=-M_PI; angle<M_PI; angle+=delta)
    {
        vectors::expIDFT(temp, FSCoeffic, EIRRotageAngle + angle);

        // Minimize l1-norm
        double sum = 0.0;
        for(int i=0; i<sim::oneYear(); i++)
        {
            double v = fabs(temp[i] - sim[i]);
            sum += v*v;
        }

        sum = sqrtf(sum);
        if(sum < min)
        {
            min = sum;
            minAngle = angle;
            // cout << angle << " " << min << " " << sum << endl;
        }

        // Or minimize peaks offset
        // int m1 = argmax(temp);
        // int m2 = argmax(sim);
        // int offset = abs(m1-m2);

        // if(offset < min)
        // {
        //     min = offset;
        //     minAngle = angle;
        // }

    }
    return minAngle;
}

class EmergenceRateEstimator
{
public:
    virtual ~EmergenceRateEstimator() = default;

    // Interface: returns true if further fitting iterations are needed.
    virtual bool estimate(AnophelesModel& m, const vector<double> &laggedKappa, double meanAvail) = 0;
};

class EmergenceRateAdaptiveFitter final : public EmergenceRateEstimator
{
public:
    EmergenceRateAdaptiveFitter(AnophelesModel &m) : scaleFactor(1.0), rotated(false), scaled(false)
    {
        // usually around 20 days; no real analysis for effect of changing EIPDuration or mosqRestDuration
        shiftAngle = m.EIRRotateAngle - (m.mosq.EIPDuration + 10) / 365. * 2. *M_PI; 
    }

    bool estimate(AnophelesModel &m, const vector<double> &laggedKappa, double meanAvail) override
    {
        std::vector<double> avgAnnualS_v(sim::oneYear(), 0.0);
        for (SimTime i = sim::fromYearsI(4); i < sim::fromYearsI(5); i = i + sim::oneDay())
        {
            avgAnnualS_v[mod_nn(i, sim::oneYear())] = m.quinquennialS_v[i];
        }

        double factor = vectors::sum(m.forcedS_v) / vectors::sum(avgAnnualS_v);

        // cout << "check: " << vectors::sum(forcedS_v) << " " << vectors::sum(avgAnnualS_v) << endl;
        // cout << "Pre-calced Sv, dynamic Sv:\t"<<sumAnnualForcedS_v<<'\t'<<vectors::sum(annualS_v)<<endl;
        if (!(factor > 1e-6 && factor < 1e6))
        {
            if (factor > 1e6 && vectors::sum(m.quinquennialS_v) < 1e-3)
            {
                throw util::base_exception("Simulated S_v is approx 0 (i.e.\
     mosquitoes are not infectious, before interventions). Simulator cannot handle this; perhaps\
     increase EIR or change the entomology model.",
                                           util::Error::VectorFitting);
            }
            if (vectors::sum(m.forcedS_v) == 0.0)
            {
                return false; // no EIR desired: nothing to do
            }
            cerr << "Input S_v for this vector:\t" << vectors::sum(m.forcedS_v) << endl;
            cerr << "Simulated S_v:\t\t\t" << vectors::sum(m.quinquennialS_v) / 5.0 << endl;
            throw TRACED_EXCEPTION("factor out of bounds", util::Error::VectorFitting);
        }

        const double LIMIT = 0.1;

        if (fabs(factor - 1.0) > LIMIT)
        {
            scaled = false;
            double factorDiff = (scaleFactor * factor - scaleFactor) * 1.0;
            scaleFactor += factorDiff;
        }
        else
            scaled = true;

        double rAngle = findAngle(m.EIRRotateAngle, m.FSCoeffic, avgAnnualS_v);
        shiftAngle += rAngle;
        rotated = true;

        // Compute forced_sv from the Fourrier Coeffs EIR
        vectors::expIDFT(m.mosqEmergeRate, m.FSCoeffic, -shiftAngle);

        // Scale the vector according to initSvFromEIR and initNv0FromSv to get the mosqEmergerate
        // scaleFactor scales the vector to correct the ratio between simulated and input EIR
        vectors::scale(m.mosqEmergeRate, scaleFactor * m.initSvFromEIR * m.initNv0FromSv);

        // initNvFromSv *= scaleFactor;     //(not currently used)

        // What factor exactly these should be scaled by isn't obvious; in any case
        // they should reach stable values quickly.
        m.scale(factor);

        return !(scaled && rotated);
    }

private:
    double scaleFactor, shiftAngle;
    bool rotated, scaled;
};

class EmergenceRateSolver final : public EmergenceRateEstimator
{
public:
    EmergenceRateSolver(AnophelesModel &m, int populationSize) : populationSize(populationSize) {}

    bool estimate(AnophelesModel &m, const vector<double> &laggedKappa, double meanAvail) override
    {
        std::vector<double> &mosqEmergeRateVector = m.mosqEmergeRate;
        int thetap = sim::DAYS_IN_YEAR;
        int tau = m.mosq.restDuration;
        int thetas = m.mosq.EIPDuration;
        std::vector<double> Ni{static_cast<double>(populationSize)};
        std::vector<double> alphai{m.initAvail * meanAvail};
        double muvA = m.mosq.seekingDeathRate;
        double thetad = m.mosq.seekingDuration;
        std::vector<double> PBi{m.mosq.probBiting};
        std::vector<double> PCi{m.mosq.probFindRestSite};
        std::vector<double> PDi{m.mosq.probResting};
        std::vector<double> PEi{m.mosq.probOvipositing};
        std::vector<double> &svInit = m.forcedS_v;
        std::vector<double> NvOguess(m.mosqEmergeRate.begin(), m.mosqEmergeRate.end());

        // Append active NHH host types (malaria-non-susceptible)
        for (const auto &[name, nhh] : m.nhhInstances) {
            // Interpret nhh.avail_i as total availability α_i * N_i
            Ni.push_back(1.0);
            alphai.push_back(nhh.avail_i);

            PBi.push_back(nhh.P_B_I);
            PCi.push_back(nhh.P_C_I);
            PDi.push_back(nhh.P_D_I);
            PEi.push_back(m.mosq.probOvipositing); // Using the same probability for nhh
        }

        int nHost = static_cast<int>(Ni.size());
        int nMalHost = 1; // only humans are malaria-susceptible here

        std::vector<double> humanKappa(thetap, 0.0);

        // Linear interpolation 73 -> daily (5-day blocks)
        for (size_t i = 0; i < laggedKappa.size(); ++i) {
            const size_t start = i * 5;
            const size_t end   = std::min(start + 5, static_cast<size_t>(thetap));

            const double v0 = laggedKappa[i];
            const double v1 = laggedKappa[(i + 1) % laggedKappa.size()];

            for (size_t d = start; d < end; ++d) {
                const double t = double(d - start) / 5.0;
                humanKappa[d] = (1.0 - t) * v0 + t * v1;
            }
        }

        std::vector<double> KviInit(static_cast<size_t>(thetap) * nHost, 0.0);

        // Fill all malaria-susceptible human groups with the humanKappa series.
        // If later each group has its own kappa, fill each column separately.
        for (int i = 0; i < nMalHost; ++i) {
            for (int k = 0; k < thetap; ++k) {
                KviInit[static_cast<size_t>(k) * nHost + i] = humanKappa[k];
            }
        }

        CalcInitMosqEmergeRate(mosqEmergeRateVector, 
            thetap, // $\theta_p$: daysInYear
            tau, // $\tau$: mosqRestDuration
            thetas, // $\theta_s$: EIPDuration
            nHost, // $n$: nHostTypes
            nMalHost, // $m$: nMalHostTypes
            Ni, // $N_i$: popSize (length n)
            alphai, // $\alpha_i$: hostAvailabilityRate	(length n)
            muvA, // $\mu_{vA}: mosqSeekingDeathRate
            thetad, // $\theta_d$: mosqSeekingDuration
            PBi, // $P_{B_i}$: mosqProbBiting (length n)
            PCi, // $P_{C_i}$: mosqProbFindRestSite (length n)
            PDi, // $P_{D_i}$: mosqProbResting (length n)
            PEi, // $P_{E_i}$: mosqProbOvipositing
            KviInit, // Kvi (size n * thetap)
            svInit // Sv (length n)
        );
        
        return false;
    }

private:
    int populationSize;
};

}
}
}
#endif
