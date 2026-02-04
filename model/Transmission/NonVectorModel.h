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

#ifndef Hmod_NonVectorModel
#define Hmod_NonVectorModel

#include "TransmissionModel.h"
#include "util/vectors.h"
#include "Host/Human.h"
#include "mon/info.h"

namespace scnXml
{
class NonVector;
}

namespace OM
{
namespace Transmission
{

//! Base transmission model, as used in Phase A
class NonVectorModel : public TransmissionModel
{
public:
    NonVectorModel(vector<double> initEIR, int interventionMode, const scnXml::Entomology &entoData, const scnXml::NonVector &nonVectorData, int eipDuration)
        : TransmissionModel(std::move(initEIR), interventionMode, 1 /*this model doesn't support multiple genotypes*/)
        , nSpore(sim::fromDays(eipDuration))
    {
        laggedKappa.resize(sim::inSteps(nSpore) + 1, 0.0);

        vector<int> nDays(sim::stepsPerYear(), 0);
        // The minimum EIR allowed in the array. The product of the average EIR and a constant.
        double minEIR = min_EIR_mult * averageEIR(nonVectorData);

        const scnXml::NonVector::EIRDailySequence &daily = nonVectorData.getEIRDaily();
        if (daily.size() < static_cast<size_t>(sim::oneYear()))
            throw util::xml_scenario_error("insufficient EIRDaily data for a year");

        for (SimTime mpcday = sim::zero(), endDay = sim::fromDays(daily.size()); mpcday < endDay; mpcday = mpcday + sim::oneDay())
        {
            double EIRdaily = std::max(static_cast<double>(daily[mpcday]), minEIR);

            // Index 0 of initialisationEIR refers to the EIR affecting the
            // first day(s) of the year. Correspondingly, the first 1 or 5 values
            // of EIRDaily affect this (1- or 5-day) time-step.
            size_t i = mod_nn(sim::inSteps(mpcday), sim::stepsPerYear());

            nDays[i] += 1;
            initialisationEIR[i] += EIRdaily;
        }

        // Calculate total annual EIR
        // divide by number of records assigned to each interval (usually one per day)
        for (size_t indTS = 0; indTS < sim::stepsPerYear(); indTS += 1)
        {
            initialisationEIR[indTS] *= sim::oneTS() / static_cast<double>(nDays[indTS]);
            annualEIR += initialisationEIR[indTS];
        }

        initialKappa.assign(sim::inSteps(sim::fromYearsI(nYearsWarmupData)), 0.0);
    }

    virtual ~NonVectorModel() {}

    virtual void init2(const vector<Host::Human> &population) { simulationMode = forcedEIR; } // no set-up needed; just indicate we're ready to roll

    virtual void scaleEIR(double factor)
    {
        util::vectors::scale(initialisationEIR, factor);
        annualEIR = util::vectors::sum(initialisationEIR);
    }

    virtual SimTime initIterate()
    {
        simulationMode = interventionMode;
        if (simulationMode != dynamicEIR) { return sim::zero(); }

        // initialKappa is used in calculateEIR
        assert(initialKappa.size() >= sim::stepsPerYear());
        assert(mod_nn(initialKappa.size(), sim::stepsPerYear()) == 0);
        for (size_t i = sim::stepsPerYear(); i < initialKappa.size(); ++i)
        {
            initialKappa[mod_nn(i, sim::stepsPerYear())] += initialKappa[i];
        }
        double factor = static_cast<double>(sim::stepsPerYear()) / static_cast<double>(initialKappa.size());
        initialKappa.resize(sim::stepsPerYear());
        for (size_t i = 0; i < initialKappa.size(); ++i)
        {
            initialKappa[i] *= factor;
            // error check:
            if (!(initialKappa[i] > 0.0)) // if not positive
                throw TRACED_EXCEPTION("initialKappa is invalid", util::Error::InitialKappa);
        }

        return sim::zero(); // nothing to do
    }

    virtual void changeEIRIntervention(const scnXml::NonVector &nonVectorData)
    {
        // Note: requires sim::intervTime() >= sim::zero(), but this can only be
        // called in intervention period anyway.
        simulationMode = transientEIRknown;

        if (nSpore != sim::fromDays(nonVectorData.getEipDuration()))
            throw util::xml_scenario_error("change-of-EIR intervention cannot change EIP duration");

        const scnXml::NonVector::EIRDailySequence &daily = nonVectorData.getEIRDaily();
        vector<int> nDays(sim::inSteps(sim::fromDays(daily.size() - 1)) + 1, 0);
        interventionEIR.assign(nDays.size(), 0.0);
        size_t required_days = static_cast<size_t>((sim::endDate() - sim::startDate()) + 1);
        if (daily.size() < required_days)
        {
            cerr << "Days: " << daily.size() << "\nIntervals: " << nDays.size() << "\nRequired: " << required_days << endl;
            throw util::xml_scenario_error("Insufficient intervention phase EIR values provided");
        }
        // The minimum EIR allowed in the array. The product of the average EIR and a constant.
        double minEIR = min_EIR_mult * averageEIR(nonVectorData);
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

        // I've no idea what this should be, so until someone asks it can be NaN.
        // It was -9.99 and later 0.0. It could of course be recalculated from interventionEIR.
        annualEIR = numeric_limits<double>::quiet_NaN();
    }

    virtual void uninfectVectors()
    {
        if (simulationMode != dynamicEIR) cerr << "Warning: uninfectVectors is not efficacious with forced EIR" << endl;
        // reset history of human infectivity, which scales dynamic EIR:
        laggedKappa.assign(laggedKappa.size(), 0.0);
    }

    virtual double updateKappa(const vector<Host::Human> &population)
    {
        double currentKappa = TransmissionModel::updateKappa(population);
        if (simulationMode == forcedEIR) { initialKappa[sim::moduloSteps(sim::ts1(), initialKappa.size())] = currentKappa; }
        return currentKappa;
    }

    virtual void calculateEIR(Host::Human &human, double ageYears, vector<double> &EIR_i, vector<double> &EIR_l) const
    {
        EIR_i.resize(1);
        EIR_l.resize(1); // no support for per-genotype tracking in this model (possible, but we're lazy)
        // where the full model, with estimates of human mosquito transmission is in use, use this:
        if (simulationMode == forcedEIR) { EIR_l[0] = initialisationEIR[sim::moduloYearSteps(sim::ts0())]; }
        else if (simulationMode == transientEIRknown)
        {
            // where the EIR for the intervention phase is known, obtain this from
            // the interventionEIR array
            EIR_l[0] = interventionEIR[sim::inSteps(sim::intervTime())];
        }
        else if (simulationMode == dynamicEIR)
        {
            EIR_l[0] = initialisationEIR[sim::moduloYearSteps(sim::ts0())];
            if (sim::intervTime() >= sim::zero())
            {
                // we modulate the initialization based on the human infectiousness time steps ago in the
                // simulation relative to infectiousness at the same time-of-year, pre-intervention.
                // nspore gives the sporozoite development delay.
                size_t t = sim::inSteps(sim::ts1() - nSpore);
                EIR_l[0] *= laggedKappa[mod_nn(t, laggedKappa.size())] / initialKappa[mod_nn(t, sim::stepsPerYear())];
            }
        }
        else
        {
            throw util::xml_scenario_error("Invalid simulation mode");
        }
#ifndef NDEBUG
        if (!(std::isfinite)(EIR_l[0]))
        {
            size_t t = sim::inSteps(sim::ts1() - nSpore);
            ostringstream msg;
            msg << "Error: non-vect eir is: " << EIR_l[0] << "\nlaggedKappa:\t" << laggedKappa[mod_nn(t, laggedKappa.size())]
                << "\ninitialKappa:\t" << initialKappa[mod_nn(t, sim::stepsPerYear())] << endl;
            throw TRACED_EXCEPTION(msg.str(), util::Error::InitialKappa);
        }
#endif
        EIR_l[0] *= human.perHostTransmission.relativeAvailabilityHetAge(ageYears);

        auto ag = human.monitoringAgeGroup.i();
        auto cs = human.getCohortSet();
        mon::reportStatMACGF(mon::MVF_INOCS, ag, cs, 0, EIR_l[0]);
    }

private:
    virtual void checkpoint(istream &stream)
    {
        TransmissionModel::checkpoint(stream);
        nSpore &stream;
        interventionEIR &stream;
        initialKappa &stream;
    }

    virtual void checkpoint(ostream &stream)
    {
        TransmissionModel::checkpoint(stream);
        nSpore &stream;
        interventionEIR &stream;
        initialKappa &stream;
    }

    //! multiplier used to calculate a positive EIR value where the measured value is zero
    /*
      0.01 was old pv(30) Now a constant. min_EIR_mult multiplies the average EIR to obtain a value used for the EIR during periods when
      it is too low to be measureable. The value of 0.01 was old pv(30) Now a constant. 0.01 was old pv(30) Now a constant.
    */
    const double min_EIR_mult = 0.01;

    /** @brief Variables set by constructor.
     *
     * There shouldn't be any need to checkpoint these, at least before
     * interventions take effect. */
    //@{
    //! Variance of Infection Rate according to fielddata
    const double totalInfectionrateVariance = 1.0;

    const int nYearsWarmupData = 5;

    //! The duration of sporogony in time steps
    // doesn't need checkpointing
    SimTime nSpore = sim::never();
    //@}

    /** EIR per time interval during the intervention period. Value at index
     * sim::intervTime().inSteps() used each time-step.
     *
     * Units: inoculations per adult per time step */
    vector<double> interventionEIR;

    /** When simulationMode == dynamicEIR, this is the annual cycle of kappa
     * from the warmup phase and has length 1 year (in time steps).
     *
     * When simulationMode == equilibriumMode, this may be multiple years long and
     * is used to collect values of kappa (human infectiousness).
     *
     * In either case, sim::ts0().moduloSteps(initialKappa.size()) is the index
     * for the current infectiousness during updates. */
    vector<double> initialKappa;
};
} // namespace Transmission
} // namespace OM
#endif
