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

#ifndef Hmod_TransmissionModel
#define Hmod_TransmissionModel

#include <string.h>
#include <fstream>

#include "Global.h"
#include "schema/interventions.h"
#include "util/errors.h"

#include "Transmission/PerHost.h"
#include "Population.h"
#include "Host/WithinHost/WHInterface.h"
#include "Host/WithinHost/Genotypes.h"
#include "mon/Continuous.h"
#include "mon/info.h"
#include "util/StreamValidator.h"
#include "util/CommandLine.h"
#include "util/vectors.h"
#include "util/ModelOptions.h"

#include <cmath>
#include <cfloat>
#include <gsl/gsl_vector.h>

using namespace std::placeholders;

namespace OM
{
class Summary;
class Population;
namespace Host
{
class Human;
}
namespace Transmission
{

inline double averageEIR(const scnXml::NonVector &nonVectorData)
{
    // Calculates the arithmetic mean of the whole daily EIR vector read from the .XML file
    double valaverageEIR = 0.0;
    size_t i = 0;
    for (const scnXml::NonVector::EIRDailySequence &daily = nonVectorData.getEIRDaily(); i < daily.size(); ++i)
    {
        valaverageEIR += (double)daily[i];
    }
    if (i == 0) throw util::xml_scenario_error("no EIRDaily values given"); // pedantic check
    return valaverageEIR / i;
}

inline double averageEIR(const scnXml::DailyValues &vectorData)
{
    // Calculates the arithmetic mean of the whole daily EIR vector read from the .XML file
    double valaverageEIR = 0.0;
    size_t i = 0;
    for (const scnXml::DailyValues::ValueSequence &daily = vectorData.getValue(); i < daily.size(); ++i)
    {
        valaverageEIR += (double)daily[i];
    }
    if (i == 0) throw util::xml_scenario_error("no daily values given"); // pedantic check
    return valaverageEIR / i;
}

/** Variable describing current simulation mode. */
enum SimulationMode
{
    /** Initial mode. Indicates that initialisation still needs to happen (i.e.
     * it is an error if this mode is still set when getEIR is called). */
    preInit,
    /** Equilibrium mode (i.e. just apply a fixed EIR)
     *
     * This is used for the warm-up period and if we want to separate direct
     * effect of an intervention from indirect effects via transmission
     * intensity. The seasonal pattern and intensity of the EIR do not change
     * over years.
     *
     * In this mode vector calculations aren't run. */
    forcedEIR,
    /** Transient EIR known
     *
     * This is used to simulate an intervention that changes EIR, and where we
     * have measurements of the EIR over time during the intervention period.
     */
    transientEIRknown,
    /** EIR changes
     *
     * The simulation is driven by the EIR which changes dynamically during the
     * intervention phase as a function of the characteristics of the
     * interventions.
     *
     * Dependending on whether the Vector or NonVector model is in use, this
     * EIR may be calculated from a mosquito emergence rate or be an input EIR
     * scaled by the relative infectiousness of the humans. */
    dynamicEIR,
};

/// Abstract base class, defines behaviour of transmission models
class TransmissionModel
{
protected:
    /// Reads all entomological parameters from the input datafile.
    /// @param entoData input configuration for model
    /// @param nGenotypes number of genotypes transmission model is using
    TransmissionModel(vector<double> initialisationEIR, int interventionMode, size_t nGenotypes)
        : initialisationEIR(std::move(initialisationEIR))
        , simulationMode(forcedEIR)
        , interventionMode(interventionMode)
        , laggedKappa(1, 0.0)
        , // if using non-vector model, it will resize this
        annualEIR(0.0)
        , _annualAverageKappa(numeric_limits<double>::signaling_NaN())
        , _sumAnnualKappa(0.0)
        , tsAdultEIR(0.0)
        , surveyInputEIR(0.0)
        , surveySimulatedEIR(0.0)
        , adultAge(PerHost::adultAge())
        , numTransmittingHumans(0)
    {
        // Set VACCINE_GENOTYPE option
        opt_vaccine_genotype = util::ModelOptions::option (util::VACCINE_GENOTYPE);

        using mon::Continuous;
        Continuous.registerCallback("input EIR", "\tinput EIR", std::bind(&TransmissionModel::ctsCbInputEIR, this, _1));
        Continuous.registerCallback("simulated EIR", "\tsimulated EIR", std::bind(&TransmissionModel::ctsCbSimulatedEIR, this, _1));
        Continuous.registerCallback("human infectiousness", "\thuman infectiousness", std::bind(&TransmissionModel::ctsCbKappa, this, _1));
        Continuous.registerCallback("num transmitting humans", "\tnum transmitting humans", std::bind(&TransmissionModel::ctsCbNumTransmittingHumans, this, _1));
    }

public:
    //! Deallocate memory for TransmissionModel parameters and clean up
    virtual ~TransmissionModel() {}

    /** Extra initialisation when not loading from a checkpoint, requiring
     * information from the human population structure. */
    virtual void init2(const vector<Host::Human> &population) = 0;

    /// Checkpointing
    template <class S>
    void operator&(S &stream)
    {
        checkpoint(stream);
    }
    //@}

    /** Set some summary items.
     *
     * Overriding functions should call this base version too. */
    virtual void summarize()
    {
        mon::reportStatMF(mon::MVF_NUM_TRANSMIT, laggedKappa[sim::moduloSteps(sim::now(), laggedKappa.size())]);
        mon::reportStatMF(mon::MVF_ANN_AVG_K, _annualAverageKappa);

        if (!mon::isReported()) return; // cannot use counters below when not reporting

        double duration = sim::inSteps(sim::now() - lastSurveyTime);
        if (duration > 0.0)
        {
            mon::reportStatMF(mon::MVF_INPUT_EIR, surveyInputEIR / duration);
            mon::reportStatMF(mon::MVF_SIM_EIR, surveySimulatedEIR / duration);
            mon::reportStatMF(mon::MVF_SIM_EIR_INTRODUCED, surveySimulatedEIR_i / duration);
            mon::reportStatMF(mon::MVF_SIM_EIR_INDIGENOUS, surveySimulatedEIR_l / duration);
        }

        surveyInputEIR = 0.0;
        surveySimulatedEIR = 0.0;
        surveySimulatedEIR_i = 0.0;
        surveySimulatedEIR_l = 0.0;
        lastSurveyTime = sim::now();
    }

    /** Scale the EIR used by the model.
     *
     * EIR is scaled in memory (so will affect this simulation).
     * XML data is not touched. */
    virtual void scaleEIR(double factor) = 0;

    virtual SimTime initIterate() = 0;

    /** Needs to be called each step of the simulation before Human::update().
     *
     * when the vector model is used this updates mosquito populations. */
    virtual void vectorUpdate(const vector<Host::Human> &population){};

    virtual void changeEIRIntervention(const scnXml::NonVector &) = 0;

    /** Does per-time-step updates and returns the EIR (inoculation rate per host
     * per time step). Should be called exactly once per time-step (at least,
     * during the intervention period when ITNs may be in use).
     *
     * Non-vector:
     * During the pre-intervention phase, the EIR is forced, using values from
     * the XML file. During the main simulation phase, it may be calculated or
     * obtained from data in the XML file.
     *
     * Vector:
     * During vector initialisation phase, EIR is forced based on the EIR given
     * in the XML file as a Fourier Series. After endVectorInitPeriod() is called
     * the simulation switches to using dynamic EIR. advanceStep _must_ be
     * called before this function in order to return the correct value.
     *
     * @param human A reference to the human who's EIR is being calculated.
     *    The human's "per host transmission" potentially needs updating.
     * @param age Age of the human in time units
     * @param ageYears Age of the human in years
     * @param EIR Out-vector of EIR per parasite genotype. The length is also set
     *    by the called function. Where genotype tracking is not supported (e.g.
     *    the non-vector model), the length is set to one.
     * @returns the sum of EIR across genotypes
     */
    double getEIR(Host::Human &human, SimTime age, double ageYears, vector<double> &EIR_i, vector<double> &EIR_l)
    {
        /* For the NonVector model, the EIR should just be multiplied by the
         * availability. For the Vector model, the availability is also required
         * for internal calculations, but again the EIR should be multiplied by the
         * availability. */
        calculateEIR(human, ageYears, EIR_i, EIR_l);

        double sum_EIR_i = util::vectors::sum(EIR_i);
        double sum_EIR_l = util::vectors::sum(EIR_l);

        double allEIR = sum_EIR_i + sum_EIR_l;
        if (age >= adultAge)
        {
            tsAdultEntoInocs_i += sum_EIR_i;
            tsAdultEntoInocs_l += sum_EIR_l;
            tsAdultEntoInocs += allEIR;
            tsNumAdults += 1;
        }
        return allEIR;
    }

    /** Remove all current infections to mosquitoes, such that without re-
     * infection, humans will then be exposed to zero EIR. */
    virtual void uninfectVectors() = 0;

    /** Needs to be called each time-step after Human::update() to update summary
     * statististics related to transmission. Also returns kappa (the average
     * human infectiousness weighted by availability to mosquitoes). */
    virtual double updateKappa(const vector<Host::Human> &population)
    {
        // We calculate kappa for output and the non-vector model.
        double sumWt_kappa = 0.0;
        double sumWeight = 0.0;
        numTransmittingHumans = 0;

        for (const Host::Human &human : population)
        {
            // NOTE: calculate availability relative to age at end of time step;
            // not my preference but consistent with TransmissionModel::getEIR().
            const double avail = human.perHostTransmission.relativeAvailabilityHetAge(sim::inYears(human.age(sim::ts1())));
            sumWeight += avail;

            vector<double> probTransGenotype_i(WithinHost::Genotypes::N());
            vector<double> probTransGenotype_l(WithinHost::Genotypes::N());
            const double pTransmit = human.withinHostModel->probTransmissionToMosquito(probTransGenotype_i, probTransGenotype_l);

            double riskTrans = 0.0;

            // Only to be consistent with old simulation runs when set to false
            // Setting this option to true will only affect reporting
            if (opt_vaccine_genotype == false)
                riskTrans = avail * pTransmit * human.vaccine.getFactor(interventions::Vaccine::TBV);
            else
            {
                for (size_t g = 0; g < WithinHost::Genotypes::N(); ++g)
                    riskTrans += (probTransGenotype_i[g] + probTransGenotype_l[g]) * human.vaccine.getFactor(interventions::Vaccine::TBV, g);
                riskTrans *= avail;
            }

            if (riskTrans > 0.0) ++numTransmittingHumans;

            sumWt_kappa += riskTrans;
        }

        size_t lKMod = sim::moduloSteps(sim::ts1(), laggedKappa.size()); // now
        if (population.size() == 0)
        {                             // this is valid
            laggedKappa[lKMod] = 0.0; // no humans: no infectiousness
        }
        else
        {
            if (!(sumWeight > DBL_MIN * 10.0))
            { // if approx. eq. 0, negative or an NaN
                ostringstream msg;
                msg << "sumWeight is invalid: " << sumWeight << ", " << sumWt_kappa << ", " << population.size();
                throw TRACED_EXCEPTION(msg.str(), util::Error::SumWeight);
            }
            laggedKappa[lKMod] = sumWt_kappa / sumWeight;
        }

        size_t tmod = sim::moduloYearSteps(sim::ts0());

        // Calculate time-weighted average of kappa
        _sumAnnualKappa += laggedKappa[lKMod] * initialisationEIR[tmod];
        if (tmod == sim::stepsPerYear() - 1)
        {
            _annualAverageKappa = _sumAnnualKappa / annualEIR; // inf or NaN when annualEIR is 0
            _sumAnnualKappa = 0.0;
        }

        return laggedKappa[lKMod]; // kappa now
    }

    virtual void surveyEIR()
    {
        tsAdultEIR = tsAdultEntoInocs / tsNumAdults;
        tsAdultEntoInocs = 0.0;

        tsAdultEIR_i = tsAdultEntoInocs_i / tsNumAdults;
        tsAdultEntoInocs_i = 0.0;

        tsAdultEIR_l = tsAdultEntoInocs_l / tsNumAdults;
        tsAdultEntoInocs_l = 0.0;

        tsNumAdults = 0;

        surveyInputEIR += initialisationEIR[sim::moduloYearSteps(sim::ts0())];

        surveySimulatedEIR += tsAdultEIR;
        surveySimulatedEIR_i += tsAdultEIR_i;
        surveySimulatedEIR_l += tsAdultEIR_l;
    }

protected:
    /** Calculates the EIR individuals are exposed to.
     *
     * Call once per time-step: updates ITNs in vector model.
     *
     * @param host Transmission data for the human to calculate EIR for.
     * @param ageGroupData Age group of this host for availablility data.
     * @param EIR Out-vector. Set to the age- and heterogeneity-specific EIR an
     *    individual human is exposed to, per parasite genotype, in units of
     *    inoculations per day. Length set by callee. 
     *    _i for imported infections and _l for local infections */
    virtual void calculateEIR(Host::Human &human, double ageYears, vector<double> &EIR_i, vector<double> &EIR_l) const = 0;

    virtual void checkpoint(istream &stream)
    {
        simulationMode &stream;
        interventionMode &stream;
        initialisationEIR &stream;
        laggedKappa &stream;
        annualEIR &stream;
        _annualAverageKappa &stream;
        _sumAnnualKappa &stream;
        tsAdultEIR &stream;
        surveyInputEIR &stream;
        surveySimulatedEIR &stream;
        surveySimulatedEIR_i &stream;
        surveySimulatedEIR_l &stream;
        lastSurveyTime &stream;
        adultAge &stream;
        numTransmittingHumans &stream;
    }

    virtual void checkpoint(ostream &stream)
    {
        simulationMode &stream;
        interventionMode &stream;
        initialisationEIR &stream;
        laggedKappa &stream;
        annualEIR &stream;
        _annualAverageKappa &stream;
        _sumAnnualKappa &stream;
        tsAdultEIR &stream;
        surveyInputEIR &stream;
        surveySimulatedEIR_i &stream;
        surveySimulatedEIR_l &stream;
        surveySimulatedEIR &stream;
        lastSurveyTime &stream;
        adultAge &stream;
        numTransmittingHumans &stream;
    }

private:
    // The times here should be for the last updated index of arrays:
    void ctsCbInputEIR(ostream &stream)
    {
        int prevStep = (sim::now() - sim::oneTS()) / sim::oneTS();
        // Note: prevStep may be negative, hence util::mod not mod_nn:
        stream << '\t' << initialisationEIR[util::mod(prevStep, sim::stepsPerYear())];
    }
    void ctsCbSimulatedEIR(ostream &stream) { stream << '\t' << tsAdultEIR; }
    void ctsCbKappa(ostream &stream)
    {
        // The latest time-step's kappa:
        stream << '\t' << laggedKappa[sim::moduloSteps(sim::now(), laggedKappa.size())];
    }
    void ctsCbNumTransmittingHumans(ostream &stream) { stream << '\t' << numTransmittingHumans; }

public:
    /** Entomological inoculation rate for adults during the
     * pre-intervention phase.
     *
     * Length: time-steps per year
     *
     * Index sim::now_mod_steps_per_year() corresponds to the EIR
     * acting on the current time-step: i.e. total inoculations since the
     * previous time-step.
     * Since time-step 0 is not calculated, initialisationEIR[0] is actually the
     * last value used (to calculate the state at the start of the second year).
     *
     * Units: infectious bites per adult per time step
     *
     * Not checkpointed; doesn't need to be except when a changeEIR intervention
     * occurs. */
    vector<double> initialisationEIR;

    /** The type of EIR calculation. Checkpointed. */
    int simulationMode;
    /** New simulation mode during intervention period. Not checkpointed. */
    int interventionMode;

    /** The probability of infection of a mosquito at each bite.
     * It is calculated as the average infectiousness per human.
     *
     * The value in index sim::ts1().moduloSteps(initialKappa.size()) is the
     * kappa from this time step (i.e. the infectiousness of humans at the end of
     * this step). Length depends on entomological incubation period from
     * non-vector model.
     *
     * Checkpointed. */
    vector<double> laggedKappa;

    /** Total annual infectious bites per adult.
     *
     * Checkpointed. */
    double annualEIR;

private:
    /*! annAvgKappa is the overall proportion of mosquitoes that get infected
     * allowing for the different densities in different seasons (approximating
     * relative mosquito density with the EIR).
     *
     * Checkpointed. */
    double _annualAverageKappa;

    /*! Used to calculate annAvgKappa.
     *
     * Checkpointed. */
    double _sumAnnualKappa;

    /// Adult-only EIR over the last update
    double tsAdultEIR, tsAdultEIR_i, tsAdultEIR_l;

    /** Per-time-step input EIR summed over inter-survey period.
     * Units: infectious bites/adult/inter-survey period. */
    double surveyInputEIR, surveyInputEIR_i, surveyInputEIR_l;

    /** Per-time-step simulated EIR summed over inter-survey period.
     * Units: infectious bites/adult/inter-survey period. */
    double surveySimulatedEIR, surveySimulatedEIR_i, surveySimulatedEIR_l;

    /** Time of last survey. */
    SimTime lastSurveyTime = sim::never();

    /// age at which an individual is considered an adult
    SimTime adultAge = sim::never();

    /// For "num transmitting humans" cts output.
    int numTransmittingHumans;

    // Reporting data. Doesn't need checkpointing due to reset every time-step.
    // accumulator for time step EIR of adults
    double tsAdultEntoInocs = 0.0, tsAdultEntoInocs_i = 0.0, tsAdultEntoInocs_l = 0.0;
    int tsNumAdults = 0; // accumulator for time step adults requesting EIR

    bool opt_vaccine_genotype = false;
};

} // namespace Transmission
} // namespace OM
#endif
