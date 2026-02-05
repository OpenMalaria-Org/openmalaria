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

#ifndef Hmod_WithinHost_Interface
#define Hmod_WithinHost_Interface

#include "Global.h"
#include "util/random.h"
#include "Host/WithinHost/Diagnostic.h"
#include "Host/WithinHost/Pathogenesis/State.h"
#include "Host/WithinHost/Infection/Infection.h"

#include "Parameters.h"

using namespace std;

namespace scnXml{
    class Scenario;
    class Model;
    class TreatmentOption;
}
class UnittestUtil;

namespace OM {
namespace Host {
    class Human;
}
namespace WithinHost {

using util::LocalRng;

/**
 * Type used to select a treatment option.
 * 
 * Pass by value; it just hides an integer.
 * 
 * Note: this struct and the Treatments class offer a sub-set of the
 * functionality offered by CMDecisionTree, and thus is technically redundant.
 */
struct TreatmentId{
    inline bool operator==( const TreatmentId that ){ return id == that.id; }
    inline bool operator!=( const TreatmentId that ){ return id != that.id; }
    
    /// Default constructor: construct to an initial value. Don't pass this value to WHInterface::treatment()!
    TreatmentId() : id(numeric_limits<uint32_t>::max()) {}
private:
    explicit TreatmentId(uint32_t id) : id(id) {}
    uint32_t id;
    friend class Treatments;
};

/**
 * Interface to the within-host models. These models encapsulate the infections
 * and related immunity factors of a single human, starting with infection
 * (i.e. assuming successful innoculation), including some drug action code,
 * and outputting parasite densities.
 */
class WHInterface {
public:
    /// @brief Static methods
    //@{
    /// Initialise static parameters
    static void init(const OM::Parameters& parameters, const scnXml::Scenario& scenario);
    
    /** Configure a new treatment option, and return the code used to select
     * that option later. */
    static TreatmentId addTreatment( const scnXml::TreatmentOption& desc );

    /// Create an instance using the appropriate model
    static unique_ptr<WHInterface> createWithinHostModel( LocalRng& rng, double comorbidityFactor );
    //@}

    /// @brief Constructors, destructors and checkpointing functions
    //@{
    WHInterface(): numInfs(0) {}
    virtual ~WHInterface() = default;
    
    WHInterface(WHInterface&&) = default;
    WHInterface& operator=(WHInterface&&) = default;

    /// Checkpointing
    template<class S>
    void operator& (S& stream) {
        checkpoint (stream);
    }
    //@}

    /** Return the infectiousness of this human to biting mosquitoes.
     * This step is independent of parasite genetics.
     * 
     * @returns the probability of this human infecting a feeding mosquito.
     * 
     * Calculates the probability of transmitting an infection of a given
     * genotype to a mosquito and store the result in probTransGenotype[g]
     * for a given genotype g 
     * 
     * Calulates the probability for imported infections _i and
     * for local infections _l */
    virtual double probTransmissionToMosquito(vector<double> &probTransGenotype_i, vector<double> &probTransGenotype_l)const = 0;

    /// @returns true if host has patent parasites
    virtual bool summarize(Host::Human& human) const =0;

    /// Create a new infection within this human
    virtual void importInfection(LocalRng& rng) =0;

    /**
     * Carry out the effects of some treatment option, optionally with intervention deployment.
     * 
     * This is equivalent to calling treatSimple, then deploying any included interventions.
     * 
     * This may be used by intervention deployment, thus should use time 
     * sim::nowOrTs1(). */
    virtual void treatment( Host::Human& human, TreatmentId treatId ) =0;
    
    /// Conditionally gives Primaquine as a treatment. Reports as appropriate.
    virtual void optionalPqTreatment( Host::Human& human ) =0;
    
    /** Treat a patient via the simple treatment model. Return true if any
     * blood-stage treatment is administered. Report any liver-stage treatments. */
    virtual bool treatSimple( Host::Human& human, SimTime timeLiver, SimTime timeBlood ) =0;
    
    /** Give a patient a course of drugs, via the Pk/Pd model
     * 
     * Note: doses sizes are modified according to age via the dosage
     * table given at the time this function is called.
     *
     * @param schedule Index of a treatment schedule
     * @param dosages Index of a dosage table
     * @param age Age of human in years
     */
    virtual void treatPkPd(size_t schedule, size_t dosages, double age, double delay_d) =0;

    /** Add new infections and update the parasite densities of existing
     * infections. Also update immune status.
     *
     * @param nNewInfs Number of inoculations this time-step
     * @param genotype_weights For use in selecting infection genotypes. See
     *  documentation of Genotypes::sampleGenotype().
     * @param ageInYears Age of human
     * @param bsvFactor Parasite survival factor for blood-stage vaccines
     */
    virtual void update(Host::Human &human, LocalRng& rng, int &nNewInfs_i, int &nNewInfs_l, 
        vector<double>& genotype_weights_i, vector<double>& genotype_weights_l, double ageInYears) =0;

    /** TODO: this should not need to be exposed. It is currently used by a
     * severe outcome (pDeath) model inside the EventScheduler "case
     * management" model, and case management diagnostics. */
    virtual double getTotalDensity() const =0;
    
    /** Simulate use of a diagnostic test.
     *
     * Does not report for costing purposes.
     * 
     * Is used both during time step updates and during monitoring.
     * 
     * @returns true when the diagnostic is positive
     */
    virtual bool diagnosticResult( LocalRng& rng, const Diagnostic& diagnostic ) const =0;
    
    /** Use the pathogenesis model to determine, based on infection status
     * and random draw, this person't morbidity.
     * 
     * This function is called after update() every time step.
     * 
     * @param human A reference to the human, used when reporting
     * @param ageYears Age of human host in years
     * @param isDoomed True if the human is already doomed to die, used for reporting
     */
    virtual Pathogenesis::StatePair determineMorbidity( Host::Human& human, double ageYears, bool isDoomed ) =0;

    /// Special intervention: clears all immunity
    virtual void clearImmunity() =0;
    
    // TODO(monitoring): these shouldn't have to be exposed (perhaps use summarize to report the data):
    virtual double getCumulative_h() const =0;
    virtual double getCumulative_Y() const =0;

    virtual InfectionOrigin getInfectionOrigin() const =0;

    /** The maximum number of infections a human can have. The only real reason
     * for this limit is to prevent incase bad input from causing the number of
     * infections to baloon stupidly.
     *
     * Exact constraint is: _MOI <= MAX_INFECTIONS. */
    static const int MAX_INFECTIONS = 21;

protected:    
    virtual void checkpoint (istream& stream);
    virtual void checkpoint (ostream& stream);

    /// Multiplicity of infection
    int numInfs;

    friend class ::UnittestUtil;
};

/// True if any by-genotype reporting is enabled (read-only)
extern bool reportInfectionsByGenotype;

}
}
#endif
