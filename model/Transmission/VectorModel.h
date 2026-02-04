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

#ifndef Hmod_VectorModel
#define Hmod_VectorModel

#include "Global.h"
#include "Transmission/TransmissionModel.h"
#include "Transmission/Anopheles/AnophelesModel.h"
#include "Transmission/Anopheles/AnophelesModelFitter.h"
#include "Host/Human.h"

namespace scnXml
{
class Vector;
}

namespace OM
{
class Population;
namespace Transmission
{

/** Transmission models, Chitnis et al.
 *
 * This class contains code for species-independent components. Per-species
 * code is in the Anopheles directory and namespace. */
class VectorModel : public TransmissionModel
{
public:
    /// Get the map of species names to indicies.
    const map<string, size_t> &getSpeciesIndexMap() { return speciesIndex; }

    VectorModel(vector<double> initEIR, int interventionMode, vector<std::unique_ptr<Anopheles::AnophelesModel>> speciesList,
                vector<std::unique_ptr<Anopheles::AnophelesModelFitter>> speciesFittersList, map<string, size_t> speciesIndexList,
                int populationSize);

    virtual ~VectorModel();

    /** Extra initialisation when not loading from a checkpoint, requiring
     * information from the human population structure. */
    virtual void init2(const vector<Host::Human> &population);

    virtual void initVectorInterv(const scnXml::Description::AnophelesSequence &list, size_t instance, const string &name);
    virtual void initVectorTrap(const scnXml::VectorTrap::DescriptionSequence list, size_t instance,
                                const scnXml::VectorTrap::NameOptional name);

    virtual void scaleEIR(double factor);

    virtual SimTime initIterate();

    virtual void vectorUpdate(const vector<Host::Human> &population);

    virtual void calculateEIR(Host::Human &human, double ageYears, vector<double> &EIR_i, vector<double> &EIR_l) const;

    virtual void deployVectorPopInterv(size_t instance);
    virtual void deployVectorTrap(size_t instance, double popSize, SimTime lifespan);

    virtual void changeEIRIntervention(const scnXml::NonVector &nonVectorData);

    virtual void uninfectVectors();

    virtual void summarize();

protected:
    virtual void checkpoint(istream &stream);
    virtual void checkpoint(ostream &stream);

private:
    void ctsCbN_v0(ostream &stream);
    void ctsCbP_A(ostream &stream);
    void ctsCbP_Amu(ostream &stream);
    void ctsCbP_A1(ostream &stream);
    void ctsCbP_Ah(ostream &stream);
    void ctsCbP_df(ostream &stream);
    void ctsCbP_dif(ostream &stream);
    void ctsCbN_v(ostream &stream);
    void ctsCbO_v(ostream &stream);
    void ctsCbS_v(ostream &stream);
    void ctsCbAlpha(const vector<Host::Human> &population, ostream &stream);
    void ctsCbP_B(const vector<Host::Human> &population, ostream &stream);
    void ctsCbP_CD(const vector<Host::Human> &population, ostream &stream);
    void ctsNetInsecticideContent(const vector<Host::Human> &population, ostream &stream);
    void ctsIRSInsecticideContent(const vector<Host::Human> &population, ostream &stream);
    void ctsIRSEffects(const vector<Host::Human> &population, ostream &stream);
    void ctsCbResAvailability(ostream &stream);
    void ctsCbResRequirements(ostream &stream);

public:
    /// RNG used by the transmission model
    LocalRng m_rng;

    /// Number of iterations performed during initialization.
    ///
    /// Special cases: 0 during initial one-human-lifespan warmup, <0 after
    /// initialisation.
    int initIterations;

    /** Per anopheles species data.
     *
     * Array will be recreated by constructor, but some members of AnophelesModel
     * need to be checkpointed. */
    vector<std::unique_ptr<Anopheles::AnophelesModel>> species;
    vector<std::unique_ptr<Anopheles::AnophelesModelFitter>> speciesFitters;

    /** A map of anopheles species/variant name to an index in species.
     *
     * When the main ento data is read from XML, each anopheles section is read
     * into one index of the species array. The name and this index are added
     * here.
     *
     * Other data read from XML should look up the name here and use the index
     * found. Doesn't need checkpointing. */
    map<string, size_t> speciesIndex;

    bool opt_vaccine_genotype = false;
    
    friend class PerHost;
    friend class AnophelesModelSuite;
};

} // namespace Transmission
} // namespace OM
#endif
