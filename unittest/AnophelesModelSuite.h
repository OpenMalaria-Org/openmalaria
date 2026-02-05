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

#ifndef Hmod_AnophelesModelSuite
#define Hmod_AnophelesModelSuite

#include "Global.h"
#include "configured/TestPaths.h"	// from config; but must be included from the build dir
#include <cxxtest/TestSuite.h>
#include "ExtraAsserts.h"
#include <yaml-cpp/yaml.h>
#include <fstream>
#include <algorithm>
#include <list>
#include <xsd/cxx/tree/exceptions.hxx>

#include "Host/Human.h"
#include "Transmission/VectorModel.h"
#include "Transmission/Anopheles/AnophelesModel.h"
#include "util/vectors.h"


/** Read an array into a vector. */
template<class T>
void operator>> (const YAML::Node& node, vector<T>& vec) {
  vec.resize (node.size());
  size_t i = 0;
  for (YAML::Iterator it = node.begin(); it != node.end(); ++it, ++i) {
    vec[i] = it->Read<double>();
  }
}
/** As above, but with function syntax. */
template<class T>
vector<T> yaml2Std (const YAML::Node& node) {
  vector<T> ret;
  node >> ret;
  return ret;
}

/** Read a YAML node into a WeibullDecayedValue.
 *
 * Should be an operator= method of WeibullDecayedValue, but that would make
 * the main simulator depend on yaml-cpp. */
WeibullDecayedValue yaml2WeibullDecayedValue (const YAML::Node& node) {
  double initial, hl, k = 1.0;
  node["initial"] >> initial;
  node["Halflife"] >> hl;
  const YAML::Node *Wbk = node.FindValue ("Weibullk");
  if (Wbk != NULL)
    *Wbk >> k;
  WeibullDecayedValue ret;
  ret.setParameters (initial, hl, k);
  return ret;
}


// The assert we use in calculateEIR unittests. It tests vtm and species params against node.
#   define TS_ASSERT_SPECIES_APPROX(node) AnophelesModelSuite::doAssertSpecies(__FILE__,__LINE__, (node))


/** Unit test for AnophelesModel code.
 *
 * Possible tests to write:
 * 
AnophelesModel unit testing:
  initialise:
    read parameters from XML
    calculate EIR array from fourier series, rotate and add to return param
    read emergence params from XML and scale
  initFeedingCycleProbs:
    calculates P_A, P_df, P_dif
    return total pop avail
  initMainSimulation:
    runs VectorEmergence code
    validates emerge params
    prints out new scenario file
 */
class AnophelesModelSuite : public CxxTest::TestSuite
{
public:
  AnophelesModelSuite() {
    Global::clResourcePath = UnittestSourceDir;
    ifstream file(Global::lookupResource ("AnophelesModelSuite.yaml").c_str());
    
    YAML::Parser parser(file);
    parser.GetNextDocument(doc);
    
    file.close();
  }
  
  /** Initialise species.
   *
   * Rather than directly initialise the elements we want to use, try to set up
   * the whole simulation. I thought it'd be easier...
   *
   * We load the emergence rates in from the scenario. They are currently wrong
   * (although still usable for alteration-notifying unittests).
   * Better if we used the cut-down model of VectorEmergenceSuite to generate.
   *
   * Note that using doc["key"] to access map elements isn't speed-optimal when
   * many key-value pairs exist, but for a small number it should be fine. */
  void setUp () {
    try {
    // "unused" node, but it still checks element exists (confirming correct data file)
    /*const YAML::Node& node =*/ doc["AnophelesModelSuite"];
    
    createDocument (UnittestScenario);
    Global::initGlobal();
    simulation = new Simulation();
    
    // Normally done by Simulation::start():
    simulation->simulationTime = 0;
    simulation->_population->estimateRemovalRates();
    simulation->_population->setupPyramid(false);
    simulation->simulationTime = simulationTime = 1;
    
    // No eir data is present in scenario.xml, so it should directly initialise all model parameters
    vtm = dynamic_cast<VectorModel*> (simulation->_population->_transmissionModel);
    assert (vtm != NULL);
    assert (vtm->numSpecies == 1);
    species = &vtm->species[0];
    population = &simulation->_population->_population;
    
    } catch (const ::xsd::cxx::tree::exception<char>& e) {
      cerr << "\nXSD Exception: " << e.what() << '\n' << e << endl;
      throw;
    } catch (const exception& e) {
      cerr << "\nException: " << e.what() << endl;
      throw;
    } catch (...) {
      cerr << "\nUnknown exception" << endl;
      throw;
    }
  }
  
  void tearDown () {
    delete simulation;
  }
  
  void testCalcInverseDFTExp () {
    const YAML::Node& node = doc["calcInverseDFTExp"];
    vector<double> out(10), fc, expected;
    node["fourierCoefficients"] >> fc;
    node["output"] >> expected;
    
    species->calcInverseDFTExp (out, fc);
    TS_ASSERT_VECTOR_APPROX (out, expected);
  }
  
  void testCalculateEIR () {
    // setUp() should have done all the hard work initialising stuff
    // current limitation: humans have no infectiousness
    //NOTE: could we load the model from a checkpoint at end of initialisation phase?
    
    // Run advancePeriod and getEIR for each human. Don't care about eir returned by getEIR;
    // just test vtm->eirPerDayOfYear and species parameters.
    vtm->timeStepNumEntoInnocs = 0;
    vtm->advancePeriod (*population, simulationTime);
    for (list<Host::Human>::iterator it = population->begin(); it != population->end(); ++it)
      vtm->getEIR (simulationTime, it->perHostTransmission, it->getAgeInYears());
    TS_ASSERT_SPECIES_APPROX (doc["calculateEIR"]["output"]);
  }
  
  void testCalculateEirDeterrency () {
    const YAML::Node& node = doc["calculateEirDeterrency"];
    species->humanBase.ITNDeterrency = yaml2WeibullDecayedValue(node["Deterrency"]);
    for (list<Host::Human>::iterator it = population->begin(); it != population->end(); ++it)
      it->setupITN();
    vtm->timeStepNumEntoInnocs = 0;
    vtm->advancePeriod (*population, simulationTime);
    for (list<Host::Human>::iterator it = population->begin(); it != population->end(); ++it)
      vtm->getEIR (simulationTime, it->perHostTransmission, it->getAgeInYears());
    TS_ASSERT_SPECIES_APPROX (node["output"]);
  }
  
  void testCalculateEirPreprandialKilling () {
    const YAML::Node& node = doc["calculateEirPreprandialKilling"];
    species->humanBase.ITNPreprandialKillingEffect = yaml2WeibullDecayedValue(node["PreprandialKilling"]);
    for (list<Host::Human>::iterator it = population->begin(); it != population->end(); ++it)
      it->setupITN();
    vtm->timeStepNumEntoInnocs = 0;
    vtm->advancePeriod (*population, simulationTime);
    for (list<Host::Human>::iterator it = population->begin(); it != population->end(); ++it)
      vtm->getEIR (simulationTime, it->perHostTransmission, it->getAgeInYears());
    TS_ASSERT_SPECIES_APPROX (node["output"]);
  }
  
  void testCalculateEirPostprandialKilling () {
    const YAML::Node& node = doc["calculateEirPostprandialKilling"];
    species->humanBase.ITNPostprandialKillingEffect = yaml2WeibullDecayedValue(node["PostprandialKilling"]);
    for (list<Host::Human>::iterator it = population->begin(); it != population->end(); ++it)
      it->setupITN();
    vtm->timeStepNumEntoInnocs = 0;
    vtm->advancePeriod (*population, simulationTime);
    for (list<Host::Human>::iterator it = population->begin(); it != population->end(); ++it)
      vtm->getEIR (simulationTime, it->perHostTransmission, it->getAgeInYears());
    TS_ASSERT_SPECIES_APPROX (node["output"]);
  }
  
  void testCalculateEirRestKilling () {
    const YAML::Node& node = doc["calculateEirRestKilling"];
    species->humanBase.IRSKillingEffect = yaml2WeibullDecayedValue(node["RestKilling"]);
    for (list<Host::Human>::iterator it = population->begin(); it != population->end(); ++it)
      it->setupIRS();
    vtm->timeStepNumEntoInnocs = 0;
    vtm->advancePeriod (*population, simulationTime);
    for (list<Host::Human>::iterator it = population->begin(); it != population->end(); ++it)
      vtm->getEIR (simulationTime, it->perHostTransmission, it->getAgeInYears());
    TS_ASSERT_SPECIES_APPROX (node["output"]);
  }
  
  void testCalculateEirLarviciding () {
    const YAML::Node& node = doc["calculateEirLarviciding"];
    node["Larviciding"]["effectiveness"] >> species->larvicidingIneffectiveness;
    species->larvicidingIneffectiveness = 1 - species->larvicidingIneffectiveness;
    species->larvicidingEndStep = 1000;
    vtm->timeStepNumEntoInnocs = 0;
    vtm->advancePeriod (*population, simulationTime);
    for (list<Host::Human>::iterator it = population->begin(); it != population->end(); ++it)
      vtm->getEIR (simulationTime, it->perHostTransmission, it->getAgeInYears());
    TS_ASSERT_SPECIES_APPROX (node["output"]);
  }
  
  void doAssertSpecies (const char *f, unsigned l, const YAML::Node& node) {
    ETS_ASSERT_EQUALS ((int)vtm->timeStepNumEntoInnocs, simulation->_population->_populationSize);
    
    double resultEIR = 0.0;
    for (vector<double>::const_iterator it = vtm->timeStepEntoInnocs.begin(); it != vtm->timeStepEntoInnocs.end(); ++it)
      resultEIR += *it;
    
    double avEIR;
    node["averageEIR"] >> avEIR;
    ETS_ASSERT_APPROX (resultEIR / vtm->timeStepNumEntoInnocs, avEIR);
    _ETS_ASSERT_VECTOR_APPROX (f,l, species->P_A, yaml2Std<double> (node["P_A"]));
    _ETS_ASSERT_VECTOR_APPROX (f,l, species->P_df, yaml2Std<double> (node["P_df"]));
    _ETS_ASSERT_VECTOR_APPROX (f,l, species->P_dif, yaml2Std<double> (node["P_dif"]));
    _ETS_ASSERT_VECTOR_APPROX (f,l, species->N_v, yaml2Std<double> (node["N_v"]));
    _ETS_ASSERT_VECTOR_APPROX (f,l, species->O_v, yaml2Std<double> (node["O_v"]));
    _ETS_ASSERT_VECTOR_APPROX (f,l, species->S_v, yaml2Std<double> (node["S_v"]));
    
    if (CxxTest::tracker().testFailed()) {
      cerr << "Unittest failed; new output:\n";
      // Print values so they can easily be copied into expected output:
      cerr << setprecision(10);
      cerr << "averageEIR: " << resultEIR / vtm->timeStepNumEntoInnocs;
      cerr << "\nP_A:\t" << species->P_A;
      cerr << "\nP_df:\t" << species->P_df;
      cerr << "\nP_dif:\t" << species->P_dif;
      cerr << "\nN_v:\t" << species->N_v;
      cerr << "\nO_v:\t" << species->O_v;
      cerr << "\nS_v:\t" << species->S_v << endl;;
    }
  }
  
private:
  YAML::Node doc;
  int simulationTime;
  Simulation *simulation;
  VectorModel *vtm;
  AnophelesModel *species;
  list<Host::Human> *population;
};

#endif
