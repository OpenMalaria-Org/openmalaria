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

#include "util/errors.h"
#include "util/CommandLine.h"
#include "util/ModelOptions.h"
#include "util/ModelNameProvider.h"
#include "util/StreamValidator.h"
#include "util/DocumentLoader.h"
#include "util/XMLChecker.h"

#include "mon/Continuous.h"
#include "mon/management.h"

#include "interventions/InterventionManager.h"
#include "Clinical/ClinicalModel.h"

#include "Host/NeonatalMortality.h"
#include "checkpoint.h"

#include "schema/scenario.h"

#include <cerrno>

namespace OM {
    using mon::Continuous;
    using interventions::InterventionManager;
    using Transmission::TransmissionModel;
}

using namespace OM;

void print_progress(int lastPercent, SimTime &estEndTime)
{
    int percent = (sim::now() * 100) / estEndTime;
    if( percent != lastPercent ) {   // avoid huge amounts of output for performance/log-file size reasons
        lastPercent = percent;
        cerr << "\r" << percent << "%\t" << flush;
    }
}

void print_errno()
{
    if( errno != 0 )
    {
       char err[256];
       sprintf(err, "t = %d Please report! Error: ", int(sim::now()));
       std::perror(err);
       errno = 0;
    }
}

// Internal simulation loop
void run(Population &population, TransmissionModel &transmission, SimTime humanWarmupLength, SimTime &endTime, SimTime &estEndTime, bool surveyOnlyNewEp, string phase)
{
    static int lastPercent = -1;

    if (util::CommandLine::option(util::CommandLine::VERBOSE)) cout << "Starting " << phase << "..." << endl;

    while (sim::now() < endTime)
    {
        if (util::CommandLine::option(util::CommandLine::VERBOSE) && sim::intervDate() > 0)
            cout << "Time step: " << sim::now() / sim::oneTS() << ", internal days: " << sim::now() << " | " << estEndTime << ", Intervention Date: " << sim::intervDate() << endl;

        // Monitoring. sim::now() gives time of end of last step,
        // and is when reporting happens in our time-series.
        Continuous.update( population );
        if( sim::intervDate() == mon::nextSurveyDate() ){
            for(Host::Human &human : population.humans)
                Host::summarize(human, surveyOnlyNewEp);
            transmission.summarize();
            mon::concludeSurvey();
        }
        
        // Deploy interventions, at time sim::now().
        InterventionManager::deploy( population.humans, transmission );
        
        // Time step updates. Time steps are mid-day to mid-day.
        // sim::ts0() gives the date at the start of the step, sim::ts1() the date at the end.
        sim::start_update();

        // This should be called before humans contract new infections in the simulation step.
        // This needs the whole population (it is an approximation before all humans are updated).
        transmission.vectorUpdate(population.humans);
        
        // NOTE: no neonatal mortalities will occur in the first 20 years of warmup
        // (until humans old enough to be pregnate get updated and can be infected).
        Host::NeonatalMortality::update (population.humans);
        
        for (Host::Human& human : population.humans)
        {
            if (human.getDOB() + sim::maxHumanAge() >= humanWarmupLength) // this is last time of possible update
                Host::update(human, transmission);
        }
       
        population.update();
        
        // Doesn't matter whether non-updated humans are included (value isn't used
        // before all humans are updated).
        transmission.updateKappa(population.humans);
        transmission.surveyEIR();

        sim::end_update();

        if (util::CommandLine::option(util::CommandLine::PROGRESS))
            print_progress(lastPercent, estEndTime);
        print_errno();
    }

    if (util::CommandLine::option(util::CommandLine::VERBOSE)) cout << "Finishing " << phase << "..." << endl;
}

/// main() — loads scenario XML and runs simulation
int main(int argc, char* argv[])
{
    int exitStatus = EXIT_SUCCESS;
    bool startedFromCheckpoint;

    string scenarioFile;
    string checkpointFileName;
    SimTime estEndTime, endTime;
    
    try {
        util::set_gsl_handler();
        
        scenarioFile = util::CommandLine::parse (argc, argv);
        unique_ptr<scnXml::Scenario> scenario = util::loadScenario(scenarioFile);

        util::XMLChecker().PerformPostValidationChecks(*scenario);

        // 1) elements with no dependencies on other elements initialised here:
        WithinHost::Genotypes::init( *scenario );
        util::master_RNG.seed( scenario->getModel().getComputationParameters().getIseed(), 0 ); // Init RNG with Iseed
        util::ModelNameProvider modelNameProvider(scenario->getModel());

        // 2) elements depending on only elements initialised in (1):
        sim::init(*scenario, modelNameProvider); // Also reads survey dates.

        // 3) elements depending on only elements initialised in (2).
        Parameters parameters( scenario->getModel().getParameters(), modelNameProvider ); // Depends on ModelNameProvider.
        util::ModelOptions::init( scenario->getModel().getModelOptions(), modelNameProvider ); // Depends on ModelNameProvider.

        // 4) elements depending on only elements initialised in (3).
        WithinHost::diagnostics::init( parameters, *scenario ); // Depends on Parameters.

        // 5) elements depending on only elements initialised in (4).
        mon::initReporting( *scenario ); // Reporting init depends on diagnostics and monitoring

        // Init models used by humans
        Transmission::PerHost::init( scenario->getModel().getHuman().getAvailabilityToMosquitoes() );
        Host::InfectionIncidenceModel::init( parameters );
        WithinHost::WHInterface::init( parameters, *scenario );
        Clinical::ClinicalModel::init( parameters, *scenario );
        Host::NeonatalMortality::init( scenario->getModel().getClinical() );
        AgeStructure::init( scenario->getDemography() );

        // 3) elements depending on other elements; dependencies on (1) are not mentioned:
        // Transmission model initialisation depends on Transmission::PerHost and
        // genotypes (both from Human, from Population::init()) and
        // mon::AgeGroup (from Surveys.init()):
        // Note: PerHost dependency can be postponed; it is only used to set adultAge
        size_t popSize = scenario->getDemography().getPopSize();

        std::unique_ptr<Population> population = std::unique_ptr<Population>(new Population(popSize));
        std::unique_ptr<TransmissionModel> transmission = std::unique_ptr<TransmissionModel>(Transmission::createTransmissionModel(scenario->getEntomology(), popSize));

        registerContinousPopulationCallbacks();

        // Depends on transmission model (for species indexes):
        // MDA1D may depend on health system (too complex to verify)
        interventions::InterventionManager::init(scenario->getInterventions(), *population, *transmission );
        Clinical::ClinicalModel::setHS( scenario->getHealthSystem() ); // Depends on interventions, PK/PD (from humanPop)
        mon::initCohorts( scenario->getMonitoring() ); // Depends on interventions

        bool surveyOnlyNewEp = scenario->getMonitoring().getSurveyOptions().getOnlyNewEpisode();

        sim::s_t0 = sim::zero();
        sim::s_t1 = sim::zero();
        
        // Make sure warmup period is at least as long as a human lifespan, as the
        // length required by vector warmup, and is a whole number of years.
        SimTime humanWarmupLength = sim::maxHumanAge();
        if(transmission->interventionMode != Transmission::forcedEIR)
            humanWarmupLength = max(humanWarmupLength, sim::fromYearsI(55)); // Data is summed over 5 years; add an extra 50 for stabilization.

        humanWarmupLength = sim::fromYearsI( static_cast<int>(ceil(sim::inYears(humanWarmupLength))) );

        // ———  End of static data initialisation  ———
        checkpointFileName = util::CommandLine::getCheckpointName();

        if(checkpointFileName == "")
            checkpointFileName = "checkpoint";

        if(util::CommandLine::option(util::CommandLine::CHECKPOINT))
        {
            ifstream checkpointFile(checkpointFileName,ios::in);
            startedFromCheckpoint = checkpointFile.is_open();
            if(startedFromCheckpoint == false)
                errno = 0; // Cleanup errno if file doesn't exist
        }
        else
            startedFromCheckpoint = false;
        
        estEndTime = humanWarmupLength + (sim::endDate() - sim::startDate()) + sim::oneTS();
        assert( estEndTime + sim::never() < sim::zero() );

        if (startedFromCheckpoint)
        {
            Continuous.init(scenario->getMonitoring(), true);
            readCheckpoint(checkpointFileName, endTime, estEndTime, *population, *transmission);

            /** Calculate ento availability percentiles **/
            Transmission::PerHostAnophParams::calcAvailabilityPercentiles();
        }
        else
        {
            Continuous.init(scenario->getMonitoring(), false);
            population->createInitialHumans();
            transmission->init2(population->humans);
            
            /** Calculate ento availability percentiles **/
            Transmission::PerHostAnophParams::calcAvailabilityPercentiles();

            /** Warm-up phase: 
             * Run the simulation using the equilibrium inoculation rates over one
             * complete lifespan (sim::maxHumanAge()) to reach immunological
             * equilibrium in all age classes. Don't report any events. */
            endTime = humanWarmupLength;
            run(*population, *transmission, humanWarmupLength, endTime, estEndTime, surveyOnlyNewEp, "Warmup");

            /** Transmission init phase:
             * Fit the emergence rate to the input EIR */
            SimTime iterate = transmission->initIterate();
            while(iterate > sim::zero())
            {
                endTime = endTime + iterate;
                // adjust estimation of final time step: end of current period + length of main phase
                estEndTime = endTime + (sim::endDate() - sim::startDate()) + sim::oneTS();
                run(*population, *transmission, humanWarmupLength, endTime, estEndTime, surveyOnlyNewEp, "EIR Calibration");
                iterate = transmission->initIterate();
            }

            /** Main phase:
             * This procedure starts with the current state of the simulation 
             * It continues updating assuming:
             * (i)         the default (exponential) demographic model
             * (ii)        the entomological input defined by the EIRs in intEIR()
             * (iii)       the intervention packages defined in Intervention()
             * (iv)        the survey times defined in Survey() */
            // reset endTime and estEndTime to their exact value after initIterate()
            estEndTime = endTime = endTime + (sim::endDate() - sim::startDate()) + sim::oneTS();
            sim::s_interv = sim::zero();
            Host::InfectionIncidenceModel::preMainSimInit();
            Clinical::InfantMortality::preMainSimInit();
            WithinHost::Genotypes::preMainSimInit();
            population->resetRecentBirths();
            transmission->summarize(); // Only to reset TransmissionModel::inoculationsPerAgeGroup
            mon::initMainSim();

            if(util::CommandLine::option (util::CommandLine::CHECKPOINT))
            {
                writeCheckpoint(startedFromCheckpoint, checkpointFileName, endTime, estEndTime, *population, *transmission);
                if( util::CommandLine::option (util::CommandLine::CHECKPOINT_STOP) )
                    throw util::cmd_exception ("Checkpoint test: checkpoint written", util::Error::None);
            }
        }

        // Main phase loop
        run(*population, *transmission, humanWarmupLength, endTime, estEndTime, surveyOnlyNewEp, "Intervention period");
       
        cerr << '\r' << flush;  // clean last line of progress-output
        
        for(Host::Human &human : population->humans)
            human.clinicalModel->flushReports();

        mon::writeSurveyData();
        
    # ifdef OM_STREAM_VALIDATOR
        util::StreamValidator.saveStream();
    # endif
        
        // simulation's destructor runs
    } catch (const OM::util::cmd_exception& e) {
        if( e.getCode() == 0 ){
            // this is not an error, but exiting due to command line
            cerr << e.what() << "; exiting..." << endl;
        }else{
            cerr << "Command-line error: "<<e.what();
            exitStatus = e.getCode();
        }
    } catch (const ::xsd::cxx::tree::exception<char>& e) {
        cerr << "XSD error: " << e.what() << '\n' << e << endl;
        exitStatus = OM::util::Error::XSD;
    } catch (const OM::util::checkpoint_error& e) {
        cerr << "Checkpoint error: " << e.what() << endl;
        cerr << e << flush;
        exitStatus = e.getCode();
    } catch (const OM::util::traced_exception& e) {
        cerr << "Code error: " << e.what() << endl;
        cerr << e << flush;
        cerr << "This is likely an error in the C++ code. Please report!" << endl;
        exitStatus = e.getCode();
    } catch (const OM::util::xml_scenario_error& e) {
        cerr << "Error: " << e.what() << endl;
        cerr << "In: " << scenarioFile << endl;
        exitStatus = e.getCode();
    } catch (const OM::util::base_exception& e) {
        cerr << "Error: " << e.message() << endl;
        exitStatus = e.getCode();
    } catch (const exception& e) {
        cerr << "Error: " << e.what() << endl;
        exitStatus = EXIT_FAILURE;
    } catch (...) {
        cerr << "Unknown error" << endl;
        exitStatus = EXIT_FAILURE;
    }
    
    // If we get to here, we already know an error occurred.
    if( errno != 0 )
        std::perror( "OpenMalaria" );
    
    return exitStatus;
}
