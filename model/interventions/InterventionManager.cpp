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

#include "interventions/InterventionManager.h"
#include "Population.h"
#include "util/CommandLine.h"
#include "util/UnitParse.h"
#include "interventions/GVI.h"
#include "interventions/IRS.h"
#include "interventions/ITN.h"
#include "interventions/Vaccine.h"
#include "interventions/HumanInterventionComponents.h"
#include "interventions/Deployments.h"
#include "Host/WithinHost/Diagnostic.h"
#include "Transmission/VectorModel.h"

namespace OM
{
namespace interventions
{

// ———  InterventionManager  ———

// static memory:

std::map<std::string, ComponentId> InterventionManager::identifierMap;
vector<unique_ptr<HumanInterventionComponent>> InterventionManager::humanComponents;
vector<ContinuousHumanDeployment> InterventionManager::continuous;
vector<unique_ptr<TimedDeployment>> InterventionManager::timed;
uint32_t InterventionManager::nextTimed;
OM::Host::ImportedInfections InterventionManager::importedInfections;

// declared in HumanComponents.h:
vector<ComponentId> removeAtIds[SubPopRemove::NUM];

// static functions:

void InterventionManager::init(const scnXml::Interventions &intervElt, const Population &population,
                               Transmission::TransmissionModel &transmission)
{
    nextTimed = 0;

    if (intervElt.getChangeHS().present())
    {
        const scnXml::ChangeHS &chs = intervElt.getChangeHS().get();
        if (chs.getTimedDeployment().size() > 0)
        {
            // timed deployments:
            for (auto it = chs.getTimedDeployment().begin(); it != chs.getTimedDeployment().end(); ++it)
            {
                try
                {
                    SimTime date = UnitParse::readDate(it->getTime(), UnitParse::STEPS /*STEPS is only for backwards compatibility*/);
                    timed.push_back(unique_ptr<TimedDeployment>(new TimedChangeHSDeployment(date, *it)));
                }
                catch (const util::format_error &e)
                {
                    throw util::xml_scenario_error(string("interventions/changeHS/timedDeployment/time: ").append(e.message()));
                }
            }
        }
    }
    if (intervElt.getChangeEIR().present())
    {
        const scnXml::ChangeEIR &eir = intervElt.getChangeEIR().get();
        if (eir.getTimedDeployment().size() > 0)
        {
            // timed deployments:
            for (auto it = eir.getTimedDeployment().begin(); it != eir.getTimedDeployment().end(); ++it)
            {
                try
                {
                    SimTime date = UnitParse::readDate(it->getTime(), UnitParse::STEPS /*STEPS is only for backwards compatibility*/);
                    timed.push_back(unique_ptr<TimedDeployment>(new TimedChangeEIRDeployment(date, *it)));
                }
                catch (const util::format_error &e)
                {
                    throw util::xml_scenario_error(string("interventions/changeEIR/timedDeployment/time: ").append(e.message()));
                }
            }
        }
    }
    // species_index_map is not available with the non-vector model or
    // non-dynamic mode, so setting it (lazily) also checks sim mode:
    map<string, size_t> species_index_map;

    {
        Transmission::VectorModel *vectorModel = dynamic_cast<Transmission::VectorModel *>(&transmission);
        if (vectorModel != nullptr) species_index_map = vectorModel->getSpeciesIndexMap();
    }

    if (intervElt.getHuman().present())
    {
        const scnXml::HumanInterventions &human = intervElt.getHuman().get();

        // 1. Read components
        for (auto it = human.getComponent().begin(), end = human.getComponent().end(); it != end; ++it)
        {
            const scnXml::HumanInterventionComponent &component = *it;
            if (identifierMap.count(component.getId()) > 0)
            {
                ostringstream msg;
                msg << "The id attribute of interventions.human.component elements must be unique; found \"" << component.getId()
                    << "\" twice.";
                throw util::xml_scenario_error(msg.str());
            }
            ComponentId id(humanComponents.size()); // i.e. index of next item
            identifierMap.insert(make_pair(component.getId(), id));

            SimTime expireAfter = sim::future();
            if (component.getSubPopRemoval().present())
            {
                const scnXml::SubPopRemoval &removeOpts = component.getSubPopRemoval().get();
                if (removeOpts.getOnFirstBout()) { removeAtIds[SubPopRemove::ON_FIRST_BOUT].push_back(id); }
                if (removeOpts.getOnFirstInfection()) { removeAtIds[SubPopRemove::ON_FIRST_INFECTION].push_back(id); }
                if (removeOpts.getOnFirstTreatment()) { removeAtIds[SubPopRemove::ON_FIRST_TREATMENT].push_back(id); }
                if (removeOpts.getAfterYears().present()) { expireAfter = sim::fromYearsN(removeOpts.getAfterYears().get()); }
            }

            HumanInterventionComponent *hiComponent;
            if (component.getScreen().present()) { hiComponent = new ScreenComponent(id, component.getScreen().get()); }
            else if (component.getTreatSimple().present())
            {
                hiComponent = new TreatSimpleComponent(id, component.getTreatSimple().get());
            }
            else if (component.getTreatPKPD().present())
            {
                hiComponent = new TreatPKPDComponent(id, component.getTreatPKPD().get());
            }
            else if (component.getDecisionTree().present())
            {
                hiComponent = new DecisionTreeComponent(id, component.getDecisionTree().get());
            }
            else if (component.getPEV().present())
            {
                hiComponent = new VaccineComponent(id, component.getPEV().get(), Vaccine::PEV);
            }
            else if (component.getBSV().present())
            {
                hiComponent = new VaccineComponent(id, component.getBSV().get(), Vaccine::BSV);
            }
            else if (component.getTBV().present())
            {
                hiComponent = new VaccineComponent(id, component.getTBV().get(), Vaccine::TBV);
            }
            else if (component.getITN().present())
            {
                hiComponent = new ITNComponent(id, component.getITN().get(), species_index_map);
            }
            else if (component.getIRS().present())
            {
                hiComponent = new IRSComponent(id, component.getIRS().get(), species_index_map);
            }
            else if (component.getGVI().present())
            {
                hiComponent = new GVIComponent(id, component.getGVI().get(), species_index_map);
            }
            else if (component.getRecruitmentOnly().present())
            {
                hiComponent = new RecruitmentOnlyComponent(id);
            }
            else if (component.getClearImmunity().present())
            {
                hiComponent = new ClearImmunityComponent(id);
            }
            else
            {
                throw util::xml_scenario_error("expected intervention.human.component element to have a "
                                               "child, didn't find it (perhaps I need updating)");
            }
            hiComponent->setExpireAfter(expireAfter);
            humanComponents.push_back(unique_ptr<HumanInterventionComponent>(hiComponent));
        }

        // 2. Read the list of deployments
        for (auto it = human.getDeployment().begin(), end = human.getDeployment().end(); it != end; ++it)
        {
            const scnXml::Deployment &elt = *it;
            // 2.a intervention components
            auto intervention = shared_ptr<HumanIntervention>(new HumanIntervention(elt.getComponent(), elt.getCondition()));

            // 2.b intervention deployments
            for (auto ctsIt = elt.getContinuous().begin(); ctsIt != elt.getContinuous().end(); ++ctsIt)
            {
                ComponentId subPop = ComponentId::wholePop();
                bool complement = false;
                if (ctsIt->getRestrictToSubPop().present())
                {
                    const string &subPopStr = ctsIt->getRestrictToSubPop().get().getId();
                    subPop = getComponentId(subPopStr);
                    complement = ctsIt->getRestrictToSubPop().get().getComplement();
                }
                const scnXml::ContinuousList::DeploySequence &ctsSeq = ctsIt->getDeploy();
                for (auto it2 = ctsSeq.begin(), end2 = ctsSeq.end(); it2 != end2; ++it2)
                {
                    try
                    {
                        SimTime begin = sim::startDate();
                        if (it2->getBegin().present())
                        {
                            begin =
                                UnitParse::readDate(it2->getBegin().get(), UnitParse::STEPS /*STEPS is only for backwards compatibility*/);
                        }
                        SimTime end = sim::future();
                        if (it2->getEnd().present())
                        {
                            end = UnitParse::readDate(it2->getEnd().get(), UnitParse::STEPS /*STEPS is only for backwards compatibility*/);
                        }
                        continuous.push_back(ContinuousHumanDeployment(begin, end, *it2, intervention, subPop, complement));
                    }
                    catch (const util::format_error &e)
                    {
                        throw util::xml_scenario_error(string("interventions/human/deployment/continuous/deploy: ").append(e.message()));
                    }
                }
            }
            for (auto timedIt = elt.getTimed().begin(); timedIt != elt.getTimed().end(); ++timedIt)
            {
                ComponentId subPop = ComponentId::wholePop();
                bool complement = false;
                if (timedIt->getRestrictToSubPop().present())
                {
                    const string &subPopStr = timedIt->getRestrictToSubPop().get().getId();
                    subPop = getComponentId(subPopStr);
                    complement = timedIt->getRestrictToSubPop().get().getComplement();
                }
                try
                {
                    multimap<SimTime, const scnXml::MassDeployment *> deployTimes;
                    for (auto it2 = timedIt->getDeploy().begin(), end2 = timedIt->getDeploy().end(); it2 != end2; ++it2)
                    {
                        SimTime date = UnitParse::readDate(it2->getTime(), UnitParse::STEPS /*STEPS is only for backwards compatibility*/);

                        if (it2->getRepeatStep().present() != it2->getRepeatEnd().present())
                        {
                            throw util::xml_scenario_error("deploy: use of repeatStep or repeatEnd without other");
                        }
                        if (it2->getRepeatStep().present())
                        {
                            SimTime step = UnitParse::readDuration(it2->getRepeatStep().get(), UnitParse::NONE);
                            if (step < sim::oneTS()) { throw util::xml_scenario_error("deploy: repeatStep must be >= 1"); }
                            SimTime end = UnitParse::readDate(it2->getRepeatEnd().get(), UnitParse::NONE);
                            while (date < end)
                            {
                                deployTimes.insert(make_pair(date, &*it2));
                                date = date + step;
                            }
                        }
                        else
                        {
                            deployTimes.insert(make_pair(date, &*it2));
                        }
                    }
                    SimTime lastDate = sim::never();
                    for (auto deploy = deployTimes.begin(), dpEnd = deployTimes.end(); deploy != dpEnd; ++deploy)
                    {
                        if (deploy->first == lastDate)
                        {
                            ostringstream msg;
                            msg << "Timed deployment of components ";
                            bool first = true;
                            for (auto cp = elt.getComponent().begin(), cpEnd = elt.getComponent().end(); cp != cpEnd; ++cp)
                            {
                                if (!first) { msg << ", "; }
                                else
                                {
                                    first = false;
                                }
                                msg << cp->getId();
                            }
                            msg << " has multiple deplyoments at date " << lastDate << " (step " << sim::inSteps(lastDate - sim::startDate())
                                << ')';
                            throw util::xml_scenario_error(msg.str());
                        }
                        lastDate = deploy->first;
                    }

                    if (timedIt->getCumulativeCoverage().present())
                    {
                        const scnXml::CumulativeCoverage &cumCov = timedIt->getCumulativeCoverage().get();
                        ComponentId cumCovComponent = getComponentId(cumCov.getComponent());
                        for (auto deploy = deployTimes.begin(), end = deployTimes.end(); deploy != end; ++deploy)
                        {
                            timed.push_back(unique_ptr<TimedDeployment>(new TimedCumulativeHumanDeployment(
                                deploy->first, *deploy->second, intervention, subPop, complement, cumCovComponent)));
                        }
                    }
                    else
                    {
                        for (auto deploy = deployTimes.begin(), end = deployTimes.end(); deploy != end; ++deploy)
                        {
                            timed.push_back(unique_ptr<TimedDeployment>(
                                new TimedHumanDeployment(deploy->first, *deploy->second, intervention, subPop, complement)));
                        }
                    }
                }
                catch (const util::format_error &e)
                {
                    throw util::xml_scenario_error(string("interventions/human/deployment/timed/deploy/time: ").append(e.message()));
                }
            }
        }
    }
    if (intervElt.getImportedInfections().present())
    {
        const scnXml::ImportedInfections &ii = intervElt.getImportedInfections().get();
        importedInfections.init(ii);
    }
    // Must come after vaccines are initialised:
    if (intervElt.getInsertR_0Case().present())
    {
        throw util::xml_scenario_error("R_0 code is disabled to reduce "
                                       "maintenance. If you need it, either use an old OpenMalaria "
                                       "version (pre 32) or request it be reenabled.");
        // code disabled; search for R_0 to find more blocks like this:
#if 0
        const scnXml::InsertR_0Case& elt = intervElt.getInsertR_0Case().get();
        if( elt.getTimedDeployment().size() > 0 ){
            Vaccine::verifyEnabledForR_0();
            // timed deployments:
            for( auto it = elt.getTimedDeployment().begin(); it != elt.getTimedDeployment().end(); ++it ){
                SimTime date = UnitParse::readDate(it->getTime(), UnitParse::STEPS /*STEPS is only for backwards compatibility*/);
                timed.push_back( new TimedR_0Deployment( date ) );
            }
        }
#endif
    }
    if (intervElt.getUninfectVectors().present())
    {
        const scnXml::UninfectVectors &elt = intervElt.getUninfectVectors().get();
        if (elt.getTimedDeployment().size() > 0)
        {
            // timed deployments:
            for (auto it = elt.getTimedDeployment().begin(); it != elt.getTimedDeployment().end(); ++it)
            {
                SimTime date = UnitParse::readDate(it->getTime(), UnitParse::STEPS /*STEPS is only for backwards compatibility*/);
                timed.push_back(unique_ptr<TimedDeployment>(new TimedUninfectVectorsDeployment(date)));
            }
        }
    }
    if (intervElt.getVectorPop().present())
    {
        Transmission::VectorModel *vectorModel = dynamic_cast<Transmission::VectorModel *>(&transmission);
        if (vectorModel)
        {
            typedef scnXml::VectorPop::InterventionSequence SeqT;
            const SeqT &seq = intervElt.getVectorPop().get().getIntervention();
            size_t instance = 0;
            for (auto it = seq.begin(), end = seq.end(); it != end; ++it)
            {
                const scnXml::VectorIntervention &elt = *it;
                if (elt.getTimed().present())
                {
                    vectorModel->initVectorInterv(elt.getDescription().getAnopheles(), instance, elt.getName());

                    const scnXml::TimedBaseList::DeploySequence &seq = elt.getTimed().get().getDeploy();
                    for (auto it = seq.begin(); it != seq.end(); ++it)
                    {
                        SimTime date = UnitParse::readDate(it->getTime(), UnitParse::STEPS /*STEPS is only for backwards compatibility*/);
                        timed.push_back(unique_ptr<TimedDeployment>(new TimedVectorDeployment(date, instance)));
                    }
                    instance++;
                }
            }
        }
        else
            throw util::xml_scenario_error("vector model interventions can not be used with the non-vector model");
    }
    if (intervElt.getAddNonHumanHosts().present())
    {
        Transmission::VectorModel *vectorModel = dynamic_cast<Transmission::VectorModel *>(&transmission);
        if (vectorModel)
        {
            typedef scnXml::AddNonHumanHosts::NonHumanHostsSequence SeqT;
            const SeqT &seq = intervElt.getAddNonHumanHosts().get().getNonHumanHosts();
            size_t instance = 0;
            for (auto it = seq.begin(), end = seq.end(); it != end; ++it)
            {
                const scnXml::NonHumanHosts2 &elt = *it;
                if (elt.getTimed().present())
                {
                    for (const scnXml::Deploy2 &deploy : elt.getTimed().get().getDeploy())
                    {
                        SimTime date =
                            UnitParse::readDate(deploy.getTime(), UnitParse::STEPS /*STEPS is only for backwards compatibility*/);
                        SimTime lifespan = UnitParse::readDuration(deploy.getLifespan(), UnitParse::NONE);
                        timed.push_back(unique_ptr<TimedDeployment>(new TimedAddNonHumanHostsDeployment(date, elt.getName(), lifespan, elt.getDescription().getAnopheles(), transmission)));
                    }
                    instance++;
                }
            }
        }
        else
            throw util::xml_scenario_error("vector model interventions can not be used with the non-vector model");
    }
    if (intervElt.getNonHumanHostsModifications().present())
    {
        Transmission::VectorModel *vectorModel = dynamic_cast<Transmission::VectorModel *>(&transmission);
        if (vectorModel)
        {
            typedef scnXml::NonHumanHostsModifications::InterventionSequence SeqT;
            const SeqT &seq = intervElt.getNonHumanHostsModifications().get().getIntervention();
            size_t instance = 0;
            for (auto it = seq.begin(), end = seq.end(); it != end; ++it)
            {
                const scnXml::NonHumanHostsIntervention &elt = *it;
                if (elt.getTimed().present())
                {
                    const scnXml::DecayFunction &decay = elt.getDecay();
                    const scnXml::TimedBaseList::DeploySequence &seq = elt.getTimed().get().getDeploy();
                    for (auto it = seq.begin(); it != seq.end(); ++it)
                    {
                        SimTime date = UnitParse::readDate(it->getTime(), UnitParse::STEPS /*STEPS is only for backwards compatibility*/);
                        timed.push_back(
                            unique_ptr<TimedDeployment>(new TimedNonHumanHostsDeployment(date, instance, elt.getNonHumanHostsName(), elt.getDescription().getAnopheles(), decay, transmission)));
                    }
                    instance++;
                }
            }
        }
        else
            throw util::xml_scenario_error("vector model interventions can not be used with the non-vector model");
    }
    if (intervElt.getVectorTrap().present())
    {
        Transmission::VectorModel *vectorModel = dynamic_cast<Transmission::VectorModel *>(&transmission);
        if (vectorModel)
        {
            size_t instance = 0;
            for (const scnXml::VectorTrap &trap : intervElt.getVectorTrap().get().getIntervention())
            {
                vectorModel->initVectorTrap(trap.getDescription(), instance, trap.getName());
                if (trap.getTimed().present())
                {
                    for (const scnXml::Deploy1 &deploy : trap.getTimed().get().getDeploy())
                    {
                        SimTime date = UnitParse::readDate(deploy.getTime(), UnitParse::STEPS);
                        double ratio = deploy.getRatioToHumans();
                        SimTime lifespan = UnitParse::readDuration(deploy.getLifespan(), UnitParse::NONE);
                        timed.push_back(unique_ptr<TimedDeployment>(new TimedTrapDeployment(date, instance, ratio, lifespan)));
                    }
                }
                instance += 1;
            }
        }
        else
            throw util::xml_scenario_error("vector model interventions can not be used with the non-vector model");
    }

    // lists must be sorted, increasing
    // For reproducability, we need to use stable_sort, not sort.
    stable_sort(continuous.begin(), continuous.end(), byDeployTime);
    stable_sort(timed.begin(), timed.end(), byDeployTime);

    // make sure the list ends with something always in the future, so we don't
    // have to check nextTimed is within range:
    timed.push_back(unique_ptr<TimedDeployment>(new DummyTimedDeployment()));

    if (util::CommandLine::option(util::CommandLine::PRINT_INTERVENTIONS))
    {
        cout << "Continuous deployments:" << endl << "begin\tend\tage\tsub pop\tcompl\tcoverag\tcomponents" << endl;
        for (auto it = continuous.begin(); it != continuous.end(); ++it)
        {
            it->print_details(std::cout);
            cout << endl;
        }
        cout << "Timed deployments:" << endl << "time\tmin age\tmax age\tsub pop\tcompl\tcoverag\tcomponents" << endl;
        for (auto &deploy : timed)
        {
            deploy->print_details(std::cout);
            cout << endl;
        }
        cout << "Human components:" << endl;
        for (auto &component : humanComponents)
        {
            component->print_details(cout);
            cout << endl;
        }
    }
}

ComponentId InterventionManager::getComponentId(const string textId)
{
    auto it = identifierMap.find(textId);
    if (it == identifierMap.end())
    {
        ostringstream msg;
        msg << "unable to find an intervention component with id \"" << textId
            << "\" (wrong name, no definition or used before definition?)";
        throw util::xml_scenario_error(msg.str());
    }
    return it->second;
}

void InterventionManager::loadFromCheckpoint(vector<Host::Human> &population, Transmission::TransmissionModel &transmission)
{
    SimTime date = sim::intervDate();
    // We need to re-deploy changeHS and changeEIR interventions, but nothing
    // else. nextTimed should be zero so we can go through all past interventions.
    // Only redeploy those which happened before this time step.
    assert(nextTimed == 0);
    while (timed[nextTimed]->date < date)
    {
        TimedDeployment *deployment = &*timed[nextTimed];
        if (dynamic_cast<TimedChangeHSDeployment *>(deployment) != 0 || dynamic_cast<TimedChangeEIRDeployment *>(deployment) != 0)
        {
            // Note: neither changeHS nor changeEIR interventions care what the
            // current time step is when they are deployed, so we don't need to
            // tell them the deployment date.
            deployment->deploy(population, transmission);
        }
        nextTimed += 1;
    }
}

void InterventionManager::deploy(vector<Host::Human> &population, Transmission::TransmissionModel &transmission)
{
    if (sim::intervTime() < sim::zero()) return;

    // deploy imported infections (not strictly speaking an intervention)
    importedInfections.import(population);

    // deploy timed interventions
    SimTime now = sim::intervDate();
    while (timed[nextTimed]->date <= now)
    {
        if (util::CommandLine::option(util::CommandLine::VERBOSE))
        {
            cout << "deploy::" << now << "::Timed";
            timed[nextTimed]->print_details(cout);
            cout << endl;
        }

        timed[nextTimed]->deploy(population, transmission);
        nextTimed += 1;
    }

    // deploy continuous interventions
    for (Host::Human &human : population)
    {
        uint32_t nextCtsDist = human.nextCtsDist;
        // deploy continuous interventions
        while (nextCtsDist < continuous.size())
        {
            if (util::CommandLine::option(util::CommandLine::VERBOSE))
            {
                cout << "deploy::" << now << "::Continuous";
                continuous[nextCtsDist].print_details(cout);
                cout << endl;
            }
            if (!continuous[nextCtsDist].filterAndDeploy(human)) break; // deployment (and all remaining) happens in the future
            nextCtsDist = ++human.nextCtsDist;
        }
    }
}

} // namespace interventions
} // namespace OM
