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

#include "Clinical/CMDecisionTree.h"
#include "Clinical/ClinicalModel.h"
#include "Host/WithinHost/WHInterface.h"
#include "Host/WithinHost/Diagnostic.h"
#include "PkPd/LSTMTreatments.h"
#include "util/random.h"
#include "util/errors.h"
#include "util/UnitParse.h"
#include "interventions/Interfaces.h"
#include "interventions/InterventionManager.h"

#include <limits>
#include <list>
#include <vector>
#include <map>

namespace OM { namespace Clinical {

using WithinHost::Diagnostic;
using WithinHost::diagnostics;
using namespace OM::util;
using std::map;
using std::vector;


// ———  special 'multiple' node  ———

/**
 * Branch out to multiple descendants.
 */
class CMDTMultiple : public CMDecisionTree {
public:
    static const CMDecisionTree& create( const ::scnXml::DTMultiple& node, bool isUC );
    
protected:
    virtual bool operator==( const CMDecisionTree& that ) const{
        if( this == &that ) return true; // short cut: same object thus equivalent
        const CMDTMultiple* p = dynamic_cast<const CMDTMultiple*>( &that );
        if( p == 0 ) return false;      // different type of node
        if( children.size() != p->children.size() ) return false;
        for( size_t i = 0; i < children.size(); ++i ){
            if( children[i] != p->children[i] ) return false;
        }
        return true;    // no tests failed; must be the same
    }
    
    virtual CMDTOut exec( CMHostData hostData ) const{
        CMDTOut result(false);
        for( auto it = children.begin(), end = children.end(); it != end; ++it ) {
            CMDTOut r2 = (*it)->exec( hostData );
            result.treated = result.treated || r2.treated;
        }
        return result;
    }
    
private:
    CMDTMultiple( /*size_t capacity*/ ){
//         children.reserve( capacity );
    }
    
    typedef vector<const CMDecisionTree*> Children_t;
    Children_t children;
};

// ———  branching nodes  ———

/**
 * Should the patient recieve first or second line treatment?
 */
class CMDTCaseType : public CMDecisionTree {
public:
    static const CMDecisionTree& create( const ::scnXml::DTCaseType& node, bool isUC );
    
protected:
    virtual bool operator==( const CMDecisionTree& that ) const{
        if( this == &that ) return true; // short cut: same object thus equivalent
        const CMDTCaseType* p = dynamic_cast<const CMDTCaseType*>( &that );
        if( p == 0 ) return false;      // different type of node
        if( firstLine != p->firstLine ) return false;
        if( secondLine != p->secondLine ) return false;
        return true;    // no tests failed; must be the same
    }
    
    virtual CMDTOut exec( CMHostData hostData ) const{
        // Uses of this in complicated cases should trigger an exception during initialisation.
        assert( (hostData.pgState & Episode::SICK) && !(hostData.pgState & Episode::COMPLICATED) );
        
        if( hostData.pgState & Episode::SECOND_CASE ) return secondLine.exec( hostData );
        else return firstLine.exec( hostData );
    }
    
private:
    CMDTCaseType( const CMDecisionTree& firstLine,
                  const CMDecisionTree& secondLine ) :
                  firstLine(firstLine), secondLine(secondLine)
                  {}
    
    const CMDecisionTree& firstLine;
    const CMDecisionTree& secondLine;
};

/**
 * Should the patient recieve a different treatment based on its infection type?
 */
class CMDTInfectionOrigin : public CMDecisionTree {
public:
    static const CMDecisionTree& create( const ::scnXml::DTInfectionOrigin& node, bool isUC );
    
protected:
    virtual bool operator==( const CMDecisionTree& that ) const{
        if( this == &that ) return true; // short cut: same object thus equivalent
        const CMDTInfectionOrigin* p = dynamic_cast<const CMDTInfectionOrigin*>( &that );
        if( p == 0 ) return false;      // different type of node
        if( imported != p->imported ) return false;
        if( introduced != p->introduced ) return false;
        if( indigenous != p->indigenous ) return false;
        return true;    // no tests failed; must be the same
    }
    
    virtual CMDTOut exec( CMHostData hostData ) const{
        WithinHost::InfectionOrigin InfectionOrigin = hostData.withinHost().getInfectionOrigin();
        if(InfectionOrigin == WithinHost::InfectionOrigin::Imported)
            return imported.exec( hostData );
        else if (InfectionOrigin == WithinHost::InfectionOrigin::Introduced)
            return introduced.exec( hostData );
        else
            return indigenous.exec( hostData );
    }
    
private:
    CMDTInfectionOrigin( const CMDecisionTree& imported,
                  const CMDecisionTree& introduced,
                  const CMDecisionTree& indigenous ) :
                  imported(imported), introduced(introduced), indigenous(indigenous)
                  {}
    
    const CMDecisionTree& imported, &introduced, &indigenous;
};

/**
 * Use a diagnostic, and call another decision tree based on the outcome.
 */
class CMDTDiagnostic : public CMDecisionTree {
public:
    static const CMDecisionTree& create( const ::scnXml::DTDiagnostic& node, bool isUC );
    
protected:
    virtual bool operator==( const CMDecisionTree& that ) const{
        if( this == &that ) return true; // short cut: same object thus equivalent
        const CMDTDiagnostic* p = dynamic_cast<const CMDTDiagnostic*>( &that );
        if( p == 0 ) return false;      // different type of node
        if( diagnostic != p->diagnostic ) return false;
        if( positive != p->positive ) return false;
        if( negative != p->negative ) return false;
        return true;    // no tests failed; must be the same
    }
    
    virtual CMDTOut exec( CMHostData hostData ) const{
        CMDTOut result;
        
        if( hostData.withinHost().diagnosticResult( hostData.human.rng, diagnostic ) )
            result = positive.exec( hostData );
        else
            result = negative.exec( hostData );
        result.screened = true;

        return result;
    }
    
private:
    CMDTDiagnostic( const Diagnostic& diagnostic,
        const CMDecisionTree& positive,
        const CMDecisionTree& negative ) :
        diagnostic(diagnostic),
        positive(positive), negative(negative)
    {}
    
    const Diagnostic& diagnostic;
    const CMDecisionTree& positive;
    const CMDecisionTree& negative;
};

class CMDTUncomplicated : public CMDecisionTree {
public:
    static const CMDecisionTree& create( const ::scnXml::DTUncomplicated& node, bool isUC );
    
protected:
    virtual bool operator==( const CMDecisionTree& that ) const{
        if( this == &that ) return true; // short cut: same object thus equivalent
        const CMDTUncomplicated* p = dynamic_cast<const CMDTUncomplicated*>( &that );
        if( p == 0 ) return false;      // different type of node
        // if( diagnostic != p->diagnostic ) return false;
        if( positive != p->positive ) return false;
        if( negative != p->negative ) return false;
        return true;    // no tests failed; must be the same
    }
    
    virtual CMDTOut exec( CMHostData hostData ) const{
        const Clinical::Episode &latest = hostData.human.clinicalModel->getLatestReport();

        if ( ((hostData.pgState & Episode::SICK) && !(hostData.pgState & Episode::COMPLICATED)) || (hostData.pgState & Episode::MALARIA) )
        {
            if (latest.time + memory >= sim::nowOrTs0())
                return positive.exec( hostData );
            else
                return negative.exec( hostData );
        }
        else
            return negative.exec( hostData );
    }
    
private:
    CMDTUncomplicated(
        const SimTime memory,
        const CMDecisionTree& positive,
        const CMDecisionTree& negative ) :
        memory(memory), positive(positive), negative(negative)
    {
        if(memory > ClinicalModel::hsMemory())
            throw util::xml_scenario_error( "<uncomplicated> memory parameter must be less than or equal to the healthsystem memory (hsmemory parameter)" );
    }
    
    const SimTime memory;
    const CMDecisionTree& positive;
    const CMDecisionTree& negative;
};

class CMDTSevere : public CMDecisionTree {
public:
    static const CMDecisionTree& create( const ::scnXml::DTSevere& node, bool isUC );
    
protected:
    virtual bool operator==( const CMDecisionTree& that ) const{
        if( this == &that ) return true; // short cut: same object thus equivalent
        const CMDTSevere* p = dynamic_cast<const CMDTSevere*>( &that );
        if( p == 0 ) return false;      // different type of node
        // if( diagnostic != p->diagnostic ) return false;
        if( positive != p->positive ) return false;
        if( negative != p->negative ) return false;
        return true;    // no tests failed; must be the same
    }
    
    virtual CMDTOut exec( CMHostData hostData ) const{
        CMDTOut result;

        if (hostData.pgState & Episode::COMPLICATED)
            result = positive.exec( hostData );
        else
            result = negative.exec( hostData );

        // result.screened = true;
        return result;
    }
    
private:
    CMDTSevere(
        const CMDecisionTree& positive,
        const CMDecisionTree& negative ) :
        positive(positive), negative(negative)
    {}
    
    const CMDecisionTree& positive;
    const CMDecisionTree& negative;
};

/**
 * Choose a branch randomly.
 * 
 * TODO(optimisation): can we re-order the trees such that this can be abstracted out or
 * used only once?
 */
class CMDTRandom : public CMDecisionTree {
public:
    static const CMDecisionTree& create( const ::scnXml::DTRandom& node, bool isUC );
    
protected:
    virtual bool operator==( const CMDecisionTree& that ) const{
        if( this == &that ) return true; // short cut: same object thus equivalent
        const CMDTRandom* p = dynamic_cast<const CMDTRandom*>( &that );
        if( p == 0 ) return false;      // different type of node
        if( branches.size() != p->branches.size() ) return false;
        for( auto it1 = branches.begin(), it2 = p->branches.begin();
            it1 != branches.end(); ++it1, ++it2 )
        {
            if( it1->first != it2->first ) return false;
            if( it1->second != it2->second ) return false;
        }
        return true;    // no tests failed; must be the same
    }
    
    virtual CMDTOut exec( CMHostData hostData ) const{
        auto it = branches.upper_bound( hostData.human.rng.uniform_01() );
        assert( it != branches.end() );
        return it->second->exec( hostData );
    }
    
private:
    CMDTRandom(){}
    
    // keys are cumulative probabilities; last entry should equal 1
    typedef map<double,const CMDecisionTree*> Branches_t;
    Branches_t branches;
};

/**
 * Choose what to do based on the patient's age.
 */
class CMDTAge : public CMDecisionTree {
public:
    static const CMDecisionTree& create( const ::scnXml::DTAge& node, bool isUC );
    
protected:
    virtual bool operator==( const CMDecisionTree& that ) const{
        if( this == &that ) return true; // short cut: same object thus equivalent
        const CMDTAge* p = dynamic_cast<const CMDTAge*>( &that );
        if( p == 0 ) return false;      // different type of node
        if( branches.size() != p->branches.size() ) return false;
        for( auto it1 = branches.begin(), it2 = p->branches.begin();
            it1 != branches.end(); ++it1, ++it2 )
        {
            if( it1->first != it2->first ) return false;
            if( it1->second != it2->second ) return false;
        }
        return true;    // no tests failed; must be the same
    }
    
    virtual CMDTOut exec( CMHostData hostData ) const{
        // age is that of human at start of time step (i.e. may be as low as 0)
        auto it = branches.upper_bound( hostData.ageYears );
        if( it == branches.end() )
            throw TRACED_EXCEPTION( "bad age-based decision tree switch", util::Error::PkPd );
        return it->second->exec( hostData );
    }
    
private:
    CMDTAge() {}
    
    // keys are upper bounds of age categories
    typedef map<double,const CMDecisionTree*> Branches_t;
    Branches_t branches;
};


// ———  action nodes  ———

/** Do nothing. **/
class CMDTNoTreatment : public CMDecisionTree {
protected:
    virtual bool operator==( const CMDecisionTree& that ) const{
        if( this == &that ) return true; // short cut: same object thus equivalent
        const CMDTNoTreatment* p = dynamic_cast<const CMDTNoTreatment*>( &that );
        return p != 0;  // same type: is equivalent
    }
    
    virtual CMDTOut exec( CMHostData hostData ) const{
        return CMDTOut(false);
    }
};

class CMDTReport : public CMDecisionTree {
public:
    CMDTReport( const scnXml::DecisionTree::ReportSequence& seq ){
        outIds.reserve( seq.size() );
        for( const scnXml::DTReport& reportElt : seq ){
            outIds.push_back( reportElt.getOutputNumber() );
        }
        assert( outIds.size() > 0 );        // CMDTTreatPKPD should not be used in this case
    }
protected:
    virtual bool operator==( const CMDecisionTree& that ) const{
     if( this == &that ) return true; // short cut: same object thus equivalent
        const CMDTReport* p = dynamic_cast<const CMDTReport*>( &that );
        if( p == 0 ) return false;      // different type of node
        if( outIds.size() != p->outIds.size() ) return false;
        for( auto it1 = outIds.begin(), it2 = p->outIds.begin();
            it1 != outIds.end(); ++it1, ++it2 )
        {
            if( *it1 != *it2 ) return false;
        }
        return true;    // no tests failed; must be the same
    }
    
    virtual CMDTOut exec( CMHostData hostData ) const{
        for( const size_t outId : outIds ){
            mon::reportEventMHI_CMDT( mon::MCD_CMDT_REPORT, hostData.human, 1, outId);
        }
        return CMDTOut(false);
    }
private:
    vector<int> outIds;
};

/** Report treament without affecting parasites. **/
class CMDTTreatFailure : public CMDecisionTree {
protected:
    virtual bool operator==( const CMDecisionTree& that ) const{
        if( this == &that ) return true; // short cut: same object thus equivalent
        const CMDTTreatFailure* p = dynamic_cast<const CMDTTreatFailure*>( &that );
        return p != 0;  // same type: is equivalent
    }
    
    virtual CMDTOut exec( CMHostData hostData ) const{
        return CMDTOut(true /*report treatment*/);
    }
};

/**
 * Deliver one or more treatments via the PK/PD model.
 */
class CMDTTreatPKPD : public CMDecisionTree {
public:
    CMDTTreatPKPD( const scnXml::DecisionTree::TreatPKPDSequence& seq ){
        treatments.reserve( seq.size() );
        for( const scnXml::DTTreatPKPD& treatElt : seq ){
            treatments.push_back( TreatInfo(
                treatElt.getSchedule(),
                treatElt.getDosage(),
                treatElt.getDelay_h()
            ) );
        }
        assert( treatments.size() > 0 );        // CMDTTreatPKPD should not be used in this case
    }
    
protected:
    virtual bool operator==( const CMDecisionTree& that ) const{
        if( this == &that ) return true; // short cut: same object thus equivalent
        const CMDTTreatPKPD* p = dynamic_cast<const CMDTTreatPKPD*>( &that );
        if( p == 0 ) return false;      // different type of node
        if( treatments.size() != p->treatments.size() ) return false;
        for( auto it1 = treatments.begin(), it2 = p->treatments.begin();
            it1 != treatments.end(); ++it1, ++it2 )
        {
            if( *it1 != *it2 ) return false;
        }
        return true;    // no tests failed; must be the same
    }
    
    virtual CMDTOut exec( CMHostData hostData ) const{
        for( const TreatInfo& treatment : treatments ){
            hostData.withinHost().treatPkPd( treatment.schedule, treatment.dosage, hostData.ageYears, 0.0 );
        }
        return CMDTOut(true);
    }
    
private:
    struct TreatInfo{
        TreatInfo( const string& s, const string& d, double h ) :
            schedule(PkPd::LSTMTreatments::findSchedule(s)),
            dosage(PkPd::LSTMTreatments::findDosages(d)),
            delay_h(h) {}
        inline bool operator!=( const TreatInfo& that )const{
            return schedule != that.schedule ||
                dosage != that.dosage ||
                delay_h != that.delay_h;
        }
        size_t schedule;        // index of the schedule
        size_t dosage;          // name of the dosage table
        double delay_h;         // delay in hours
    };
    vector<TreatInfo> treatments;
};

/**
 * Deliver one or more simple treatments
 */
class CMDTTreatSimple : public CMDecisionTree {
public:
    CMDTTreatSimple( const scnXml::DecisionTree::TreatSimpleSequence& elt ) :
        timeLiver(sim::zero()), timeBlood(sim::zero())
    {
        for (const auto &treatElt : elt) {
            //NOTE: this code is currently identical to that in SimpleTreatComponent
            try{
                SimTime durL = UnitParse::readShortDuration( treatElt.getDurationLiver(), UnitParse::NONE ),
                    durB = UnitParse::readShortDuration( treatElt.getDurationBlood(), UnitParse::NONE );
                SimTime neg1 = -sim::oneTS();
                if( durL < neg1 || durB < neg1 ){
                    throw util::xml_scenario_error( "treatSimple: cannot have durationBlood or durationLiver less than -1" );
                }
                timeLiver.push_back(durL);
                timeBlood.push_back(durB);
            }catch( const util::format_error& e ){
                throw util::xml_scenario_error( string("treatSimple: ").append(e.message()) );
            }
        }
    }
    
protected:
    virtual bool operator==( const CMDecisionTree& that ) const{
     if( this == &that ) return true; // short cut: same object thus equivalent
        const CMDTTreatSimple* p = dynamic_cast<const CMDTTreatSimple*>( &that );
        if( p == 0 ) return false;      // different type of node
        if( timeLiver.size() != p->timeLiver.size() ) return false;
        if( timeBlood.size() != p->timeBlood.size() ) return false;
        for( auto it1 = timeLiver.begin(), it2 = p->timeLiver.begin();
            it1 != timeLiver.end(); ++it1, ++it2 )
        {
            if( *it1 != *it2 ) return false;
        }
        return true;    // no tests failed; must be the same
    }

    
    virtual CMDTOut exec( CMHostData hostData ) const{
        bool bsTreatment = false;
        for(size_t i=0; i<timeLiver.size(); i++)
        {
            bsTreatment = hostData.withinHost().treatSimple( hostData.human, timeLiver[i], timeBlood[i] );
        }
        return CMDTOut(bsTreatment);
    }
    
private:
    vector<SimTime> timeLiver, timeBlood;
};

/**
 * Deploy one or more interventions.
 */
/*NOTE: we use interventions::HumanIntervention to read and sort the list of
 * interventions to be deployed, though we don't use pointers it finds. */
class CMDTDeploy : public CMDecisionTree, interventions::HumanIntervention {
public:
    CMDTDeploy( const scnXml::DecisionTree::DeploySequence& seq ) :
        HumanIntervention(seq) {}
    
protected:
    virtual bool operator==( const CMDecisionTree& that ) const{
        if( this == &that ) return true; // short cut: same object thus equivalent
        const CMDTDeploy* p = dynamic_cast<const CMDTDeploy*>( &that );
        if( p == 0 ) return false;      // different type of node
        if( components.size() != p->components.size() ) return false;
        for( auto it1 = components.begin(), it2 = p->components.begin();
            it1 != components.end(); ++it1, ++it2 )
        {
            // here we can compare pointers since components are always de-duplicated
            if( *it1 != *it2 ) return false;
        }
        return true;    // no tests failed; must be the same
    }
    
    virtual CMDTOut exec( CMHostData hostData ) const{
        deploy( hostData.human,
                  mon::Deploy::TREAT,
                  interventions::VaccineLimits(/*default initialise: no limits*/) );
        
        //NOTE: it's not intuitively obvious what value should be returned here
        // in the case of intervention deployment. This at least means that
        // repeat seekers get second-line treatment.
        return CMDTOut(false);
    }
};

class CMDTCohort : public CMDecisionTree {
public:
    static const CMDecisionTree& create( const ::scnXml::DTCohort& node, bool isUC );
    
protected:
    virtual bool operator==( const CMDecisionTree& that ) const{
        if( this == &that ) return true; // short cut: same object thus equivalent
        const CMDTCohort* p = dynamic_cast<const CMDTCohort*>( &that );
        if( p == 0 ) return false;      // different type of node
        // if( diagnostic != p->diagnostic ) return false;
        if( positive != p->positive ) return false;
        if( negative != p->negative ) return false;
        return true;    // no tests failed; must be the same
    }
    
    virtual CMDTOut exec( CMHostData hostData ) const{
        CMDTOut result;

        // Rely on the health system memory to not count the same episode twice
        if(hostData.human.isInSubPop(component))
            result = positive.exec( hostData );
        else
            result = negative.exec( hostData );

        // result.screened = true;
        return result;
    }
    
private:
    CMDTCohort(
        const std::string component,
        const CMDecisionTree& positive,
        const CMDecisionTree& negative ) :
        component(interventions::InterventionManager::getComponentId(component)), positive(positive), negative(negative)
    {
    }
    
    interventions::ComponentId component;
    const CMDecisionTree& positive;
    const CMDecisionTree& negative;
};

// ———  static functions  ———

// Memory management: lists all decisions and frees memory at program exit.
// We store pointers into this list, so elements must not move.
vector<unique_ptr<CMDecisionTree>> decision_library;

// Saves a decision to decision_library, making it const.
// Also optimises away duplicates
const CMDecisionTree& save_decision( CMDecisionTree* decision ){
    // We search the library for a duplicate, and delete this one if there is
    // a duplicate. Note that this is not implemented efficiently, but a little
    // wasted time at start up is hardly a concern.
    for( auto& d : decision_library ){
        if( *d == *decision ){
            return *d;
        }
    }
    
    // No match: add to the library.
    decision_library.push_back( unique_ptr<CMDecisionTree>(decision) );
    return *decision_library.back();
}

const CMDecisionTree& CMDecisionTree::create( const scnXml::DecisionTree& node, bool isUC){
    if( node.getMultiple().present() ) return CMDTMultiple::create( node.getMultiple().get(), isUC);
    // branching nodes
    if( node.getCaseType().present() ) return CMDTCaseType::create( node.getCaseType().get(), isUC);
    if( node.getInfectionOrigin().present() ) return CMDTInfectionOrigin::create( node.getInfectionOrigin().get(), isUC);
    if( node.getDiagnostic().present() ) return CMDTDiagnostic::create( node.getDiagnostic().get(), isUC);
    if( node.getUncomplicated().present() ) return CMDTUncomplicated::create( node.getUncomplicated().get(), isUC);
    if( node.getSevere().present() ) return CMDTSevere::create( node.getSevere().get(), true);
    if( node.getRandom().present() ) return CMDTRandom::create( node.getRandom().get(), isUC);
    if( node.getAge().present() ) return CMDTAge::create( node.getAge().get(), isUC);
    if( node.getCohort().present() ) return CMDTCohort::create( node.getCohort().get(), isUC);

    // action nodes
    if( node.getNoTreatment().present() ) return save_decision( new CMDTNoTreatment() );
    if( node.getReport().size() ) return save_decision( new CMDTReport(node.getReport()) );
    if( node.getTreatFailure().present() ) return save_decision( new CMDTTreatFailure() );
    if( node.getTreatPKPD().size() ) return save_decision(
        new CMDTTreatPKPD( node.getTreatPKPD() ) );
    if( node.getTreatSimple().size() ) return save_decision(
        new CMDTTreatSimple( node.getTreatSimple() ) );
    if( node.getDeploy().size() ) return save_decision(
        new CMDTDeploy( node.getDeploy() ) );
    throw xml_scenario_error( "unterminated decision tree" );
}

const CMDecisionTree& CMDTMultiple::create( const scnXml::DTMultiple& node, bool isUC ){
    CMDTMultiple* self = new CMDTMultiple();
    for( const scnXml::DTCaseType& sn : node.getCaseType() ){
        self->children.push_back( &CMDTCaseType::create(sn, isUC) );
    }
    for( const scnXml::DTInfectionOrigin& sn : node.getInfectionOrigin() ){
        self->children.push_back( &CMDTInfectionOrigin::create(sn, isUC) );
    }
    for( const scnXml::DTDiagnostic& sn : node.getDiagnostic() ){
        self->children.push_back( &CMDTDiagnostic::create(sn, isUC) );
    }
    for( const scnXml::DTUncomplicated& sn : node.getUncomplicated() ){
        self->children.push_back( &CMDTUncomplicated::create(sn, isUC) );
    }
    for( const scnXml::DTSevere& sn : node.getSevere() ){
        self->children.push_back( &CMDTSevere::create(sn, isUC) );
    }
    for( const scnXml::DTRandom& sn : node.getRandom() ){
        self->children.push_back( &CMDTRandom::create(sn, isUC) );
    }
    for( const scnXml::DTAge& sn : node.getAge() ){
        self->children.push_back( &CMDTAge::create(sn, isUC) );
    }
    for( const scnXml::DTCohort& sn : node.getCohort() ){
        self->children.push_back( &CMDTCohort::create(sn, isUC) );
    }
    if( node.getTreatPKPD().size() ){
        self->children.push_back( &save_decision(new CMDTTreatPKPD(node.getTreatPKPD())) );
    }
    if( node.getTreatSimple().size() ){
        self->children.push_back( &save_decision(new CMDTTreatSimple(node.getTreatSimple())) );
    }
    if( node.getReport().size() ){
        self->children.push_back( &save_decision(new CMDTReport(node.getReport())) );
    }
    if( node.getDeploy().size() ){
        self->children.push_back( &save_decision(new CMDTDeploy(node.getDeploy())) );
    }
    return save_decision( self );
}

const CMDecisionTree& CMDTCaseType::create( const scnXml::DTCaseType& node, bool isUC ){
    if( !isUC ){
        throw util::xml_scenario_error( "decision tree: caseType can only be used for uncomplicated cases" );
    }
    return save_decision( new CMDTCaseType(
        CMDecisionTree::create( node.getFirstLine(), isUC ),
        CMDecisionTree::create( node.getSecondLine(), isUC )
    ) );
}

const CMDecisionTree& CMDTInfectionOrigin::create( const scnXml::DTInfectionOrigin& node, bool isUC ){
    return save_decision( new CMDTInfectionOrigin(
        CMDecisionTree::create( node.getImported(), isUC ),
        CMDecisionTree::create( node.getIntroduced(), isUC ),
        CMDecisionTree::create( node.getIndigenous(), isUC )
    ) );
}

const CMDecisionTree& CMDTDiagnostic::create( const scnXml::DTDiagnostic& node, bool isUC ){
    return save_decision( new CMDTDiagnostic(
        diagnostics::get( node.getDiagnostic() ),
        CMDecisionTree::create( node.getPositive(), isUC ),
        CMDecisionTree::create( node.getNegative(), isUC )
    ) );
}

const CMDecisionTree& CMDTUncomplicated::create( const scnXml::DTUncomplicated& node, bool isUC ){
    return save_decision( new CMDTUncomplicated(
        UnitParse::readShortDuration(node.getMemory(), UnitParse::STEPS),
        CMDecisionTree::create( node.getPositive(), isUC ),
        CMDecisionTree::create( node.getNegative(), isUC )
    ) );
}

const CMDecisionTree& CMDTSevere::create( const scnXml::DTSevere& node, bool isUC ){
    return save_decision( new CMDTSevere(
        CMDecisionTree::create( node.getPositive(), isUC ),
        CMDecisionTree::create( node.getNegative(), isUC )
    ) );
}

const CMDecisionTree& CMDTCohort::create( const scnXml::DTCohort& node, bool isUC ){
    return save_decision( new CMDTCohort(
        node.getComponent(),
        CMDecisionTree::create( node.getPositive(), isUC ),
        CMDecisionTree::create( node.getNegative(), isUC )
    ) );
}


const CMDecisionTree& CMDTRandom::create(
    const scnXml::DTRandom& node, bool isUC )
{
    CMDTRandom* result = new CMDTRandom();
    
    double cum_p = 0.0;
    for( const scnXml::Outcome& outcome : node.getOutcome() ){
        cum_p += outcome.getP();
        result->branches.insert( make_pair(cum_p, &CMDecisionTree::create( outcome, isUC )) );
    }
    
    // Test cum_p is approx. 1.0 in case the input tree is wrong. We require no
    // less than one to make sure generated random numbers are not greater than
    // the last option.
    if (cum_p < 1.0 || cum_p > 1.001){
        throw util::xml_scenario_error (
            "decision tree (random node): expected probability sum to be 1.0 but found " + to_string(cum_p));
    }
    
    return save_decision( result );
}

const CMDecisionTree& CMDTAge::create(const scnXml::DTAge& node, bool isUC){
    CMDTAge* result = new CMDTAge();
    
    double lastAge = numeric_limits<double>::quiet_NaN();
    const CMDecisionTree* lastNode = nullptr;
    for( const scnXml::Age& age : node.getAge() ){
        if( lastAge != lastAge ){
            if( age.getLb() != 0.0 )
                throw util::xml_scenario_error( "decision tree age switch must have first lower bound equal 0" );
        }else{
            if( age.getLb() <= lastAge ){
                throw util::xml_scenario_error( "decision tree age switch must list age groups in increasing order" );
            }
            double lb = age.getLb();
            result->branches.insert( make_pair(lb, lastNode) );
        }
        lastNode = &CMDecisionTree::create(age, isUC);
        lastAge = age.getLb();
    }
    double noLb = numeric_limits<double>::infinity();
    result->branches.insert( make_pair(noLb, lastNode) );
    
    return save_decision( result );
}

} }
