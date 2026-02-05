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

// Unittest for the EventScheduler case management

#ifndef Hmod_ESCaseManagementSuite
#define Hmod_ESCaseManagementSuite

#include <cxxtest/TestSuite.h>
#include "Clinical/ESCaseManagement.h"
#include "UnittestUtil.h"
#include "ExtraAsserts.h"
#include <list>

using namespace OM::Clinical;

#if 0
TODO: replace this test suite
class ESCaseManagementSuite : public CxxTest::TestSuite
{
public:
    void setUp () {
	vector<string> vals;
	vals += "extra","poor";
	dMap.dvMap.add_decision_values( "modQty", vals );
	vals.clear();
	vals += "0","5";
	dMap.dvMap.add_decision_values( "modD1", vals );
	vals.clear();
	vals += "B2";
	dMap.dvMap.add_decision_values( "modD2", vals );
	
	// We need to add treatment, hospitalised and test decisions to make the thing work, though we won't be using them:
	xsd::cxx::tree::sequence<scnXml::HSESDecision, false> decisionSeq;
	decisionSeq.push_back(
	    scnXml::HSESDecision(
		"treatment_1",	// tree
		"treatment",	// decision
		"",	// depends
		"treatment_1"	// values
	    )
	);
	decisionSeq.push_back(
	    scnXml::HSESDecision(
		"age(0-5): immediate age(5-inf): none",	// tree
		"hospitalisation",	// decision
		"age",	// depends
		"immediate,delayed,none"	// values
	    )
	);
	decisionSeq.push_back(
	    scnXml::HSESDecision(
		"RDT",	// tree
		"test",	// decision
		"",	// depends
		"none,microscopy,RDT"	// values
	    )
	);
	// dummy-decisions, required to determine modifiers (though not used)
	decisionSeq.push_back(
	    scnXml::HSESDecision(
		"extra",	// tree
		"modQty",	// decision
		"",	// depends
		"poor,extra"	// values
	    )
	);
	decisionSeq.push_back(
	    scnXml::HSESDecision(
		"0",	// tree
		"modD1",	// decision
		"",	// depends
		"0,5"	// values
	    )
	);
	decisionSeq.push_back(
	    scnXml::HSESDecision(
		"B2",	// tree
		"modD2",	// decision
		"",	// depends
		"B2"	// values
	    )
	);
	scnXml::HSESDecisions decisions;
	decisions.setDecision( decisionSeq );
	
	// Base schedule
	xsd::cxx::tree::sequence<scnXml::Medicate, false> medicateSeq;
	medicateSeq.push_back( scnXml::Medicate( "A", 1000.0, 0 ) );
	medicateSeq.push_back( scnXml::Medicate( "B", 3000.0, 0 ) );
	medicateSeq.push_back( scnXml::Medicate( "B", 3000.0, 12 ) );
	scnXml::HSESTreatmentSchedule treatSched;
	treatSched.setMedicate( medicateSeq );
	
	// Modifiers
	xsd::cxx::tree::sequence<scnXml::HSESTreatmentModifierEffect, false> modQtySeq;
	modQtySeq.push_back( scnXml::HSESTreatmentModifierEffect( "extra", "A(2),B(1.3)" ) );
	modQtySeq.push_back( scnXml::HSESTreatmentModifierEffect( "poor", " A(0.5) , B( 0.2 ) " ) );
	scnXml::HSESTreatmentModifier modQty( "modQty" );
	modQty.setMultiplyQty( modQtySeq );
	
	xsd::cxx::tree::sequence<scnXml::HSESTreatmentModifierEffect, false> modD1Seq;
	modD1Seq.push_back( scnXml::HSESTreatmentModifierEffect( "0", "A(0),B(0)" ) );
	modD1Seq.push_back( scnXml::HSESTreatmentModifierEffect( "5", "A(5),B(5)" ) );
	scnXml::HSESTreatmentModifier modD1( "modD1" );
	modD1.setDelay( modD1Seq );
	
	xsd::cxx::tree::sequence<scnXml::HSESTreatmentModifierEffect, false> modD2Seq;
	modD2Seq.push_back( scnXml::HSESTreatmentModifierEffect( "B2", "B(2 ),A(0) " ) );	// note: backwards
	scnXml::HSESTreatmentModifier modD2( "modD2" );
	modD2.setDelay( modD2Seq );
	
	xsd::cxx::tree::sequence<scnXml::HSESTreatmentModifier, false> modifierSeq;
	modifierSeq.push_back( modQty );
	modifierSeq.push_back( modD1 );
	modifierSeq.push_back( modD2 );
	
	// Treatment
	scnXml::HSESTreatment treatment_1( treatSched, "treatment_1" );
	treatment_1.setModifier( modifierSeq );
	
	scnXml::HSESTreatments treatments;
	xsd::cxx::tree::sequence<scnXml::HSESTreatment, false> treatmentSeq( 1, treatment_1 );
	treatments.setTreatment( treatmentSeq );
	
	// Final CaseManagement element
	::scnXml::HSESCaseManagement xmlCM( decisions, treatments );
	
	// use complicated tree, because it doesn't add so many unwanted decisions
	dMap.initialize( xmlCM, ESDecisionMap::Complicated, true );
    }
    
    // Note: use ETS_..., not TS_..., when assertion should throw if false (prevent dangerous operations)
    
    void testTreatments() {
	// When we give all decisions, we should get the expected medications:
	ESDecisionValue treatment_1 = dMap.dvMap.get( "treatment", "treatment_1" );	// has 3 treatments; A at time 0 and B at times 0,12
	treatment_1 |= dMap.dvMap.get( "modQty", "poor" );	// reduce quantities
	treatment_1 |= dMap.dvMap.get( "modD1", "5" );	// delay by 5 hours
	treatment_1 |= dMap.dvMap.get( "modD2", "B2" );	// delay B by 2 hours
	
	const ESTreatmentSchedule& sched = dMap.getSchedule( treatment_1 );
	
	list<MedicateData> medQueue;
	sched.apply( medQueue );
	ETS_ASSERT_EQUALS( medQueue.size(), 3u );
	
	const MedicateData& md1 = medQueue.front();
	TS_ASSERT_EQUALS( md1.abbrev, "A" );
	TS_ASSERT_EQUALS( md1.qty, 500.0 );
	TS_ASSERT_DELTA( md1.time, 5.0/24.0, 1.0e-10 );
        medQueue.pop_front();
        const MedicateData& md2 = medQueue.front();
        TS_ASSERT_EQUALS( md2.abbrev, "B" );
        TS_ASSERT_EQUALS( md2.qty, 600.0 );
        TS_ASSERT_DELTA( md2.time, 7.0/24.0, 1.0e-10 );
	medQueue.pop_front();
	const MedicateData& md3 = medQueue.front();
	TS_ASSERT_EQUALS( md3.abbrev, "B" );
	TS_ASSERT_EQUALS( md3.qty, 600.0 );
	TS_ASSERT_DELTA( md3.time, 19.0/24.0, 1.0e-10 );
    }
    
    void testExecution() {
	// Again, test output treatment, but this time evaluating decision trees.
	// Aim: check wanted decisions are _not_ optimised out and the unwanted
	// "result" decision is optimised out.
	
	UnittestUtil::initTime( 1 );
	UnittestUtil::EmpiricalWHM_setup();
	WHFalciparum* whm = dynamic_cast<WHFalciparum*>( WHInterface::createWithinHostModel( 1.0 ) );
        ETS_ASSERT( whm != 0 );
	UnittestUtil::setTotalParasiteDensity( *whm, numeric_limits< double >::infinity() );	// infinite, which means P(true outcome) should be 1.0 with an RDT test
	CMHostData hd( 16., *whm, Episode::NONE );
	ESDecisionValue outcome = dMap.determine( hd );
	delete whm;
	
	// If result test was executed, result should be true (or false). But it
	// should have been optimised out (since unused), leaving result none.
	ESDecisionValue result_mask = dMap.dvMap.getDecisionMask( "result" );
	ESDecisionValue result_none = dMap.dvMap.get( "result", "none" );
	TS_ASSERT_EQUALS( outcome & result_mask, result_none );
	
	const ESTreatmentSchedule& sched = dMap.getSchedule( outcome );
	
	list<MedicateData> medQueue;
	sched.apply( medQueue );
	ETS_ASSERT_EQUALS( medQueue.size(), 3u );
	
	const MedicateData& md1 = medQueue.front();
	TS_ASSERT_EQUALS( md1.abbrev, "A" );
	TS_ASSERT_EQUALS( md1.qty, 2000.0 );
	TS_ASSERT_DELTA( md1.time, 0.0/24.0, 1.0e-10 );
	medQueue.pop_front();
	const MedicateData& md2 = medQueue.front();
	TS_ASSERT_EQUALS( md2.abbrev, "B" );
	TS_ASSERT_EQUALS( md2.qty, 3900.0 );
	TS_ASSERT_DELTA( md2.time, 2.0/24.0, 1.0e-10 );
	medQueue.pop_front();
	const MedicateData& md3 = medQueue.front();
	TS_ASSERT_EQUALS( md3.abbrev, "B" );
	TS_ASSERT_EQUALS( md3.qty, 3900.0 );
	TS_ASSERT_DELTA( md3.time, 14.0/24.0, 1.0e-10 );
	
	TS_ASSERT_EQUALS( dMap.hospitalisation( outcome ), CMAuxOutput::NONE );
	TS_ASSERT_EQUALS( dMap.RDT_used( outcome ), true );
    }
    
private:
    ESDecisionMap dMap;
};
#endif

#endif
