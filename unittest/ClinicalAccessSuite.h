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

#ifndef Hmod_ClinicalAccessSuite
#define Hmod_ClinicalAccessSuite

#include <cxxtest/TestSuite.h>
#include "Clinical/DecisionTree5Day.h"
#include "Clinical/ImmediateOutcomes.h"
#include "UnittestUtil.h"

class ClinicalAccessSuite : public CxxTest::TestSuite {
    class Access : public OM::Clinical::CM5DayCommon {
    public:
        static double firstLine() { return accessUCAny[FirstLine]; }
        static double secondLine() { return accessUCAny[SecondLine]; }
    };

    template<class ClinicalModel, class Description>
    void checkAccess(Description description) {
        for (double probability : {0.0, 0.5, 1.0}) {
            description.getPSeekOfficialCareUncomplicated2().setValue(probability);
            ClinicalModel::setHealthSystem(description);
            TS_ASSERT_DELTA(Access::firstLine(), 0.5, 1e-12);
            TS_ASSERT_EQUALS(Access::secondLine(), probability);
        }
        for (double probability : {-0.1, 1.1}) {
            description.getPSeekOfficialCareUncomplicated2().setValue(probability);
            TS_ASSERT_THROWS(ClinicalModel::setHealthSystem(description), const OM::util::xml_scenario_error&);
        }
    }

public:
    void setUp() {
        UnittestUtil::initTime(5);
        UnittestUtil::DescriptiveInfection_init();
    }

    void testDecisionTreeSecondLineAccess() {
        scnXml::DecisionTree tree;
        tree.setNoTreatment(scnXml::DTNoTreatment());
        scnXml::HSDT5Day description(
            scnXml::DoubleValue(0.3), scnXml::DoubleValue(0.2),
            scnXml::DoubleValue(0.0), scnXml::DoubleValue(0.5),
            tree, tree, scnXml::DoubleValue(1.0), scnXml::TreatmentOption());
        checkAccess<OM::Clinical::DecisionTree5Day>(description);
    }

    void testImmediateOutcomesSecondLineAccess() {
        scnXml::TreatmentDetails effective(scnXml::DoubleValue(1.0));
        effective.setACT(scnXml::DoubleValue(1.0));
        scnXml::TreatmentActions actions;
        actions.setACT(scnXml::TreatmentOption());
        scnXml::HSImmediateOutcomes description(
            scnXml::DrugRegimen("ACT", "ACT", "ACT"),
            effective, effective, effective, actions,
            scnXml::DoubleValue(0.3), scnXml::DoubleValue(0.2),
            scnXml::DoubleValue(0.0), scnXml::DoubleValue(0.5));
        checkAccess<OM::Clinical::ImmediateOutcomes>(description);
    }
};

#endif
