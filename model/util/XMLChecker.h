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

#ifndef Hmod_util_XMLChecker
#define Hmod_util_XMLChecker

#include "Global.h"
#include <schema/scenario.h>

namespace OM
{
    namespace util
    {
        class XMLChecker
        {
        public:
            XMLChecker(){}
            ~XMLChecker(){}

            /*
            * Should be called after the input XML has been loaded and before all subsequent
            * initialisation. I.e. well before the simulation starts running.
            *
            * Performs custom checks on the input XML which are not enforced in the schema itself.
            *
            * It would be better to handle as many such issues as possible in the schema itself
            * However, some required checks are not possible to do using the version of the XSD
            * spec supported by the XML validation library used by OpenMalaria.
            *
            * Throws iff a check fails.
            *
            * The purpose of these XML checks to:
            *  - identify certain problems in the input XML as early as possible, and
            *  - enable the user to obtain a more informative error message.
            */
            void PerformPostValidationChecks(const scnXml::Scenario& scenario)
            {
                CheckModelOptionsAndParams(scenario);
            }
        private:

            /*
            * Verifies that, if no model name is written in input XML, then both parameters and
            * model options and written explicitly.
            */
            void CheckModelOptionsAndParams(const scnXml::Scenario& scenario)
            {
                // For each relevant element, determine whether it's specified explicitly in XML.
                const bool modelName = scenario.getModel().getModelName().present();
                const bool parameters = scenario.getModel().getParameters().present();
                const bool modelOptions = scenario.getModel().getModelOptions().present();
                if (!modelName && (!parameters || !modelOptions))
                {
                    throw util::xml_scenario_error(
                        "If a model name is not specified then both <ModelOptions> and <parameters> must be specified");
                }
            }
        };
    }
}

#endif
