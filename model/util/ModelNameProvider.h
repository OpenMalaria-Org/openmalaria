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

#ifndef Hmod_util_ModelNameProvider
#define Hmod_util_ModelNameProvider

#include <schema/scenario.h>
#include "util/errors.h"

namespace OM
{
    namespace util
    {
        /*
        * It is intended that this holds the entire collection of named models.
        */
        enum class ModelNames
        {
            none, // Represents the case where input XML contains no model name.
            base
        };

        /*
        * Stores model name (if any) specified in XML at initialization time
        * then exposes it as an enum class.  This is exists so that clients
        * don't have to hardcode strings referring (or possibly, erroneously,
        * not referring) to model names.
        *
        * If new named models are added, support for them should be added to
        * this class.
        */
        class ModelNameProvider
        {
        public:
            ModelNameProvider(const scnXml::Model& model)
            {
                const bool useNamedModel = model.getModelName().present();
                if (useNamedModel)
                {
                    const std::string nameFromXML = model.getModelName().get().getName();
                    if (nameFromXML == "base")
                    {
                        modelInUse = ModelNames::base;
                    }
                    else
                    {
                        throw util::xml_scenario_error("Unrecognized model name: " + nameFromXML);
                    }
                }
            }

            ~ModelNameProvider(){}

            ModelNames GetModelName() { return modelInUse; }

        private:
            ModelNames modelInUse = ModelNames::none;
        };
    }
}

#endif
