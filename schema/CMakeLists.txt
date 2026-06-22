# CMake configuration for all code in this dir and below

# This file is part of OpenMalaria.
# 
# Copyright (C) 2005-2025 Swiss Tropical and Public Health Institute
# Copyright (C) 2005-2015 Liverpool School Of Tropical Medicine
# Copyright (C) 2020-2025 University of Basel
# Copyright (C) 2025 The Kids Research Institute Australia
#
# OpenMalaria is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or (at
# your option) any later version.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

set (SCHEMA_NAMES
  scenario
  demography
  monitoring
  interventions
  healthSystem
  entomology
  pharmacology
  vivax
  util
)
set (SCHEMA_CPP "")
set (SCHEMA_H "")
set (SCHEMA_XSD "")
foreach (XSD_NAME ${SCHEMA_NAMES})
  list (APPEND SCHEMA_CPP "${CMAKE_CURRENT_BINARY_DIR}/${XSD_NAME}.cpp")
  list (APPEND SCHEMA_H "${CMAKE_CURRENT_BINARY_DIR}/${XSD_NAME}.h")
  list (APPEND SCHEMA_XSD "${XSD_NAME}.xsd")
  
  set (XSD_FILE ${CMAKE_CURRENT_SOURCE_DIR}/${XSD_NAME}.xsd)
  add_custom_command (
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${XSD_NAME}.cpp ${CMAKE_CURRENT_BINARY_DIR}/${XSD_NAME}.h
    COMMAND ${XSD_EXECUTABLE} cxx-tree
        --std c++11
        --type-naming ucc --function-naming java
        --namespace-map http://openmalaria.org/schema/scenario_49=scnXml
#         --generate-serialization
        --generate-doxygen
        --generate-intellisense
        --hxx-suffix .h --cxx-suffix .cpp
        ${XSD_FILE}
    MAIN_DEPENDENCY ${XSD_FILE}
    COMMENT "Compiling ${XSD_FILE}"
    VERBATIM
  )
endforeach (XSD_NAME)

set( INLINED_XSD ${CMAKE_CURRENT_BINARY_DIR}/scenario_current.xsd )
add_custom_command (OUTPUT ${INLINED_XSD}
  DEPENDS ${SCHEMA_XSD}
  COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/util/inlineSchema.py
  scenario.xsd ${INLINED_XSD}
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  COMMENT "Inlining scenario.xsd and dependencies into ${INLINED_XSD}"
  VERBATIM
)
add_custom_target (inlined_xsd ALL DEPENDS ${INLINED_XSD})

include_directories (SYSTEM
  ${XSD_INCLUDE_DIRS}
  ${XERCESC_INCLUDE_DIRS}
)
include_directories (
  ${CMAKE_BINARY_DIR}
)

add_library (schema
  ${SCHEMA_CPP}
  ${SCHEMA_H}
)
