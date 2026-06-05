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

# Configure this file
file(STRINGS ${CMAKE_SOURCE_DIR}/version.txt OM_VERSION)
# make sure version.h gets rebuilt everytime cmake is run
unset(${CMAKE_CURRENT_BINARY_DIR}/util/version.h CACHE)
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/util/version.h.in
    ${CMAKE_CURRENT_BINARY_DIR}/util/version.h
    @ONLY
)

# Don't use aux_source_directory on . because we don't want to compile openMalaria.cpp in to the lib.
set (Model_CPP
  sim.cpp
  Population.cpp
  PopulationAgeStructure.cpp
  checkpoint.cpp
  
  Host/Human.cpp
  Host/InfectionIncidenceModel.cpp
  Host/NeonatalMortality.cpp
  
  Clinical/ClinicalModel.cpp
  Clinical/EventScheduler.cpp
  Clinical/ESCaseManagement.cpp
  Clinical/CMDecisionTree.cpp
  Clinical/ImmediateOutcomes.cpp
  Clinical/Episode.cpp
  Clinical/DecisionTree5Day.cpp
  Clinical/CM5DayCommon.cpp
  
  PkPd/LSTMModel.cpp
  PkPd/Drug/LSTMDrug.cpp
  PkPd/Drug/LSTMDrugOneComp.cpp
  PkPd/Drug/LSTMDrugThreeComp.cpp
  PkPd/Drug/LSTMDrugConversion.cpp
  PkPd/Drug/LSTMDrugType.cpp
  PkPd/LSTMTreatments.cpp
  
  Transmission/VectorModel.cpp
  Transmission/PerHost.cpp
  Transmission/Anopheles/AnophelesModel.cpp
  Transmission/Anopheles/CalcNvO.cpp
  
  Host/WithinHost/WHInterface.cpp
  Host/WithinHost/WHFalciparum.cpp
  Host/WithinHost/WHVivax.cpp
  Host/WithinHost/CommonWithinHost.cpp
  Host/WithinHost/DescriptiveWithinHost.cpp
  Host/WithinHost/Genotypes.cpp
  Host/WithinHost/Infection/CommonInfection.cpp
  Host/WithinHost/Infection/DescriptiveInfection.cpp
  Host/WithinHost/Infection/DummyInfection.cpp
  Host/WithinHost/Infection/EmpiricalInfection.cpp
  Host/WithinHost/Infection/MolineauxInfection.cpp
  Host/WithinHost/Infection/PennyInfection.cpp
  Host/WithinHost/Pathogenesis/PathogenesisModel.cpp
  Host/WithinHost/Pathogenesis/Submodels.cpp
  Host/WithinHost/Diagnostic.cpp
  Host/WithinHost/Treatments.cpp
  
  mon/mon.cpp
  mon/misc.cpp
  mon/Continuous.cpp
  
  util/errors.cpp
  util/checkpoint.cpp
  util/ModelOptions.cpp
  util/CommandLine.cpp
  util/AgeGroupInterpolation.cpp
  util/sampler.cpp
  util/SpeciesIndexChecker.cpp
  util/DocumentLoader.cpp
  util/random.cpp
  util/UnitParse.cpp
  util/DecayFunction.cpp
  
  interventions/InterventionManager.cpp
  interventions/ITN.cpp
  interventions/IRS.cpp
  interventions/GVI.cpp
  interventions/Vaccine.cpp
)

if (${OM_STREAM_VALIDATOR})
  list(APPEND Model_CPP util/StreamValidator.cpp)
endif (${OM_STREAM_VALIDATOR})
# Headers - only included so they show up in IDEs:
# This misses loads of headers. Fix if you care.
file (GLOB_RECURSE Model_H "${CMAKE_SOURCE_DIR}/model/*.h")


include_directories (SYSTEM
  ${XSD_INCLUDE_DIRS}
  ${XERCESC_INCLUDE_DIRS}
  ${GSL_INCLUDE_DIRS}
  ${Z_INCLUDE_DIRS}
  # ${SPIRIT_INCLUDE_DIRS}
  
  ${CMAKE_BINARY_DIR}/model
  ${CMAKE_SOURCE_DIR}/contrib
  # ${CMAKE_SOURCE_DIR}/contrib/floating_point_utilities_v3
)
include_directories (
  ${CMAKE_SOURCE_DIR}/model ${CMAKE_BINARY_DIR}
)

add_library (model
  ${Model_CPP}
  ${Model_H}
  ${CMAKE_CURRENT_BINARY_DIR}/util/version.h
)

if (MSVC)
  set_target_properties (model PROPERTIES
    LINK_FLAGS "${OM_LINK_FLAGS}"
    COMPILE_FLAGS "${OM_COMPILE_FLAGS}"
  )
endif (MSVC)
