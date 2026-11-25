# CMake configuration for openmalaria's unit tests
# Copyright © 2005-2009 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
# Licence: GNU General Public Licence version 2 or later (see COPYING)

set (OM_CXXTEST_INTERPRETER ${PYTHON_EXECUTABLE})
set (OM_CXXTEST_SCRIPT ${CMAKE_SOURCE_DIR}/contrib/cxxtest/bin/cxxtestgen)

# We want abort-on-fail for safety checks (like bad memory accesses) to work
set (OM_CXXTEST_OPTIONS --have-std --have-eh --abort-on-fail)


# For yaml-cpp: I did use this to read in data for the AnophelesModel unittest.
# It's another not-the-easiest dependency to install, so maybe don't use it?
#find_library (YAML_CPP_LIB NAMES yaml-cpp)
#if (NOT YAML_CPP_LIB)
#  message (FATAL_ERROR "Unable to find yaml-cpp library. (You can avoid this dependency by disabling OM_CXXTEST_ENABLE)")
#endif (NOT YAML_CPP_LIB)

set (CMAKE_CONFIGURED 1)
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/TestPaths.h ${CMAKE_CURRENT_BINARY_DIR}/configured/TestPaths.h @ONLY)
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/scenario.xml ${CMAKE_CURRENT_BINARY_DIR}/configured/scenario.xml @ONLY)
include_directories( SYSTEM
    ${CMAKE_SOURCE_DIR}/contrib/cxxtest
    ${CMAKE_SOURCE_DIR}/contrib
)
include_directories (
  ${CMAKE_SOURCE_DIR}/model
  ${CMAKE_SOURCE_DIR}/unittest
  ${CMAKE_BINARY_DIR}/unittest
)

set (OM_CXXTEST_HEADERS
  #AnophelesModelSuite.h needs updating for vector model changes
  ExtraAsserts.h	# must appear after at least some of the above
  LSTMPkPdSuite.h
  CheckpointSuite.h
  DummyInfectionSuite.h
  EmpiricalInfectionSuite.h
  InfectionImmunitySuite.h
  #DescriptiveInfectionSuite.h
  #ESCaseManagementSuite.h
  CMDecisionTreeSuite.h
  AgeGroupInterpolationSuite.h
  DecayFunctionSuite.h
  SamplerSuite.h
  PennyInfectionSuite.h
  MolineauxInfectionSuite.h
  #MosqLifeCycleSuite.h
  UtilVectorsSuite.h
  PkPdComplianceSuite.h
  ChaChaSuite.h
  XoshiroSuite.h
  ClinicalAccessSuite.h
)

add_custom_command (OUTPUT tests.cpp
    COMMAND ${PYTHON_EXECUTABLE} ${OM_CXXTEST_SCRIPT} ${OM_CXXTEST_OPTIONS} --runner=ParenPrinter -o ${CMAKE_CURRENT_BINARY_DIR}/tests.cpp ${OM_CXXTEST_HEADERS}
    DEPENDS ${OM_CXXTEST_HEADERS}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Generating unittest code with cxxtestgen"
    VERBATIM
)

# These files are needed to run various infection model tests:
configure_file (${CMAKE_SOURCE_DIR}/test/autoRegressionParameters.csv ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/PennyCirDens.txt ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/PennySeqDens.txt ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/MolineauxCirDens.txt ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/MolineauxStatsPairwise ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/MolineauxStatsPairwiseRG ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/MolineauxStatsOrig ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/MolineauxStatsOrigRG ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
# configure_file (${CMAKE_CURRENT_SOURCE_DIR}/MolineauxStats1MG ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
# configure_file (${CMAKE_CURRENT_SOURCE_DIR}/MolineauxStats1MGRG ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
# configure_file (${CMAKE_CURRENT_SOURCE_DIR}/MolineauxStatsMDG ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
# configure_file (${CMAKE_CURRENT_SOURCE_DIR}/MolineauxStatsMDGRG ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
# configure_file (${CMAKE_CURRENT_SOURCE_DIR}/MolineauxStats1MGMDG ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
# configure_file (${CMAKE_CURRENT_SOURCE_DIR}/MolineauxStats1MGMDGRG ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
# and these by MosqLifeCycleSuite:
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/MLCS_EmergenceRateConstantRA.txt ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/MLCS_AdultMosqConstantRA.txt ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/MLCS_OvipositingAdultsConstantRA.txt ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/MLCS_EmergenceRatePeriodicRA.txt ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/MLCS_AdultMosqPeriodicRA.txt ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/MLCS_OvipositingAdultsPeriodicRA.txt ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/MLCS_PeriodicResourceAvailability.txt ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

add_executable (unittest
  WHMock.cpp
  tests.cpp
  ${OM_CXXTEST_HEADERS} # for IDEs
)
target_link_libraries (unittest
  model
  schema
  contrib
  ${OM_CXXTEST_GUI_LIB}
  ${GSL_LIBRARIES}
  ${XERCESC_LIBRARIES}
  ${Z_LIBRARIES}
  ${PTHREAD_LIBRARIES}
  ${OM_STD_LIBS}
)

if (MSVC)
  set_target_properties (unittest PROPERTIES
    LINK_FLAGS "${OM_LINK_FLAGS}"
    COMPILE_FLAGS "${OM_COMPILE_FLAGS}"
  )
endif (MSVC)

add_test (unittest unittest)

mark_as_advanced (
  OM_CXXTEST_OPTIONS
  OM_CXXTEST_GUI_LIB
)
