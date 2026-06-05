# CMake tests configuration for OpenMalaria
# Copyright © 2005-2009 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
# Licence: GNU General Public Licence version 2 or later (see COPYING)

# Testing via run.sh or run.py script:
configure_file (
  ${CMAKE_CURRENT_SOURCE_DIR}/run.py
  ${CMAKE_CURRENT_BINARY_DIR}/run.py
  @ONLY
)

# working tests (with checkpointing):
set (OM_BOXTEST_NAMES
  #disabled because of IPT:
  # 2 3 IPT VT
  1 4 5 6 9 10 11 12
  2ITNs
  Cohort
  EffectiveDrug
  Empirical
  # ESCMTest (cannot be updated to new format due to dropped features and weirdness)
  ESTS
  Genotypes
  IRS30
  KK_20150612
  LifeNet1
  LifeNet2
  Molineaux
  MolPairwise
  MSAT
  MSAT_HRP2
  # disabled due to test failure (TODO: since when?)
  # MultiVecPopInterventions
  NamawalaArabiensis
  NoInterv
  NoMPDLarviciding
  Penny
  # model deactivated: R_0
  Rach5IC
  SimpleMPDLarviciding
  SimpleMPDTest
  SubPopRemoval
  TrapTest
  TriggeredMSAT
  VecFullTest
  # disabled because the vector lifecycle model isn't yet done:
  # VecLCMonthly VecLCTest
  VecMonthly
  VecTest
  Vivax
  NonHumanHosts
  PEV
  BSV
  TBV
  DecisionTree5DayDielmo
  HetVecDailyEIR
  AvailabilityFilter
  GaussianCopulaGamma
  GaussianCopulaLognormal
  InfectionOrigin
  ModelNameNoOverrides
  ModelNameParamOverrides
  ModelNameModelOptionOverrides
  ModelNameManyOverrides
  ExactNv0
)
# tests with broken checkpointing:
set (OM_BOXTEST_NC_NAMES)
# Disabled due to "in-progress" work: (none)

# Use python script.
foreach (TEST_NAME ${OM_BOXTEST_NAMES})
    add_test (${TEST_NAME} ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/run.py ${TEST_NAME} -- --checkpoint-stop)
endforeach (TEST_NAME)
foreach (TEST_NAME ${OM_BOXTEST_NC_NAMES})
    add_test (${TEST_NAME} ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/run.py -- ${TEST_NAME})
endforeach (TEST_NAME)
