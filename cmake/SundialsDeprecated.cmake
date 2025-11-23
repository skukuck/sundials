# ---------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2025-2026, Lawrence Livermore National Security,
# University of Maryland Baltimore County, and the SUNDIALS contributors.
# Copyright (c) 2013-2025, Lawrence Livermore National Security
# and Southern Methodist University.
# Copyright (c) 2002-2013, Lawrence Livermore National Security.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ---------------------------------------------------------------

# ---------------------------------------------------------------
# Option for testing CI deprecated options
# ---------------------------------------------------------------

if(DEFINED ENV{SUNDIALS_ENABLE_UNSET_DEPRECATED})
  set(_sundials_enable_unset_deprecated_default "$ENV{SUNDIALS_ENABLE_UNSET_DEPRECATED}")
else()
  set(_sundials_enable_unset_deprecated_default OFF)
endif()

sundials_option(SUNDIALS_ENABLE_UNSET_DEPRECATED BOOL
                "Unset deprecated CMake options" ${_sundials_enable_unset_deprecated_default} ADVANCED)

if(SUNDIALS_ENABLE_UNSET_DEPRECATED)
  message(WARNING "Unsetting deprecated SUNDIALS options.")
endif()

# ---------------------------------------------------------------
# Deprecated options that can not use the DEPRECATES_NAMES option
# to sundials_option
# ---------------------------------------------------------------

# Deprecated with SUNDIALS 6.4.0
if(DEFINED SUPERLUDIST_LIBRARY_DIR)
  message(DEPRECATION "The CMake option SUPERLUDIST_LIBRARY_DIR is deprecated. "
                      "Use SUPERLUDIST_DIR instead.")
  set(SUPERLUDIST_DIR
      "${SUPERLUDIST_LIBRARY_DIR}/../"
      CACHE BOOL "SuperLU_DIST root directory" FORCE)
  if(SUNDIALS_ENABLE_UNSET_DEPRECATED)
    unset(SUPERLUDIST_LIBRARY_DIR CACHE)
  endif()
endif()

# Deprecated CUDA_ARCH option
if(DEFINED CUDA_ARCH)
  message(DEPRECATION "The CMake option CUDA_ARCH is deprecated. "
                      "Use CMAKE_CUDA_ARCHITECTURES instead.")
  # convert sm_** to just **
  string(REGEX MATCH "[0-9]+" arch_name "${CUDA_ARCH}")
  set(CMAKE_CUDA_ARCHITECTURES
      ${arch_name}
      CACHE STRING "CUDA Architectures" FORCE)
  if(SUNDIALS_ENABLE_UNSET_DEPRECATED)
    unset(CUDA_ARCH)
  endif()
endif()

#
# Deprecated Testing Options
#

if(SUNDIALS_TEST_PROFILE)
  message(DEPRECATION "The CMake option SUNDIALS_TEST_PROFILE is deprecated. "
                      "Use SUNDIALS_TEST_ENABLE_PROFILING instead.")
  set(SUNDIALS_TEST_ENABLE_PROFILING
      ${SUNDIALS_TEST_PROFILE}
      CACHE BOOL "Profile tests" FORCE)
  unset(SUNDIALS_TEST_PROFILE)
endif()

if(SUNDIALS_TEST_NODIFF)
  message(DEPRECATION "The CMake option SUNDIALS_TEST_NODIFF is deprecated. "
                      "Use SUNDIALS_TEST_ENABLE_DIFF_OUTPUT instead.")
  if(SUNDIALS_TEST_NODIFF)
    set(_new_value OFF)
  else()
    set(_new_value ON)
  endif()
  set(SUNDIALS_TEST_ENABLE_DIFF_OUTPUT
      ${_new_value}
      CACHE BOOL "Compare test output with saved answer files" FORCE)
  if(SUNDIALS_ENABLE_UNSET_DEPRECATED)
    unset(SUNDIALS_TEST_NODIFF)
  endif()
  unset(_new_value)
endif()

if(SUNDIALS_TEST_DEVTESTS)
  message(DEPRECATION "The CMake option SUNDIALS_TEST_DEVTESTS is deprecated. "
                      "Use SUNDIALS_TEST_ENABLE_DEV_TESTS instead.")
  set(SUNDIALS_TEST_ENABLE_DEV_TESTS
      ${SUNDIALS_TEST_DEVTESTS}
      CACHE BOOL "Include development tests" FORCE)
  unset(SUNDIALS_TEST_DEVTESTS)
endif()

if(SUNDIALS_TEST_UNITTESTS)
  message(DEPRECATION "The CMake option SUNDIALS_TEST_UNITTESTS is deprecated. "
                      "Use SUNDIALS_TEST_ENABLE_UNIT_TESTS instead.")
  set(SUNDIALS_TEST_ENABLE_UNIT_TESTS
      ${SUNDIALS_TEST_UNITTESTS}
      CACHE BOOL "Include units tests" FORCE)
  unset(SUNDIALS_TEST_UNITTESTS)
endif()

if(SUNDIALS_CALIPER_OUTPUT_DIR)
  message(
    DEPRECATION
      "The CMake option SUNDIALS_CALIPER_OUTPUT_DIR is deprecated. "
      "Use SUNDIALS_TEST_CALIPER_OUTPUT_DIR and SUNDIALS_BENCHMARK_CALIPER_OUTPUT_DIR instead."
  )
  set(SUNDIALS_TEST_CALIPER_OUTPUT_DIR
      ${SUNDIALS_TEST_CALIPER_OUTPUT_DIR}
      CACHE PATH "Location to write test caliper files" FORCE)
  set(SUNDIALS_BENCHMARK_CALIPER_OUTPUT_DIR
      ${SUNDIALS_BENCHMARK_CALIPER_OUTPUT_DIR}
      CACHE PATH "Location to write benchmark caliper files" FORCE)
  if(SUNDIALS_ENABLE_UNSET_DEPRECATED)
    unset(SUNDIALS_CALIPER_OUTPUT_DIR)
  endif()
endif()
