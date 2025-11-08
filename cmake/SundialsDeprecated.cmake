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

#
# Deprecated TPL options
#

if(DEFINED SUPERLUDIST_ENABLE)
  message(DEPRECATION "The CMake option SUPERLUDIST_ENABLE is deprecated. "
                      "Use ENABLE_SUPERLUDIST instead.")
  set(ENABLE_SUPERLUDIST
      ${SUPERLUDIST_ENABLE}
      CACHE BOOL "Enable SuperLU_DIST support" FORCE)
  unset(SUPERLUDIST_ENABLE CACHE)
endif()

# Deprecated with SUNDIALS 6.4.0
if(DEFINED SUPERLUDIST_LIBRARY_DIR)
  message(DEPRECATION "The CMake option SUPERLUDIST_LIBRARY_DIR is deprecated. "
                      "Use SUPERLUDIST_DIR instead.")
  set(SUPERLUDIST_DIR
      "${SUPERLUDIST_LIBRARY_DIR}/../"
      CACHE BOOL "SuperLU_DIST root directory" FORCE)
  unset(SUPERLUDIST_LIBRARY_DIR CACHE)
endif()
if(DEFINED SUPERLUDIST_INCLUDE_DIR)
  message(DEPRECATION "The CMake option SUPERLUDIST_INCLUDE_DIR is deprecated. "
                      "Use SUPERLUDIST_INCLUDE_DIRS instead.")
  set(SUPERLUDIST_INCLUDE_DIRS
      "${SUPERLUDIST_INCLUDE_DIR}"
      CACHE BOOL "SuperLU_DIST include directoroes" FORCE)
  unset(SUPERLUDIST_INCLUDE_DIR CACHE)
endif()

if(DEFINED SUPERLUMT_ENABLE)
  message(DEPRECATION "The CMake option SUPERLUMT_ENABLE is deprecated. "
                      "Use ENABLE_SUPERLUMT instead.")
  set(ENABLE_SUPERLUMT
      ${SUPERLUMT_ENABLE}
      CACHE BOOL "Enable SuperLU_MT support" FORCE)
  unset(SUPERLUMT_ENABLE CACHE)
endif()

if(DEFINED KLU_ENABLE)
  message(DEPRECATION "The CMake option KLU_ENABLE is deprecated. "
                      "Use ENABLE_KLU instead.")
  set(ENABLE_KLU
      ${KLU_ENABLE}
      CACHE BOOL "Enable KLU support" FORCE)
  unset(KLU_ENABLE CACHE)
endif()

if(DEFINED HYPRE_ENABLE)
  message(DEPRECATION "The CMake option HYPRE_ENABLE is deprecated. "
                      "Use ENABLE_HYPRE instead.")
  set(ENABLE_HYPRE
      ${HYPRE_ENABLE}
      CACHE BOOL "Enable HYPRE support" FORCE)
  unset(HYPRE_ENABLE CACHE)
endif()

if(DEFINED PETSC_ENABLE)
  message(DEPRECATION "The CMake option PETSC_ENABLE is deprecated. "
                      "Use ENABLE_PETSC instead.")
  set(ENABLE_PETSC
      ${PETSC_ENABLE}
      CACHE BOOL "Enable PETSC support" FORCE)
  unset(PETSC_ENABLE CACHE)
endif()

if(DEFINED Trilinos_ENABLE)
  message(DEPRECATION "The CMake option Trilinos_ENABLE is deprecated. "
                      "Use ENABLE_TRILINOS instead.")
  set(ENABLE_TRILINOS
      ${Trilinos_ENABLE}
      CACHE BOOL "Enable Trilinos support" FORCE)
  unset(Trilinos_ENABLE CACHE)
endif()

if(DEFINED RAJA_ENABLE)
  message(DEPRECATION "The CMake option RAJA_ENABLE is deprecated. "
                      "Use ENABLE_RAJA instead.")
  set(ENABLE_RAJA
      ${RAJA_ENABLE}
      CACHE BOOL "Enable RAJA support" FORCE)
  unset(RAJA_ENABLE CACHE)
endif()

#
# Deprecated CUDA_ARCH option
#

if(DEFINED CUDA_ARCH)
  message(DEPRECATION "The CMake option CUDA_ARCH is deprecated. "
                      "Use CMAKE_CUDA_ARCHITECTURES instead.")
  # convert sm_** to just **
  string(REGEX MATCH "[0-9]+" arch_name "${CUDA_ARCH}")
  set(CMAKE_CUDA_ARCHITECTURES
      ${arch_name}
      CACHE STRING "CUDA Architectures" FORCE)
  unset(CUDA_ARCH)
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
  unset(SUNDIALS_TEST_NODIFF)
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
  unset(SUNDIALS_CALIPER_OUTPUT_DIR)
endif()
