/* -----------------------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ UMBC
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2025-2026, Lawrence Livermore National Security,
 * University of Maryland Baltimore County, and the SUNDIALS contributors.
 * Copyright (c) 2013-2025, Lawrence Livermore National Security
 * and Southern Methodist University.
 * Copyright (c) 2002-2013, Lawrence Livermore National Security.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
 * Utility routines for various "ark_test_stageinfo" unit tests.
 * ---------------------------------------------------------------------------*/

#ifndef STAGEINFO_HPP_
#define STAGEINFO_HPP_

#include "arkode/arkode.h"
#include "sundials/sundials_math.h"

int preprocess_step(sunrealtype t, N_Vector y, void* user_data)
{
  sunrealtype tn, tcur;
  if (ARKodeGetLastTime(arkode_mem, &tn) != ARK_SUCCESS)
  {
    std::cerr << "Error in ARKodeGetLastTime" << std::endl;
    return -1;
  }
  if (ARKodeGetCurrentTime(arkode_mem, &tcur) != ARK_SUCCESS)
  {
    std::cerr << "Error in ARKodeGetCurrentTime" << std::endl;
    return -1;
  }
  std::cout << "    [Pre-step processing at t = " << std::setprecision(2) << t
            << " (tn = " << tn << " , tcur = " << tcur << "),"
            << std::setprecision(10)
            << "||y||_2 = " << SUNRsqrt(N_VDotProd(y, y)) << "]" << std::endl
            << std::flush;
  return 0;
}

int postprocess_step(sunrealtype t, N_Vector y, void* user_data)
{
  sunrealtype tn, tcur;
  if (ARKodeGetLastTime(arkode_mem, &tn) != ARK_SUCCESS)
  {
    std::cerr << "Error in ARKodeGetLastTime" << std::endl;
    return -1;
  }
  if (ARKodeGetCurrentTime(arkode_mem, &tcur) != ARK_SUCCESS)
  {
    std::cerr << "Error in ARKodeGetCurrentTime" << std::endl;
    return -1;
  }
  std::cout << "    [Post-step processing at t = " << std::setprecision(2) << t
            << " (tn = " << tn << " , tcur = " << tcur << "),"
            << std::setprecision(10)
            << "||y||_2 = " << SUNRsqrt(N_VDotProd(y, y)) << "]" << std::endl
            << std::flush;
  return 0;
}

int postprocess_step_fail(sunrealtype t, N_Vector y, void* user_data)
{
  sunrealtype tn, tcur;
  if (ARKodeGetLastTime(arkode_mem, &tn) != ARK_SUCCESS)
  {
    std::cerr << "Error in ARKodeGetLastTime" << std::endl;
    return -1;
  }
  if (ARKodeGetCurrentTime(arkode_mem, &tcur) != ARK_SUCCESS)
  {
    std::cerr << "Error in ARKodeGetCurrentTime" << std::endl;
    return -1;
  }
  std::cout << "    [Post-step failure processing at t = "
            << std::setprecision(2) << t << " (tn = " << tn
            << " , tcur = " << tcur << ")," << std::setprecision(10)
            << "||y||_2 = " << SUNRsqrt(N_VDotProd(y, y)) << "]" << std::endl
            << std::flush;
  return 0;
}

int preprocess_stage(sunrealtype t, N_Vector y, void* user_data)
{
  int stage, max_stages;
  sunrealtype tn, tcur;
  if (ARKodeGetLastTime(arkode_mem, &tn) != ARK_SUCCESS)
  {
    std::cerr << "Error in ARKodeGetLastTime" << std::endl;
    return -1;
  }
  if (ARKodeGetCurrentTime(arkode_mem, &tcur) != ARK_SUCCESS)
  {
    std::cerr << "Error in ARKodeGetCurrentTime" << std::endl;
    return -1;
  }
  if (ARKodeGetStageIndex(arkode_mem, &stage, &max_stages) != ARK_SUCCESS)
  {
    std::cerr << "Error in ARKodeGetStageIndex" << std::endl;
    return -1;
  }
  std::cout << "    [Pre-RHS processing (stage " << stage << " of " << max_stages
            << ") at t = " << std::setprecision(2) << t << " (tn = " << tn
            << " , tcur = " << tcur << "), " << std::setprecision(10)
            << "||y||_2 = " << SUNRsqrt(N_VDotProd(y, y)) << "]" << std::endl
            << std::flush;
  return 0;
}

int postprocess_stage(sunrealtype t, N_Vector y, void* user_data)
{
  int stage, max_stages;
  sunrealtype tn, tcur;
  if (ARKodeGetLastTime(arkode_mem, &tn) != ARK_SUCCESS)
  {
    std::cerr << "Error in ARKodeGetLastTime" << std::endl;
    return -1;
  }
  if (ARKodeGetCurrentTime(arkode_mem, &tcur) != ARK_SUCCESS)
  {
    std::cerr << "Error in ARKodeGetCurrentTime" << std::endl;
    return -1;
  }
  if (ARKodeGetStageIndex(arkode_mem, &stage, &max_stages) != ARK_SUCCESS)
  {
    std::cerr << "Error in ARKodeGetStageIndex" << std::endl;
    return -1;
  }
  std::cout << "    [Post-stage processing (stage " << stage << " of "
            << max_stages << ") at t = " << std::setprecision(2) << t
            << " (tn = " << tn << " , tcur = " << tcur << "), "
            << std::setprecision(10)
            << "||y||_2 = " << SUNRsqrt(N_VDotProd(y, y)) << "]" << std::endl
            << std::flush;
  return 0;
}

#endif // STAGEINFO_HPP_
