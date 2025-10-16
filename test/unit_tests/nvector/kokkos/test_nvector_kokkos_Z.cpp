/* -----------------------------------------------------------------------------
 * Programmer(s): Daniel McGreer and Cody J. Balos @ LLNL
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2025, Lawrence Livermore National Security,
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
 * This is the testing routine for the NVector implementation using Kokkos.
 * ---------------------------------------------------------------------------*/

#include <nvector/nvector_kokkos.hpp>
#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>

#include "test_nvector_complex.h"

#if defined(USE_CUDA)
using ExecSpace = Kokkos::Cuda;
#elif defined(USE_HIP)
#if KOKKOS_VERSION / 10000 > 3
using ExecSpace = Kokkos::HIP;
#else
using ExecSpace = Kokkos::Experimental::HIP;
#endif
#elif defined(USE_OPENMP)
using ExecSpace = Kokkos::OpenMP;
#else
using ExecSpace = Kokkos::Serial;
#endif

using VecType  = sundials::kokkos::Vector<ExecSpace>;
using SizeType = VecType::size_type;

/* ----------------------------------------------------------------------
 * Main NVector Testing Routine
 * --------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
  int fails{0};        /* counter for test failures */
  sunindextype length; /* vector length             */
  int print_timing;    /* turn timing on/off        */

  Test_Init_Z(SUN_COMM_NULL);

  /* check input and set vector length */
  if (argc < 3)
  {
    printf("ERROR: TWO (2) Inputs required: vector length, print timing \n");
    return (-1);
  }

  length = (sunindextype)atol(argv[1]);
  if (length <= 0)
  {
    printf("ERROR: length of vector must be a positive integer \n");
    return (-1);
  }

  print_timing = atoi(argv[2]);
  SetTiming_Z(print_timing, 0);

  printf("Testing KOKKOS N_Vector \n");

  printf("Vector length %ld \n\n", (long int)length);

  Kokkos::initialize(argc, argv);
  {
    VecType X{static_cast<SizeType>(length), sunctx};

    /* Check vector ID */
    fails += Test_N_VGetVectorID_Z(X, SUNDIALS_NVEC_KOKKOS, 0);

    /* Test clone functions */
    fails += Test_N_VClone_Z(X, length, 0);
    fails += Test_N_VCloneVectorArray_Z(5, X, length, 0);

    /* Check vector length */
    fails += Test_N_VGetLength_Z(X, 0);

    /* Check vector communicator */
    fails += Test_N_VGetCommunicator_Z(X, SUN_COMM_NULL, 0);

    /* Clone additional vectors for testing */
    VecType Y{X};
    VecType Z{X};

    /* Standard vector operation tests */
    printf("\nTesting standard vector operations:\n\n");

    fails += Test_N_VAbs_Z(X, Z, length, 0);
    fails += Test_N_VAddConst_Z(X, Z, length, 0);
    fails += Test_N_VCompare_Z(X, Z, length, 0);
    fails += Test_N_VConst_Z(X, length, 0);
    fails += Test_N_VConstrMask_Z(X, Y, Z, length, 0);
    fails += Test_N_VDiv_Z(X, Y, Z, length, 0);
    fails += Test_N_VDotProd_Z(X, Y, length, 0);
    fails += Test_N_VInv_Z(X, Z, length, 0);
    fails += Test_N_VInvTest_Z(X, Z, length, 0);
    fails += Test_N_VL1Norm_Z(X, length, 0);
    fails += Test_N_VLinearSum_Z(X, Y, Z, length, 0);
    fails += Test_N_VMaxNorm_Z(X, length, 0);
    fails += Test_N_VMin_Z(X, length, 0);
    fails += Test_N_VMinQuotient_Z(X, Y, length, 0);
    fails += Test_N_VProd_Z(X, Y, Z, length, 0);
    fails += Test_N_VScale_Z(X, Z, length, 0);
    fails += Test_N_VWL2Norm_Z(X, Y, length, 0);
    fails += Test_N_VWrmsNorm_Z(X, Y, length, 0);
    fails += Test_N_VWrmsNormMask_Z(X, Y, Z, length, 0);

    /* Fused and vector array operations tests (disabled) */
    printf("\nTesting fused and vector array operations (disabled):\n\n");

    /* create vector and test vector array operations */
    VecType U{X};

    /* fused operations */
    fails += Test_N_VLinearCombination_Z(U, length, 0);
    fails += Test_N_VScaleAddMulti_Z(U, length, 0);
    fails += Test_N_VDotProdMulti_Z(U, length, 0);

    /* vector array operations */
    fails += Test_N_VLinearSumVectorArray_Z(U, length, 0);
    fails += Test_N_VScaleVectorArray_Z(U, length, 0);
    fails += Test_N_VConstVectorArray_Z(U, length, 0);
    fails += Test_N_VWrmsNormVectorArray_Z(U, length, 0);
    fails += Test_N_VWrmsNormMaskVectorArray_Z(U, length, 0);
    fails += Test_N_VScaleAddMultiVectorArray_Z(U, length, 0);
    fails += Test_N_VLinearCombinationVectorArray_Z(U, length, 0);

    /* local reduction operations */
    printf("\nTesting local reduction operations:\n\n");

    fails += Test_N_VDotProdLocal_Z(X, Y, length, 0);
    fails += Test_N_VMaxNormLocal_Z(X, length, 0);
    fails += Test_N_VMinLocal_Z(X, length, 0);
    fails += Test_N_VL1NormLocal_Z(X, length, 0);
    fails += Test_N_VWSqrSumLocal_Z(X, Y, length, 0);
    fails += Test_N_VWSqrSumMaskLocal_Z(X, Y, Z, length, 0);
    fails += Test_N_VInvTestLocal_Z(X, Z, length, 0);
    fails += Test_N_VConstrMaskLocal_Z(X, Y, Z, length, 0);
    fails += Test_N_VMinQuotientLocal_Z(X, Y, length, 0);
  }
  Kokkos::finalize();

  /* Print result */
  if (fails) { printf("FAIL: NVector module failed %i tests \n\n", fails); }
  else { printf("SUCCESS: NVector module passed all tests \n\n"); }

  Test_Finalize_Z();

  return (fails);
}

/* ----------------------------------------------------------------------
 * Implementation specific utility functions for vector tests
 * --------------------------------------------------------------------*/

int check_ans(sunrealtype ans, N_Vector X, sunindextype local_length)
{
  return check_ans_Z(ans, X, local_length);
}

int check_ans_Z(sunscalartype ans, N_Vector X, sunindextype local_length)
{
  int failure{0};
  auto Xvec{static_cast<VecType*>(X->content)};
  auto Xdata{Xvec->HostView()};

  sundials::kokkos::CopyFromDevice<VecType>(*Xvec);
  for (sunindextype i = 0; i < local_length; i++)
  {
    failure += SUNCompare(Xdata[i], ans);
  }

  return (failure > ZERO) ? (1) : (0);
}

sunbooleantype has_data_Z(N_Vector X)
{
  /* check if vector data is non-null */
  return SUNTRUE;
}

void set_element_Z(N_Vector X, sunindextype i, sunscalartype val)
{
  /* set i-th element of data array */
  set_element_range_Z(X, i, i, val);
}

void set_element_range_Z(N_Vector X, sunindextype is, sunindextype ie,
                         sunscalartype val)
{
  auto Xvec{static_cast<VecType*>(X->content)};
  auto Xdata{Xvec->HostView()};

  /* set elements [is,ie] of the data array */
  sundials::kokkos::CopyFromDevice<VecType>(X);
  for (sunindextype i = is; i <= ie; i++) { Xdata[i] = val; }
  sundials::kokkos::CopyToDevice<VecType>(X);
}

sunscalartype get_element_Z(N_Vector X, sunindextype i)
{
  /* get i-th element of data array */
  auto Xvec{static_cast<VecType*>(X->content)};
  auto Xdata{Xvec->HostView()};
  sundials::kokkos::CopyFromDevice<VecType>(X);
  return Xdata[i];
}

double max_time_Z(N_Vector X, double time)
{
  /* not running in parallel, just return input time */
  return time;
}

void sync_device_Z(N_Vector x)
{
  /* sync with GPU */
  Kokkos::fence();
  return;
}
