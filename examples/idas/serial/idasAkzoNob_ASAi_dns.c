/* -----------------------------------------------------------------
 * Programmer(s): Radu Serban and Cosmin Petra @ LLNL
 * -----------------------------------------------------------------
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
 * -----------------------------------------------------------------
 * Adjoint sensitivity example problem
 *
 * This IVP is a stiff system of 6 non-linear DAEs of index 1. The
 * problem originates from Akzo Nobel Central research in Arnhern,
 * The Netherlands, and describes a chemical process in which 2
 * species are mixed, while carbon dioxide is continuously added.
 * See http://pitagora.dm.uniba.it/~testset/report/chemakzo.pdf
 *
 * IDAS also computes the sensitivities of the integral
 *   G = int_t0^tf y1 dt
 * with respect to the initial values of the first components of y
 * (the differential components). These sensitivities are the first
 * five components of the solution of the adjoint system, at t = 0.
 * -----------------------------------------------------------------*/

#include <idas/idas.h> /* prototypes for IDA fcts., consts.    */
#include <math.h>
#include <nvector/nvector_serial.h> /* access to serial N_Vector            */
#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_math.h> /* def. of SUNRpowerI                   */
#include <sundials/sundials_types.h> /* defs. of sunrealtype, sunindextype      */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define ESYM "Le"
#define FSYM "Lf"
#else
#define ESYM "e"
#define FSYM "f"
#endif

/* Accessor macros */


/* Problem Constants */
#define NEQ 6
#define T0  SUN_RCONST(0.0)

#define TF SUN_RCONST(180.0) /* Final time. */

#define RTOL  SUN_RCONST(1.0e-08)
#define ATOL  SUN_RCONST(1.0e-10)
#define RTOLB SUN_RCONST(1.0e-06)
#define ATOLB SUN_RCONST(1.0e-08)
#define RTOLQ SUN_RCONST(1.0e-10)
#define ATOLQ SUN_RCONST(1.0e-12)

#define ZERO SUN_RCONST(0.0)
#define HALF SUN_RCONST(0.5)
#define ONE  SUN_RCONST(1.0)
#define TWO  SUN_RCONST(2.0)

#define STEPS 150

typedef struct
{
  sunrealtype k1, k2, k3, k4;
  sunrealtype K, klA, Ks, pCO2, H;
}* UserData;

static int res(sunrealtype t, N_Vector yy, N_Vector yd, N_Vector resval,
               void* userdata);

static int resB(sunrealtype tt, N_Vector yy, N_Vector yp, N_Vector yyB,
                N_Vector ypB, N_Vector rrB, void* user_dataB);

static int rhsQ(sunrealtype t, N_Vector yy, N_Vector yp, N_Vector qdot,
                void* user_data);

static void PrintOutput(sunrealtype tfinal, N_Vector yB, N_Vector ypB);
static int check_retval(void* returnvalue, const char* funcname, int opt);

/* Main program */
int main(void)
{
  UserData data;
  void* mem;
  N_Vector yy, yp, rr, q;
  N_Vector yB, ypB;
  int ncheck, retval;
  sunrealtype time;
  long int nst, nstB;
  int indexB;
  SUNMatrix A, AB;
  SUNLinearSolver LS, LSB;
  SUNContext ctx;

  /* Consistent IC for  y, y'. */
  const sunrealtype y00 = SUN_RCONST(0.444);
  const sunrealtype y01 = SUN_RCONST(0.00123);
  const sunrealtype y02 = SUN_RCONST(0.0);
  const sunrealtype y03 = SUN_RCONST(0.007);
  const sunrealtype y04 = SUN_RCONST(0.0);

  mem = NULL;
  yy = yp = NULL;
  A = AB = NULL;
  LS = LSB = NULL;

  printf("\nAdjoint Sensitivity Example for Akzo-Nobel Chemical Kinetics\n");
  printf("-------------------------------------------------------------\n");
  printf("Sensitivity of G = int_t0^tf (y1) dt with respect to IC.\n");
  printf("-------------------------------------------------------------\n\n");

  /* Create the SUNDIALS context object for this simulation */
  retval = SUNContext_Create(SUN_COMM_NULL, &ctx);
  if (check_retval(&retval, "SUNContext_Create", 1)) { return 1; }

  /* Allocate user data. */
  data = (UserData)malloc(sizeof(*data));

  /* Fill user's data with the appropriate values for coefficients. */
  data->k1   = SUN_RCONST(18.7);
  data->k2   = SUN_RCONST(0.58);
  data->k3   = SUN_RCONST(0.09);
  data->k4   = SUN_RCONST(0.42);
  data->K    = SUN_RCONST(34.4);
  data->klA  = SUN_RCONST(3.3);
  data->Ks   = SUN_RCONST(115.83);
  data->pCO2 = SUN_RCONST(0.9);
  data->H    = SUN_RCONST(737.0);

  /* Allocate N-vectors. */
  yy = N_VNew_Serial(NEQ, ctx);
  if (check_retval((void*)yy, "N_VNew_Serial", 0)) { return (1); }
  sunrealtype* yy_data = N_VGetArrayPointer(yy);
  yp = N_VClone(yy);
  if (check_retval((void*)yp, "N_VNew_Serial", 0)) { return (1); }

  /* Set IC */
  yy_data[0] = y00;
  yy_data[1] = y01;
  yy_data[2] = y02;
  yy_data[3] = y03;
  yy_data[4] = y04;
  yy_data[5] = data->Ks * y00 * y03;

  /* Get y' = - res(t0, y, 0) */
  N_VConst(ZERO, yp);

  rr = N_VClone(yy);
  res(T0, yy, yp, rr, data);
  N_VScale(-ONE, rr, yp);
  N_VDestroy(rr);

  /* Create and initialize q0 for quadratures. */
  q = N_VNew_Serial(1, ctx);
  if (check_retval((void*)q, "N_VNew_Serial", 0)) { return (1); }
  sunrealtype* q_data = N_VGetArrayPointer(q);
  q_data[0] = ZERO;

  /* Call IDACreate and IDAInit to initialize IDA memory */
  mem = IDACreate(ctx);
  if (check_retval((void*)mem, "IDACreate", 0)) { return (1); }

  retval = IDAInit(mem, res, T0, yy, yp);
  if (check_retval(&retval, "IDAInit", 1)) { return (1); }

  /* Set tolerances. */
  retval = IDASStolerances(mem, RTOL, ATOL);
  if (check_retval(&retval, "IDASStolerances", 1)) { return (1); }

  /* Attach user data. */
  retval = IDASetUserData(mem, data);
  if (check_retval(&retval, "IDASetUserData", 1)) { return (1); }

  /* Create dense SUNMatrix for use in linear solves */
  A = SUNDenseMatrix(NEQ, NEQ, ctx);
  if (check_retval((void*)A, "SUNDenseMatrix", 0)) { return (1); }

  /* Create dense SUNLinearSolver object */
  LS = SUNLinSol_Dense(yy, A, ctx);
  if (check_retval((void*)LS, "SUNLinSol_Dense", 0)) { return (1); }

  /* Attach the matrix and linear solver */
  retval = IDASetLinearSolver(mem, LS, A);
  if (check_retval(&retval, "IDASetLinearSolver", 1)) { return (1); }

  /* Initialize QUADRATURE(S). */
  retval = IDAQuadInit(mem, rhsQ, q);
  if (check_retval(&retval, "IDAQuadInit", 1)) { return (1); }

  /* Set tolerances and error control for quadratures. */
  retval = IDAQuadSStolerances(mem, RTOLQ, ATOLQ);
  if (check_retval(&retval, "IDAQuadSStolerances", 1)) { return (1); }

  retval = IDASetQuadErrCon(mem, SUNTRUE);
  if (check_retval(&retval, "IDASetQuadErrCon", 1)) { return (1); }

  /* Prepare ADJOINT. */
  retval = IDAAdjInit(mem, STEPS, IDA_HERMITE);
  if (check_retval(&retval, "IDAAdjInit", 1)) { return (1); }

  /* FORWARD run. */
  printf("Forward integration ... ");
  retval = IDASolveF(mem, TF, &time, yy, yp, IDA_NORMAL, &ncheck);
  if (check_retval(&retval, "IDASolveF", 1)) { return (1); }

  retval = IDAGetNumSteps(mem, &nst);
  if (check_retval(&retval, "IDAGetNumSteps", 1)) { return (1); }

  printf("done ( nst = %ld )\n", nst);

  retval = IDAGetQuad(mem, &time, q);
  if (check_retval(&retval, "IDAGetQuad", 1)) { return (1); }

  printf("G:          %24.16" FSYM " \n", q_data[0]);
  printf("--------------------------------------------------------\n\n");

  /* BACKWARD run */

  /* Initialize yB */
  yB = N_VClone(yy);
  if (check_retval((void*)yB, "N_VNew_Serial", 0)) { return (1); }
  N_VConst(ZERO, yB);

  ypB = N_VClone(yB);
  if (check_retval((void*)ypB, "N_VNew_Serial", 0)) { return (1); }
  sunrealtype* ypB_data = N_VGetArrayPointer(ypB);
  N_VConst(ZERO, ypB);
  ypB_data[0] = -ONE;

  retval = IDACreateB(mem, &indexB);
  if (check_retval(&retval, "IDACreateB", 1)) { return (1); }

  retval = IDAInitB(mem, indexB, resB, TF, yB, ypB);
  if (check_retval(&retval, "IDAInitB", 1)) { return (1); }

  retval = IDASStolerancesB(mem, indexB, RTOLB, ATOLB);
  if (check_retval(&retval, "IDASStolerancesB", 1)) { return (1); }

  retval = IDASetUserDataB(mem, indexB, data);
  if (check_retval(&retval, "IDASetUserDataB", 1)) { return (1); }

  retval = IDASetMaxNumStepsB(mem, indexB, 1000);

  /* Create dense SUNMatrix for use in linear solves */
  AB = SUNDenseMatrix(NEQ, NEQ, ctx);
  if (check_retval((void*)AB, "SUNDenseMatrix", 0)) { return (1); }

  /* Create dense SUNLinearSolver object */
  LSB = SUNLinSol_Dense(yB, AB, ctx);
  if (check_retval((void*)LSB, "SUNLinSol_Dense", 0)) { return (1); }

  /* Attach the matrix and linear solver */
  retval = IDASetLinearSolverB(mem, indexB, LSB, AB);
  if (check_retval(&retval, "IDASetLinearSolverB", 1)) { return (1); }

  printf("Backward integration ... ");

  retval = IDASolveB(mem, T0, IDA_NORMAL);
  if (check_retval(&retval, "IDASolveB", 1)) { return (1); }

  IDAGetNumSteps(IDAGetAdjIDABmem(mem, indexB), &nstB);
  printf("done ( nst = %ld )\n", nstB);

  retval = IDAGetB(mem, indexB, &time, yB, ypB);
  if (check_retval(&retval, "IDAGetB", 1)) { return (1); }

  PrintOutput(time, yB, ypB);

  IDAFree(&mem);
  SUNLinSolFree(LS);
  SUNMatDestroy(A);
  SUNLinSolFree(LSB);
  SUNMatDestroy(AB);
  N_VDestroy(yy);
  N_VDestroy(yp);
  N_VDestroy(yB);
  N_VDestroy(ypB);
  N_VDestroy(q);
  free(data);
  SUNContext_Free(&ctx);

  return (0);
}

static int res(sunrealtype t, N_Vector yy, N_Vector yd, N_Vector resval,
               void* userdata)
{
  UserData data;
  sunrealtype k1, k2, k3, k4;
  sunrealtype K, klA, Ks, pCO2, H;

  sunrealtype y0, y1, y2, y3, y4, y5;
  sunrealtype yd0, yd1, yd2, yd3, yd4;

  sunrealtype r1, r2, r3, r4, r5, Fin;
  sunrealtype* yy_data = N_VGetArrayPointer(yy);
  sunrealtype* yd_data = N_VGetArrayPointer(yd);
  sunrealtype* resval_data = N_VGetArrayPointer(resval);

  data = (UserData)userdata;
  k1   = data->k1;
  k2   = data->k2;
  k3   = data->k3;
  k4   = data->k4;
  K    = data->K;
  klA  = data->klA;
  Ks   = data->Ks;
  pCO2 = data->pCO2;
  H    = data->H;

  y0 = yy_data[0];
  y1 = yy_data[1];
  y2 = yy_data[2];
  y3 = yy_data[3];
  y4 = yy_data[4];
  y5 = yy_data[5];

  yd0 = yd_data[0];
  yd1 = yd_data[1];
  yd2 = yd_data[2];
  yd3 = yd_data[3];
  yd4 = yd_data[4];

  r1  = k1 * SUNRpowerI(y0, 4) * sqrt(y1);
  r2  = k2 * y2 * y3;
  r3  = k2 / K * y0 * y4;
  r4  = k3 * y0 * y3 * y3;
  r5  = k4 * y5 * y5 * sqrt(y1);
  Fin = klA * (pCO2 / H - y1);

  resval_data[0] = yd0 + TWO * r1 - r2 + r3 + r4;
  resval_data[1] = yd1 + HALF * r1 + r4 + HALF * r5 - Fin;
  resval_data[2] = yd2 - r1 + r2 - r3;
  resval_data[3] = yd3 + r2 - r3 + TWO * r4;
  resval_data[4] = yd4 - r2 + r3 - r5;
  resval_data[5] = Ks * y0 * y3 - y5;

  return (0);
}

/*
 * rhsQ routine. Computes quadrature(t,y).
 */

static int rhsQ(sunrealtype t, N_Vector yy, N_Vector yp, N_Vector qdot,
                void* user_data)
{
  sunrealtype* yy_data = N_VGetArrayPointer(yy);
  sunrealtype* qdot_data = N_VGetArrayPointer(qdot);
  qdot_data[0] = yy_data[0];

  return (0);
}

#define QUARTER SUN_RCONST(0.25)
#define FOUR    SUN_RCONST(4.0)
#define EIGHT   SUN_RCONST(8.0)

/*
 * resB routine. Residual for adjoint system.
 */
static int resB(sunrealtype tt, N_Vector yy, N_Vector yp, N_Vector yyB,
                N_Vector ypB, N_Vector rrB, void* user_dataB)
{
  UserData data;

  sunrealtype y0, y1, y2, y3, y4, y5;

  sunrealtype yB0, yB1, yB2, yB3, yB4, yB5;
  sunrealtype ypB0, ypB1, ypB2, ypB3, ypB4;

  sunrealtype k1, k2, k3, k4;
  sunrealtype K, klA, Ks;

  sunrealtype y2tohalf, y1to3, k2overK, tmp1, tmp2;
  sunrealtype* yyB_data = N_VGetArrayPointer(yyB);
  sunrealtype* yy_data = N_VGetArrayPointer(yy);
  sunrealtype* ypB_data = N_VGetArrayPointer(ypB);
  sunrealtype* rrB_data = N_VGetArrayPointer(rrB);

  data = (UserData)user_dataB;
  k1   = data->k1;
  k2   = data->k2;
  k3   = data->k3;
  k4   = data->k4;
  K    = data->K;
  klA  = data->klA;
  Ks   = data->Ks;

  y0 = yy_data[0];
  y1 = yy_data[1];
  y2 = yy_data[2];
  y3 = yy_data[3];
  y4 = yy_data[4];
  y5 = yy_data[5];

  yB0 = yyB_data[0];
  yB1 = yyB_data[1];
  yB2 = yyB_data[2];
  yB3 = yyB_data[3];
  yB4 = yyB_data[4];
  yB5 = yyB_data[5];

  ypB0 = ypB_data[0];
  ypB1 = ypB_data[1];
  ypB2 = ypB_data[2];
  ypB3 = ypB_data[3];
  ypB4 = ypB_data[4];

  y2tohalf = sqrt(y1);
  y1to3    = y0 * y0 * y0;
  k2overK  = k2 / K;

  tmp1        = k1 * y1to3 * y2tohalf;
  tmp2        = k3 * y3 * y3;
  rrB_data[0] = 1 + ypB0 - (EIGHT * tmp1 + k2overK * y4 + tmp2) * yB0 -
                (TWO * tmp1 + tmp2) * yB1 + (FOUR * tmp1 + k2overK * y4) * yB2 +
                k2overK * y4 * (yB3 - yB4) - TWO * tmp2 * yB3 + Ks * y3 * yB5;

  tmp1        = k1 * y0 * y1to3 * (y2tohalf / y1);
  tmp2        = k4 * y5 * y5 * (y2tohalf / y1);
  rrB_data[1] = ypB1 - tmp1 * yB0 -
                (QUARTER * tmp1 + QUARTER * tmp2 + klA) * yB1 +
                HALF * tmp1 * yB2 + HALF * tmp2 * yB4;

  rrB_data[2] = ypB2 + k2 * y3 * (yB0 - yB2 - yB3 + yB4);

  tmp1 = k3 * y0 * y3;
  tmp2 = k2 * y2;
  rrB_data[3] = ypB3 + (tmp2 - TWO * tmp1) * yB0 - TWO * tmp1 * yB1 - tmp2 * yB2 -
                (tmp2 + FOUR * tmp1) * yB3 + tmp2 * yB4 + Ks * y0 * yB5;

  rrB_data[4] = ypB4 - k2overK * y0 * (yB0 - yB2 - yB3 + yB4);

  rrB_data[5] = k4 * y5 * y2tohalf * (2 * yB4 - yB1) - yB5;

  return 0;
}

/*
 * Print results after backward integration
 */
static void PrintOutput(sunrealtype tfinal, N_Vector yB, N_Vector ypB)
{
  sunrealtype* yB_data = N_VGetArrayPointer(yB);
  printf("dG/dy0: \t%12.4" ESYM "\n\t\t%12.4" ESYM "\n\t\t%12.4" ESYM
         "\n\t\t%12.4" ESYM "\n\t\t%12.4" ESYM "\n",
         yB_data[0], yB_data[1], yB_data[2], yB_data[3], yB_data[4]);
  printf("--------------------------------------------------------\n\n");
}

/*
 * Check function return value.
 *    opt == 0 means SUNDIALS function allocates memory so check if
 *             returned NULL pointer
 *    opt == 1 means SUNDIALS function returns an integer value so check if
 *             retval < 0
 *    opt == 2 means function allocates memory so check if returned
 *             NULL pointer
 */

static int check_retval(void* returnvalue, const char* funcname, int opt)
{
  int* retval;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && returnvalue == NULL)
  {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return (1);
  }

  /* Check if retval < 0 */
  else if (opt == 1)
  {
    retval = (int*)returnvalue;
    if (*retval < 0)
    {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with retval = %d\n\n",
              funcname, *retval);
      return (1);
    }
  }

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && returnvalue == NULL)
  {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return (1);
  }

  return (0);
}
