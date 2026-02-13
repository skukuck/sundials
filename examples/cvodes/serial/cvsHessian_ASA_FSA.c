/* -----------------------------------------------------------------
 * Programmer(s): Radu Serban @ LLNL
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
 *
 * Hessian through adjoint sensitivity example problem.
 *
 *        [ - p1 * y1^2 - y3 ]           [ 1 ]
 *   y' = [    - y2          ]    y(0) = [ 1 ]
 *        [ -p2^2 * y2 * y3  ]           [ 1 ]
 *
 *   p1 = 1.0
 *   p2 = 2.0
 *
 *           2
 *          /
 *   G(p) = |  0.5 * ( y1^2 + y2^2 + y3^2 ) dt
 *          /
 *          0
 *
 * Compute the gradient (ASA) and Hessian (FSA over ASA) of G(p).
 *
 * See D.B. Ozyurt and P.I. Barton, SISC 26(5) 1725-1743, 2005.
 *
 * -----------------------------------------------------------------*/

#include <cvodes/cvodes.h>          /* prototypes for CVODES fcts., consts. */
#include <nvector/nvector_serial.h> /* access to serial N_Vector            */
#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_math.h> /* definition of SUNRabs, SUNRexp, etc. */
#include <sunlinsol/sunlinsol_dense.h> /* access to band SUNLinearSolver       */
#include <sunmatrix/sunmatrix_dense.h> /* access to band SUNMatrix             */

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define ESYM "Le"
#define GSYM "Lg"
#else
#define ESYM "e"
#define GSYM "g"
#endif

#define ZERO SUN_RCONST(0.0)
#define ONE  SUN_RCONST(1.0)

typedef struct
{
  sunrealtype p1, p2;
}* UserData;

static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
static int fQ(sunrealtype t, N_Vector y, N_Vector qdot, void* user_data);
static int fS(int Ns, sunrealtype t, N_Vector y, N_Vector ydot, N_Vector* yS,
              N_Vector* ySdot, void* user_data, N_Vector tmp1, N_Vector tmp2);
static int fQS(int Ns, sunrealtype t, N_Vector y, N_Vector* yS, N_Vector yQdot,
               N_Vector* yQSdot, void* user_data, N_Vector tmp, N_Vector tmpQ);

static int fB1(sunrealtype t, N_Vector y, N_Vector* yS, N_Vector yB,
               N_Vector yBdot, void* user_dataB);
static int fQB1(sunrealtype t, N_Vector y, N_Vector* yS, N_Vector yB,
                N_Vector qBdot, void* user_dataB);

static int fB2(sunrealtype t, N_Vector y, N_Vector* yS, N_Vector yB,
               N_Vector yBdot, void* user_dataB);
static int fQB2(sunrealtype t, N_Vector y, N_Vector* yS, N_Vector yB,
                N_Vector qBdot, void* user_dataB);

int PrintFwdStats(void* cvode_mem);
int PrintBckStats(void* cvode_mem, int idx);

/* Private function to check function return values */

static int check_retval(void* returnvalue, const char* funcname, int opt);

/*
 *--------------------------------------------------------------------
 * MAIN PROGRAM
 *--------------------------------------------------------------------
 */

int main(int argc, char* argv[])
{
  SUNContext sunctx;
  UserData data;

  SUNMatrix A, AB1, AB2;
  SUNLinearSolver LS, LSB1, LSB2;
  void* cvode_mem;

  sunindextype Neq, Np2;
  int Np;

  sunrealtype t0, tf;

  sunrealtype reltol;
  sunrealtype abstol, abstolQ, abstolB, abstolQB;

  N_Vector y, yQ;
  N_Vector *yS, *yQS;
  N_Vector yB1, yB2;
  N_Vector yQB1, yQB2;

  int steps, ncheck;
  int indexB1, indexB2;

  int retval;
  sunrealtype time;

  sunrealtype dp;
  sunrealtype G, Gp, Gm;
  sunrealtype grdG_fwd[2], grdG_bck[2], grdG_cntr[2];
  sunrealtype H11, H22;

  data = NULL;
  y = yQ = NULL;
  yB1 = yB2 = NULL;
  yQB1 = yQB2 = NULL;
  A = AB1 = AB2 = NULL;
  LS = LSB1 = LSB2 = NULL;
  cvode_mem        = NULL;

  /* User data structure */

  data     = (UserData)malloc(sizeof *data);
  data->p1 = SUN_RCONST(1.0);
  data->p2 = SUN_RCONST(2.0);

  /* Problem size, integration interval, and tolerances */

  Neq = 3;
  Np  = 2;
  Np2 = 2 * Np;

  t0 = 0.0;
  tf = 2.0;

  reltol = 1.0e-8;

  abstol  = 1.0e-8;
  abstolQ = 1.0e-8;

  abstolB  = 1.0e-8;
  abstolQB = 1.0e-8;

  /* Create the SUNDIALS simulation context that all SUNDIALS objects require */
  retval = SUNContext_Create(SUN_COMM_NULL, &sunctx);
  if (check_retval(&retval, "SUNContext_Create", 1)) { return (1); }

  /* Initializations for forward problem */

  y = N_VNew_Serial(Neq, sunctx);
  if (check_retval((void*)y, "N_VNew_Serial", 0)) { return (1); }
  sunrealtype* y_data = N_VGetArrayPointer(y);
  N_VConst(ONE, y);

  yQ = N_VNew_Serial(1, sunctx);
  if (check_retval((void*)yQ, "N_VNew_Serial", 0)) { return (1); }
  sunrealtype* yQ_data = N_VGetArrayPointer(yQ);
  N_VConst(ZERO, yQ);

  yS = N_VCloneVectorArray(Np, y);
  if (check_retval((void*)yS, "N_VCloneVectorArray", 0)) { return (1); }
  N_VConst(ZERO, yS[0]);
  N_VConst(ZERO, yS[1]);

  yQS = N_VCloneVectorArray(Np, yQ);
  if (check_retval((void*)yQS, "N_VCloneVectorArray", 0)) { return (1); }
  N_VConst(ZERO, yQS[0]);
  N_VConst(ZERO, yQS[1]);

  /* Create and initialize forward problem */

  cvode_mem = CVodeCreate(CV_BDF, sunctx);
  if (check_retval((void*)cvode_mem, "CVodeCreate", 0)) { return (1); }

  retval = CVodeInit(cvode_mem, f, t0, y);
  if (check_retval(&retval, "CVodeInit", 1)) { return (1); }

  retval = CVodeSStolerances(cvode_mem, reltol, abstol);
  if (check_retval(&retval, "CVodeSStolerances", 1)) { return (1); }

  retval = CVodeSetUserData(cvode_mem, data);
  if (check_retval(&retval, "CVodeSetUserData", 1)) { return (1); }

  /* Create a dense SUNMatrix */
  A = SUNDenseMatrix(Neq, Neq, sunctx);
  if (check_retval((void*)A, "SUNDenseMatrix", 0)) { return (1); }

  /* Create banded SUNLinearSolver for the forward problem */
  LS = SUNLinSol_Dense(y, A, sunctx);
  if (check_retval((void*)LS, "SUNLinSol_Dense", 0)) { return (1); }

  /* Attach the matrix and linear solver */
  retval = CVodeSetLinearSolver(cvode_mem, LS, A);
  if (check_retval(&retval, "CVodeSetLinearSolver", 1)) { return (1); }

  retval = CVodeQuadInit(cvode_mem, fQ, yQ);
  if (check_retval(&retval, "CVodeQuadInit", 1)) { return (1); }

  retval = CVodeQuadSStolerances(cvode_mem, reltol, abstolQ);
  if (check_retval(&retval, "CVodeQuadSStolerances", 1)) { return (1); }

  retval = CVodeSetQuadErrCon(cvode_mem, SUNTRUE);
  if (check_retval(&retval, "CVodeSetQuadErrCon", 1)) { return (1); }

  retval = CVodeSensInit(cvode_mem, Np, CV_SIMULTANEOUS, fS, yS);
  if (check_retval(&retval, "CVodeSensInit", 1)) { return (1); }

  retval = CVodeSensEEtolerances(cvode_mem);
  if (check_retval(&retval, "CVodeSensEEtolerances", 1)) { return (1); }

  retval = CVodeSetSensErrCon(cvode_mem, SUNTRUE);
  if (check_retval(&retval, "CVodeSetSensErrCon", 1)) { return (1); }

  retval = CVodeQuadSensInit(cvode_mem, fQS, yQS);
  if (check_retval(&retval, "CVodeQuadSensInit", 1)) { return (1); }

  retval = CVodeQuadSensEEtolerances(cvode_mem);
  if (check_retval(&retval, "CVodeQuadSensEEtolerances", 1)) { return (1); }

  retval = CVodeSetQuadSensErrCon(cvode_mem, SUNTRUE);
  if (check_retval(&retval, "CVodeSetQuadSensErrCon", 1)) { return (1); }

  /* Initialize ASA */

  steps  = 100;
  retval = CVodeAdjInit(cvode_mem, steps, CV_POLYNOMIAL);
  if (check_retval(&retval, "CVodeAdjInit", 1)) { return (1); }

  /* Forward integration */

  printf("-------------------\n");
  printf("Forward integration\n");
  printf("-------------------\n\n");

  retval = CVodeF(cvode_mem, tf, y, &time, CV_NORMAL, &ncheck);
  if (check_retval(&retval, "CVodeF", 1)) { return (1); }

  retval = CVodeGetQuad(cvode_mem, &time, yQ);
  if (check_retval(&retval, "CVodeGetQuad", 1)) { return (1); }

  G = yQ_data[0];

  retval = CVodeGetSens(cvode_mem, &time, yS);
  if (check_retval(&retval, "CVodeGetSens", 1)) { return (1); }

  retval = CVodeGetQuadSens(cvode_mem, &time, yQS);
  if (check_retval(&retval, "CVodeGetQuadSens", 1)) { return (1); }

  sunrealtype* yS0_data = N_VGetArrayPointer(yS[0]);
  sunrealtype* yS1_data = N_VGetArrayPointer(yS[1]);
  sunrealtype* yQS0_data = N_VGetArrayPointer(yQS[0]);
  sunrealtype* yQS1_data = N_VGetArrayPointer(yQS[1]);

  printf("ncheck = %d\n", ncheck);
  printf("\n");
  printf("     y:    %12.4" ESYM " %12.4" ESYM " %12.4" ESYM, y_data[0], y_data[1], y_data[2]);
  printf("     G:    %12.4" ESYM "\n", yQ_data[0]);
  printf("\n");
  printf("     yS1:  %12.4" ESYM " %12.4" ESYM " %12.4" ESYM "\n", yS0_data[0], yS0_data[1],
         yS0_data[2]);
  printf("     yS2:  %12.4" ESYM " %12.4" ESYM " %12.4" ESYM "\n", yS1_data[0], yS1_data[1],
         yS1_data[2]);
  printf("\n");
  printf("   dG/dp:  %12.4" ESYM " %12.4" ESYM "\n", yQS0_data[0], yQS1_data[0]);
  printf("\n");

  printf("Final Statistics for forward pb.\n");
  printf("--------------------------------\n");
  retval = PrintFwdStats(cvode_mem);
  if (check_retval(&retval, "PrintFwdStats", 1)) { return (1); }

  /* Initializations for backward problems */

  yB1 = N_VNew_Serial(2 * Neq, sunctx);
  if (check_retval((void*)yB1, "N_VNew_Serial", 0)) { return (1); }
  N_VConst(ZERO, yB1);

  yQB1 = N_VNew_Serial(Np2, sunctx);
  if (check_retval((void*)yQB1, "N_VNew_Serial", 0)) { return (1); }
  sunrealtype* yQB1_data = N_VGetArrayPointer(yQB1);
  N_VConst(ZERO, yQB1);

  yB2 = N_VNew_Serial(2 * Neq, sunctx);
  if (check_retval((void*)yB2, "N_VNew_Serial", 0)) { return (1); }
  N_VConst(ZERO, yB2);

  yQB2 = N_VNew_Serial(Np2, sunctx);
  if (check_retval((void*)yQB2, "N_VNew_Serial", 0)) { return (1); }
  sunrealtype* yQB2_data = N_VGetArrayPointer(yQB2);
  N_VConst(ZERO, yQB2);

  /* Create and initialize backward problems (one for each column of the Hessian) */

  /* -------------------------
     First backward problem
     -------------------------*/

  retval = CVodeCreateB(cvode_mem, CV_BDF, &indexB1);
  if (check_retval(&retval, "CVodeCreateB", 1)) { return (1); }

  retval = CVodeInitBS(cvode_mem, indexB1, fB1, tf, yB1);
  if (check_retval(&retval, "CVodeInitBS", 1)) { return (1); }

  retval = CVodeSStolerancesB(cvode_mem, indexB1, reltol, abstolB);
  if (check_retval(&retval, "CVodeSStolerancesB", 1)) { return (1); }

  retval = CVodeSetUserDataB(cvode_mem, indexB1, data);
  if (check_retval(&retval, "CVodeSetUserDataB", 1)) { return (1); }

  retval = CVodeQuadInitBS(cvode_mem, indexB1, fQB1, yQB1);
  if (check_retval(&retval, "CVodeQuadInitBS", 1)) { return (1); }

  retval = CVodeQuadSStolerancesB(cvode_mem, indexB1, reltol, abstolQB);
  if (check_retval(&retval, "CVodeQuadSStolerancesB", 1)) { return (1); }

  retval = CVodeSetQuadErrConB(cvode_mem, indexB1, SUNTRUE);
  if (check_retval(&retval, "CVodeSetQuadErrConB", 1)) { return (1); }

  /* Create a dense SUNMatrix */
  AB1 = SUNDenseMatrix(2 * Neq, 2 * Neq, sunctx);
  if (check_retval((void*)A, "SUNDenseMatrix", 0)) { return (1); }

  /* Create dense SUNLinearSolver for the forward problem */
  LSB1 = SUNLinSol_Dense(yB1, AB1, sunctx);
  if (check_retval((void*)LSB1, "SUNLinSol_Dense", 0)) { return (1); }

  /* Attach the matrix and linear solver */
  retval = CVodeSetLinearSolverB(cvode_mem, indexB1, LSB1, AB1);
  if (check_retval(&retval, "CVodeSetLinearSolverB", 1)) { return (1); }

  /* -------------------------
     Second backward problem
     -------------------------*/

  retval = CVodeCreateB(cvode_mem, CV_BDF, &indexB2);
  if (check_retval(&retval, "CVodeCreateB", 1)) { return (1); }

  retval = CVodeInitBS(cvode_mem, indexB2, fB2, tf, yB2);
  if (check_retval(&retval, "CVodeInitBS", 1)) { return (1); }

  retval = CVodeSStolerancesB(cvode_mem, indexB2, reltol, abstolB);
  if (check_retval(&retval, "CVodeSStolerancesB", 1)) { return (1); }

  retval = CVodeSetUserDataB(cvode_mem, indexB2, data);
  if (check_retval(&retval, "CVodeSetUserDataB", 1)) { return (1); }

  retval = CVodeQuadInitBS(cvode_mem, indexB2, fQB2, yQB2);
  if (check_retval(&retval, "CVodeQuadInitBS", 1)) { return (1); }

  retval = CVodeQuadSStolerancesB(cvode_mem, indexB2, reltol, abstolQB);
  if (check_retval(&retval, "CVodeQuadSStolerancesB", 1)) { return (1); }

  retval = CVodeSetQuadErrConB(cvode_mem, indexB2, SUNTRUE);
  if (check_retval(&retval, "CVodeSetQuadErrConB", 1)) { return (1); }

  /* Create a dense SUNMatrix */
  AB2 = SUNDenseMatrix(2 * Neq, 2 * Neq, sunctx);
  if (check_retval((void*)AB2, "SUNDenseMatrix", 0)) { return (1); }

  /* Create dense SUNLinearSolver for the forward problem */
  LSB2 = SUNLinSol_Dense(yB2, AB2, sunctx);
  if (check_retval((void*)LSB2, "SUNLinSol_Dense", 0)) { return (1); }

  /* Attach the matrix and linear solver */
  retval = CVodeSetLinearSolverB(cvode_mem, indexB2, LSB2, AB2);
  if (check_retval(&retval, "CVodeSetLinearSolverB", 1)) { return (1); }

  /* Backward integration */

  printf("---------------------------------------------\n");
  printf("Backward integration ... (2 adjoint problems)\n");
  printf("---------------------------------------------\n\n");

  retval = CVodeB(cvode_mem, t0, CV_NORMAL);
  if (check_retval(&retval, "CVodeB", 1)) { return (1); }

  retval = CVodeGetB(cvode_mem, indexB1, &time, yB1);
  if (check_retval(&retval, "CVodeGetB", 1)) { return (1); }

  retval = CVodeGetQuadB(cvode_mem, indexB1, &time, yQB1);
  if (check_retval(&retval, "CVodeGetQuadB", 1)) { return (1); }

  retval = CVodeGetB(cvode_mem, indexB2, &time, yB2);
  if (check_retval(&retval, "CVodeGetB", 1)) { return (1); }

  retval = CVodeGetQuadB(cvode_mem, indexB2, &time, yQB2);
  if (check_retval(&retval, "CVodeGetQuadB", 1)) { return (1); }

  printf("   dG/dp:  %12.4" ESYM " %12.4" ESYM "   (from backward pb. 1)\n", -yQB1_data[0],
         -yQB1_data[1]);
  printf("           %12.4" ESYM " %12.4" ESYM "   (from backward pb. 2)\n", -yQB2_data[0],
         -yQB2_data[1]);
  printf("\n");
  printf("   H = d2G/dp2:\n");
  printf("        (1)            (2)\n");
  printf("  %12.4" ESYM "   %12.4" ESYM "\n", -yQB1_data[2], -yQB2_data[2]);
  printf("  %12.4" ESYM "   %12.4" ESYM "\n", -yQB1_data[3], -yQB2_data[3]);
  printf("\n");

  printf("Final Statistics for backward pb. 1\n");
  printf("-----------------------------------\n");
  retval = PrintBckStats(cvode_mem, indexB1);
  if (check_retval(&retval, "PrintBckStats", 1)) { return (1); }

  printf("Final Statistics for backward pb. 2\n");
  printf("-----------------------------------\n");
  retval = PrintBckStats(cvode_mem, indexB2);
  if (check_retval(&retval, "PrintBckStats", 1)) { return (1); }

  /* Free memory */

  CVodeFree(&cvode_mem);
  SUNLinSolFree(LS);
  SUNMatDestroy(A);
  SUNLinSolFree(LSB1);
  SUNMatDestroy(AB1);
  SUNLinSolFree(LSB2);
  SUNMatDestroy(AB2);

  /* Finite difference tests */

  dp = SUN_RCONST(1.0e-2);

  printf("-----------------------\n");
  printf("Finite Difference tests\n");
  printf("-----------------------\n\n");

  printf("del_p = %" GSYM "\n\n", dp);

  cvode_mem = CVodeCreate(CV_BDF, sunctx);

  N_VConst(ONE, y);
  N_VConst(ZERO, yQ);

  retval = CVodeInit(cvode_mem, f, t0, y);
  if (check_retval(&retval, "CVodeInit", 1)) { return (1); }

  retval = CVodeSStolerances(cvode_mem, reltol, abstol);
  if (check_retval(&retval, "CVodeSStolerances", 1)) { return (1); }

  retval = CVodeSetUserData(cvode_mem, data);
  if (check_retval(&retval, "CVodeSetUserData", 1)) { return (1); }

  /* Create a dense SUNMatrix */
  A = SUNDenseMatrix(Neq, Neq, sunctx);
  if (check_retval((void*)A, "SUNDenseMatrix", 0)) { return (1); }

  /* Create dense SUNLinearSolver for the forward problem */
  LS = SUNLinSol_Dense(y, A, sunctx);
  if (check_retval((void*)LS, "SUNLinSol_Dense", 0)) { return (1); }

  /* Attach the matrix and linear solver */
  retval = CVodeSetLinearSolver(cvode_mem, LS, A);
  if (check_retval(&retval, "CVodeSetLinearSolver", 1)) { return (1); }

  retval = CVodeQuadInit(cvode_mem, fQ, yQ);
  if (check_retval(&retval, "CVodeQuadInit", 1)) { return (1); }

  retval = CVodeQuadSStolerances(cvode_mem, reltol, abstolQ);
  if (check_retval(&retval, "CVodeQuadSStolerances", 1)) { return (1); }

  retval = CVodeSetQuadErrCon(cvode_mem, SUNTRUE);
  if (check_retval(&retval, "CVodeSetQuadErrCon", 1)) { return (1); }

  data->p1 += dp;

  retval = CVode(cvode_mem, tf, y, &time, CV_NORMAL);
  if (check_retval(&retval, "CVode", 1)) { return (1); }

  retval = CVodeGetQuad(cvode_mem, &time, yQ);
  if (check_retval(&retval, "CVodeGetQuad", 1)) { return (1); }

  Gp = yQ_data[0];

  printf("p1+  y:   %12.4" ESYM " %12.4" ESYM " %12.4" ESYM, y_data[0], y_data[1], y_data[2]);
  printf("     G:   %12.4" ESYM "\n", yQ_data[0]);
  data->p1 -= 2.0 * dp;

  N_VConst(ONE, y);
  N_VConst(ZERO, yQ);

  CVodeReInit(cvode_mem, t0, y);
  CVodeQuadReInit(cvode_mem, yQ);

  retval = CVode(cvode_mem, tf, y, &time, CV_NORMAL);
  if (check_retval(&retval, "CVode", 1)) { return (1); }

  retval = CVodeGetQuad(cvode_mem, &time, yQ);
  if (check_retval(&retval, "CVodeGetQuad", 1)) { return (1); }

  Gm = yQ_data[0];
  printf("p1-  y:   %12.4" ESYM " %12.4" ESYM " %12.4" ESYM, y_data[0], y_data[1], y_data[2]);
  printf("     G:   %12.4" ESYM "\n", yQ_data[0]);
  data->p1 += dp;

  grdG_fwd[0]  = (Gp - G) / dp;
  grdG_bck[0]  = (G - Gm) / dp;
  grdG_cntr[0] = (Gp - Gm) / (2.0 * dp);
  H11          = (Gp - 2.0 * G + Gm) / (dp * dp);

  data->p2 += dp;

  N_VConst(ONE, y);
  N_VConst(ZERO, yQ);

  CVodeReInit(cvode_mem, t0, y);
  CVodeQuadReInit(cvode_mem, yQ);

  retval = CVode(cvode_mem, tf, y, &time, CV_NORMAL);
  if (check_retval(&retval, "CVode", 1)) { return (1); }

  retval = CVodeGetQuad(cvode_mem, &time, yQ);
  if (check_retval(&retval, "CVodeGetQuad", 1)) { return (1); }

  Gp = yQ_data[0];
  printf("p2+  y:   %12.4" ESYM " %12.4" ESYM " %12.4" ESYM, y_data[0], y_data[1], y_data[2]);
  printf("     G:   %12.4" ESYM "\n", yQ_data[0]);
  data->p2 -= 2.0 * dp;

  N_VConst(ONE, y);
  N_VConst(ZERO, yQ);

  CVodeReInit(cvode_mem, t0, y);
  CVodeQuadReInit(cvode_mem, yQ);

  retval = CVode(cvode_mem, tf, y, &time, CV_NORMAL);
  if (check_retval(&retval, "CVode", 1)) { return (1); }

  retval = CVodeGetQuad(cvode_mem, &time, yQ);
  if (check_retval(&retval, "CVodeGetQuad", 1)) { return (1); }

  Gm = yQ_data[0];
  printf("p2-  y:   %12.4" ESYM " %12.4" ESYM " %12.4" ESYM, y_data[0], y_data[1], y_data[2]);
  printf("     G:   %12.4" ESYM "\n", yQ_data[0]);
  data->p2 += dp;

  grdG_fwd[1]  = (Gp - G) / dp;
  grdG_bck[1]  = (G - Gm) / dp;
  grdG_cntr[1] = (Gp - Gm) / (2.0 * dp);
  H22          = (Gp - 2.0 * G + Gm) / (dp * dp);

  printf("\n");

  printf("   dG/dp:  %12.4" ESYM " %12.4" ESYM "   (fwd FD)\n", grdG_fwd[0], grdG_fwd[1]);
  printf("           %12.4" ESYM " %12.4" ESYM "   (bck FD)\n", grdG_bck[0], grdG_bck[1]);
  printf("           %12.4" ESYM " %12.4" ESYM "   (cntr FD)\n", grdG_cntr[0], grdG_cntr[1]);
  printf("\n");
  printf("  H(1,1):  %12.4" ESYM "\n", H11);
  printf("  H(2,2):  %12.4" ESYM "\n", H22);

  /* Free memory */

  CVodeFree(&cvode_mem);
  SUNLinSolFree(LS);
  SUNMatDestroy(A);

  N_VDestroy(y);
  N_VDestroy(yQ);

  N_VDestroyVectorArray(yS, Np);
  N_VDestroyVectorArray(yQS, Np);

  N_VDestroy(yB1);
  N_VDestroy(yQB1);
  N_VDestroy(yB2);
  N_VDestroy(yQB2);

  SUNContext_Free(&sunctx);

  free(data);

  return (0);
}

/*
 *--------------------------------------------------------------------
 * FUNCTIONS CALLED BY CVODES
 *--------------------------------------------------------------------
 */

static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  sunrealtype y0, y1, y2;
  UserData data;
  sunrealtype p1, p2;
  sunrealtype* ydot_data = N_VGetArrayPointer(ydot);
  sunrealtype* y_data = N_VGetArrayPointer(y);

  data = (UserData)user_data;
  p1   = data->p1;
  p2   = data->p2;

  y0 = y_data[0];
  y1 = y_data[1];
  y2 = y_data[2];

  ydot_data[0] = -p1 * y0 * y0 - y2;
  ydot_data[1] = -y1;
  ydot_data[2] = -p2 * p2 * y1 * y2;

  return (0);
}

static int fQ(sunrealtype t, N_Vector y, N_Vector qdot, void* user_data)
{
  sunrealtype y0, y1, y2;
  sunrealtype* y_data = N_VGetArrayPointer(y);
  sunrealtype* qdot_data = N_VGetArrayPointer(qdot);

  y0 = y_data[0];
  y1 = y_data[1];
  y2 = y_data[2];

  qdot_data[0] = 0.5 * (y0 * y0 + y1 * y1 + y2 * y2);

  return (0);
}

static int fS(int Ns, sunrealtype t, N_Vector y, N_Vector ydot, N_Vector* yS,
              N_Vector* ySdot, void* user_data, N_Vector tmp1, N_Vector tmp2)
{
  UserData data;
  sunrealtype y0, y1, y2;
  sunrealtype s0, s1, s2;
  sunrealtype fys0, fys1, fys2;
  sunrealtype p1, p2;
  sunrealtype* y_data = N_VGetArrayPointer(y);

  data = (UserData)user_data;
  p1   = data->p1;
  p2   = data->p2;

  y0 = y_data[0];
  y1 = y_data[1];
  y2 = y_data[2];

  sunrealtype* yS0_data = N_VGetArrayPointer(yS[0]);
  sunrealtype* ySdot0_data = N_VGetArrayPointer(ySdot[0]);
  sunrealtype* ySdot1_data = N_VGetArrayPointer(ySdot[1]);

  /* 1st sensitivity RHS */

  s0 = yS0_data[0];
  s1 = yS0_data[1];
  s2 = yS0_data[2];

  fys0 = -2.0 * p1 * y0 * s0 - s2;
  fys1 = -s1;
  fys2 = -p2 * p2 * y2 * s1 - p2 * p2 * y1 * s2;

  ySdot0_data[0] = fys0 - y0 * y0;
  ySdot0_data[1] = fys1;
  ySdot0_data[2] = fys2;

  /* 2nd sensitivity RHS */

  sunrealtype* yS1_data = N_VGetArrayPointer(yS[1]);
  s0 = yS1_data[0];
  s1 = yS1_data[1];
  s2 = yS1_data[2];

  fys0 = -2.0 * p1 * y0 * s0 - s2;
  fys1 = -s1;
  fys2 = -p2 * p2 * y2 * s1 - p2 * p2 * y1 * s2;

  ySdot1_data[0] = fys0;
  ySdot1_data[1] = fys1;
  ySdot1_data[2] = fys2 - 2.0 * p2 * y1 * y2;

  return (0);
}

static int fQS(int Ns, sunrealtype t, N_Vector y, N_Vector* yS, N_Vector yQdot,
               N_Vector* yQSdot, void* user_data, N_Vector tmp, N_Vector tmpQ)
{
  sunrealtype y0, y1, y2;
  sunrealtype s0, s1, s2;
  sunrealtype* y_data = N_VGetArrayPointer(y);

  y0 = y_data[0];
  y1 = y_data[1];
  y2 = y_data[2];

  sunrealtype* yS0_data = N_VGetArrayPointer(yS[0]);
  sunrealtype* yQSdot0_data = N_VGetArrayPointer(yQSdot[0]);
  sunrealtype* yQSdot1_data = N_VGetArrayPointer(yQSdot[1]);

  /* 1st sensitivity RHS */

  s0 = yS0_data[0];
  s1 = yS0_data[1];
  s2 = yS0_data[2];

  yQSdot0_data[0] = y0 * s0 + y1 * s1 + y2 * s2;

  /* 2nd sensitivity RHS */

  sunrealtype* yS1_data = N_VGetArrayPointer(yS[1]);
  s0 = yS1_data[0];
  s1 = yS1_data[1];
  s2 = yS1_data[2];

  yQSdot1_data[0] = y0 * s0 + y1 * s1 + y2 * s2;

  return (0);
}

static int fB1(sunrealtype t, N_Vector y, N_Vector* yS, N_Vector yB,
               N_Vector yBdot, void* user_dataB)
{
  UserData data;
  sunrealtype p1, p2;
  sunrealtype y0, y1, y2; /* solution */
  sunrealtype s0, s1, s2; /* sensitivity 1 */
  sunrealtype l0, l1, l2; /* lambda */
  sunrealtype m0, m1, m2; /* mu */
  sunrealtype* yBdot_data = N_VGetArrayPointer(yBdot);
  sunrealtype* yB_data = N_VGetArrayPointer(yB);
  sunrealtype* y_data = N_VGetArrayPointer(y);

  data = (UserData)user_dataB;
  p1   = data->p1;
  p2   = data->p2;

  y0 = y_data[0];
  y1 = y_data[1];
  y2 = y_data[2];

  sunrealtype* yS0_data = N_VGetArrayPointer(yS[0]);
  s0 = yS0_data[0];
  s1 = yS0_data[1];
  s2 = yS0_data[2];

  l0 = yB_data[0];
  l1 = yB_data[1];
  l2 = yB_data[2];

  m0 = yB_data[3];
  m1 = yB_data[4];
  m2 = yB_data[5];

  yBdot_data[0] = 2.0 * p1 * y0 * l0 - y0;
  yBdot_data[1] = l1 + p2 * p2 * y2 * l2 - y1;
  yBdot_data[2] = l0 + p2 * p2 * y1 * l2 - y2;

  yBdot_data[3] = 2.0 * p1 * y0 * m0 + l0 * 2.0 * (y0 + p1 * s0) - s0;
  yBdot_data[4] = m1 + p2 * p2 * y2 * m2 + l2 * p2 * p2 * s2 - s1;
  yBdot_data[5] = m0 + p2 * p2 * y1 * m2 + l2 * p2 * p2 * s1 - s2;

  return (0);
}

static int fQB1(sunrealtype t, N_Vector y, N_Vector* yS, N_Vector yB,
                N_Vector qBdot, void* user_dataB)
{
  UserData data;
  sunrealtype p2;
  sunrealtype y0, y1, y2; /* solution */
  sunrealtype s0, s1, s2; /* sensitivity 1 */
  sunrealtype l0, l2;     /* lambda */
  sunrealtype m0, m2;     /* mu */
  sunrealtype* yB_data = N_VGetArrayPointer(yB);
  sunrealtype* y_data = N_VGetArrayPointer(y);
  sunrealtype* qBdot_data = N_VGetArrayPointer(qBdot);

  data = (UserData)user_dataB;

  p2 = data->p2;

  y0 = y_data[0];
  y1 = y_data[1];
  y2 = y_data[2];

  sunrealtype* yS0_data = N_VGetArrayPointer(yS[0]);
  s0 = yS0_data[0];
  s1 = yS0_data[1];
  s2 = yS0_data[2];

  l0 = yB_data[0];
  l2 = yB_data[2];

  m0 = yB_data[3];
  m2 = yB_data[5];

  qBdot_data[0] = -y0 * y0 * l0;
  qBdot_data[1] = -2.0 * p2 * y1 * y2 * l2;

  qBdot_data[2] = -y0 * y0 * m0 - l0 * 2.0 * y0 * s0;
  qBdot_data[3] = -2.0 * p2 * y1 * y2 * m2 -
                  l2 * 2.0 * (p2 * y2 * s1 + p2 * y1 * s2);

  return (0);
}

static int fB2(sunrealtype t, N_Vector y, N_Vector* yS, N_Vector yB,
               N_Vector yBdot, void* user_dataB)
{
  UserData data;
  sunrealtype p1, p2;
  sunrealtype y0, y1, y2; /* solution */
  sunrealtype s0, s1, s2; /* sensitivity 2 */
  sunrealtype l0, l1, l2; /* lambda */
  sunrealtype m0, m1, m2; /* mu */
  sunrealtype* yBdot_data = N_VGetArrayPointer(yBdot);
  sunrealtype* yB_data = N_VGetArrayPointer(yB);
  sunrealtype* y_data = N_VGetArrayPointer(y);

  data = (UserData)user_dataB;
  p1   = data->p1;
  p2   = data->p2;

  y0 = y_data[0];
  y1 = y_data[1];
  y2 = y_data[2];

  sunrealtype* yS1_data = N_VGetArrayPointer(yS[1]);
  s0 = yS1_data[0];
  s1 = yS1_data[1];
  s2 = yS1_data[2];

  l0 = yB_data[0];
  l1 = yB_data[1];
  l2 = yB_data[2];

  m0 = yB_data[3];
  m1 = yB_data[4];
  m2 = yB_data[5];

  yBdot_data[0] = 2.0 * p1 * y0 * l0 - y0;
  yBdot_data[1] = l1 + p2 * p2 * y2 * l2 - y1;
  yBdot_data[2] = l0 + p2 * p2 * y1 * l2 - y2;

  yBdot_data[3] = 2.0 * p1 * y0 * m0 + l0 * 2.0 * p1 * s0 - s0;
  yBdot_data[4] = m1 + p2 * p2 * y2 * m2 + l2 * (2.0 * p2 * y2 + p2 * p2 * s2) -
                  s1;
  yBdot_data[5] = m0 + p2 * p2 * y1 * m2 + l2 * (2.0 * p2 * y1 + p2 * p2 * s1) -
                  s2;

  return (0);
}

static int fQB2(sunrealtype t, N_Vector y, N_Vector* yS, N_Vector yB,
                N_Vector qBdot, void* user_dataB)
{
  UserData data;
  sunrealtype p2;
  sunrealtype y0, y1, y2; /* solution */
  sunrealtype s0, s1, s2; /* sensitivity 2 */
  sunrealtype l0, l2;     /* lambda */
  sunrealtype m0, m2;     /* mu */
  sunrealtype* yB_data = N_VGetArrayPointer(yB);
  sunrealtype* y_data = N_VGetArrayPointer(y);
  sunrealtype* qBdot_data = N_VGetArrayPointer(qBdot);

  data = (UserData)user_dataB;

  p2 = data->p2;

  y0 = y_data[0];
  y1 = y_data[1];
  y2 = y_data[2];

  sunrealtype* yS1_data = N_VGetArrayPointer(yS[1]);
  s0 = yS1_data[0];
  s1 = yS1_data[1];
  s2 = yS1_data[2];

  l0 = yB_data[0];
  l2 = yB_data[2];

  m0 = yB_data[3];
  m2 = yB_data[5];

  qBdot_data[0] = -y0 * y0 * l0;
  qBdot_data[1] = -2.0 * p2 * y1 * y2 * l2;

  qBdot_data[2] = -y0 * y0 * m0 - l0 * 2.0 * y0 * s0;
  qBdot_data[3] = -2.0 * p2 * y1 * y2 * m2 -
                  l2 * 2.0 * (p2 * y2 * s1 + p2 * y1 * s2 + y1 * y2);

  return (0);
}

/*
 *--------------------------------------------------------------------
 * PRIVATE FUNCTIONS
 *--------------------------------------------------------------------
 */

int PrintFwdStats(void* cvode_mem)
{
  long int nst, nfe, nsetups, nni, ncfn, netf;
  long int nfQe, netfQ;
  long int nfSe, nfeS, nsetupsS, netfS;
  long int nfQSe, netfQS;

  int qlast, qcur;
  sunrealtype h0u, hlast, hcur, tcur;

  int retval;

  retval = CVodeGetIntegratorStats(cvode_mem, &nst, &nfe, &nsetups, &netf,
                                   &qlast, &qcur, &h0u, &hlast, &hcur, &tcur);

  retval = CVodeGetNonlinSolvStats(cvode_mem, &nni, &ncfn);

  retval = CVodeGetQuadStats(cvode_mem, &nfQe, &netfQ);

  retval = CVodeGetSensStats(cvode_mem, &nfSe, &nfeS, &netfS, &nsetupsS);

  retval = CVodeGetQuadSensStats(cvode_mem, &nfQSe, &netfQS);

  printf(" Number steps: %5ld\n\n", nst);
  printf(" Function evaluations:\n");
  printf("  f:        %5ld\n  fQ:       %5ld\n  fS:       %5ld\n  fQS:      "
         "%5ld\n",
         nfe, nfQe, nfSe, nfQSe);
  printf(" Error test failures:\n");
  printf("  netf:     %5ld\n  netfQ:    %5ld\n  netfS:    %5ld\n  netfQS:   "
         "%5ld\n",
         netf, netfQ, netfS, netfQS);
  printf(" Linear solver setups:\n");
  printf("  nsetups:  %5ld\n  nsetupsS: %5ld\n", nsetups, nsetupsS);
  printf(" Nonlinear iterations:\n");
  printf("  nni:      %5ld\n", nni);
  printf(" Convergence failures:\n");
  printf("  ncfn:     %5ld\n", ncfn);

  printf("\n");

  return (retval);
}

int PrintBckStats(void* cvode_mem, int idx)
{
  void* cvode_mem_bck;

  long int nst, nfe, nsetups, nni, ncfn, netf;
  long int nfQe, netfQ;

  int qlast, qcur;
  sunrealtype h0u, hlast, hcur, tcur;

  int retval;

  cvode_mem_bck = CVodeGetAdjCVodeBmem(cvode_mem, idx);

  retval = CVodeGetIntegratorStats(cvode_mem_bck, &nst, &nfe, &nsetups, &netf,
                                   &qlast, &qcur, &h0u, &hlast, &hcur, &tcur);

  retval = CVodeGetNonlinSolvStats(cvode_mem_bck, &nni, &ncfn);

  retval = CVodeGetQuadStats(cvode_mem_bck, &nfQe, &netfQ);

  printf(" Number steps: %5ld\n\n", nst);
  printf(" Function evaluations:\n");
  printf("  f:        %5ld\n  fQ:       %5ld\n", nfe, nfQe);
  printf(" Error test failures:\n");
  printf("  netf:     %5ld\n  netfQ:    %5ld\n", netf, netfQ);
  printf(" Linear solver setups:\n");
  printf("  nsetups:  %5ld\n", nsetups);
  printf(" Nonlinear iterations:\n");
  printf("  nni:      %5ld\n", nni);
  printf(" Convergence failures:\n");
  printf("  ncfn:     %5ld\n", ncfn);

  printf("\n");

  return (retval);
}

/*
 * Check function return value...
 *   opt == 0 means SUNDIALS function allocates memory so check if
 *            returned NULL pointer
 *   opt == 1 means SUNDIALS function returns an integer value so check if
 *            retval < 0
 *   opt == 2 means function allocates memory so check if returned
 *            NULL pointer
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
