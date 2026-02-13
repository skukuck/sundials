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
 *
 * Hessian using adjoint sensitivity example problem.
 *
 * This simple example problem for IDAS, due to Robertson,
 * is from chemical kinetics, and consists of the following three
 * equations:
 *
 *   y1' + p1 * y1 - p2 * y2 * y3             = 0
 *   y2' - p1 * y1 + p2 * y2 * y3 + p3 * y2^2 = 0
 *   y1 + y2 + y3 - 1                         = 0
 *
 *        [1]        [-p1]
 *   y(0)=[0]  y'(0)=[ p1]   p1 = 0.04   p2 = 1e4   p3 = 1e07
 *        [0]        [ 0 ]
 *
 *       80
 *      /
 *  G = | 0.5 * (y1^2 + y2^2 + y3^2) dt
 *      /
 *      0
 * Compute the gradient (using FSA and ASA) and Hessian (FSA over ASA)
 * of G with respect to parameters p1 and p2.
 *
 * Reference: D.B. Ozyurt and P.I. Barton, SISC 26(5) 1725-1743, 2005.
 *
 * Error handling was suppressed for code readability reasons.
 * -----------------------------------------------------------------*/

#include <idas/idas.h> /* prototypes for IDA fcts., consts.    */
#include <math.h>
#include <nvector/nvector_serial.h> /* access to serial N_Vector            */
#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_math.h> /* defs. of SUNRabs, SUNRexp, etc.      */
#include <sundials/sundials_types.h> /* defs. of sunrealtype, sunindextype      */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */

/* Accessor macros */

/* Problem Constants */
#define NEQ 3 /* number of equations                  */
#define NP  2 /* number of sensitivities              */

#define T0 SUN_RCONST(0.0)  /* Initial time. */
#define TF SUN_RCONST(80.0) /* Final time. */

/* Tolerances */
#define RTOL  SUN_RCONST(1e-08) /* scalar relative tolerance            */
#define ATOL  SUN_RCONST(1e-10) /* vector absolute tolerance components */
#define RTOLA SUN_RCONST(1e-08) /* for adjoint integration              */
#define ATOLA SUN_RCONST(1e-08) /* for adjoint integration              */

/* Parameters */
#define P1 SUN_RCONST(0.04)
#define P2 SUN_RCONST(1.0e4)
#define P3 SUN_RCONST(3.0e7)

/* Predefined consts */
#define HALF SUN_RCONST(0.5)
#define ZERO SUN_RCONST(0.0)
#define ONE  SUN_RCONST(1.0)
#define TWO  SUN_RCONST(2.0)

/* User defined struct */
typedef struct
{
  sunrealtype p[3];
}* UserData;

/* residual for forward problem */
static int res(sunrealtype t, N_Vector yy, N_Vector yp, N_Vector resval,
               void* user_data);

static int resS(int Ns, sunrealtype t, N_Vector yy, N_Vector yp,
                N_Vector resval, N_Vector* yyS, N_Vector* ypS, N_Vector* resvalS,
                void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

static int rhsQ(sunrealtype t, N_Vector yy, N_Vector yp, N_Vector qdot,
                void* user_data);

static int rhsQS(int Ns, sunrealtype t, N_Vector yy, N_Vector yp, N_Vector* yyS,
                 N_Vector* ypS, N_Vector rrQ, N_Vector* rhsvalQS,
                 void* user_data, N_Vector yytmp, N_Vector yptmp, N_Vector tmpQS);

static int resBS1(sunrealtype tt, N_Vector yy, N_Vector yp, N_Vector* yyS,
                  N_Vector* ypS, N_Vector yyB, N_Vector ypB, N_Vector resvalBQ,
                  void* user_dataB);

static int rhsQBS1(sunrealtype tt, N_Vector yy, N_Vector yp, N_Vector* yyS,
                   N_Vector* ypS, N_Vector yyB, N_Vector ypB, N_Vector rhsBQS,
                   void* user_dataB);

static int resBS2(sunrealtype tt, N_Vector yy, N_Vector yp, N_Vector* yyS,
                  N_Vector* ypS, N_Vector yyB, N_Vector ypB, N_Vector resvalBQ,
                  void* user_dataB);

static int rhsQBS2(sunrealtype tt, N_Vector yy, N_Vector yp, N_Vector* yyS,
                   N_Vector* ypS, N_Vector yyB, N_Vector ypB, N_Vector rhsBQS,
                   void* user_dataB);

static int check_retval(void* returnvalue, const char* funcname, int opt);

int main(int argc, char* argv[])
{
  N_Vector yy, yp, q, *yyS, *ypS, *qS;
  N_Vector yyB1, ypB1, qB1, yyB2, ypB2, qB2;
  void* ida_mem;
  UserData data;
  sunrealtype time, ti, tf;
  int retval, nckp, indexB1, indexB2, is;
  sunrealtype G, Gm, Gp, dp1, dp2, grdG_fwd[2], grdG_bck[2], grdG_cntr[2], H11,
    H22;
  sunrealtype rtolFD, atolFD;
  sunrealtype* yy_data;
  sunrealtype* yp_data;
  SUNMatrix A, AB1, AB2;
  SUNLinearSolver LS, LSB1, LSB2;

  /* Create the SUNDIALS context object for this simulation. */
  SUNContext ctx = NULL;
  SUNContext_Create(SUN_COMM_NULL, &ctx);

  /* Print problem description */
  printf("\nAdjoint Sensitivity Example for Chemical Kinetics\n");
  printf("---------------------------------------------------------\n");
  printf("DAE: dy1/dt + p1*y1 - p2*y2*y3 = 0\n");
  printf("     dy2/dt - p1*y1 + p2*y2*y3 + p3*(y2)^2 = 0\n");
  printf("               y1  +  y2  +  y3 = 0\n\n");
  printf("Find dG/dp and d^2G/dp^2, where p=[p1,p2] for\n");
  printf("     G = int_t0^tB0 g(t,p,y) dt\n");
  printf("     g(t,p,y) = y3\n\n\n");

  /* Allocate and initialize user data. */
  data       = (UserData)malloc(sizeof(*data));
  data->p[0] = P1;
  data->p[1] = P2;
  data->p[2] = P3;

  /* Consistent IC */
  yy         = N_VNew_Serial(NEQ, ctx);
  yp         = N_VClone(yy);
  yy_data    = N_VGetArrayPointer(yy);
  yp_data    = N_VGetArrayPointer(yp);
  yy_data[0] = ONE;
  yy_data[1] = ZERO;
  yy_data[2] = ZERO;
  yp_data[0] = -P1;
  yp_data[1] = P1;
  yp_data[2] = 0;

  q                   = N_VNew_Serial(1, ctx);
  sunrealtype* q_data = N_VGetArrayPointer(q);
  N_VConst(ZERO, q);

  yyS = N_VCloneVectorArray(NP, yy);
  ypS = N_VCloneVectorArray(NP, yp);
  N_VConst(ZERO, yyS[0]);
  N_VConst(ZERO, yyS[1]);
  N_VConst(ZERO, ypS[0]);
  N_VConst(ZERO, ypS[1]);

  qS = N_VCloneVectorArray(NP, q);
  for (is = 0; is < NP; is++) { N_VConst(ZERO, qS[is]); }

  ida_mem = IDACreate(ctx);

  ti     = T0;
  retval = IDAInit(ida_mem, res, ti, yy, yp);

  /* Forward problem's setup. */
  retval = IDASStolerances(ida_mem, RTOL, ATOL);

  /* Create dense SUNMatrix for use in linear solves */
  A = SUNDenseMatrix(NEQ, NEQ, ctx);
  if (check_retval((void*)A, "SUNDenseMatrix", 0)) { return (1); }

  /* Create dense SUNLinearSolver object */
  LS = SUNLinSol_Dense(yy, A, ctx);
  if (check_retval((void*)LS, "SUNLinSol_Dense", 0)) { return (1); }

  /* Attach the matrix and linear solver */
  retval = IDASetLinearSolver(ida_mem, LS, A);
  if (check_retval(&retval, "IDASetLinearSolver", 1)) { return (1); }

  retval = IDASetUserData(ida_mem, data);
  retval = IDASetMaxNumSteps(ida_mem, 1500);

  /* Quadrature's setup. */
  retval = IDAQuadInit(ida_mem, rhsQ, q);
  retval = IDAQuadSStolerances(ida_mem, RTOL, ATOL);
  retval = IDASetQuadErrCon(ida_mem, SUNTRUE);

  /* Sensitivity's setup. */
  retval = IDASensInit(ida_mem, NP, IDA_SIMULTANEOUS, resS, yyS, ypS);
  retval = IDASensEEtolerances(ida_mem);
  retval = IDASetSensErrCon(ida_mem, SUNTRUE);

  /* Setup of quadrature's sensitivities */
  retval = IDAQuadSensInit(ida_mem, rhsQS, qS);
  retval = IDAQuadSensEEtolerances(ida_mem);
  retval = IDASetQuadSensErrCon(ida_mem, SUNTRUE);

  /* Initialize ASA. */
  retval = IDAAdjInit(ida_mem, 100, IDA_HERMITE);

  printf("---------------------------------------------------------\n");
  printf("Forward integration\n");
  printf("---------------------------------------------------------\n\n");

  tf     = TF;
  retval = IDASolveF(ida_mem, tf, &time, yy, yp, IDA_NORMAL, &nckp);

  IDAGetQuad(ida_mem, &time, q);
  G = q_data[0];
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("     G:    %12.4Le\n", G);
#else
  printf("     G:    %12.4e\n", G);
#endif

  /* Sensitivities are needed for IC of backward problems. */
  IDAGetSensDky(ida_mem, tf, 0, yyS);
  IDAGetSensDky(ida_mem, tf, 1, ypS);

  IDAGetQuadSens(ida_mem, &time, qS);
  sunrealtype* qS0_data  = N_VGetArrayPointer(qS[0]);
  sunrealtype* qS1_data  = N_VGetArrayPointer(qS[1]);
  sunrealtype* yyS0_data = N_VGetArrayPointer(yyS[0]);
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("   dG/dp:  %12.4Le %12.4Le\n", qS0_data[0], qS1_data[0]);
#else
  printf("   dG/dp:  %12.4e %12.4e\n", qS0_data[0], qS1_data[0]);
#endif
  printf("\n");
  /******************************
  * BACKWARD PROBLEM #1
  *******************************/

  /* Consistent IC. */
  yyB1                   = N_VNew_Serial(2 * NEQ, ctx);
  sunrealtype* yyB1_data = N_VGetArrayPointer(yyB1);
  ypB1                   = N_VClone(yyB1);
  sunrealtype* ypB1_data = N_VGetArrayPointer(ypB1);

  N_VConst(ZERO, yyB1);
  yyB1_data[2] = yy_data[2];
  yyB1_data[5] = yyS0_data[2];

  N_VConst(ZERO, ypB1);
  ypB1_data[0] = yy_data[2] - yy_data[0];
  ypB1_data[1] = yy_data[2] - yy_data[1];
  ypB1_data[3] = yyS0_data[2] - yyS0_data[0];
  ypB1_data[4] = yyS0_data[2] - yyS0_data[1];

  qB1                   = N_VNew_Serial(2 * NP, ctx);
  sunrealtype* qB1_data = N_VGetArrayPointer(qB1);
  N_VConst(ZERO, qB1);

  retval = IDACreateB(ida_mem, &indexB1);
  retval = IDAInitBS(ida_mem, indexB1, resBS1, tf, yyB1, ypB1);
  retval = IDASStolerancesB(ida_mem, indexB1, RTOLA, ATOLA);
  retval = IDASetUserDataB(ida_mem, indexB1, data);
  retval = IDASetMaxNumStepsB(ida_mem, indexB1, 5000);

  /* Create dense SUNMatrix for use in linear solves */
  AB1 = SUNDenseMatrix(2 * NEQ, 2 * NEQ, ctx);
  if (check_retval((void*)AB1, "SUNDenseMatrix", 0)) { return (1); }

  /* Create dense SUNLinearSolver object */
  LSB1 = SUNLinSol_Dense(yyB1, AB1, ctx);
  if (check_retval((void*)LSB1, "SUNLinSol_Dense", 0)) { return (1); }

  /* Attach the matrix and linear solver */
  retval = IDASetLinearSolverB(ida_mem, indexB1, LSB1, AB1);
  if (check_retval(&retval, "IDASetLinearSolverB", 1)) { return (1); }

  retval = IDAQuadInitBS(ida_mem, indexB1, rhsQBS1, qB1);

  /******************************
  * BACKWARD PROBLEM #2
  *******************************/

  /* Consistent IC. */
  yyB2                   = N_VNew_Serial(2 * NEQ, ctx);
  sunrealtype* yyB2_data = N_VGetArrayPointer(yyB2);
  ypB2                   = N_VNew_Serial(2 * NEQ, ctx);
  sunrealtype* ypB2_data = N_VGetArrayPointer(ypB2);
  sunrealtype* yyS1_data = N_VGetArrayPointer(yyS[1]);

  N_VConst(ZERO, yyB2);
  yyB2_data[2] = yy_data[2];
  yyB2_data[5] = yyS1_data[2];

  N_VConst(ZERO, ypB2);
  ypB2_data[0] = yy_data[2] - yy_data[0];
  ypB2_data[1] = yy_data[2] - yy_data[1];
  ypB2_data[3] = yyS1_data[2] - yyS1_data[0];
  ypB2_data[4] = yyS1_data[2] - yyS1_data[1];

  qB2                   = N_VNew_Serial(2 * NP, ctx);
  sunrealtype* qB2_data = N_VGetArrayPointer(qB2);
  N_VConst(ZERO, qB2);

  retval = IDACreateB(ida_mem, &indexB2);
  retval = IDAInitBS(ida_mem, indexB2, resBS2, tf, yyB2, ypB2);
  retval = IDASStolerancesB(ida_mem, indexB2, RTOLA, ATOLA);
  retval = IDASetUserDataB(ida_mem, indexB2, data);
  retval = IDASetMaxNumStepsB(ida_mem, indexB2, 2500);

  /* Create dense SUNMatrix for use in linear solves */
  AB2 = SUNDenseMatrix(2 * NEQ, 2 * NEQ, ctx);
  if (check_retval((void*)AB2, "SUNDenseMatrix", 0)) { return (1); }

  /* Create dense SUNLinearSolver object */
  LSB2 = SUNLinSol_Dense(yyB2, AB2, ctx);
  if (check_retval((void*)LSB2, "SUNLinSol_Dense", 0)) { return (1); }

  /* Attach the matrix and linear solver */
  retval = IDASetLinearSolverB(ida_mem, indexB2, LSB2, AB2);
  if (check_retval(&retval, "IDASetLinearSolverB", 1)) { return (1); }

  retval = IDAQuadInitBS(ida_mem, indexB2, rhsQBS2, qB2);

  /* Integrate backward problems. */
  printf("---------------------------------------------------------\n");
  printf("Backward integration \n");
  printf("---------------------------------------------------------\n\n");

  retval = IDASolveB(ida_mem, ti, IDA_NORMAL);

  retval = IDAGetB(ida_mem, indexB1, &time, yyB1, ypB1);
  /*
     retval = IDAGetNumSteps(IDAGetAdjIDABmem(ida_mem, indexB1), &nst);
     printf("at time=%g \tpb 1 Num steps:%d\n", time, nst);
     retval = IDAGetNumSteps(IDAGetAdjIDABmem(ida_mem, indexB2), &nst);
     printf("at time=%g \tpb 2 Num steps:%d\n\n", time, nst);
  */

  retval = IDAGetQuadB(ida_mem, indexB1, &time, qB1);
  retval = IDAGetQuadB(ida_mem, indexB2, &time, qB2);
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("   dG/dp:  %12.4Le %12.4Le   (from backward pb. 1)\n", qB1_data[0],
         qB1_data[1]);
  printf("   dG/dp:  %12.4Le %12.4Le   (from backward pb. 2)\n", qB2_data[0],
         qB2_data[1]);
#else
  printf("   dG/dp:  %12.4e %12.4e   (from backward pb. 1)\n", qB1_data[0],
         qB1_data[1]);
  printf("   dG/dp:  %12.4e %12.4e   (from backward pb. 2)\n", qB2_data[0],
         qB2_data[1]);
#endif

  printf("\n");
  printf("   H = d2G/dp2:\n");
  printf("        (1)            (2)\n");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("  %12.4Le  %12.4Le\n", qB1_data[2], qB2_data[2]);
  printf("  %12.4Le  %12.4Le\n", qB1_data[3], qB2_data[3]);
#else
  printf("  %12.4e  %12.4e\n", qB1_data[2], qB2_data[2]);
  printf("  %12.4e  %12.4e\n", qB1_data[3], qB2_data[3]);
#endif

  IDAFree(&ida_mem);
  SUNLinSolFree(LS);
  SUNMatDestroy(A);
  SUNLinSolFree(LSB1);
  SUNMatDestroy(AB1);
  SUNLinSolFree(LSB2);
  SUNMatDestroy(AB2);

  /*********************************
  * Use Finite Differences to verify
  **********************************/

  /* Perturbations are of different magnitudes as p1 and p2 are. */
  dp1 = SUN_RCONST(1.0e-3);
  dp2 = SUN_RCONST(2.5e+2);

  printf("\n");
  printf("---------------------------------------------------------\n");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("Finite Differences ( dp1=%6.1Le and dp2 = %6.1Le )\n", dp1, dp2);
#else
  printf("Finite Differences ( dp1=%6.1e and dp2 = %6.1e )\n", dp1, dp2);
#endif
  printf("---------------------------------------------------------\n\n");

  ida_mem = IDACreate(ctx);

  /********************
  * Forward FD for p1
  ********************/
  data->p[0] += dp1;

  yy_data[0] = ONE;
  yy_data[1] = ZERO;
  yy_data[2] = ZERO;
  yp_data[0] = -data->p[0];
  yp_data[1] = -yp_data[0];
  yp_data[2] = 0;
  N_VConst(ZERO, q);
  ti = T0;
  tf = TF;

  retval = IDAInit(ida_mem, res, ti, yy, yp);

  rtolFD = SUN_RCONST(1.0e-12);
  atolFD = SUN_RCONST(1.0e-14);

  retval = IDASStolerances(ida_mem, rtolFD, atolFD);

  /* Create dense SUNMatrix for use in linear solves */
  A = SUNDenseMatrix(NEQ, NEQ, ctx);
  if (check_retval((void*)A, "SUNDenseMatrix", 0)) { return (1); }

  /* Create dense SUNLinearSolver object */
  LS = SUNLinSol_Dense(yy, A, ctx);
  if (check_retval((void*)LS, "SUNLinSol_Dense", 0)) { return (1); }

  /* Attach the matrix and linear solver */
  retval = IDASetLinearSolver(ida_mem, LS, A);
  if (check_retval(&retval, "IDASetLinearSolver", 1)) { return (1); }

  retval = IDASetUserData(ida_mem, data);
  retval = IDASetMaxNumSteps(ida_mem, 10000);

  retval = IDAQuadInit(ida_mem, rhsQ, q);
  retval = IDAQuadSStolerances(ida_mem, rtolFD, atolFD);
  retval = IDASetQuadErrCon(ida_mem, SUNTRUE);

  retval = IDASolve(ida_mem, tf, &time, yy, yp, IDA_NORMAL);
  retval = IDAGetQuad(ida_mem, &time, q);
  Gp     = q_data[0];

  /********************
  * Backward FD for p1
  ********************/
  data->p[0] -= 2 * dp1;

  yy_data[0] = ONE;
  yy_data[1] = ZERO;
  yy_data[2] = ZERO;
  yp_data[0] = -data->p[0];
  yp_data[1] = -yp_data[0];
  yp_data[2] = 0;
  N_VConst(ZERO, q);

  retval = IDAReInit(ida_mem, ti, yy, yp);
  retval = IDAQuadReInit(ida_mem, q);

  retval = IDASolve(ida_mem, tf, &time, yy, yp, IDA_NORMAL);
  retval = IDAGetQuad(ida_mem, &time, q);
  Gm     = q_data[0];

  /* Compute FD for p1. */
  grdG_fwd[0]  = (Gp - G) / dp1;
  grdG_bck[0]  = (G - Gm) / dp1;
  grdG_cntr[0] = (Gp - Gm) / (2.0 * dp1);
  H11          = (Gp - 2.0 * G + Gm) / (dp1 * dp1);

  /********************
  * Forward FD for p2
  ********************/
  /*restore p1*/
  data->p[0] += dp1;
  data->p[1] += dp2;

  yy_data[0] = ONE;
  yy_data[1] = ZERO;
  yy_data[2] = ZERO;
  yp_data[0] = -data->p[0];
  yp_data[1] = -yp_data[0];
  yp_data[2] = 0;
  N_VConst(ZERO, q);

  retval = IDAReInit(ida_mem, ti, yy, yp);
  retval = IDAQuadReInit(ida_mem, q);

  retval = IDASolve(ida_mem, tf, &time, yy, yp, IDA_NORMAL);
  retval = IDAGetQuad(ida_mem, &time, q);
  Gp     = q_data[0];

  /********************
  * Backward FD for p2
  ********************/
  data->p[1] -= 2 * dp2;

  yy_data[0] = ONE;
  yy_data[1] = ZERO;
  yy_data[2] = ZERO;
  yp_data[0] = -data->p[0];
  yp_data[1] = -yp_data[0];
  yp_data[2] = 0;
  N_VConst(ZERO, q);

  retval = IDAReInit(ida_mem, ti, yy, yp);
  retval = IDAQuadReInit(ida_mem, q);

  retval = IDASolve(ida_mem, tf, &time, yy, yp, IDA_NORMAL);
  retval = IDAGetQuad(ida_mem, &time, q);
  Gm     = q_data[0];

  /* Compute FD for p2. */
  grdG_fwd[1]  = (Gp - G) / dp2;
  grdG_bck[1]  = (G - Gm) / dp2;
  grdG_cntr[1] = (Gp - Gm) / (2.0 * dp2);
  H22          = (Gp - 2.0 * G + Gm) / (dp2 * dp2);

  printf("\n");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("   dG/dp:  %12.4Le  %12.4Le   (fwd FD)\n", grdG_fwd[0], grdG_fwd[1]);
  printf("           %12.4Le  %12.4Le   (bck FD)\n", grdG_bck[0], grdG_bck[1]);
  printf("           %12.4Le  %12.4Le   (cntr FD)\n", grdG_cntr[0], grdG_cntr[1]);
  printf("\n");
  printf("  H(1,1):  %12.4Le\n", H11);
  printf("  H(2,2):  %12.4Le\n", H22);
#else
  printf("   dG/dp:  %12.4e  %12.4e   (fwd FD)\n", grdG_fwd[0], grdG_fwd[1]);
  printf("           %12.4e  %12.4e   (bck FD)\n", grdG_bck[0], grdG_bck[1]);
  printf("           %12.4e  %12.4e   (cntr FD)\n", grdG_cntr[0], grdG_cntr[1]);
  printf("\n");
  printf("  H(1,1):  %12.4e\n", H11);
  printf("  H(2,2):  %12.4e\n", H22);
#endif

  IDAFree(&ida_mem);
  SUNLinSolFree(LS);
  SUNMatDestroy(A);

  N_VDestroy(yyB1);
  N_VDestroy(ypB1);
  N_VDestroy(qB1);

  N_VDestroy(yyB2);
  N_VDestroy(ypB2);
  N_VDestroy(qB2);

  N_VDestroy(yy);
  N_VDestroy(yp);
  N_VDestroy(q);
  N_VDestroyVectorArray(yyS, NP);
  N_VDestroyVectorArray(ypS, NP);
  N_VDestroyVectorArray(qS, NP);

  free(data);

  SUNContext_Free(&ctx);
  return 0;
}

static int res(sunrealtype tres, N_Vector yy, N_Vector yp, N_Vector rr,
               void* user_data)
{
  sunrealtype y1, y2, y3, yp1, yp2, *rval;
  UserData data;
  sunrealtype p1, p2, p3;
  sunrealtype* yy_data = N_VGetArrayPointer(yy);
  sunrealtype* yp_data = N_VGetArrayPointer(yp);

  y1   = yy_data[0];
  y2   = yy_data[1];
  y3   = yy_data[2];
  yp1  = yp_data[0];
  yp2  = yp_data[1];
  rval = N_VGetArrayPointer(rr);

  data = (UserData)user_data;
  p1   = data->p[0];
  p2   = data->p[1];
  p3   = data->p[2];

  rval[0] = p1 * y1 - p2 * y2 * y3;
  rval[1] = -rval[0] + p3 * y2 * y2 + yp2;
  rval[0] += yp1;
  rval[2] = y1 + y2 + y3 - 1;

  return (0);
}

static int resS(int Ns, sunrealtype t, N_Vector yy, N_Vector yp,
                N_Vector resval, N_Vector* yyS, N_Vector* ypS, N_Vector* resvalS,
                void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  UserData data;
  sunrealtype p1, p2, p3;
  sunrealtype y1, y2, y3;
  sunrealtype s1, s2, s3;
  sunrealtype sd1, sd2;
  sunrealtype rs1, rs2, rs3;
  int is;
  sunrealtype* yy_data = N_VGetArrayPointer(yy);

  data = (UserData)user_data;
  p1   = data->p[0];
  p2   = data->p[1];
  p3   = data->p[2];

  y1 = yy_data[0];
  y2 = yy_data[1];
  y3 = yy_data[2];

  for (is = 0; is < NP; is++)
  {
    sunrealtype* yySis_data     = N_VGetArrayPointer(yyS[is]);
    sunrealtype* ypSis_data     = N_VGetArrayPointer(ypS[is]);
    sunrealtype* resvalSis_data = N_VGetArrayPointer(resvalS[is]);

    s1 = yySis_data[0];
    s2 = yySis_data[1];
    s3 = yySis_data[2];

    sd1 = ypSis_data[0];
    sd2 = ypSis_data[1];

    rs1 = sd1 + p1 * s1 - p2 * y3 * s2 - p2 * y2 * s3;
    rs2 = sd2 - p1 * s1 + p2 * y3 * s2 + p2 * y2 * s3 + TWO * p3 * y2 * s2;
    rs3 = s1 + s2 + s3;

    switch (is)
    {
    case 0:
      rs1 += y1;
      rs2 -= y1;
      break;
    case 1:
      rs1 -= y2 * y3;
      rs2 += y2 * y3;
      break;
    }

    resvalSis_data[0] = rs1;
    resvalSis_data[1] = rs2;
    resvalSis_data[2] = rs3;
  }

  return (0);
}

static int rhsQ(sunrealtype t, N_Vector yy, N_Vector yp, N_Vector qdot,
                void* user_data)
{
  sunrealtype y1, y2, y3;
  sunrealtype* yy_data   = N_VGetArrayPointer(yy);
  sunrealtype* qdot_data = N_VGetArrayPointer(qdot);

  y1           = yy_data[0];
  y2           = yy_data[1];
  y3           = yy_data[2];
  qdot_data[0] = HALF * (y1 * y1 + y2 * y2 + y3 * y3);

  return (0);
}

static int rhsQS(int Ns, sunrealtype t, N_Vector yy, N_Vector yp, N_Vector* yyS,
                 N_Vector* ypS, N_Vector rrQ, N_Vector* rhsvalQS,
                 void* user_data, N_Vector yytmp, N_Vector yptmp, N_Vector tmpQS)
{
  sunrealtype y1, y2, y3;
  sunrealtype s1, s2, s3;
  sunrealtype* yy_data        = N_VGetArrayPointer(yy);
  sunrealtype* yyS0_data      = N_VGetArrayPointer(yyS[0]);
  sunrealtype* yyS1_data      = N_VGetArrayPointer(yyS[1]);
  sunrealtype* rhsvalQS0_data = N_VGetArrayPointer(rhsvalQS[0]);
  sunrealtype* rhsvalQS1_data = N_VGetArrayPointer(rhsvalQS[1]);

  y1 = yy_data[0];
  y2 = yy_data[1];
  y3 = yy_data[2];

  /* 1st sensitivity RHS */
  s1                = yyS0_data[0];
  s2                = yyS0_data[1];
  s3                = yyS0_data[2];
  rhsvalQS0_data[0] = y1 * s1 + y2 * s2 + y3 * s3;

  /* 2nd sensitivity RHS */
  s1                = yyS1_data[0];
  s2                = yyS1_data[1];
  s3                = yyS1_data[2];
  rhsvalQS1_data[0] = y1 * s1 + y2 * s2 + y3 * s3;

  return (0);
}

/* Residuals for adjoint model. */
static int resBS1(sunrealtype tt, N_Vector yy, N_Vector yp, N_Vector* yyS,
                  N_Vector* ypS, N_Vector yyB, N_Vector ypB, N_Vector rrBS,
                  void* user_dataB)

{
  UserData data;
  sunrealtype y1, y2, y3;
  sunrealtype p1, p2, p3;
  sunrealtype l1, l2, l3, m1, m2, m3;
  sunrealtype lp1, lp2, mp1, mp2;
  sunrealtype s1, s2, s3;
  sunrealtype l21;
  sunrealtype* yyB_data  = N_VGetArrayPointer(yyB);
  sunrealtype* yy_data   = N_VGetArrayPointer(yy);
  sunrealtype* ypB_data  = N_VGetArrayPointer(ypB);
  sunrealtype* rrBS_data = N_VGetArrayPointer(rrBS);

  data = (UserData)user_dataB;

  /* The parameters. */
  p1 = data->p[0];
  p2 = data->p[1];
  p3 = data->p[2];

  /* The y vector. */
  y1 = yy_data[0];
  y2 = yy_data[1];
  y3 = yy_data[2];

  /* The lambda vector. */
  l1 = yyB_data[0];
  l2 = yyB_data[1];
  l3 = yyB_data[2];
  /* The mu vector. */
  m1 = yyB_data[3];
  m2 = yyB_data[4];
  m3 = yyB_data[5];

  /* The lambda dot vector. */
  lp1 = ypB_data[0];
  lp2 = ypB_data[1];
  /* The mu dot vector. */
  mp1 = ypB_data[3];
  mp2 = ypB_data[4];

  /* The sensitivity with respect to p1 */
  sunrealtype* yyS0_data = N_VGetArrayPointer(yyS[0]);
  s1                     = yyS0_data[0];
  s2                     = yyS0_data[1];
  s3                     = yyS0_data[2];

  /* Temporary variables */
  l21 = l2 - l1;

  rrBS_data[0] = lp1 + p1 * l21 - l3 + y1;
  rrBS_data[1] = lp2 - p2 * y3 * l21 - TWO * p3 * y2 * l2 - l3 + y2;
  rrBS_data[2] = -p2 * y2 * l21 - l3 + y3;

  rrBS_data[3] = mp1 + p1 * (-m1 + m2) - m3 + l21 + s1;
  rrBS_data[4] = mp2 + p2 * y3 * m1 - (p2 * y3 + TWO * p3 * y2) * m2 - m3 +
                 p2 * s3 * l1 - (TWO * p3 * s2 + p2 * s3) * l2 + s2;
  rrBS_data[5] = p2 * y2 * (m1 - m2) - m3 - p2 * s2 * l21 + s3;

  return (0);
}

static int rhsQBS1(sunrealtype tt, N_Vector yy, N_Vector yp, N_Vector* yyS,
                   N_Vector* ypS, N_Vector yyB, N_Vector ypB, N_Vector rhsBQS,
                   void* user_dataB)
{
  sunrealtype y1, y2, y3;
  sunrealtype l1, l2, m1, m2;
  sunrealtype s1, s2, s3;
  sunrealtype l21;
  sunrealtype* yyB_data    = N_VGetArrayPointer(yyB);
  sunrealtype* yy_data     = N_VGetArrayPointer(yy);
  sunrealtype* rhsBQS_data = N_VGetArrayPointer(rhsBQS);

  /* The y vector */
  y1 = yy_data[0];
  y2 = yy_data[1];
  y3 = yy_data[2];

  /* The lambda vector. */
  l1 = yyB_data[0];
  l2 = yyB_data[1];

  /* The mu vector. */
  m1 = yyB_data[3];
  m2 = yyB_data[4];

  /* The sensitivity with respect to p1 */
  sunrealtype* yyS0_data = N_VGetArrayPointer(yyS[0]);
  s1                     = yyS0_data[0];
  s2                     = yyS0_data[1];
  s3                     = yyS0_data[2];

  /* Temporary variables */
  l21 = l2 - l1;

  rhsBQS_data[0] = -y1 * l21;
  rhsBQS_data[1] = y2 * y3 * l21;

  rhsBQS_data[2] = y1 * (m1 - m2) - s1 * l21;
  rhsBQS_data[3] = y2 * y3 * (m2 - m1) + (y3 * s2 + y2 * s3) * l21;

  return (0);
}

static int resBS2(sunrealtype tt, N_Vector yy, N_Vector yp, N_Vector* yyS,
                  N_Vector* ypS, N_Vector yyB, N_Vector ypB, N_Vector rrBS,
                  void* user_dataB)

{
  UserData data;
  sunrealtype y1, y2, y3;
  sunrealtype p1, p2, p3;
  sunrealtype l1, l2, l3, m1, m2, m3;
  sunrealtype lp1, lp2, mp1, mp2;
  sunrealtype s1, s2, s3;
  sunrealtype l21;
  sunrealtype* yyB_data  = N_VGetArrayPointer(yyB);
  sunrealtype* yy_data   = N_VGetArrayPointer(yy);
  sunrealtype* ypB_data  = N_VGetArrayPointer(ypB);
  sunrealtype* rrBS_data = N_VGetArrayPointer(rrBS);

  data = (UserData)user_dataB;

  /* The parameters. */
  p1 = data->p[0];
  p2 = data->p[1];
  p3 = data->p[2];

  /* The y vector. */
  y1 = yy_data[0];
  y2 = yy_data[1];
  y3 = yy_data[2];

  /* The lambda vector. */
  l1 = yyB_data[0];
  l2 = yyB_data[1];
  l3 = yyB_data[2];
  /* The mu vector. */
  m1 = yyB_data[3];
  m2 = yyB_data[4];
  m3 = yyB_data[5];

  /* The lambda dot vector. */
  lp1 = ypB_data[0];
  lp2 = ypB_data[1];

  /* The mu dot vector. */
  mp1 = ypB_data[3];
  mp2 = ypB_data[4];

  /* The sensitivity with respect to p2 */
  sunrealtype* yyS1_data = N_VGetArrayPointer(yyS[1]);
  s1                     = yyS1_data[0];
  s2                     = yyS1_data[1];
  s3                     = yyS1_data[2];

  /* Temporary variables */
  l21 = l2 - l1;

  rrBS_data[0] = lp1 + p1 * l21 - l3 + y1;
  rrBS_data[1] = lp2 - p2 * y3 * l21 - TWO * p3 * y2 * l2 - l3 + y2;
  rrBS_data[2] = -p2 * y2 * l21 - l3 + y3;

  rrBS_data[3] = mp1 + p1 * (-m1 + m2) - m3 + s1;
  rrBS_data[4] = mp2 + p2 * y3 * m1 - (p2 * y3 + TWO * p3 * y2) * m2 - m3 +
                 (y3 + p2 * s3) * l1 - (y3 + TWO * p3 * s2 + p2 * s3) * l2 + s2;
  rrBS_data[5] = p2 * y2 * (m1 - m2) - m3 - (y2 + p2 * s2) * l21 + s3;

  return (0);
}

static int rhsQBS2(sunrealtype tt, N_Vector yy, N_Vector yp, N_Vector* yyS,
                   N_Vector* ypS, N_Vector yyB, N_Vector ypB, N_Vector rhsBQS,
                   void* user_dataB)
{
  sunrealtype y1, y2, y3;
  sunrealtype l1, l2, m1, m2;
  sunrealtype s1, s2, s3;
  sunrealtype l21;
  sunrealtype* yyB_data    = N_VGetArrayPointer(yyB);
  sunrealtype* yy_data     = N_VGetArrayPointer(yy);
  sunrealtype* rhsBQS_data = N_VGetArrayPointer(rhsBQS);

  /* The y vector */
  y1 = yy_data[0];
  y2 = yy_data[1];
  y3 = yy_data[2];

  /* The lambda vector. */
  l1 = yyB_data[0];
  l2 = yyB_data[1];

  /* The mu vector. */
  m1 = yyB_data[3];
  m2 = yyB_data[4];

  /* The sensitivity with respect to p2 */
  sunrealtype* yyS1_data = N_VGetArrayPointer(yyS[1]);
  s1                     = yyS1_data[0];
  s2                     = yyS1_data[1];
  s3                     = yyS1_data[2];

  /* Temporary variables */
  l21 = l2 - l1;

  rhsBQS_data[0] = -y1 * l21;
  rhsBQS_data[1] = y2 * y3 * l21;

  rhsBQS_data[2] = y1 * (m1 - m2) - s1 * l21;
  rhsBQS_data[3] = y2 * y3 * (m2 - m1) + (y3 * s2 + y2 * s3) * l21;

  return (0);
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
