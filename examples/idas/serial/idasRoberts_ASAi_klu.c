/* -----------------------------------------------------------------
 * Programmer(s): Ting Yan @ UMBC
 *      Based on idasRoberts_ASAi_dns.c and modified to use KLU
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
 * Adjoint sensitivity example problem.
 *
 * This simple example problem for IDAS, due to Robertson,
 * is from chemical kinetics, and consists of the following three
 * equations:
 *
 *      dy1/dt + p1*y1 - p2*y2*y3            = 0
 *      dy2/dt - p1*y1 + p2*y2*y3 + p3*y2**2 = 0
 *                 y1  +  y2  +  y3  -  1    = 0
 *
 * on the interval from t = 0.0 to t = 4.e10, with initial
 * conditions: y1 = 1, y2 = y3 = 0.The reaction rates are: p1=0.04,
 * p2=1e4, and p3=3e7
 *
 * It uses a scalar relative tolerance and a vector absolute
 * tolerance.
 *
 * IDAS can also compute sensitivities with respect to
 * the problem parameters p1, p2, and p3 of the following quantity:
 *   G = int_t0^t1 g(t,p,y) dt
 * where
 *   g(t,p,y) = y3
 *
 * The gradient dG/dp is obtained as:
 *   dG/dp = int_t0^t1 (g_p - lambda^T F_p ) dt -
 *           lambda^T*F_y'*y_p | _t0^t1
 *         = int_t0^t1 (lambda^T*F_p) dt
 * where lambda and are solutions of the adjoint system:
 *   d(lambda^T * F_y' )/dt -lambda^T F_y = -g_y
 *
 * During the backward integration, IDAS also evaluates G as
 *   G = - phi(t0)
 * where
 *   d(phi)/dt = g(t,y,p)
 *   phi(t1) = 0
 * -----------------------------------------------------------------*/

#include <idas/idas.h>              /* prototypes for IDA fcts., consts.    */
#include <nvector/nvector_serial.h> /* access to serial N_Vector            */
#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_math.h> /* defs. of SUNRabs, SUNRexp, etc.      */
#include <sundials/sundials_types.h> /* defs. of sunrealtype, sunindextype      */
#include <sunlinsol/sunlinsol_klu.h> /* access to KLU linear solver          */
#include <sunmatrix/sunmatrix_sparse.h> /* access to sparse SUNMatrix           */

/* Accessor macros */

/* Problem Constants */

#define NEQ 3 /* number of equations                  */

#define RTOL SUN_RCONST(1e-06) /* scalar relative tolerance            */

#define ATOL1 SUN_RCONST(1e-08) /* vector absolute tolerance components */
#define ATOL2 SUN_RCONST(1e-12)
#define ATOL3 SUN_RCONST(1e-08)

#define ATOLA SUN_RCONST(1e-08) /* absolute tolerance for adjoint vars. */
#define ATOLQ SUN_RCONST(1e-06) /* absolute tolerance for quadratures   */

#define T0   SUN_RCONST(0.0)  /* initial time                         */
#define TOUT SUN_RCONST(4e10) /* final time                           */

#define TB1 SUN_RCONST(50.0) /* starting point for adjoint problem   */
#define TB2 TOUT             /* starting point for adjoint problem   */

#define T1B SUN_RCONST(49.0) /* for IDACalcICB                       */

#define STEPS 100 /* number of steps between check points */

#define NP 3 /* number of problem parameters         */

#define ONE  SUN_RCONST(1.0)
#define ZERO SUN_RCONST(0.0)

/* Type : UserData */

typedef struct
{
  sunrealtype p[3];
}* UserData;

/* Prototypes of user-supplied functions */

static int res(sunrealtype t, N_Vector yy, N_Vector yp, N_Vector resval,
               void* user_data);
static int Jac(sunrealtype t, sunrealtype cj, N_Vector yy, N_Vector yp,
               N_Vector resvec, SUNMatrix JJ, void* user_data, N_Vector tmp1,
               N_Vector tmp2, N_Vector tmp3);

static int rhsQ(sunrealtype t, N_Vector yy, N_Vector yp, N_Vector qdot,
                void* user_data);
static int ewt(N_Vector y, N_Vector w, void* user_data);

static int resB(sunrealtype tt, N_Vector yy, N_Vector yp, N_Vector yyB,
                N_Vector ypB, N_Vector rrB, void* user_dataB);

static int JacB(sunrealtype tt, sunrealtype cjB, N_Vector yy, N_Vector yp,
                N_Vector yyB, N_Vector ypB, N_Vector rrB, SUNMatrix JB,
                void* user_data, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B);

static int rhsQB(sunrealtype tt, N_Vector yy, N_Vector yp, N_Vector yyB,
                 N_Vector ypB, N_Vector rrQB, void* user_dataB);

/* Prototypes of private functions */
static void PrintOutput(sunrealtype tfinal, N_Vector yB, N_Vector ypB,
                        N_Vector qB);
static int check_retval(void* returnvalue, const char* funcname, int opt);

/*
 *--------------------------------------------------------------------
 * MAIN PROGRAM
 *--------------------------------------------------------------------
 */

int main(int argc, char* argv[])
{
  SUNContext ctx;

  UserData data;

  void* ida_mem;
  SUNMatrix A, AB;
  SUNLinearSolver LS, LSB;

  sunrealtype reltolQ, abstolQ;
  N_Vector yy, yp, q;
  N_Vector yyTB1, ypTB1;
  N_Vector id;

  int steps;

  int indexB;

  sunrealtype reltolB, abstolB, abstolQB;
  N_Vector yB, ypB, qB;
  sunrealtype time;
  int retval, nnz, ncheck;

  IDAadjCheckPointRec* ckpnt;

  long int nst, nstB;
  sunrealtype* qB_data = N_VGetArrayPointer(qB);
  sunrealtype* id_data = N_VGetArrayPointer(id);

  data    = NULL;
  ckpnt   = NULL;
  ida_mem = NULL;
  yy = yp = yB = qB = NULL;
  A = AB = NULL;
  LS = LSB = NULL;

  /* Print problem description */
  printf("\nAdjoint Sensitivity Example for Chemical Kinetics\n");
  printf("-------------------------------------------------\n\n");
  printf("DAE: dy1/dt + p1*y1 - p2*y2*y3 = 0\n");
  printf("     dy2/dt - p1*y1 + p2*y2*y3 + p3*(y2)^2 = 0\n");
  printf("               y1  +  y2  +  y3 = 0\n\n");
  printf("Find dG/dp for\n");
  printf("     G = int_t0^tB0 g(t,p,y) dt\n");
  printf("     g(t,p,y) = y3\n\n\n");

  /* Create the SUNDIALS context object for this simulation */
  retval = SUNContext_Create(SUN_COMM_NULL, &ctx);
  if (check_retval(&retval, "SUNContext_Create", 1)) { return 1; }

  /* User data structure */
  data = (UserData)malloc(sizeof *data);
  if (check_retval((void*)data, "malloc", 2)) { return (1); }
  data->p[0] = SUN_RCONST(0.04);
  data->p[1] = SUN_RCONST(1.0e4);
  data->p[2] = SUN_RCONST(3.0e7);

  /* Initialize y */
  yy = N_VNew_Serial(NEQ, ctx);
  if (check_retval((void*)yy, "N_VNew_Serial", 0)) { return (1); }
  sunrealtype* yy_data = N_VGetArrayPointer(yy);
  yy_data[0]           = ONE;
  yy_data[1]           = ZERO;
  yy_data[2]           = ZERO;

  /* Initialize yprime */
  yp = N_VClone(yy);
  if (check_retval((void*)yp, "N_VNew_Serial", 0)) { return (1); }
  sunrealtype* yp_data = N_VGetArrayPointer(yp);
  yp_data[0]           = SUN_RCONST(-0.04);
  yp_data[1]           = SUN_RCONST(0.04);
  yp_data[2]           = ZERO;

  /* Initialize q */
  q = N_VNew_Serial(1, ctx);
  if (check_retval((void*)q, "N_VNew_Serial", 0)) { return (1); }
  sunrealtype* q_data = N_VGetArrayPointer(q);
  q_data[0]           = ZERO;

  /* Set the scalar relative and absolute tolerances reltolQ and abstolQ */
  reltolQ = RTOL;
  abstolQ = ATOLQ;

  /* Create and allocate IDAS memory for forward run */
  printf("Create and allocate IDAS memory for forward runs\n");

  ida_mem = IDACreate(ctx);
  if (check_retval((void*)ida_mem, "IDACreate", 0)) { return (1); }

  retval = IDAInit(ida_mem, res, T0, yy, yp);
  if (check_retval(&retval, "IDAInit", 1)) { return (1); }

  retval = IDAWFtolerances(ida_mem, ewt);
  if (check_retval(&retval, "IDAWFtolerances", 1)) { return (1); }

  retval = IDASetUserData(ida_mem, data);
  if (check_retval(&retval, "IDASetUserData", 1)) { return (1); }

  /* Create sparse SUNMatrix for use in linear solves */
  nnz = NEQ * NEQ;
  A   = SUNSparseMatrix(NEQ, NEQ, nnz, SUN_CSC_MAT, ctx);
  if (check_retval((void*)A, "SUNSparseMatrix", 0)) { return (1); }

  /* Create KLU SUNLinearSolver object (one thread) */
  LS = SUNLinSol_KLU(yy, A, ctx);
  if (check_retval((void*)LS, "SUNLinSol_KLU", 0)) { return (1); }

  /* Attach the matrix and linear solver */
  retval = IDASetLinearSolver(ida_mem, LS, A);
  if (check_retval(&retval, "IDASetLinearSolver", 1)) { return (1); }

  /* Set the user-supplied Jacobian routine */
  retval = IDASetJacFn(ida_mem, Jac);
  if (check_retval(&retval, "IDASetJacFn", 1)) { return (1); }

  /* Setup quadrature integration */
  retval = IDAQuadInit(ida_mem, rhsQ, q);
  if (check_retval(&retval, "IDAQuadInit", 1)) { return (1); }

  retval = IDAQuadSStolerances(ida_mem, reltolQ, abstolQ);
  if (check_retval(&retval, "IDAQuadSStolerances", 1)) { return (1); }

  retval = IDASetQuadErrCon(ida_mem, SUNTRUE);
  if (check_retval(&retval, "IDASetQuadErrCon", 1)) { return (1); }

  /* Call IDASetMaxNumSteps to set the maximum number of steps the
   * solver will take in an attempt to reach the next output time
   * during forward integration. */
  retval = IDASetMaxNumSteps(ida_mem, 2500);
  if (check_retval(&retval, "IDASetMaxNumSteps", 1)) { return (1); }

  /* Allocate global memory */

  steps  = STEPS;
  retval = IDAAdjInit(ida_mem, steps, IDA_HERMITE);
  /*retval = IDAAdjInit(ida_mem, steps, IDA_POLYNOMIAL);*/
  if (check_retval(&retval, "IDAAdjInit", 1)) { return (1); }

  /* Perform forward run */
  printf("Forward integration ... ");

  /* Integrate till TB1 and get the solution (y, y') at that time. */
  retval = IDASolveF(ida_mem, TB1, &time, yy, yp, IDA_NORMAL, &ncheck);
  if (check_retval(&retval, "IDASolveF", 1)) { return (1); }

  yyTB1 = N_VClone(yy);
  ypTB1 = N_VClone(yp);
  /* Save the states at t=TB1. */
  N_VScale(ONE, yy, yyTB1);
  N_VScale(ONE, yp, ypTB1);

  /* Continue integrating till TOUT is reached. */
  retval = IDASolveF(ida_mem, TOUT, &time, yy, yp, IDA_NORMAL, &ncheck);
  if (check_retval(&retval, "IDASolveF", 1)) { return (1); }

  retval = IDAGetNumSteps(ida_mem, &nst);
  if (check_retval(&retval, "IDAGetNumSteps", 1)) { return (1); }

  printf("done ( nst = %ld )\n", nst);

  retval = IDAGetQuad(ida_mem, &time, q);
  if (check_retval(&retval, "IDAGetQuad", 1)) { return (1); }

  printf("--------------------------------------------------------\n");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("G:          %12.4Le \n", q_data[0]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("G:          %12.4e \n", q_data[0]);
#else
  printf("G:          %12.4e \n", q_data[0]);
#endif
  printf("--------------------------------------------------------\n\n");

  /* Test check point linked list
     (uncomment next block to print check point information) */

  /*
  {
    int i;

    printf("\nList of Check Points (ncheck = %d)\n\n", ncheck);
    ckpnt = (IDAadjCheckPointRec *) malloc ( (ncheck+1)*sizeof(IDAadjCheckPointRec));
    IDAGetAdjCheckPointsInfo(ida_mem, ckpnt);
    for (i=0;i<=ncheck;i++) {
      printf("Address:       %p\n",ckpnt[i].my_addr);
      printf("Next:          %p\n",ckpnt[i].next_addr);
      printf("Time interval: %le  %le\n",ckpnt[i].t0, ckpnt[i].t1);
      printf("Step number:   %ld\n",ckpnt[i].nstep);
      printf("Order:         %d\n",ckpnt[i].order);
      printf("Step size:     %le\n",ckpnt[i].step);
      printf("\n");
    }

  }
  */

  /* Create BACKWARD problem. */

  /* Allocate yB (i.e. lambda_0). */
  yB = N_VClone(yy);
  if (check_retval((void*)yB, "N_VNew_Serial", 0)) { return (1); }
  sunrealtype* yB_data = N_VGetArrayPointer(yB);

  /* Consistently initialize yB. */
  yB_data[0] = ZERO;
  yB_data[1] = ZERO;
  yB_data[2] = ONE;

  /* Allocate ypB (i.e. lambda'_0). */
  ypB = N_VClone(yy);
  if (check_retval((void*)ypB, "N_VNew_Serial", 0)) { return (1); }
  sunrealtype* ypB_data = N_VGetArrayPointer(ypB);

  /* Consistently initialize ypB. */
  ypB_data[0] = ONE;
  ypB_data[1] = ONE;
  ypB_data[2] = ZERO;

  /* Set the scalar relative tolerance reltolB */
  reltolB = RTOL;

  /* Set the scalar absolute tolerance abstolB */
  abstolB = ATOLA;

  /* Set the scalar absolute tolerance abstolQB */
  abstolQB = ATOLQ;

  /* Create and allocate IDAS memory for backward run */
  printf("Create and allocate IDAS memory for backward run\n");

  retval = IDACreateB(ida_mem, &indexB);
  if (check_retval(&retval, "IDACreateB", 1)) { return (1); }

  retval = IDAInitB(ida_mem, indexB, resB, TB2, yB, ypB);
  if (check_retval(&retval, "IDAInitB", 1)) { return (1); }

  retval = IDASStolerancesB(ida_mem, indexB, reltolB, abstolB);
  if (check_retval(&retval, "IDASStolerancesB", 1)) { return (1); }

  retval = IDASetUserDataB(ida_mem, indexB, data);
  if (check_retval(&retval, "IDASetUserDataB", 1)) { return (1); }

  retval = IDASetMaxNumStepsB(ida_mem, indexB, 1000);
  if (check_retval(&retval, "IDASetMaxNumStepsB", 1)) { return (1); }

  /* Create sparse SUNMatrix for use in linear solves */
  AB = SUNSparseMatrix(NEQ, NEQ, nnz, SUN_CSC_MAT, ctx);
  if (check_retval((void*)AB, "SUNSparseMatrix", 0)) { return (1); }

  /* Create KLU SUNLinearSolver object (one thread) */
  LSB = SUNLinSol_KLU(yB, AB, ctx);
  if (check_retval((void*)LSB, "SUNLinSol_KLU", 0)) { return (1); }

  /* Attach the matrix and linear solver */
  retval = IDASetLinearSolverB(ida_mem, indexB, LSB, AB);
  if (check_retval(&retval, "IDASetLinearSolverB", 1)) { return (1); }

  /* Set the user-supplied Jacobian routine */
  retval = IDASetJacFnB(ida_mem, indexB, JacB);
  if (check_retval(&retval, "IDASetJacFnB", 1)) { return (1); }

  /* Quadrature for backward problem. */

  /* Initialize qB */
  qB = N_VNew_Serial(NP, ctx);
  if (check_retval((void*)qB, "N_VNew", 0)) { return (1); }
  qB_data[0] = ZERO;
  qB_data[1] = ZERO;
  qB_data[2] = ZERO;

  retval = IDAQuadInitB(ida_mem, indexB, rhsQB, qB);
  if (check_retval(&retval, "IDAQuadInitB", 1)) { return (1); }

  retval = IDAQuadSStolerancesB(ida_mem, indexB, reltolB, abstolQB);
  if (check_retval(&retval, "IDAQuadSStolerancesB", 1)) { return (1); }

  /* Include quadratures in error control. */
  retval = IDASetQuadErrConB(ida_mem, indexB, SUNTRUE);
  if (check_retval(&retval, "IDASetQuadErrConB", 1)) { return (1); }

  /* Backward Integration */
  printf("Backward integration ... ");

  retval = IDASolveB(ida_mem, T0, IDA_NORMAL);
  if (check_retval(&retval, "IDASolveB", 1)) { return (1); }

  IDAGetNumSteps(IDAGetAdjIDABmem(ida_mem, indexB), &nstB);
  printf("done ( nst = %ld )\n", nstB);

  retval = IDAGetB(ida_mem, indexB, &time, yB, ypB);
  if (check_retval(&retval, "IDAGetB", 1)) { return (1); }

  retval = IDAGetQuadB(ida_mem, indexB, &time, qB);
  if (check_retval(&retval, "IDAGetB", 1)) { return (1); }

  PrintOutput(TB2, yB, ypB, qB);

  /* Reinitialize backward phase and start from a different time (TB1). */
  printf("Re-initialize IDAS memory for backward run\n");

  /* Both algebraic part from y and the entire y' are computed by IDACalcIC. */
  yB_data[0] = ZERO;
  yB_data[1] = ZERO;
  yB_data[2] = SUN_RCONST(0.50); /* not consistent */

  /* Rough guess for ypB. */
  ypB_data[0] = SUN_RCONST(0.80);
  ypB_data[1] = SUN_RCONST(0.75);
  ypB_data[2] = ZERO;

  /* Initialize qB */
  qB_data[0] = ZERO;
  qB_data[1] = ZERO;
  qB_data[2] = ZERO;

  retval = IDAReInitB(ida_mem, indexB, TB1, yB, ypB);
  if (check_retval(&retval, "IDAReInitB", 1)) { return (1); }

  /* Also reinitialize quadratures. */
  retval = IDAQuadReInitB(ida_mem, indexB, qB);
  if (check_retval(&retval, "IDAQuadReInitB", 1)) { return (1); }

  /* Use IDACalcICB to compute consistent initial conditions
     for this backward problem. */

  id         = N_VClone(yy);
  id_data[0] = 1.0;
  id_data[1] = 1.0;
  id_data[2] = 0.0;

  /* Specify which variables are differential (1) and which algebraic (0).*/
  retval = IDASetIdB(ida_mem, indexB, id);
  if (check_retval(&retval, "IDASetId", 1)) { return (1); }

  retval = IDACalcICB(ida_mem, indexB, T1B, yyTB1, ypTB1);
  if (check_retval(&retval, "IDACalcICB", 1)) { return (1); }

  /* Get the consistent IC found by IDAS. */
  retval = IDAGetConsistentICB(ida_mem, indexB, yB, ypB);
  if (check_retval(&retval, "IDAGetConsistentICB", 1)) { return (1); }

  printf("Backward integration ... ");

  retval = IDASolveB(ida_mem, T0, IDA_NORMAL);
  if (check_retval(&retval, "IDASolveB", 1)) { return (1); }

  IDAGetNumSteps(IDAGetAdjIDABmem(ida_mem, indexB), &nstB);
  printf("done ( nst = %ld )\n", nstB);

  retval = IDAGetB(ida_mem, indexB, &time, yB, ypB);
  if (check_retval(&retval, "IDAGetB", 1)) { return (1); }

  retval = IDAGetQuadB(ida_mem, indexB, &time, qB);
  if (check_retval(&retval, "IDAGetQuadB", 1)) { return (1); }

  PrintOutput(TB1, yB, ypB, qB);

  /* Free any memory used.*/

  printf("Free memory\n\n");

  IDAFree(&ida_mem);
  SUNLinSolFree(LS);
  SUNMatDestroy(A);
  SUNLinSolFree(LSB);
  SUNMatDestroy(AB);
  N_VDestroy(yy);
  N_VDestroy(yp);
  N_VDestroy(q);
  N_VDestroy(yB);
  N_VDestroy(ypB);
  N_VDestroy(qB);
  N_VDestroy(id);
  N_VDestroy(yyTB1);
  N_VDestroy(ypTB1);

  if (ckpnt != NULL) { free(ckpnt); }
  free(data);

  SUNContext_Free(&ctx);

  return (0);
}

/*
 *--------------------------------------------------------------------
 * FUNCTIONS CALLED BY IDAS
 *--------------------------------------------------------------------
 */

/*
 * f routine. Compute f(t,y).
*/

static int res(sunrealtype t, N_Vector yy, N_Vector yp, N_Vector resval,
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
  rval = N_VGetArrayPointer(resval);

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

/*
 * Jacobian routine. Compute J(t,y).
*/

static int Jac(sunrealtype t, sunrealtype cj, N_Vector yy, N_Vector yp,
               N_Vector resvec, SUNMatrix JJ, void* user_data, N_Vector tmp1,
               N_Vector tmp2, N_Vector tmp3)
{
  sunrealtype* yval;
  sunindextype* colptrs = SUNSparseMatrix_IndexPointers(JJ);
  sunindextype* rowvals = SUNSparseMatrix_IndexValues(JJ);
  sunrealtype* data     = SUNSparseMatrix_Data(JJ);

  UserData userdata;
  sunrealtype p1, p2, p3;

  yval = N_VGetArrayPointer(yy);

  userdata = (UserData)user_data;
  p1       = userdata->p[0];
  p2       = userdata->p[1];
  p3       = userdata->p[2];

  SUNMatZero(JJ);

  colptrs[0] = 0;
  colptrs[1] = 3;
  colptrs[2] = 6;
  colptrs[3] = 9;

  /* column 0 */
  data[0]    = p1 + cj;
  rowvals[0] = 0;
  data[1]    = -p1;
  rowvals[1] = 1;
  data[2]    = ONE;
  rowvals[2] = 2;

  /* column 1 */
  data[3]    = -p2 * yval[2];
  rowvals[3] = 0;
  data[4]    = p2 * yval[2] + 2 * p3 * yval[1] + cj;
  rowvals[4] = 1;
  data[5]    = ONE;
  rowvals[5] = 2;

  /* column 2 */
  data[6]    = -p2 * yval[1];
  rowvals[6] = 0;
  data[7]    = p2 * yval[1];
  rowvals[7] = 1;
  data[8]    = ONE;
  rowvals[8] = 2;

  return (0);
}

/*
 * rhsQ routine. Compute fQ(t,y).
*/

static int rhsQ(sunrealtype t, N_Vector yy, N_Vector yp, N_Vector qdot,
                void* user_data)
{
  sunrealtype* yy_data   = N_VGetArrayPointer(yy);
  sunrealtype* qdot_data = N_VGetArrayPointer(qdot);
  qdot_data[0]           = yy_data[2];
  return (0);
}

/*
 * EwtSet function. Computes the error weights at the current solution.
 */

static int ewt(N_Vector y, N_Vector w, void* user_data)
{
  int i;
  sunrealtype yy, ww, rtol, atol[3];
  sunrealtype* y_data = N_VGetArrayPointer(y);
  sunrealtype* w_data = N_VGetArrayPointer(w);

  rtol    = RTOL;
  atol[0] = ATOL1;
  atol[1] = ATOL2;
  atol[2] = ATOL3;

  for (i = 0; i < 3; i++)
  {
    yy = y_data[i];
    ww = rtol * SUNRabs(yy) + atol[i];
    if (ww <= 0.0) { return (-1); }
    w_data[i] = 1.0 / ww;
  }

  return (0);
}

/*
 * resB routine.
*/

static int resB(sunrealtype tt, N_Vector yy, N_Vector yp, N_Vector yyB,
                N_Vector ypB, N_Vector rrB, void* user_dataB)
{
  UserData data;
  sunrealtype y2, y3;
  sunrealtype p1, p2, p3;
  sunrealtype l1, l2, l3;
  sunrealtype lp1, lp2;
  sunrealtype l21;
  sunrealtype* yyB_data = N_VGetArrayPointer(yyB);
  sunrealtype* yy_data  = N_VGetArrayPointer(yy);
  sunrealtype* ypB_data = N_VGetArrayPointer(ypB);
  sunrealtype* rrB_data = N_VGetArrayPointer(rrB);

  data = (UserData)user_dataB;

  /* The p vector */
  p1 = data->p[0];
  p2 = data->p[1];
  p3 = data->p[2];

  /* The y  vector */
  y2 = yy_data[1];
  y3 = yy_data[2];

  /* The lambda vector */
  l1 = yyB_data[0];
  l2 = yyB_data[1];
  l3 = yyB_data[2];

  /* The lambda dot vector */
  lp1 = ypB_data[0];
  lp2 = ypB_data[1];

  /* Temporary variables */
  l21 = l2 - l1;

  /* Load residual. */
  rrB_data[0] = lp1 + p1 * l21 - l3;
  rrB_data[1] = lp2 - p2 * y3 * l21 - SUN_RCONST(2.0) * p3 * y2 * l2 - l3;
  rrB_data[2] = -p2 * y2 * l21 - l3 + SUN_RCONST(1.0);

  return (0);
}

/*Jacobian for backward problem. */
static int JacB(sunrealtype tt, sunrealtype cjB, N_Vector yy, N_Vector yp,
                N_Vector yyB, N_Vector ypB, N_Vector rrB, SUNMatrix JB,
                void* user_data, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B)
{
  sunrealtype* yvalB;
  sunindextype* colptrsB = SUNSparseMatrix_IndexPointers(JB);
  sunindextype* rowvalsB = SUNSparseMatrix_IndexValues(JB);
  sunrealtype* dataB     = SUNSparseMatrix_Data(JB);

  UserData userdata;
  sunrealtype p1, p2, p3;

  yvalB = N_VGetArrayPointer(yy);

  userdata = (UserData)user_data;
  p1       = userdata->p[0];
  p2       = userdata->p[1];
  p3       = userdata->p[2];

  SUNMatZero(JB);

  colptrsB[0] = 0;
  colptrsB[1] = 3;
  colptrsB[2] = 6;
  colptrsB[3] = 9;

  /* column 0 */
  dataB[0]    = -p1 + cjB;
  rowvalsB[0] = 0;
  dataB[1]    = p2 * yvalB[2];
  rowvalsB[1] = 1;
  dataB[2]    = p2 * yvalB[1];
  rowvalsB[2] = 2;

  /* column 1 */
  dataB[3]    = p1;
  rowvalsB[3] = 0;
  dataB[4]    = -(p2 * yvalB[2] + SUN_RCONST(2.0) * p3 * yvalB[1]) + cjB;
  rowvalsB[4] = 1;
  dataB[5]    = -p2 * yvalB[1];
  rowvalsB[5] = 2;

  /* column 2 */
  dataB[6]    = -ONE;
  rowvalsB[6] = 0;
  dataB[7]    = -ONE;
  rowvalsB[7] = 1;
  dataB[8]    = -ONE;
  rowvalsB[8] = 2;

  return (0);
}

static int rhsQB(sunrealtype tt, N_Vector yy, N_Vector yp, N_Vector yyB,
                 N_Vector ypB, N_Vector rrQB, void* user_dataB)
{
  sunrealtype y1, y2, y3;
  sunrealtype l1, l2;
  sunrealtype l21;
  sunrealtype* yyB_data  = N_VGetArrayPointer(yyB);
  sunrealtype* yy_data   = N_VGetArrayPointer(yy);
  sunrealtype* rrQB_data = N_VGetArrayPointer(rrQB);

  /* The y vector */
  y1 = yy_data[0];
  y2 = yy_data[1];
  y3 = yy_data[2];

  /* The lambda vector */
  l1 = yyB_data[0];
  l2 = yyB_data[1];

  /* Temporary variables */
  l21 = l2 - l1;

  rrQB_data[0] = y1 * l21;
  rrQB_data[1] = -y3 * y2 * l21;
  rrQB_data[2] = -y2 * y2 * l2;

  return (0);
}

/*
 *--------------------------------------------------------------------
 * PRIVATE FUNCTIONS
 *--------------------------------------------------------------------
 */

/*
 * Print results after backward integration
 */

static void PrintOutput(sunrealtype tfinal, N_Vector yB, N_Vector ypB, N_Vector qB)
{
  sunrealtype* yB_data = N_VGetArrayPointer(yB);
  sunrealtype* qB_data = N_VGetArrayPointer(qB);
  printf("--------------------------------------------------------\n");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("tB0:        %12.4Le\n", tfinal);
  printf("dG/dp:      %12.4Le %12.4Le %12.4Le\n", -qB_data[0], -qB_data[1],
         -qB_data[2]);
  printf("lambda(t0): %12.4Le %12.4Le %12.4Le\n", yB_data[0], yB_data[1],
         yB_data[2]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("tB0:        %12.4e\n", tfinal);
  printf("dG/dp:      %12.4e %12.4e %12.4e\n", -qB_data[0], -qB_data[1],
         -qB_data[2]);
  printf("lambda(t0): %12.4e %12.4e %12.4e\n", yB_data[0], yB_data[1],
         yB_data[2]);
#else
  printf("tB0:        %12.4e\n", tfinal);
  printf("dG/dp:      %12.4e %12.4e %12.4e\n", -qB_data[0], -qB_data[1],
         -qB_data[2]);
  printf("lambda(t0): %12.4e %12.4e %12.4e\n", yB_data[0], yB_data[1],
         yB_data[2]);
#endif
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
