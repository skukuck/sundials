/*-----------------------------------------------------------------
 * Programmer(s): Sylvia Amihere @ UMBC
 *---------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2025, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 *---------------------------------------------------------------
 * Example problem:
 *
 * The following test employs a manufactured true solution to
 * find the accuracy of numerical solution.  We construct a
 * large system of ODEs by replicating a 3 component system
 * n=100 times, for an ODE system with 3n components.  If each of
 * u,v,w are vectors of length n, then the overall solution
 * vector is given by X = [u,v,w], where for each "location" k,
 * the u_k, v_k and w_k satisfy the equations
 *
 *   du_k/dt = v - k*u + (1 + k*t)*exp(i*t)
 *   dv_k/dt = w - t + i*v
 *   dw_k/dt = 1 + i*(w-t)
 *
 * for t in the interval [0.0, 5.0], with initial conditions
 * Y0 = [u0,v0,w0], computed substituting t=0 in the true solution,
 *
 *   u_k(t) = t*exp(i*t) + exp(-k*t)
 *   v_k(t) = i*t*exp(i*t)
 *   w_k(t) = i*exp(i*t) + t
 *
 * We trivially convert this to the DAE system:
 *
 * 0 = du/dt - (v - k*u + (1 + k*t)*exp(i*t))
 * 0 = dv/dt - (w - t + i*v)
 * 0 = dw/dt - (1 + i*(w-t))
 *
 * This program solves the problem with a DIRK method.  Based on
 * a command-line input, it uses either a dense or banded
 * Jacobian matrix and solver, where the Jacobian is computed
 * internally by IDA using finite differences, or it will use
 * the GMRES iterative solver without preconditioning.
 *
 * Specifically, the first command-line argument following the
 * executable name should be an integer:
 *   0 => dense linear solver
 *   1 => banded linear solver (default)
 *   2 => GMRES iterative linear solver (no preconditioning)
 *   3 => GMRES iterative linear solver with IDABBDPRE
 *
 * 5 outputs are printed at equal intervals, and run statistics
 * are printed at the end.
 *-----------------------------------------------------------------*/

/* Header files */
#include <ida/ida.h> /* prototypes for IDA fcts., consts */
#include <math.h>
#include <nvector/nvector_serial.h> /* serial N_Vector types, fcts., macros */
#include <stdio.h>
#include <sundials/sundials_types.h> /* def. of type sunscalartype */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_band.h>  /* access to banded SUNLinearSolver     */
#include <sunmatrix/sunmatrix_band.h>  /* access to banded SUNMatrix           */
#include <sunlinsol/sunlinsol_spgmr.h> /* access to SPGMR SUNLinearSolver      */
#include <ida/ida_bbdpre.h>      /* prototypes for IDABBDPRE module      */

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define GSYM "Lg"
#define ESYM "Le"
#define FSYM "Lf"
#else
#define GSYM "g"
#define ESYM "e"
#define FSYM "f"
#endif

// Problem parameters
#define NX 100               // number of locations
#define NEQ (3 * NX)         // number of equations

// Macro to access variables at a specific location
#define IDX(v,x) ((v) + 3*(x))

/* DAE residual */
static int res(sunrealtype tres, N_Vector yy, N_Vector yp, N_Vector resval,
               void* user_data);
static int resloc(sunindextype Nlocal, sunrealtype t, N_Vector yy, N_Vector yp, 
                  N_Vector rr, void* user_data);

// Compute the true solution and derivative
static int Solution(sunrealtype t, N_Vector u);
static int SolutionDerivative(sunrealtype t, N_Vector up);

// Compute the numerical solution error and derivative error
static int SolutionError(sunrealtype t, N_Vector u, N_Vector e);
static int DerivativeError(sunrealtype t, N_Vector u_p, N_Vector e_p);

// Private function to check function return values
static int check_flag(void* flagvalue, const char* funcname, int opt);

// Main Program
int main(int argc, char* argv[])
{
  /* general problem parameters */
  sunrealtype T0     = SUN_RCONST(0.0);       /* initial time */
  sunrealtype Tf     = SUN_RCONST(5.0);       /* final time */
  sunrealtype dTout  = SUN_RCONST(1.0);       /* time between outputs */
  int Nt             = (int)ceil(Tf / dTout); /* number of output times */
    #if defined(SUNDIALS_SINGLE_PRECISION)
  sunrealtype reltol = SUN_RCONST(1.0e-3);                /* tolerances */
  #else
  sunrealtype reltol = SUN_RCONST(1.0e-6);
  #endif
  sunrealtype abstol = SUN_RCONST(1.0e-10);
  int maxl           = 10; /* max linear solver iterations */

  /* general problem variables */
  int flag;                  /* reusable error-checking flag */
  N_Vector y          = NULL; /* empty vector for storing solution */
  N_Vector y_p        = NULL; /* empty vector for storing solution derivative */
  N_Vector True_Sol   = NULL; /* vector for storing true solution */
  N_Vector True_Sol_p = NULL; /* vector for storing true solution derivative */
  N_Vector Error      = NULL; /* vector for storing the error */
  N_Vector Error_p    = NULL; /* vector for storing the derivative error */
  SUNMatrix A         = NULL; /* empty matrix for linear solver */
  SUNLinearSolver LS  = NULL; /* empty linear solver object */
  void* ida_mem       = NULL; /* empty IDA memory structure */
  sunscalartype rdata[3];
  sunrealtype t, tout;
  int iout;

  /* Retrieve the command-line option specifying the linear solver type:
     0 => dense, 1 => banded (default), 2 => GMRES (no preconditioning), 
     3 => GMRES (BBDPRE) */
  int linear_solver_type = 1;
  if (argc > 1)
  {
    linear_solver_type = atoi(argv[1]);
  }
  if (linear_solver_type < 0 || linear_solver_type > 3)
  {
    fprintf(
      stderr,
      "ERROR: Unrecognized linear solver type %d. Valid options are:\n"
      "  0 => dense linear solver\n"
      "  1 => banded linear solver (default)\n"
      "  2 => GMRES iterative linear solver (no preconditioning)\n"
      "  3 => GMRES iterative linear solver (BBDPRE)\n",
      linear_solver_type);
    return 1;
  }

  /* Create the SUNDIALS context object for this testing */
  SUNContext ctx;
  flag = SUNContext_Create(SUN_COMM_NULL, &ctx);
  if (check_flag(&flag, "SUNContext_Create", 1)) { return 1; }

  y = N_VNew_Serial(NEQ, ctx); /* Create serial vector for solution */
  if (check_flag((void*)y, "N_VNew_Serial", 0)) { return 1; }
  y_p                  = N_VClone(y);
  True_Sol             = N_VClone(y);
  True_Sol_p           = N_VClone(y);
  Error                = N_VClone(y);
  Error_p              = N_VClone(y);
  sunscalartype* yvals = N_VGetArrayPointer(y);
  if (check_flag(yvals, "N_VGetArrayPointer", 0)) { return 1; }
  sunscalartype* ypvals = N_VGetArrayPointer(y_p);
  if (check_flag(ypvals, "N_VGetArrayPointer", 0)) { return 1; }

   // set up the problem data (unused in this case)
  rdata[0] = SUN_CCONST(1.0, 0.0);
  rdata[1] = SUN_CCONST(1.0, 0.0);
  rdata[2] = SUN_CCONST(1.0, 0.0);

  // Set initial condition
  flag = Solution(0.0, True_Sol);
  if (check_flag(&flag, "Solution", 1)) { return 1; }
  flag = SolutionDerivative(0.0, True_Sol_p);
  if (check_flag(&flag, "SolutionDerivative", 1)) { return 1; }

  N_VScale(1.0, True_Sol, y); /* Set initial conditions */
  N_VScale(1.0, True_Sol_p, y_p);

  /* Initial problem output */
  if (linear_solver_type == 0)
  {
    printf("\nAnalytic ODE test in complex arithmetic with dense linear solver:\n");
  }
  else if (linear_solver_type == 1)
  {
    printf("\nAnalytic ODE test in complex arithmetic with banded linear solver:\n");
  }
  else
  {
    printf("\nAnalytic ODE test in complex arithmetic with GMRES iterative linear solver:\n");
    printf("    maxl = %i\n", maxl);
    if (linear_solver_type == 3)
    {
      printf("    preconditioning: IDABBDPRE\n");
    }
    else
    {
      printf("    no preconditioning\n");
    }
  }
  printf("    reltol = %.1" ESYM ",  abstol = %.1" ESYM "\n\n", reltol, abstol);

  /* Call IDACreate and IDAInit to initialize IDA memory */
  ida_mem = IDACreate(ctx);
  if (check_flag((void*)ida_mem, "IDACreate", 0)) { return 1; }

  flag = IDAInit(ida_mem, res, T0, y, y_p);
  if (check_flag(&flag, "IDAInit", 1)) { return (1); }

  /* Set routines */
  flag = IDASetUserData(ida_mem, (void*)rdata); /* Pass rdata to user functions */
  if (check_flag(&flag, "IDASetUserData", 1)) { return 1; }

  flag = IDASStolerances(ida_mem, reltol, abstol); /* Specify tolerances */
  if (check_flag(&flag, "IDASStolerances", 1)) { return 1; }

  /* Initialize matrix data structure and solver */
  if (linear_solver_type == 0)
  {
    A = SUNDenseMatrix(NEQ, NEQ, ctx);
    if (check_flag((void*)A, "SUNDenseMatrix", 0)) { return 1; }
    LS = SUNLinSol_Dense(y, A, ctx);
    if (check_flag((void*)LS, "SUNLinSol_Dense", 0)) { return 1; }
  }
  else if (linear_solver_type == 1)
  {
    A = SUNBandMatrix(NEQ, 2, 2, ctx);
    if (check_flag((void*)A, "SUNBandMatrix", 0)) { return 1; }
    LS = SUNLinSol_Band(y, A, ctx);
    if (check_flag((void*)LS, "SUNLinSol_Band", 0)) { return 1; }
  }
  else if (linear_solver_type == 2)
  {
    LS = SUNLinSol_SPGMR(y, SUN_PREC_NONE, maxl, ctx);
    if (check_flag((void*)LS, "SUNLinSol_SPGMR", 0)) { return 1; }
  }
  else
  {
    LS = SUNLinSol_SPGMR(y, SUN_PREC_RIGHT, maxl, ctx);
    if (check_flag((void*)LS, "SUNLinSol_SPGMR", 0)) { return 1; }
  }

  /* Linear solver interface */
  flag = IDASetLinearSolver(ida_mem, LS, A); /* Attach matrix and linear solver */
  if (check_flag(&flag, "IDASetLinearSolver", 1)) { return 1; }

  /* Set preconditioner if requested */
  if (linear_solver_type == 3)
  {
    flag = IDABBDPrecInit(ida_mem, NEQ, 2, 2, 2, 2, SUN_RCONST(0.0), resloc, NULL);
    if (check_flag(&flag, "IDABBDPrecInit", 1)) { return 1; }
  }

  /* Override any current settings with command-line options */
  flag = IDASetOptions(ida_mem, NULL, NULL, argc, argv);
  if (check_flag(&flag, "IDASetOptions", 1)) { return 1; }

  /* Main time-stepping loop: calls IDASolve to perform the integration, then
     prints results.  Stops when the final time has been reached */
  t    = T0;
  tout = T0 + dTout;
  printf(
    "      t              u_90                        u_90'                       v_90"
    "                        v_90'                       w_90                        "
    "w_90'\n");
  printf("   "
         "---------------------------------------------------------------------"
         "------"
         "---------------------------------------------------------------------"
         "-----------\n");
  printf(" %8.3" FSYM " | "
           "%8.5" FSYM " + "
           "%8.5" FSYM "i  |  "
           "%8.5" FSYM " + "
           "%8.5" FSYM "i  |  "
           "%8.5" FSYM " + "
           "%8.5" FSYM "i  |  "
           "%8.5" FSYM " + "
           "%8.5" FSYM "i  |  "
           "%8.5" FSYM " + "
           "%8.5" FSYM "i  |  "
           "%8.5" FSYM " + "
           "%8.5" FSYM "i  \n",
           t, SUN_REAL(yvals[IDX(0,90)]), SUN_IMAG(yvals[IDX(0,90)]), SUN_REAL(ypvals[IDX(0,90)]), 
           SUN_IMAG(ypvals[IDX(0,90)]), SUN_REAL(yvals[IDX(1,90)]), SUN_IMAG(yvals[IDX(1,90)]),
           SUN_REAL(ypvals[IDX(1,90)]), SUN_IMAG(ypvals[IDX(1,90)]), SUN_REAL(yvals[IDX(2,90)]), 
           SUN_IMAG(yvals[IDX(2,90)]), SUN_REAL(ypvals[IDX(2,90)]), SUN_IMAG(ypvals[IDX(2,90)]));

  for (iout = 0; iout < Nt; iout++)
  {
    flag = IDASolve(ida_mem, tout, &t, y, y_p, IDA_NORMAL); /* call integrator */
    if (check_flag(&flag, "IDASolve", 1)) { break; }
    printf(" %8.3" FSYM " | "
           "%8.5" FSYM " + "
           "%8.5" FSYM "i  |  "
           "%8.5" FSYM " + "
           "%8.5" FSYM "i  |  "
           "%8.5" FSYM " + "
           "%8.5" FSYM "i  |  "
           "%8.5" FSYM " + "
           "%8.5" FSYM "i  |  "
           "%8.5" FSYM " + "
           "%8.5" FSYM "i  |  "
           "%8.5" FSYM " + "
           "%8.5" FSYM "i  \n",
           t, SUN_REAL(yvals[IDX(0,90)]), SUN_IMAG(yvals[IDX(0,90)]), SUN_REAL(ypvals[IDX(0,90)]), 
           SUN_IMAG(ypvals[IDX(0,90)]), SUN_REAL(yvals[IDX(1,90)]), SUN_IMAG(yvals[IDX(1,90)]),
           SUN_REAL(ypvals[IDX(1,90)]), SUN_IMAG(ypvals[IDX(1,90)]), SUN_REAL(yvals[IDX(2,90)]), 
           SUN_IMAG(yvals[IDX(2,90)]), SUN_REAL(ypvals[IDX(2,90)]), SUN_IMAG(ypvals[IDX(2,90)]));
           

    if (flag >= 0)
    { /* successful solve: update time */
      tout += dTout;
      tout = (tout > Tf) ? Tf : tout;
    }
    else
    { /* unsuccessful solve: break */
      fprintf(stderr, "Solver failure, stopping integration\n");
      break;
    }
  }
  printf("   "
         "---------------------------------------------------------------------"
         "------"
         "---------------------------------------------------------------------"
         "-----------\n");

  SolutionError(Tf, y, Error);
  DerivativeError(Tf, y_p, Error_p);

  /* Print all final statistics */
  printf("\nFinal Solver Statistics:\n");
  flag = IDAPrintAllStats(ida_mem, stdout, SUN_OUTPUTFORMAT_TABLE);
  if (check_flag(&flag, "IDAPrintAllStats", 1)) { return 1; }

  /* Clean up and return with successful completion */
  N_VDestroy(y); /* Free vectors */
  N_VDestroy(y_p); 
  N_VDestroy(True_Sol);
  N_VDestroy(True_Sol_p);
  N_VDestroy(Error);
  N_VDestroy(Error_p);
  IDAFree(&ida_mem); /* Free integrator memory */
  SUNLinSolFree(LS);       /* Free linear solver */
  SUNMatDestroy(A);        /* Free A matrix */
  SUNContext_Free(&ctx);   /* Free context */

  return 0;
}

/*-------------------------------
 * Functions called by the solver
 *-------------------------------*/

/* f routine to compute the ODE RHS function f(t,y)
      du_k/dt = v - k*u + (1 + k*t)*exp(i*t)
      dv_k/dt = w - t + i*v
      dw_k/dt = 1 + i*(w-t)  */

static int res(sunrealtype t, N_Vector yy, N_Vector yp, N_Vector rr, void* user_data)
{
  sunscalartype* rdata =
    (sunscalartype*)user_data; /* cast user_data to sunscalartype */
  sunscalartype* yvals   = N_VGetArrayPointer(yy);
  sunscalartype* ypvals  = N_VGetArrayPointer(yp);
  sunscalartype* resvals = N_VGetArrayPointer(rr);
  for (sunindextype k=0; k<NX; k++) {
    sunscalartype a        = rdata[0]; /* access data entries (unused) */
    sunscalartype b        = rdata[1];
    sunscalartype c        = rdata[2];
    sunscalartype u   = yvals[IDX(0,k)]; /* access solution values */
    sunscalartype v   = yvals[IDX(1,k)];
    sunscalartype w   = yvals[IDX(2,k)];
    sunscalartype u_p = ypvals[IDX(0,k)]; /* access derivative values */
    sunscalartype v_p = ypvals[IDX(1,k)];
    sunscalartype w_p = ypvals[IDX(2,k)];
    const sunscalartype ONE = SUN_RCONST(1.0);
    const sunscalartype TWO = SUN_RCONST(2.0);
    const sunscalartype K = (sunscalartype)k;

    /* fill in the RHS function */
    resvals[IDX(0,k)] = u_p - v + k * u - (ONE + K * t)*SUNexp(SUN_I * t);
    resvals[IDX(1,k)] = v_p - w + t - SUN_I * v;
    resvals[IDX(2,k)] = w_p - ONE - SUN_I * (w - t);
  }

  return 0; /* Return with success */
}

/* f routine to compute the "local" portion of the ODE RHS function f(t,y)
   (for use by the IDABBDPRE module) */
static int resloc(sunindextype Nlocal, sunrealtype t, N_Vector yy, N_Vector yp, N_Vector rr, void* user_data)
{
   sunscalartype* rdata =
    (sunscalartype*)user_data; /* cast user_data to sunscalartype */
    return res(t, yy, yp, rr, user_data);
}

/*-------------------------------
 * Private helper functions
 *-------------------------------*/

/* Check function return value...
    opt == 0 means SUNDIALS function allocates memory so check if
             returned NULL pointer
    opt == 1 means SUNDIALS function returns a flag so check if
             flag >= 0
    opt == 2 means function allocates memory so check if returned
             NULL pointer
*/
static int check_flag(void* flagvalue, const char* funcname, int opt)
{
  int* errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL)
  {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return 1;
  }

  /* Check if flag < 0 */
  else if (opt == 1)
  {
    errflag = (int*)flagvalue;
    if (*errflag < 0)
    {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
              funcname, *errflag);
      return 1;
    }
  }

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL)
  {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return 1;
  }

  return 0;
}

/* Compute the exact solution
    u_k(t) = t*exp(i*t) + exp(-k*t)
    v_k(t) = i*t*exp(i*t)
    w_k(t) = i*exp(i*t) + t  */
static int Solution(sunrealtype t, N_Vector u)
{
  sunscalartype* uarray = N_VGetArrayPointer(u);
  if (check_flag((void*)uarray, "N_VGetArrayPointer", 0)) { return -1; }

  for (sunindextype k=0; k<NX; k++) {
    sunscalartype K = (sunscalartype)k;
    uarray[IDX(0,k)] = t * SUNexp(SUN_I * t) + SUNexp(-K * t);
    uarray[IDX(1,k)] = t * SUN_I * SUNexp(SUN_I * t);
    uarray[IDX(2,k)] = SUN_I * SUNexp(SUN_I * t) + t;
  }

  return 0;
}

/* Compute the derivative of the exact solution */
static int SolutionDerivative(sunrealtype t, N_Vector u_p)
{
  sunscalartype* uparray = N_VGetArrayPointer(u_p);
  if (check_flag((void*)uparray, "N_VGetArrayPointer", 0)) { return -1; }

  for (sunindextype k=0; k<NX; k++) {
    sunscalartype K = (sunscalartype)k;
    uparray[IDX(0,k)] = SUNexp(SUN_I * t) + (t * SUN_I * SUNexp(SUN_I * t)) -
                        K * SUNexp(-K * t);
    uparray[IDX(1,k)] = SUN_I * SUNexp(SUN_I * t) - t * SUNexp(SUN_I * t);
    uparray[IDX(2,k)] = SUN_RCONST(1.0) - SUNexp(SUN_I * t);
  }

  return 0;
}

/* Compute the solution error */
static int SolutionError(sunrealtype t, N_Vector u, N_Vector e)
{
  /* Compute true solution */
  int flag = Solution(t, e);
  if (flag != 0) { return -1; }

  /* Compute max-norm of the error */
  N_VLinearSum(1.0, u, -1.0, e, e);
  sunrealtype error_norm = N_VMaxNorm(e);
  printf("    Max-norm of the error is %.5" ESYM "\n", error_norm);
  return 0;
}

/* Compute the solution derivative error */
static int DerivativeError(sunrealtype t, N_Vector u_p, N_Vector e_p)
{
  /* Compute true solution */
  int flag = SolutionDerivative(t, e_p);
  if (flag != 0) { return -1; }

  /* Compute max-norm of the error */
  N_VLinearSum(1.0, u_p, -1.0, e_p, e_p);
  sunrealtype error_norm = N_VMaxNorm(e_p);
  printf("    Max-norm of the derivative error is %.5" ESYM "\n", error_norm);
  return 0;
}

/*---- end of file ----*/
