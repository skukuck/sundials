#!/bin/python
# -----------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
# -----------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2025-2025, Lawrence Livermore National Security,
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
# -----------------------------------------------------------------

import pytest
import numpy as np
from fixtures import *
from sundials4py.core import *
from sundials4py.cvodes import *
from problems import AnalyticODE


@pytest.mark.skipif(
    sunrealtype == np.float32, reason="Test not supported for sunrealtype=np.float32"
)
def test_cvodes_ivp(sunctx):
    NEQ = 1
    y = N_VNew_Serial(NEQ, sunctx)
    ls = SUNLinSol_SPGMR(y, 0, 0, sunctx)

    ode_problem = AnalyticODE()
    ode_problem.set_init_cond(y)

    cvode = CVodeCreate(CV_BDF, sunctx)

    status = CVodeInit(cvode.get(), lambda t, y, ydot, _: ode_problem.f(t, y, ydot), 0, y)
    assert status == CV_SUCCESS

    status = CVodeSStolerances(cvode.get(), SUNREALTYPE_RTOL, SUNREALTYPE_ATOL)
    assert status == CV_SUCCESS

    status = CVodeSetLinearSolver(cvode.get(), ls, None)
    assert status == CV_SUCCESS

    nrtfn = 2

    def rootfn(t, y, gout, _):
        # just a smoke test of the root finding callback
        assert len(gout) == nrtfn
        return 0

    status = CVodeRootInit(cvode.get(), nrtfn, rootfn)
    assert status == CV_SUCCESS

    tout = 10.0
    status, tret = CVode(cvode.get(), tout, y, CV_NORMAL)
    assert status == CV_SUCCESS

    sol = N_VClone(y)
    ode_problem.solution(y, sol, tret)
    assert np.allclose(N_VGetArrayPointer(sol), N_VGetArrayPointer(y), atol=1e-2)

    status, num_steps = CVodeGetNumSteps(cvode.get())
    assert status == CV_SUCCESS
    assert num_steps > 0


@pytest.mark.skipif(
    sunrealtype == np.float32, reason="Test not supported for sunrealtype=np.float32"
)
def test_cvodes_fsa(sunctx):
    # Forward Sensitivity Analysis (FSA) with respect to initial condition
    NEQ = 1
    y = N_VNew_Serial(NEQ, sunctx)
    ls = SUNLinSol_SPGMR(y, 0, 0, sunctx)

    ode_problem = AnalyticODE()
    ode_problem.set_init_cond(y)

    cvode = CVodeCreate(CV_BDF, sunctx)

    # This problem requires tighter tolerances in order to get the forward
    # sensitivity to converge to the expected solution within a reasonable tolerance
    atol = 1e-10
    rtol = 1e-10

    status = CVodeInit(cvode.get(), lambda t, y, ydot, _: ode_problem.f(t, y, ydot), 0, y)
    assert status == CV_SUCCESS

    status = CVodeSStolerances(cvode.get(), rtol, atol)
    assert status == CV_SUCCESS

    status = CVodeSetLinearSolver(cvode.get(), ls, None)
    assert status == CV_SUCCESS

    # Sensitivity setup
    Ns = 1  # number of sensitivities (wrt y0)
    ism = CV_SIMULTANEOUS
    yS0 = [N_VClone(y)]
    N_VConst(1.0, yS0[0])

    def fS(Ns, t, y, ydot, yS, ySdot, _, tmp1, tmp2):
        # Sensitivity RHS: df/dy * yS + df/dp (here, p = y0, so df/dp = 0)
        yarr = N_VGetArrayPointer(y)
        ySarr = N_VGetArrayPointer(yS[0])
        ySdotarr = N_VGetArrayPointer(ySdot[0])
        # df/dy = lambda
        lamb = ode_problem.lamb
        ySdotarr[0] = lamb * ySarr[0]
        return 0

    status = CVodeSensInit(cvode.get(), Ns, ism, fS, yS0)
    assert status == CV_SUCCESS

    status = CVodeSensSStolerances(cvode.get(), rtol, np.array([atol], dtype=sunrealtype))
    assert status == CV_SUCCESS

    status = CVodeSetMaxNumSteps(cvode.get(), 100000)
    assert status == CV_SUCCESS

    tout = 10.0
    ySout = [N_VClone(y) for _ in range(Ns)]
    status, tret = CVode(cvode.get(), tout, y, CV_NORMAL)
    assert status == CV_SUCCESS

    # Check IVP solution
    sol = N_VClone(y)
    ode_problem.solution(y, sol, tret)
    assert np.allclose(N_VGetArrayPointer(sol), N_VGetArrayPointer(y), atol=1e-2)

    # Get sensitivities
    status, tret_sens = CVodeGetSens(cvode.get(), ySout)
    assert status == CV_SUCCESS

    lamb = ode_problem.lamb
    expected = np.exp(lamb * tret)
    sens_val = N_VGetArrayPointer(ySout[0])[0]
    assert np.allclose(sens_val, expected, atol=1e-2)


@pytest.mark.skipif(
    sunrealtype == np.float32, reason="Test not supported for sunrealtype=np.float32"
)
def test_cvodes_adjoint(sunctx):
    # Adjoint Sensitivity Analysis (ASA) for the same ODE problem as FSA
    y = N_VNew_Serial(1, sunctx)
    ls = SUNLinSol_SPGMR(y, 0, 0, sunctx)

    ode_problem = AnalyticODE()
    ode_problem.set_init_cond(y)

    cvode = CVodeCreate(CV_BDF, sunctx)

    status = CVodeInit(cvode.get(), lambda t, y, ydot, _: ode_problem.f(t, y, ydot), 0, y)
    assert status == CV_SUCCESS

    status = CVodeSStolerances(cvode.get(), SUNREALTYPE_RTOL, SUNREALTYPE_ATOL)
    assert status == CV_SUCCESS

    status = CVodeSetLinearSolver(cvode.get(), ls, None)
    assert status == CV_SUCCESS

    # Adjoint setup
    steps = 1000
    interp = CV_HERMITE
    status = CVodeAdjInit(cvode.get(), steps, interp)
    assert status == CV_SUCCESS

    # Forward solve
    tout = 10.0
    status, tret, ncheck = CVodeF(cvode.get(), tout, y, CV_NORMAL)
    assert status == CV_SUCCESS

    # Define a simple functional: g = y(T)
    # The gradient dg/dy0 = dy(T)/dy0, which is exp(lambda*T)

    # Backward problem setup
    yB = N_VNew_Serial(1, sunctx)
    N_VConst(1.0, yB)
    lsB = SUNLinSol_SPGMR(yB, 0, 0, sunctx)

    status, whichB = CVodeCreateB(cvode.get(), CV_BDF)
    assert status == CV_SUCCESS

    def fB(t, y, yB, yBdot, _):
        # Adjoint RHS: -df/dy^T * lambdaB
        yarr = N_VGetArrayPointer(y)
        yBarr = N_VGetArrayPointer(yB)
        yBdotarr = N_VGetArrayPointer(yBdot)
        lamb = ode_problem.lamb
        yBdotarr[0] = -lamb * yBarr[0]
        return 0

    tB0 = tret
    status = CVodeInitB(cvode.get(), whichB, fB, tB0, yB)
    assert status == CV_SUCCESS

    status = CVodeSStolerancesB(cvode.get(), whichB, SUNREALTYPE_RTOL, SUNREALTYPE_ATOL)
    assert status == CV_SUCCESS

    status = CVodeSetLinearSolverB(cvode.get(), whichB, lsB, None)
    assert status == CV_SUCCESS

    # Set terminal condition for lambda(T) = dG/dy(T)
    yBarr = N_VGetArrayPointer(yB)
    yBarr[0] = 1.0

    # Integrate backward
    tBout = 0.0
    status = CVodeB(cvode.get(), tBout, CV_NORMAL)
    assert status >= CV_SUCCESS

    # Get lambda(0) = dg/dy0
    yB0 = N_VNew_Serial(1, sunctx)
    status, tBret = CVodeGetB(cvode.get(), whichB, yB0)
    assert status == CV_SUCCESS

    # Analytical result: exp(lambda * T)
    lamb = ode_problem.lamb
    expected = np.exp(lamb * tret)
    sens_val = N_VGetArrayPointer(yB0)[0]
    assert np.allclose(sens_val, expected, atol=1e-2)


@pytest.mark.skipif(
    sunrealtype == np.float32, reason="Test not supported for sunrealtype=np.float32"
)
def test_cvodes_adjoint_quad(sunctx):
    # Adjoint Sensitivity Analysis (ASA) for the same ODE problem as before but add quadratures
    NEQ = 1

    y = N_VNew_Serial(NEQ, sunctx)
    ls = SUNLinSol_SPGMR(y, 0, 0, sunctx)

    ode_problem = AnalyticODE()
    ode_problem.set_init_cond(y)

    # This problem requires tighter tolerances in order to get the forward
    # sensitivity to converge to the expected solution within a reasonable tolerance
    atol = 1e-10
    rtol = 1e-10

    cvode = CVodeCreate(CV_BDF, sunctx)

    status = CVodeInit(cvode.get(), lambda t, y, ydot, _: ode_problem.f(t, y, ydot), 0, y)
    assert status == CV_SUCCESS

    status = CVodeSStolerances(cvode.get(), rtol, atol)
    assert status == CV_SUCCESS

    status = CVodeSetLinearSolver(cvode.get(), ls, None)
    assert status == CV_SUCCESS

    # Forward sensitivity setup
    Ns = 1  # number of sensitivities (wrt y0)
    ism = CV_SIMULTANEOUS
    yS0 = [N_VClone(y)]
    N_VConst(1.0, yS0[0])

    def fS(Ns, t, y, ydot, yS, ySdot, _, tmp1, tmp2):
        # Sensitivity RHS: df/dy * yS + df/dp (here, p = y0, so df/dp = 0)
        yarr = N_VGetArrayPointer(y)
        ySarr = N_VGetArrayPointer(yS[0])
        ySdotarr = N_VGetArrayPointer(ySdot[0])
        # df/dy = lambda
        lamb = ode_problem.lamb
        ySdotarr[0] = lamb * ySarr[0]
        return 0

    status = CVodeSensInit(cvode.get(), Ns, ism, fS, yS0)
    assert status == CV_SUCCESS

    status = CVodeSensSStolerances(cvode.get(), rtol, np.array([atol], dtype=sunrealtype))
    assert status == CV_SUCCESS

    # Forward quadrature setup
    def fQ(t, y, qdot, _):
        # Smoke test
        qdotarr = N_VGetArrayPointer(qdot)
        qdotarr[:] = 0
        return 0

    Nq = 2
    yQ = N_VNew_Serial(Nq, sunctx)
    N_VConst(0.0, yQ)

    status = CVodeQuadInit(cvode.get(), fQ, yQ)
    assert status == CV_SUCCESS

    status = CVodeQuadSStolerances(cvode.get(), rtol, atol)
    assert status == CV_SUCCESS

    # Adjoint setup
    steps = 1000
    interp = CV_HERMITE
    status = CVodeAdjInit(cvode.get(), steps, interp)
    assert status == CV_SUCCESS

    # Forward solve
    tout = 10.0
    status, tret, ncheck = CVodeF(cvode.get(), tout, y, CV_NORMAL)
    assert status == CV_SUCCESS

    # Define a simple functional: g = y(T)
    # The gradient dg/dy0 = dy(T)/dy0, which is exp(lambda*T)

    # Backward problem setup
    Ns = 1
    yB = N_VNew_Serial(NEQ, sunctx)
    yQB = N_VNew_Serial(Ns, sunctx)
    N_VConst(1.0, yB)
    lsB = SUNLinSol_SPGMR(yB, 0, 0, sunctx)

    status, whichB = CVodeCreateB(cvode.get(), CV_BDF)
    assert status == CV_SUCCESS

    def fBS(t, y, yS, yB, yBdot, _):
        # Adjoint RHS: -df/dy^T * lambdaB
        yarr = N_VGetArrayPointer(y)
        yBarr = N_VGetArrayPointer(yB)
        yBdotarr = N_VGetArrayPointer(yBdot)
        lamb = ode_problem.lamb
        yBdotarr[0] = -lamb * yBarr[0]
        return 0

    def fQB(t, y, yS, yB, qBdot, _):
        # Smoke test
        assert len(yS) == Ns
        return 0

    tB0 = tret
    status = CVodeInitBS(cvode.get(), whichB, fBS, tB0, yB)
    assert status == CV_SUCCESS

    # Setup backward quadratures
    status = CVodeQuadInitBS(cvode.get(), whichB, fQB, yQB)

    status = CVodeSStolerancesB(cvode.get(), whichB, rtol, atol)
    assert status == CV_SUCCESS

    status = CVodeSetLinearSolverB(cvode.get(), whichB, lsB, None)
    assert status == CV_SUCCESS

    # Set terminal condition for lambda(T) = dG/dy(T)
    yBarr = N_VGetArrayPointer(yB)
    yBarr[0] = 1.0

    # Integrate backward
    tBout = 0.0
    status = CVodeB(cvode.get(), tBout, CV_NORMAL)
    assert status >= CV_SUCCESS

    # Get lambda(0) = dg/dy0
    yB0 = N_VNew_Serial(NEQ, sunctx)
    status, tBret = CVodeGetB(cvode.get(), whichB, yB0)
    assert status == CV_SUCCESS

    # Analytical result: exp(lambda * T)
    lamb = ode_problem.lamb
    expected = np.exp(lamb * tret)
    sens_val = N_VGetArrayPointer(yB0)[0]
    assert np.allclose(sens_val, expected, atol=1e-2)
