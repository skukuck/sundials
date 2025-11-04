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
from sundials4py.core import *
from sundials4py.idas import *
from problems import AnalyticDAE
from fixtures import *


@pytest.mark.skipif(
    sunrealtype == np.float32, reason="Test not supported for sunrealtype=np.float32"
)
def test_idas_ivp(sunctx):
    ode_problem = AnalyticDAE()

    solver = IDACreate(sunctx)
    yy = N_VNew_Serial(2, sunctx)
    yp = N_VNew_Serial(2, sunctx)

    # y and y' initial conditions
    ode_problem.set_init_cond(yy, yp, ode_problem.T0)

    ls = SUNLinSol_SPGMR(yy, SUN_PREC_LEFT, 0, sunctx)

    def resfn(t, yy, yp, rr, _):
        return ode_problem.res(t, yy, yp, rr)

    def psolve(t, yy, yp, rr, r, z, cj, delta, _):
        return ode_problem.psolve(t, yy, yp, rr, r, z, cj, delta)

    status = IDAInit(solver.get(), resfn, 0.0, yy, yp)
    assert status == IDA_SUCCESS

    status = IDASStolerances(solver.get(), 1e-4, 1e-4)
    assert status == IDA_SUCCESS

    status = IDASetLinearSolver(solver.get(), ls, None)
    assert status == IDA_SUCCESS

    status = IDASetPreconditioner(solver.get(), None, psolve)
    assert status == IDA_SUCCESS

    # TODO(CJB): enable rootfinding once we have nrtfn as a callback argument
    # nrtfn = 2
    # def rootfn(t, yy, yp, gout, _):
    #     # just a smoke test of the root finding callback
    #     assert len(gout) == nrtfn
    #     return 0

    # status = IDARootInit(solver.get(), nrtfn, rootfn)
    # assert status == IDA_SUCCESS

    tout = ode_problem.TF
    status, tret = IDASolve(solver.get(), tout, yy, yp, IDA_NORMAL)
    assert status == IDA_SUCCESS

    status, num_steps = IDAGetNumSteps(solver.get())
    assert status == IDA_SUCCESS
    print("Number of steps: ", num_steps)

    sol_yy = N_VClone(yy)
    sol_yp = N_VClone(yp)

    ode_problem.solution(sol_yy, sol_yp, tret)
    assert np.allclose(N_VGetArrayPointer(sol_yy), N_VGetArrayPointer(yy), rtol=1e-2)
    assert np.allclose(N_VGetArrayPointer(sol_yp), N_VGetArrayPointer(yp), rtol=1e-2)


def test_idas_fsa(sunctx):
    # Forward Sensitivity Analysis (FSA) with respect to initial condition for DAE
    ode_problem = AnalyticDAE()

    solver = IDACreate(sunctx)
    yy = N_VNew_Serial(2, sunctx)
    yp = N_VNew_Serial(2, sunctx)

    # y and y' initial conditions
    ode_problem.set_init_cond(yy, yp, ode_problem.T0)

    ls = SUNLinSol_SPGMR(yy, SUN_PREC_LEFT, 0, sunctx)

    def resfn(t, yy, yp, rr, _):
        return ode_problem.res(t, yy, yp, rr)

    def psolve(t, yy, yp, rr, r, z, cj, delta, _):
        return ode_problem.psolve(t, yy, yp, rr, r, z, cj, delta)

    status = IDAInit(solver.get(), resfn, 0.0, yy, yp)
    assert status == IDA_SUCCESS

    status = IDASStolerances(solver.get(), 1e-4, 1e-4)
    assert status == IDA_SUCCESS

    status = IDASetLinearSolver(solver.get(), ls, None)
    assert status == IDA_SUCCESS

    status = IDASetPreconditioner(solver.get(), None, psolve)
    assert status == IDA_SUCCESS

    # Sensitivity setup
    Ns = 1  # number of sensitivities (wrt yy0)
    ism = IDA_SIMULTANEOUS
    yyS0 = [N_VClone(yy) for _ in range(Ns)]
    ypS0 = [N_VClone(yp) for _ in range(Ns)]
    for v in yyS0:
        N_VConst(0.0, v)
    for v in ypS0:
        N_VConst(0.0, v)

    def resS(Ns, t, yy, yp, resval, yS, ypS, resvalS, _, tmp1, tmp2, tmp3):
        # Sensitivity residuals: d(res)/d(yy) * yyS + d(res)/d(yp) * ypS
        # For this analytic DAE, the Jacobians are known, but for smoke test, just zero
        # print(type(resvalS))
        for i in range(Ns):
            N_VConst(0.0, resvalS[i])
        return 0

    status = IDASensInit(solver.get(), Ns, ism, resS, yyS0, ypS0)
    assert status == IDA_SUCCESS

    status = IDASensSStolerances(solver.get(), 1e-4, np.array([1e-4], dtype=sunrealtype))
    assert status == IDA_SUCCESS

    tout = ode_problem.TF
    yySout = [N_VClone(yy) for _ in range(Ns)]
    ypSout = [N_VClone(yp) for _ in range(Ns)]
    status, tret = IDASolve(solver.get(), tout, yy, yp, IDA_NORMAL)
    assert status == IDA_SUCCESS

    # Get sensitivities (smoke test: just check call and shape)
    status, tret_sens = IDAGetSens(solver.get(), yySout)
    assert status == IDA_SUCCESS
    assert len(yySout) == Ns


def test_idas_adjoint(sunctx):
    # TODO(CJB): implement
    pass
