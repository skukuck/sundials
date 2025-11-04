#!/usr/bin/env python3
# -----------------------------------------------------------------
# Programmer(s): Cody J. Balos
# -----------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2025, Lawrence Livermore National Security,
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
# This is a direct port of the C example,
#   examples/kinsol/serial/kinLaplace_bnd.c
#
# This example solves a 2D elliptic PDE
#
#    d^2 u / dx^2 + d^2 u / dy^2 = u^3 - u - 2.0
#
# subject to homogeneous Dirichlet boundary conditions.
# The PDE is discretized on a uniform NX+2 by NY+2 grid with
# central differencing, and with boundary values eliminated,
# leaving a system of size NEQ = NX*NY.
# The nonlinear system is solved by KINSOL using the SUNBAND linear
# solver.
# -----------------------------------------------------------------

import numpy as np
from sundials4py.core import *
from sundials4py.kinsol import *

# Problem Constants
NX = 31
NY = 31
NEQ = NX * NY
SKIP = 3
FTOL = 1e-12
ZERO = 0.0
ONE = 1.0
TWO = 2.0


# Helper to get index in y for (i, j)
def idx(i, j):
    return (j - 1) + (i - 1) * NY


class Laplace2D:
    def __init__(self, NX, NY):
        self.NX = NX
        self.NY = NY
        self.NEQ = NX * NY
        self.dx = ONE / (NX + 1)
        self.dy = ONE / (NY + 1)
        self.hdc = ONE / (self.dx * self.dx)
        self.vdc = ONE / (self.dy * self.dy)

    def set_init_cond(self, yvec):
        y = N_VGetArrayPointer(yvec)
        # Initial guess: zero everywhere
        for i in range(1, self.NX + 1):
            for j in range(1, self.NY + 1):
                y[idx(i, j)] = ZERO
        return 0

    def func(self, uvec, fvec):
        u = N_VGetArrayPointer(uvec)
        f = N_VGetArrayPointer(fvec)
        # Reshape to 2D for vectorized operations
        u2d = np.zeros((self.NX + 2, self.NY + 2), dtype=sunrealtype)
        # Fill interior points
        for i in range(1, self.NX + 1):
            for j in range(1, self.NY + 1):
                u2d[i, j] = u[idx(i, j)]
        # Vectorized finite difference
        uij = u2d[1:-1, 1:-1]
        udn = u2d[1:-1, 0:-2]
        uup = u2d[1:-1, 2:]
        ult = u2d[0:-2, 1:-1]
        urt = u2d[2:, 1:-1]
        hdiff = self.hdc * (ult - TWO * uij + urt)
        vdiff = self.vdc * (uup - TWO * uij + udn)
        f2d = hdiff + vdiff + uij - uij * uij * uij + 2.0
        # Write back to 1D f
        for i in range(1, self.NX + 1):
            for j in range(1, self.NY + 1):
                f[idx(i, j)] = f2d[i - 1, j - 1]
        return 0

    def print_output(self, yvec):
        y = N_VGetArrayPointer(yvec)
        print("            ", end="")
        for i in range(1, self.NX + 1, SKIP):
            x = i * self.dx
            print(f"{x:<8.5f} ", end="")
        print("\n")
        for j in range(1, self.NY + 1, SKIP):
            yval = j * self.dy
            print(f"{yval:<8.5f}    ", end="")
            for i in range(1, self.NX + 1, SKIP):
                print(f"{y[idx(i, j)]:<8.5f} ", end="")
            print()


def main():
    print("\n2D elliptic PDE on unit square")
    print("   d^2 u / dx^2 + d^2 u / dy^2 = u^3 - u + 2.0")
    print(" + homogeneous Dirichlet boundary conditions\n")
    print("Solution method: Modified Newton with band linear solver")
    print(f"Problem size: {NX} x {NY} = {NEQ}\n")

    status, sunctx = SUNContext_Create(SUN_COMM_NULL)
    y = N_VNew_Serial(NEQ, sunctx)
    scale = N_VNew_Serial(NEQ, sunctx)

    # Create problem instance and set initial guess
    problem = Laplace2D(NX, NY)
    problem.set_init_cond(y)

    # Initialize and allocate memory for KINSOL
    kin = KINCreate(sunctx)
    status = KINInit(kin.get(), lambda u, f, _: problem.func(u, f), y)
    assert status == KIN_SUCCESS

    # Set function norm tolerance
    status = KINSetFuncNormTol(kin.get(), FTOL)
    assert status == KIN_SUCCESS

    # Create band matrix and linear solver
    J = SUNBandMatrix(NEQ, NX, NX, sunctx)
    LS = SUNLinSol_Band(y, J, sunctx)
    status = KINSetLinearSolver(kin.get(), LS, J)
    assert status == KIN_SUCCESS

    # Set Modified Newton parameters
    status = KINSetMaxSetupCalls(kin.get(), 100)
    assert status == KIN_SUCCESS
    status = KINSetMaxSubSetupCalls(kin.get(), 1)
    assert status == KIN_SUCCESS

    # No scaling used
    N_VConst(ONE, scale)

    # Call KINSol to solve problem
    status = KINSol(kin.get(), y, KIN_LINESEARCH, scale, scale)
    assert status == KIN_SUCCESS

    # Get scaled norm of the system function
    status, fnorm = KINGetFuncNorm(kin.get())
    assert status == KIN_SUCCESS
    print(f"\nComputed solution (||F|| = {fnorm}):\n")
    problem.print_output(y)

    # Print final statistics (faithful to C PrintFinalStats)
    status, nni = KINGetNumNonlinSolvIters(kin.get())
    assert status == KIN_SUCCESS

    status, nfe = KINGetNumFuncEvals(kin.get())
    assert status == KIN_SUCCESS

    status, nbcfails = KINGetNumBetaCondFails(kin.get())
    assert status == KIN_SUCCESS

    status, nbacktr = KINGetNumBacktrackOps(kin.get())
    assert status == KIN_SUCCESS

    status, nje = KINGetNumJacEvals(kin.get())
    assert status == KIN_SUCCESS

    status, nfeD = KINGetNumLinFuncEvals(kin.get())
    assert status == KIN_SUCCESS

    print("\nFinal Statistics.. \n")
    print(f"nni      = {nni:6d}    nfe     = {nfe:6d}")
    print(f"nbcfails = {nbcfails:6d}    nbacktr = {nbacktr:6d}")
    print(f"nje      = {nje:6d}    nfeB    = {nfeD:6d}")


def test_kinLaplace_bnd():
    main()


if __name__ == "__main__":
    main()
