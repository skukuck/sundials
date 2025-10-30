#!/usr/bin/env python3
# -----------------------------------------------------------------
# Programmer(s): Cody J. Balos
# -----------------------------------------------------------------
# Python port of the SUNDIALS ARKODE Lotka-Volterra adjoint sensitivity example
# (ark_lotka_volterra_ASA.c)
# -----------------------------------------------------------------
import numpy as np
import sundials4py.core as sun
import sundials4py.arkode as ark


class LotkaVolterraODE:
    def __init__(self, p):
        self.p = np.array(p, dtype=sun.sunrealtype)
        self.NEQ = 2
        self.NP = 4

    def set_init_cond(self, yvec):
        y = sun.N_VGetArrayPointer(yvec)
        y[0] = 1.0
        y[1] = 1.0
        return 0

    def f(self, t, yvec, ydotvec):
        p = self.p
        y = sun.N_VGetArrayPointer(yvec)
        ydot = sun.N_VGetArrayPointer(ydotvec)
        ydot[0] = p[0] * y[0] - p[1] * y[0] * y[1]
        ydot[1] = -p[2] * y[1] + p[3] * y[0] * y[1]
        return 0

    def vjp(self, vvec, Jvvec, t, yvec):
        p = self.p
        v = sun.N_VGetArrayPointer(vvec)
        Jv = sun.N_VGetArrayPointer(Jvvec)
        y = sun.N_VGetArrayPointer(yvec)
        Jv[0] = (p[0] - p[1] * y[1]) * v[0] + p[3] * y[1] * v[1]
        Jv[1] = -p[1] * y[0] * v[0] + (-p[2] + p[3] * y[0]) * v[1]
        return 0

    def parameter_vjp(self, vvec, Jvvec, t, yvec):
        v = sun.N_VGetArrayPointer(vvec)
        Jv = sun.N_VGetArrayPointer(Jvvec)
        y = sun.N_VGetArrayPointer(yvec)
        Jv[0] = y[0] * v[0]
        Jv[1] = -y[0] * y[1] * v[0]
        Jv[2] = -y[1] * v[1]
        Jv[3] = y[0] * y[1] * v[1]
        return 0

    def dgdu(self, yvec):
        y = sun.N_VGetArrayPointer(yvec)
        return np.array([-1.0 + y[0], -1.0 + y[1]], dtype=sun.sunrealtype)

    def dgdp(self, yvec):
        return np.zeros(self.NP, dtype=sun.sunrealtype)

    def adj_rhs(self, t, y, sens, sens_dot):
        l = sun.N_VGetSubvector_ManyVector(sens, 0)
        ldot = sun.N_VGetSubvector_ManyVector(sens_dot, 0)
        nu = sun.N_VGetSubvector_ManyVector(sens_dot, 1)
        self.vjp(l, ldot, t, y)
        self.parameter_vjp(l, nu, t, y)
        return 0

    def quad_rhs(self, t, yvec, muvec, qBdotvec):
        self.parameter_vjp(muvec, qBdotvec, t, yvec)
        return 0


def main():
    # Program args
    tf = 10.0
    dt = 1e-3
    order = 4
    check_freq = 2
    keep_checks = True

    # Problem parameters
    p = [1.5, 1.0, 3.0, 1.0]
    t0 = 0.0
    reltol = 1e-10
    abstol = 1e-14
    ode = LotkaVolterraODE(p)
    NEQ = ode.NEQ
    NP = ode.NP

    #
    # Create the initial conditions vector
    #
    status, sunctx = sun.SUNContext_Create(sun.SUN_COMM_NULL)
    assert status == sun.SUN_SUCCESS

    y = sun.N_VNew_Serial(NEQ, sunctx)
    ode.set_init_cond(y)

    #
    # Create the ARKODE stepper that will be used for the forward evolution.
    #
    arkode = ark.ARKodeView.Create(
        ark.ARKStepCreate(lambda t, y, ydot, _: ode.f(t, y, ydot), None, t0, y, sunctx)
    )
    status = ark.ARKodeSetOrder(arkode.get(), 4)
    assert status == ark.ARK_SUCCESS
    status = ark.ARKodeSStolerances(arkode.get(), reltol, abstol)
    assert status == ark.ARK_SUCCESS
    status = ark.ARKodeSetFixedStep(arkode.get(), dt)
    assert status == ark.ARK_SUCCESS

    # Due to roundoff in the `t` accumulation within the integrator,
    # the integrator may actually use nsteps + 1 time steps to reach tf
    status = ark.ARKodeSetMaxNumSteps(arkode.get(), int((tf - t0) / dt) + 1)
    assert status == ark.ARK_SUCCESS

    # # Enable checkpointing during the forward run
    nsteps = int(np.ceil((tf - t0) / dt))
    ncheck = nsteps * order
    mem_helper = sun.SUNMemoryHelper_Sys(sunctx)
    status, checkpoint_scheme = sun.SUNAdjointCheckpointScheme_Create_Fixed(
        sun.SUNDATAIOMODE_INMEM, mem_helper, check_freq, ncheck, keep_checks, sunctx
    )
    assert status == sun.SUN_SUCCESS
    # status = ark.ARKodeSetAdjointCheckpointScheme(arkode.get(), checkpoint_scheme)
    # assert status == ark.ARK_SUCCESS

    #
    # Compute the forward solution
    #

    print("Initial condition:")
    yarr = sun.N_VGetArrayPointer(y)
    print(yarr)

    tret = t0
    status, tret = ark.ARKodeEvolve(arkode.get(), tf, y, ark.ARK_NORMAL)
    assert status == ark.ARK_SUCCESS
    print("Forward Solution:")
    print(sun.N_VGetArrayPointer(y))
    # print("ARKODE Stats for Forward Solution:")
    # ARKodePrintAllStats(arkode.get(), None, 0)
    # print()

    #
    # Create the adjoint stepper
    #

    # Adjoint terminal condition
    uB = sun.N_VNew_Serial(NEQ, sunctx)
    arr_uB = ode.dgdu(y)
    uB_arr = sun.N_VGetArrayPointer(uB)
    uB_arr[:] = arr_uB
    qB = sun.N_VNew_Serial(NP, sunctx)
    qB_arr = sun.N_VGetArrayPointer(qB)
    qB_arr[:] = ode.dgdp(y)

    # Combine adjoint vectors into a ManyVector
    sens = [uB, qB]
    sf = sun.N_VNew_ManyVector(2, sens, sunctx)
    print("Adjoint terminal condition:")
    print(sun.N_VGetArrayPointer(uB))
    print(sun.N_VGetArrayPointer(qB))

    # # Create ARKStep adjoint stepper
    # status, adj_stepper = ark.ARKStepCreateAdjointStepper(
    #     arkode.get(),
    #     lambda t, yv, lv, ldotv, _: ode.adj_rhs(t, yv, lv, ldotv),
    #     None,
    #     tf,
    #     sf,
    #     sunctx,
    # )
    # adj_stepper = sun.SUNAdjointStepperView.Create(adj_stepper)

    # #
    # # Now compute the adjoint solution
    # #

    # status, tret = sun.SUNAdjointStepper_Evolve(adj_stepper, t0, sf)
    # assert status == ark.ARK_SUCCESS

    # print("Adjoint Solution:")
    # print(sun.N_VGetArrayPointer(uB))
    # print(sun.N_VGetArrayPointer(qB))

    # print("\nARKStep Adjoint Stats:")
    # ARKStepAdjointStepperPrintAllStats(adj_stepper.get(), None, 0)

    # # print(type(cs_view))


if __name__ == "__main__":
    main()
