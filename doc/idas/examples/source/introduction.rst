..
   -----------------------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2025-2026, Lawrence Livermore National Security,
   University of Maryland Baltimore County, and the SUNDIALS contributors.
   Copyright (c) 2013-2025, Lawrence Livermore National Security
   and Southern Methodist University.
   Copyright (c) 2002-2013, Lawrence Livermore National Security.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   -----------------------------------------------------------------------------

.. _IDAS.Examples.Intro:

Introduction
============

This report is intended to serve as a companion document to the IDAS User
Guide :cite:p:`idas_ug`. It provides details, with listings, on the example
programs supplied with the IDAS distribution package.

The IDAS distribution contains examples of the following types: serial and
parallel examples of Initial Value Problem (IVP) integration, serial and
parallel examples of forward sensitivity analysis (FSA), and serial and
parallel examples of adjoint sensitivity analysis (ASA). In addition, there are
two examples using OpenMP, and there are two examples in Fortran using the
SUNDIALS Fortran interface modules.

Example naming convention
--------------------------

With the exception of "demo"-type example files, the names of all the examples
distributed with SUNDIALS are of the form ``[slv][PbName]_[SA]_[ls]_[prec]_[p]``,
where:

- ``[slv]`` identifies the solver (for IDAS examples this is ``idas``).
- ``[PbName]`` identifies the problem.
- ``[SA]`` identifies sensitivity analysis examples. This field can be one of:
  ``FSA`` for forward sensitivity examples, ``ASAi`` for adjoint sensitivity
  examples using an integral-form model output, or ``ASAp`` for adjoint
  sensitivity examples using a pointwise model output.
- ``[ls]`` identifies the linear solver module used.
- ``[prec]`` indicates the IDAS preconditioner module used (if applicable — for
  examples using a Krylov linear solver and the ``IDABBDPRE`` module, this will
  be ``bbd``).
- ``[p]`` indicates an example using the MPI parallel vector.

Example overview
----------------

The following table lists all IDAS examples by type:

.. table:: IDAS example programs

   +-----+------------------------------------+----------------------------------+
   |     | Serial examples                    | Parallel examples                |
   +=====+====================================+==================================+
   | IVP | ``idasRoberts_dns``                | ``idasHeat2D_kry_p``             |
   |     |                                    |                                  |
   |     | ``idasRoberts_klu``                | ``idasHeat2D_kry_bbd_p``         |
   |     |                                    |                                  |
   |     | ``idasRoberts_sps``                | ``idasFoodWeb_kry_p``            |
   |     |                                    |                                  |
   |     | ``idasAkzoNob_dns``                | ``idasFoodWeb_kry_bbd_p``        |
   |     |                                    |                                  |
   |     | ``idasSlCrank_dns``                | ``idasBruss_kry_bbd_p``          |
   |     |                                    |                                  |
   |     | ``idasHeat2D_bnd``                 |                                  |
   |     |                                    |                                  |
   |     | ``idasHeat2D_kry``                 |                                  |
   |     |                                    |                                  |
   |     | ``idasFoodWeb_bnd``                |                                  |
   |     |                                    |                                  |
   |     | ``idasFoodWeb_bnd_omp``            |                                  |
   |     |                                    |                                  |
   |     | ``idasFoodWeb_kry_omp``            |                                  |
   |     |                                    |                                  |
   |     | ``idasKrylovDemo_ls``              |                                  |
   |     |                                    |                                  |
   |     | ``idasAnalytic_mels``              |                                  |
   |     |                                    |                                  |
   |     | ``idasHeat2D_kry_f2003``           |                                  |
   +-----+------------------------------------+----------------------------------+
   | FSA | ``idasRoberts_FSA_dns``            | ``idasBruss_FSA_kry_bbd_p``      |
   |     |                                    |                                  |
   |     | ``idasRoberts_FSA_klu``            | ``idasHeat2D_FSA_kry_bbd_p``     |
   |     |                                    |                                  |
   |     | ``idasRoberts_FSA_sps``            |                                  |
   |     |                                    |                                  |
   |     | ``idasSlCrank_FSA_dns``            |                                  |
   +-----+------------------------------------+----------------------------------+
   | ASA | ``idasRoberts_ASAi_dns``           | ``idasBruss_ASAp_kry_bbd_p``     |
   |     |                                    |                                  |
   |     | ``idasRoberts_ASAi_klu``           |                                  |
   |     |                                    |                                  |
   |     | ``idasRoberts_ASAi_sps``           |                                  |
   |     |                                    |                                  |
   |     | ``idasAkzoNob_ASAi_dns``           |                                  |
   |     |                                    |                                  |
   |     | ``idasAkzoNob_ASAi_dns_f2003``     |                                  |
   |     |                                    |                                  |
   |     | ``idasHessian_ASA_FSA``            |                                  |
   +-----+------------------------------------+----------------------------------+

Serial examples
---------------

Supplied in the ``examples/idas/serial`` directory are the following serial
examples (using the :ref:`serial vector <NVectors.NVSerial>`):

IVP integration examples
^^^^^^^^^^^^^^^^^^^^^^^^^

- ``idasRoberts_dns`` solves the Robertson chemical kinetics problem
  :cite:p:`Rob:66`, which consists of two differential equations and one
  algebraic constraint. It also uses the rootfinding feature of IDAS.

  The problem is solved with the :ref:`dense linear solver <SUNLinSol_Dense>`
  using a user-supplied Jacobian.

- ``idasRoberts_klu`` is the same as ``idasRoberts_dns`` but uses the :ref:`KLU
  sparse direct linear solver <SUNLinSol.KLU>`.

- ``idasRoberts_sps`` is the same as ``idasRoberts_dns`` but uses the
  :ref:`SuperLU_MT sparse direct linear solver <SUNLinSol.SuperLUMT>` (with one
  thread).

- ``idasAkzoNob_dns`` solves the Akzo-Nobel chemical kinetics problem, which
  consists of six nonlinear DAEs of index 1. The problem originates from Akzo
  Nobel Central research in Arnhem, The Netherlands, and describes a chemical
  process in which two species are mixed, while carbon dioxide is continuously
  added.

  The problem is solved with the :ref:`dense linear solver <SUNLinSol_Dense>`
  using the default difference quotient dense Jacobian approximation.

- ``idasHeat2D_bnd`` solves a 2D heat equation, semidiscretized to a DAE on the
  unit square.

  This program solves the problem with the :ref:`banded direct linear solver
  <SUNLinSol_Band>` and the default difference-quotient Jacobian
  approximation. For purposes of illustration, ``IDACalcIC`` is called to
  compute correct values at the boundary, given incorrect values as input
  initial guesses. The constraint :math:`u > 0.0` is imposed for all components.

- ``idasHeat2D_kry`` solves the same 2D heat equation problem as
  ``idasHeat2D_bnd``, with the :ref:`GMRES iterative linear solver
  <SUNLinSol.SPGMR>`. The preconditioner uses only the diagonal elements of the
  Jacobian.

- ``idasFoodWeb_bnd`` solves a system of PDEs modeling a food web problem, with
  predator-prey interaction and diffusion, on the unit square in 2D.

  The PDEs are discretized in space to a system of DAEs which are solved using
  the :ref:`banded direct linear solver <SUNLinSol_Band>` with the default
  difference-quotient Jacobian approximation.

- ``idasSlCrank_dns`` solves a system of index-2 DAEs, modeling a planar
  slider-crank mechanism.

  The problem is obtained through a stabilized index reduction
  (Gear-Gupta-Leimkuhler) starting from the index-3 DAE equations of motion
  derived using three generalized coordinates and two algebraic position
  constraints. The program also computes the time-averaged kinetic energy as a
  quadrature.

- ``idasKrylovDemo_ls`` solves the same problem as ``idasHeat2D_kry``, with
  three Krylov linear solvers: :ref:`GMRES <SUNLinSol.SPGMR>`, :ref:`BiCGSTAB
  <SUNLinSol.SPBCGS>`, and :ref:`TFQMR <SUNLinSol.SPTFQMR>`. The preconditioner
  uses only the diagonal elements of the Jacobian.

- ``idasAnalytic_mels`` solves a problem having a simple analytic solution,
  using a custom matrix-embedded linear solver.

Forward sensitivity analysis examples
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- ``idasRoberts_FSA_dns`` solves the same kinetics problem as in
  ``idasRoberts_dns``.

  IDAS also computes both its solution and solution sensitivities with respect
  to the three reaction rate constants appearing in the model. This program
  solves the problem with the :ref:`dense linear solver <SUNLinSol_Dense>`, and
  a user-supplied Jacobian routine.

- ``idasRoberts_FSA_klu`` solves the same problem as in ``idasRoberts_FSA_dns``
  but uses the :ref:`KLU sparse direct linear solver <SUNLinSol.KLU>`.

- ``idasRoberts_FSA_sps`` solves the same problem as in ``idasRoberts_FSA_dns``
  but uses the :ref:`SuperLU_MT sparse direct linear solver
  <SUNLinSol.SuperLUMT>`.

- ``idasSlCrank_FSA_dns`` solves a system of index-2 DAEs, modeling a planar
  slider-crank mechanism.

  This example computes both its solution and solution sensitivities with
  respect to the problem parameters :math:`k` (spring constant) and :math:`c`
  (damper constant), and then uses them to evaluate the gradient of the
  cumulative kinetic energy of the system.

Adjoint sensitivity analysis examples
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- ``idasRoberts_ASAi_dns`` solves the same kinetics problem as in
  ``idasRoberts_dns``.

  Here the adjoint capability of IDAS is used to compute gradients of a
  functional of the solution with respect to the three reaction rate constants
  appearing in the model. This program solves both the forward and backward
  problems with the :ref:`dense linear solver <SUNLinSol_Dense>`, and
  user-supplied Jacobian routines.

- ``idasRoberts_ASAi_klu`` solves the same problem as in
  ``idasRoberts_ASAi_dns``, but uses the :ref:`KLU sparse direct linear solver
  <SUNLinSol.KLU>`.

- ``idasRoberts_ASAi_sps`` solves the same problem as in
  ``idasRoberts_ASAi_dns``, but uses the :ref:`SuperLU_MT sparse direct linear
  solver <SUNLinSol.SuperLUMT>`.

- ``idasAkzoNob_ASAi_dns`` solves the Akzo-Nobel chemical kinetics problem.

  The adjoint capability of IDAS is used to compute gradients with respect to
  the initial conditions of the integral over time of the concentration of the
  first species.

- ``idasHessian_ASA_FSA`` is an example of using the *forward-over-adjoint*
  method for computing 2nd-order derivative information, in the form of
  Hessian-times-vector products.

MPI parallel examples
---------------------

Supplied in the ``examples/idas/parallel`` directory are the following parallel
examples (using the :ref:`MPI parallel vector <NVectors.NVParallel>`):

IVP integration examples
^^^^^^^^^^^^^^^^^^^^^^^^^

- ``idasHeat2D_kry_p`` solves the same 2D heat equation problem as
  ``idasHeat2D_kry``, with :ref:`GMRES <SUNLinSol.SPGMR>` in parallel, and with
  a user-supplied diagonal preconditioner.

- ``idasHeat2D_kry_bbd_p`` solves the same problem as ``idasHeat2D_kry_p``.

  This program uses the :ref:`GMRES iterative linear solver <SUNLinSol.SPGMR>`
  in parallel, and the :ref:`IDABBDPRE band-block-diagonal preconditioner
  <IDAS.Usage.precond.idabbdpre>` with half-bandwidths equal to 1.

- ``idasFoodWeb_kry_p`` solves the same food web problem as ``idasFoodWeb_bnd``,
  but with :ref:`GMRES <SUNLinSol.SPGMR>` and a user-supplied preconditioner.

  The preconditioner supplied to SPGMR is the block-diagonal part of the
  Jacobian with :math:`n_s \times n_s` blocks arising from the reaction terms
  only (:math:`n_s` = number of species).

- ``idasFoodWeb_kry_bbd_p`` solves the same food web problem as
  ``idasFoodWeb_kry_p``.

  This program solves the problem using :ref:`GMRES <SUNLinSol.SPGMR>` in
  parallel and the :ref:`IDABBDPRE preconditioner
  <IDAS.Usage.precond.idabbdpre>`.

- ``idasBruss_kry_bbd_p`` solves the two-species time-dependent PDE known as the
  Brusselator problem, using the :ref:`GMRES iterative linear solver
  <SUNLinSol.SPGMR>` and the :ref:`IDABBDPRE preconditioner
  <IDAS.Usage.precond.idabbdpre>`.

  The PDEs are discretized by central differencing on a 2D spatial mesh. The
  system is actually implemented on submeshes, processor by processor.

Forward sensitivity analysis examples
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- ``idasBruss_FSA_kry_bbd_p`` solves the Brusselator problem with the forward
  sensitivity capability in IDAS used to compute solution sensitivities with
  respect to two of the problem parameters, and then the gradient of a model
  output functional, written as the final time value of the spatial integral of
  the first PDE component.

- ``idasHeat2D_FSA_kry_bbd_p`` solves the same problem as
  ``idasHeat2D_kry_p``, but using the :ref:`IDABBDPRE preconditioner
  <IDAS.Usage.precond.idabbdpre>`, and with forward sensitivity enabled to compute
  the solution sensitivity with respect to two coefficients of the original PDE.

Adjoint sensitivity analysis examples
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- ``idasBruss_ASAp_kry_bbd_p`` solves the same problem as
  ``idasBruss_FSA_kry_bbd_p`` but using an adjoint sensitivity approach for
  computing the gradient of the model output functional.

OpenMP examples
---------------

Supplied in the ``examples/idas/C_openmp`` directory are the following examples,
using the :ref:`OpenMP vector <NVectors.OpenMP>`:

- ``idasFoodWeb_bnd_omp`` solves the same problem as in ``idasFoodWeb_bnd`` but
  uses the OpenMP vector.

- ``idasFoodWeb_kry_omp`` solves the same problem as in ``idasFoodWeb_kry`` but
  uses the OpenMP vector.

Fortran examples
----------------

Supplied in the ``examples/idas/F2003_serial`` directory are the following
examples, both in Fortran90 and using the SUNDIALS Fortran interface modules:

- ``idasHeat2D_kry_f2003`` solves the same problem as ``idasHeat2D_kry``.

- ``idasAkzoNob_ASAi_dns_f2003`` solves the same problem as
  ``idasAkzoNob_ASAi_dns``.

Additional notes
----------------

In the following sections, we give detailed descriptions of some (but not all)
of the sensitivity analysis examples. We do not discuss the examples for IVP
integration; for those, the interested reader should consult the IDA Examples
document :cite:p:`ida_ex`. Any IDA problem will work with IDAS with only two
modifications: (1) the main program should include the header file ``idas.h``
instead of ``ida.h``, and (2) the loader command must reference the
``libsundials_idas`` library instead of ``libsundials_ida``.

We also give our output files for each of the examples described below, but
users should be cautioned that their results may differ slightly from
these. Differences in solution values may differ within the tolerances, and
differences in cumulative counters, such as numbers of steps or Newton
iterations, may differ from one machine environment to another by as much as 10%
to 20%.

In the descriptions below, we make frequent references to the IDAS User Guide
:cite:p:`idas_ug`. All citations to specific sections are references to parts
of that user guide, unless explicitly stated otherwise.

.. note::

   The examples in the IDAS distribution were written in such a way as to
   compile and run for any combination of configuration options during the
   installation of SUNDIALS. As a consequence, they contain portions of code
   that will not typically be present in a user program. For example, all
   example programs make use of the variables ``SUNDIALS_EXTENDED_PRECISION``
   and ``SUNDIALS_DOUBLE_PRECISION`` to test if the solver libraries were built
   in extended or double precision, and use the appropriate conversion
   specifiers in ``printf`` functions. Similarly, all forward sensitivity
   examples can be run with or without sensitivity computations enabled and, in
   the former case, with various combinations of methods and error control
   strategies. This is achieved in these example through the program arguments.

For all the IDAS examples, either of the two sensitivity method options,
``IDA_SIMULTANEOUS`` or ``IDA_STAGGERED``, can be used, and sensitivities may
be included in the error test or not (``errconS`` set to ``SUNTRUE`` or
``SUNFALSE``, respectively, in the call to ``IDASetSensErrCon``).
