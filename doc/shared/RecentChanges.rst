.. For package-specific references use :ref: rather than :numref: so intersphinx
   links to the appropriate place on read the docs

**Major Features**

**New Features and Enhancements**

The functions ``CVodeGetUserDataB`` and ``IDAGetUserDataB`` were added to CVODES
and IDAS, respectively.

Added the ``ARKODE_SSP_ERK_2_1_2``, ``ARKODE_SSP_ERK_3_1_2``, ``ARKODE_SSP_ERK_4_1_2``,
``ARKODE_SSP_ERK_10_1_2``, ``ARKODE_SSP_ERK_4_2_3``, ``ARKODE_SSP_ERK_9_2_3``,
``ARKODE_SSP_ERK_10_3_4``, ``ARKODE_SSP_LSPUM_ERK_3_1_2``, and ``ARKODE_ASCHER_ERK_3_1_2``
embedded explicit Runge-Kutta Butcher tables.

Added the ``ARKODE_SSP_SDIRK_2_1_2``, ``ARKODE_SSP_DIRK_3_1_2``,
``ARKODE_SSP_LSPUM_SDIRK_3_1_2``, ``ARKODE_SSP_ESDIRK_4_2_3``,
``ARKODE_SSP_ESDIRK_9_2_3``, and ``ARKODE_ASCHER_SDIRK_3_1_2`` embedded diagonally
implicit Runge-Kutta Butcher tables.

Of these, embedded additive Runge-Kutta methods may be formed using
``ARKODE_SSP_ERK_2_1_2`` + ``ARKODE_SSP_SDIRK_2_1_2``,
``ARKODE_SSP_ERK_3_1_2`` + ``ARKODE_SSP_DIRK_3_1_2``,
``ARKODE_SSP_ERK_4_2_3`` + ``ARKODE_SSP_ESDIRK_4_2_3``,
``ARKODE_SSP_ERK_9_2_3`` + ``ARKODE_SSP_ESDIRK_9_2_3``,
``ARKODE_SSP_LSPUM_ERK_3_1_2`` + ``ARKODE_SSP_LSPUM_SDIRK_3_1_2``,
and ``ARKODE_ASCHER_ERK_3_1_2`` + ``ARKODE_ASCHER_SDIRK_3_1_2``.

Added the ``ARKODE_IMEX_MRI_GARK_ARS222`` and ``ARKODE_IMEX_MRI_GARK_GIRALDO``
embedded implicit-explicit MRI-GARK coupling tables.

**Bug Fixes**

On the initial time step with a user-supplied initial step size, ARKODE and
CVODE(S) will now return ``ARK_TOO_CLOSE`` or ``CV_TOO_CLOSE``, respectively,
when the requested output time is the same as the initial time (or within
numerical roundoff of the initial time). Before a ``TOO_CLOSE`` error would only
be returned when internally estimating the initial step size. In IDA(S), added a
``IDA_TOO_CLOSE`` return value for when the initial and output time are too
close. Previously, IDA(S) would return ``IDA_ILL_INPUT``.

Fixed a bug in ARKODE, CVODE(S), and IDA(S) where the linear solver counters
were not reinitialized until the next call to advance the system. As such,
non-zero linear solver statistics could be returned if retrieving or printing
linear solver counters between the initialization and the next call to advance
the system.

The interface to Ginkgo batched linear solvers has been updated to fix build
errors when using 64-bit index types. Note, only the batched dense matrix in
Ginkgo is currently compatible with 64-bit indexing (as of Ginkgo 1.10).

The SPRKStep module now accounts for zero coefficients in the SPRK tables,
eliminating extraneous function evaluations.

A bug preventing a user supplied :c:func:`SUNStepper_ResetCheckpointIndex`
function from being called was fixed.

The Kokkos N_Vector now properly handles unmanaged views. Previously, if a
Kokkos ``N_Vector`` was created from an unmanaged view, the view would become a
managed view and the data would be freed unexpectedly.

A bug was fixed in KINSOL where the information logging function would always be
called even when informational logging was disabled.

**Deprecation Notices**

``SUNDIALSFileOpen`` and ``SUNDIALSFileClose`` will be removed in the next major release.
Use :c:func:`SUNFileOpen` and :c:func:`SUNFileClose` instead.

The ``Convert`` methods on the ``sundials::kokkos:Vector``, ``sundials::kokkos::DenseMatrix``,
``sundials::ginkgo::Matrix``, ``sundials::ginkgo::BatchMatrix``, ``sundials::kokkos::DenseLinearSolver``,
``sundials::ginkgo::LinearSolver``, and ``sundials::ginkgo::BatchLinearSolver`` classes have
been deprecated and will be removed in the next major release. The method ``get``, should
be used instead.

The ``CSC_MAT`` and ``CSR_MAT`` macros defined in ``sunmatrix_sparse.h`` will be removed in
the next major release. Use ``SUN_CSC_MAT`` and ``SUN_CSR_MAT`` instead.

The ``N_Vector_S`` typedef to ``N_Vector*`` is deprecated and will be removed in the next major release.
