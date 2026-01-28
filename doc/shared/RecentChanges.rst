.. For package-specific references use :ref: rather than :numref: so intersphinx
   links to the appropriate place on read the docs

**Major Features**

**New Features and Enhancements**

Multiple minor updates were made to the ARKODE package.  We removed an extraneous
copy of the output vector when using ARKODE in ``ARK_ONE_STEP`` mode.  We
standardized calls to the user-supplied right-hand-side functions so that these
are provided the user-supplied solution vector passed to :c:func:`ARKodeEvolve` whenever
possible -- the notable exceptions are the Hermite temporal interpolation module,
the provided preconditioners ARKBANDPRE and ARKBBDPRE, banded or dense linear
solvers with automatically-approximated Jacobian matrices, iterative linear solvers
with automatically-approximated Jacobian-times-vector product, temporal root-finding,
discrete adjoint modules in ARKStep or ERKStep, the SPRKStep stepper, and LSRKStep's
use of the automated dominant eigenvalue estimation module.

**Bug Fixes**

**Deprecation Notices**
