.. For package-specific references use :ref: rather than :numref: so intersphinx
   links to the appropriate place on read the docs

**Major Features**

**New Features and Enhancements**

Multiple minor updates were made to the ARKODE package.  We removed an extraneous
copy of the output vector when using ARKODE in ``ARK_ONE_STEP`` mode.  We added the
function :c:func:`ARKodeGetStageIndex` that returns the index of the stage currently
being processed, and the total number of stages in the method, for users who must
compute auxiliary quantities in their IVP right-hand side functions during some
stages and not others (e.g., in all but the first or last stage). We added the
function :c:func:`ARKodeAllocateInternalData` to ARKODE to enable stage-related
data allocation before the first call to :c:func:`ARKodeEvolve` (but after all other
optional input routines have been called), to support users who measure memory
usage before beginning a simulation.  Finally, we standardized calls to the user-supplied
right-hand-side functions so that these are provided the user-supplied solution vector
passed to :c:func:`ARKodeEvolve` whenever possible -- the notable exceptions are the
Hermite temporal interpolation module, the provided preconditioners ARKBANDPRE and
ARKBBDPRE, banded or dense linear solvers with automatically-approximated Jacobian
matrices, iterative linear solvers with automatically-approximated Jacobian-times-vector
product, temporal root-finding, discrete adjoint modules in ARKStep or ERKStep, the
SPRKStep stepper, and LSRKStep's use of the automated dominant eigenvalue estimation module.

An optional N_Vector routine, :c:func:`N_VCopy`, was added, to streamline data copies between two
vectors.  For user-supplied N_Vector modules that do not provide this function, :c:func:`N_VScale`
will be used instead.


**Bug Fixes**

**Deprecation Notices**
