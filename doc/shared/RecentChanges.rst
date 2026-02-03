.. For package-specific references use :ref: rather than :numref: so intersphinx
   links to the appropriate place on read the docs

**Major Features**

**New Features and Enhancements**

  The default numbers of stages for the SSP Runge--Kutta methods :c:enumerator:`ARKODE_LSRK_SSP_S_2`
  and :c:enumerator:`ARKODE_LSRK_SSP_S_3` in LSRKStep were changed from 10 and 9, respectively, to
  their minimum allowable values of 2 and 4.  Users may revert to the previous values by calling
  :c:func:`LSRKStepSetNumSSPStages`.

**Bug Fixes**

**Deprecation Notices**
