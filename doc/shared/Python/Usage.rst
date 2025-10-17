.. ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2025, Lawrence Livermore National Security,
   University of Maryland Baltimore County, and the SUNDIALS contributors.
   Copyright (c) 2013-2025, Lawrence Livermore National Security
   and Southern Methodist University.
   Copyright (c) 2002-2013, Lawrence Livermore National Security.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _Python.Usage:

Using sundials4py
=================

At a high level, using SUNDIALS from Python via sundials4py looks a lot like using SUNDIALS from C or C++.

The few notable differences are:

View Classes and Memory Management
----------------------------------

The SUNDIALS C objects (e.g., an ``N_Vector``) need to be wrapped in their corresponding "View" class to ensure memory is cleaned up in accordance with the Python garbage collection.
While we must wrap the C objects in a View, the SUNDIALS functions must be passed the raw C object. The object is extracted from the View with the ``get`` function. Together, these
two points result in code like:

.. code-block:: python
   
   from sundials4py.core import *

   sunctx = SUNContextView.Create()
   x = NVectorView.Create(N_VNew_Serial(10, sunctx.get()))


Return-by-Pointer Parameters
----------------------------

Functions that return values via pointer arguments in the C API are mapped to Python functions that return a tuple:

- **First element:** The function's return value (typically an error code).
- **Subsequent elements:** Values that would be returned via pointer arguments in C, in the same order as the C function signature.

**Example 1: Single Return-by-Pointer Value**

C:
   .. code-block:: C

      int CVodeGetNumSteps(void *cvode_mem, long int *numsteps);

Python:
   .. code-block:: python

      retval, numsteps = CVodeGetNumSteps(cvode_mem)
      print(f"Number of steps: {numsteps}")

**Example 2: Multiple Return-by-Pointer Values**

C:
   .. code-block:: C

      int CVodeGetIntegratorStats(void *cvode_mem,
                                 long int *nsteps,
                                 long int *nfevals,
                                 long int *nlinsetups,
                                 long int *netfails);

Python:
   .. code-block:: python

      retval, nsteps, nfevals, nlinsetups, netfails = CVodeGetIntegratorStats(cvode_mem)
      print(f"Steps: {nsteps}, Function evals: {nfevals}, Linear setups: {nlinsetups}, Error test fails: {netfails}")


Arrays
------

``N_Vector`` objects in sundials4py are compatible with numpy's `ndarray`. Each ``N_Vector`` can work on a numpy arrays without copies, and you can access 
and modify the underlying data directly using :py:func:`N_VGetArrayPointer`, which returns a numpy `ndarray` view of the data.

- SUNDIALS matrix types (dense, banded, sparse) are also exposed as Python objects that provide access to their underlying data as numpy arrays (e.g., via :py:func:`SUNDenseMatrix_Data`).
- Arrays of scalars (e.g., scaling factors passed to :py:func:`N_VLinearCombination`) are also represented as numpy arrays.

**Example: Accessing and modifying an N_Vector**

.. code-block:: python

   y_nvec = NVectorView.Create(N_VNew_Serial(10, sunctx.get()))
   y = N_VGetArrayPointer(y_nvec.get())
   y[:] = np.linspace(0, 1, 10)  # Set values using numpy

**Example: Using a matrix as a numpy array**

.. code-block:: python

   mat = SUNMatrixView.Create(SUNDenseMatrix(3, 3, sunctx.get()))
   arr = SUNDenseMatrix_Data(mat.get())
   arr = np.eye(3)  # Set to identity matrix

This allows you to use numpy operations for vector and matrix data, and to pass numpy arrays to and from SUNDIALS routines efficiently and without unnecessary copies.


User-Supplied Callback Functions
--------------------------------

SUNDIALS packages and several modules/classes require user-supplied callback functions to define problem-specific behavior, 
such as the right-hand side of an ODE or a nonlinear system function. In sundials4py, you can provide these as standard Python functions or lambdas.
Somethings to note:

- The callback signatures follow the C API exactly. As such, ``N_Vector`` arguments are passed as ``N_Vector`` objects and the underlying ndarray must be extracted in the user code.
- Most callback signatures include a ``void* user_data`` argument. In Python, this argument must be present in the signature, but it should be ignored.

**Example: ODE right-hand side for ARKStep**

.. code-block:: python

   def rhs(t, y_nvector, ydot_nvector, _):
      # Compute ydot = f(t, y)
      y = N_VGetArrayPointer(y_nvector)
      ydot = N_VGetArrayPointer(ydot_nvector) 
      ydot[:] = -y
      return 0

   ark = ARKodeView.Create(ARKStepCreate(rhs, None, t0, y.get(), sunctx.get()))

**Example: Nonlinear system for KINSOL**

.. code-block:: python

   def fp_function(u_nvector, g_nvector, _):
      # Compute g = F(u)
      u = N_VGetArrayPointer(u_nvector)
      g = N_VGetArrayPointer(g_nvector)
      g[:] = u**2 - 1
      return 0

   kin = KINView.Create(KINCreate(sunctx.get()))
   KINInit(kin.get(), fp_function, u.get())


.. warning::

   The ``user_data`` argument is almost always ``None`` on the Python side, but if it is not it should be ignored to avoid causing catastrophic errors.


Examples
--------

Examples can be found in ``bindings/sundials4py/examples``.
