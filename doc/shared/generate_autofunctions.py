# -----------------------------------------------------------------------------
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
# -----------------------------------------------------------------------------
# Script that generates sphinx-autodoc autofunction directives for each
# function in the Python modules. This is necessary due to automodule not 
# yet being able to handle nanobind generated functions.
# See https://github.com/sphinx-doc/sphinx/issues/13868.
# -----------------------------------------------------------------------------

import os
import importlib
import sundials4py


def generate_autofunctions_for_submodule(module_name: str):
    module = importlib.import_module(f'sundials4py.{module_name}')
    autogen_file = os.path.join(os.path.dirname(__file__), f'./Python/sundials4py-{module_name}-functions.rst')
    with open(autogen_file, 'w') as f:
        f.write('Functions\n')
        f.write('^^^^^^^^^\n\n')
        for func_name in dir(module):
            obj = getattr(module, func_name)
            if type(obj).__name__ == 'nb_func':
                f.write(f'.. autofunction:: sundials4py.{module_name}.{func_name}\n')
                f.write('  :no-index:\n\n')
                f.write(f'  See :c:func:`{func_name}`.\n\n')


def generate_autofunctions_for_sundials4py():
    generate_autofunctions_for_submodule('core')
    generate_autofunctions_for_submodule('arkode')
    generate_autofunctions_for_submodule('cvodes')
    generate_autofunctions_for_submodule('idas')
    generate_autofunctions_for_submodule('kinsol')
