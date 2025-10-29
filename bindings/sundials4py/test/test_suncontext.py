#!/bin/python

import pytest
import numpy as np
from sundials4py.core import *


def test_suncontext_wo_comm():
    # Create without comm
    sunctx = SUNContextCreate()

    # Try calling a SUNContext_ function
    last_err = SUNContext_GetLastError(sunctx)

    assert last_err == SUN_SUCCESS


def test_with_null_comm():
    # Create a new context with a null comm
    sunctx = SUNContextCreate(SUNContext_Create(SUN_COMM_NULL))

    # Try calling a SUNContext_ function
    last_err = SUNContext_GetLastError(sunctx)

    assert last_err == SUN_SUCCESS
