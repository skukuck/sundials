#!/bin/python

import pytest
import numpy as np
from sundials4py.core import *


def test_with_null_comm():
    # Create a new context with a null comm
    err, sunctx = SUNContext_Create(SUN_COMM_NULL)
    assert err == SUN_SUCCESS

    # Try calling a SUNContext_ function
    last_err = SUNContext_GetLastError(sunctx)

    assert last_err == SUN_SUCCESS
