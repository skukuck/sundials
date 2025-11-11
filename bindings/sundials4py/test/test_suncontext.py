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


def test_push_pop_err_handlers():
    # Create a new context with a null comm
    err, sunctx = SUNContext_Create(SUN_COMM_NULL)
    assert err == SUN_SUCCESS

    called = {"err_fn1": False}

    def err_fn1(line, func_name, file_name, msg, err_code, _, sunctx):
        called["err_fn1"] = True

    status = SUNContext_PushErrHandler(sunctx, err_fn1)
    assert status == SUN_SUCCESS

    # # TODO(CJB): this will only cause an error if error checks are turned on
    # # Try and make negative length vector to trigger error handler
    # x = N_VNew_Serial(-1, sunctx)
    # assert called["err_fn1"]

if __name__ == "__main__":
    test_push_pop_err_handlers()
