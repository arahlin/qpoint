"""
The C-side self-test of the constructor and copy API.

qp_init_det, qp_init_detarr, qp_init_point, qp_init_map_from_arrays,
qp_copy_pixhash and qp_copy_memory are public C API that nothing here
calls: the ctypes layer builds those structures itself, field by field
through its mirrors. So nothing else in this suite reaches them.

src/qp_selftest.c drives them as an external C consumer would, and this is
the one call that runs it. A failure comes back as the number and
description of the first check that failed rather than a bare exit code.
"""

import pytest
import qpoint
from qpoint._libqpoint import selftest


def test_c_constructors_and_copies():
    """
    Run all of it. The C reports which check failed, so the assertion
    message names the specific broken invariant.
    """
    failure = selftest()
    assert failure is None, "C self-test: {}".format(failure)


def test_the_selftest_is_actually_wired_up():
    """
    A guard against the binding silently going missing: if the symbol
    were dropped from the build, `selftest()` would raise rather than
    return, and a test asserting only "no failure" could not tell the
    difference between passing and never running.
    """
    assert callable(selftest)
    assert hasattr(qpoint._libqpoint.libqp, "qp_selftest")
