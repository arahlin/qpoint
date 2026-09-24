"""
The C-side self-test of the constructor and copy API.

`qp_init_det`, `qp_init_detarr`, `qp_init_point`, `qp_init_map_from_arrays`,
`qp_copy_pixhash` and `qp_copy_memory` are public C API that neither Python
package calls: `qpoint/_libqpoint.py` builds those structures itself, field
by field through its ctypes mirrors, and `qpoint2` replaced them with C++
classes. So nothing in the rest of this suite reaches them, and their
coverage read zero -- not because they were known-good, but because they
were unreachable from here.

`src/qp_selftest.c` drives them the way an external C consumer would, with
its own assertions, and this is the one call that runs it. A failure comes
back as the number and description of the first check that failed rather
than as a bare exit code, because a C abort would say nothing useful.

The first thing it found was real: `qp_init_point` marked `ctime` and
`q_hwp` as malloc'd whether or not it allocated them, so `qp_free_point`
handed `free()` two pointers the constructor never assigned. Fixed on
master, below this branch.
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
