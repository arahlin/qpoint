"""
Shared pytest configuration.

The benchmarks in test_benchmarks.py are opt-in: they measure rather than
assert, and a shared runner is too noisy for the numbers to mean much. Run
them with `pytest --benchmark`, which prints a table at the end; without
the flag they are collected and skipped.
"""

import os
import subprocess
import sys
import threading
import time

import pytest


def openmp_runtime():
    """
    The OpenMP runtime the extension is linked against, or None.

    There is nothing in the C library to ask, so this reads the linkage
    of the built extension. It is worth knowing: without OpenMP the
    threading tests still pass, because one thread trivially agrees with
    itself, so a build that quietly lost its parallelism would look no
    different from one that kept it.
    """
    from qpoint._libqpoint import libqp

    path = getattr(libqp, "_name", None)
    if not path or not os.path.exists(path):
        return None
    cmd = ["ldd", path] if sys.platform.startswith("linux") else ["otool", "-L", path]
    try:
        out = subprocess.run(cmd, capture_output=True, text=True, timeout=60).stdout
    except (OSError, subprocess.SubprocessError):
        return None
    for line in out.splitlines():
        name = line.split()[0] if line.split() else ""
        # libgomp is gcc's, libomp/libiomp LLVM's; "omp" alone matches
        # libSystem on macOS
        if any(k in name.lower() for k in ("gomp", "libomp", "iomp")):
            return os.path.basename(name)
    return None


@pytest.fixture(scope="session")
def omp_runtime():
    return openmp_runtime()


def pytest_report_header(config):
    """
    Say so in the header, so a run's parallelism is never a guess.

    Both extensions are reported. They are built separately and the
    threading tests are split across them, so knowing about one says
    nothing about the other -- and the pair of tests that brackets
    `HAS_OPENMP` skips one way or the other either way, which makes the
    skip count useless for telling them apart.
    """
    parts = ["qpoint {}".format(openmp_runtime() or "serial")]
    try:
        import qpoint2._libqpoint2 as lib2

        parts.append("qpoint2 {}".format("threaded" if lib2.HAS_OPENMP else "serial"))
    except ImportError:
        pass
    return "openmp: " + ", ".join(parts)


def pytest_addoption(parser):
    parser.addoption(
        "--benchmark",
        action="store_true",
        default=False,
        help="run the timing benchmarks and print a summary table",
    )


def pytest_configure(config):
    config.addinivalue_line("markers", "benchmark: a timing measurement, not a check")
    config._benchmark_rows = []


@pytest.fixture(scope="session")
def benchmark_rows(request):
    return request.config._benchmark_rows


@pytest.fixture
def bench(request, benchmark_rows):
    """
    Time a call and record it for the summary table.

    Reports the fastest of several runs: the spread is contention, the floor
    the cost of the work.

    The work runs on a worker thread, and has to. At unlucky main-thread stack
    offsets a loop's spilled temporaries alias the heap buffers it streams,
    moving a measurement by 3x reproducibly. A worker thread gets a freshly
    mapped stack, and only waits, so nothing contends for the GIL.
    """
    if not request.config.getoption("--benchmark"):
        pytest.skip("timing benchmark; pass --benchmark to run")

    def run(label, fn, samples=None, repeat=5, unit="sample"):
        box = {}

        def measure():
            try:
                box["result"] = fn()  # warm up: first touch of every buffer
                times = []
                for _ in range(repeat):
                    start = time.perf_counter()
                    box["result"] = fn()
                    times.append(time.perf_counter() - start)
                box["times"] = times
            except BaseException as err:  # reported on the calling thread
                box["error"] = err

        worker = threading.Thread(target=measure)
        worker.start()
        worker.join()
        if "error" in box:
            raise box["error"]
        result, times = box["result"], box["times"]
        best = min(times)
        rate = samples / best if samples else None
        benchmark_rows.append(
            (label, 1e3 * best, 1e3 * sorted(times)[len(times) // 2], rate, unit)
        )
        return result

    return run


def pytest_terminal_summary(terminalreporter, exitstatus, config):
    rows = getattr(config, "_benchmark_rows", [])
    if not rows:
        return
    width = max(len(r[0]) for r in rows)
    write = terminalreporter.write_line
    write("")
    write("benchmarks (fastest of 5)")
    write(f"  {'':{width}}  {'best':>9}  {'median':>9}  rate")
    for label, best, median, rate, unit in rows:
        speed = f"{rate:,.0f} {unit}/s" if rate else ""
        write(f"  {label:{width}}  {best:>8.2f}ms  {median:>8.2f}ms  {speed}")
