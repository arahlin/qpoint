"""
Shared pytest configuration.

The benchmarks in test_benchmarks.py are opt-in: they measure rather than
assert, they take a few seconds, and a shared CI runner is too noisy for
the numbers to mean much. Run them with

    pytest --benchmark

and they print a table at the end. Without the flag they skip, so the
ordinary suite stays fast and the benchmark code still gets collected,
which is enough to catch it going stale.
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

    Reports the fastest of several runs rather than the mean: the spread
    is contention with everything else on the machine, and the floor is
    the closest thing to the cost of the work itself.

    The work runs on a worker thread, which is not a detail. The main
    thread's stack sits directly below argv and the environment, so its
    address is a deterministic function of the command line, while every
    other mapping is randomized. At unlucky stack offsets a tight loop's
    spilled temporaries alias with the heap buffers it streams, and the
    measurement changes by 3x -- reproducibly, so it reads as a real
    result rather than as noise. Adding any pytest argument moves the
    stack and moves the number with it, which is how this was found. A
    worker thread gets a freshly mapped stack instead, and the same
    command measured with and without an extra flag then agrees row for
    row. The main thread only waits, so nothing contends for the GIL.
    """
    if not request.config.getoption("--benchmark"):
        pytest.skip("timing benchmark; pass --benchmark to run")

    # which package this test was parametrized over, if it was
    spec = getattr(request.node, "callspec", None)
    impl = spec.params.get("mod") if spec else None
    impl = getattr(impl, "__name__", "")

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
            (label, impl, 1e3 * best, 1e3 * sorted(times)[len(times) // 2], rate, unit)
        )
        return result

    return run


def pytest_terminal_summary(terminalreporter, exitstatus, config):
    """
    One row per benchmark. Where the suite ran over more than one
    package, they get a column each and a ratio, which is the number
    worth looking at; where it ran over one, the median and the rate.
    """
    rows = getattr(config, "_benchmark_rows", [])
    if not rows:
        return
    labels = list(dict.fromkeys(r[0] for r in rows))
    impls = [i for i in dict.fromkeys(r[1] for r in rows) if i]
    best = {(r[0], r[1]): r for r in rows}
    width = max(len(l) for l in labels)
    write = terminalreporter.write_line
    write("")
    write("benchmarks (fastest of 5)")

    if len(impls) < 2:
        write(f"  {'':{width}}  {'best':>9}  {'median':>9}  rate")
        for label, _, b, median, rate, unit in rows:
            speed = f"{rate:,.0f} {unit}/s" if rate else ""
            write(f"  {label:{width}}  {b:>8.2f}ms  {median:>8.2f}ms  {speed}")
        return

    write(f"  {'':{width}}" + "".join(f"  {i:>10}" for i in impls) + "   ratio")
    for label in labels:
        line = f"  {label:{width}}"
        times = []
        for impl in impls:
            row = best.get((label, impl))
            times.append(row[2] if row else None)
            line += f"  {row[2]:>8.2f}ms" if row else f"  {'':>10}"
        if len(times) == 2 and all(times):
            line += f"  {times[0] / times[1]:>5.2f}x"
        write(line)
