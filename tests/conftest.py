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

import threading
import time

import pytest


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
