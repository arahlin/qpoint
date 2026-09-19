import os
import subprocess
import sys

from setuptools import setup, find_packages, Extension
from extension_helpers import add_openmp_flags_if_available
from extension_helpers._openmp_helpers import check_openmp_support


def _libomp_prefix():
    """Locate Homebrew's libomp, which Apple clang does not ship."""
    if sys.platform != "darwin":
        return None

    candidates = []
    if os.environ.get("LIBOMP_PREFIX"):
        candidates.append(os.environ["LIBOMP_PREFIX"])
    try:
        candidates.append(
            subprocess.run(
                ["brew", "--prefix", "libomp"],
                capture_output=True,
                text=True,
                check=True,
            ).stdout.strip()
        )
    except (OSError, subprocess.CalledProcessError):
        pass
    candidates += ["/opt/homebrew/opt/libomp", "/usr/local/opt/libomp"]

    for path in candidates:
        if path and os.path.exists(os.path.join(path, "include", "omp.h")):
            return path
    return None


def add_openmp(ext):
    """
    Enable OpenMP, falling back to Apple clang's spelling of the flags.

    extension_helpers probes with a bare -fopenmp, which Apple clang rejects
    outright even when libomp is installed, so on macOS it always reports no
    OpenMP support. Apple clang needs -Xpreprocessor -fopenmp and explicit
    paths to the separately installed runtime. Set LIBOMP_PREFIX to override
    the search.
    """
    if add_openmp_flags_if_available(ext):
        return True

    prefix = _libomp_prefix()
    if prefix is None:
        return False

    flags = {
        "compiler_flags": [
            "-Xpreprocessor",
            "-fopenmp",
            "-I{}/include".format(prefix),
        ],
        "linker_flags": [
            "-L{}/lib".format(prefix),
            "-Wl,-rpath,{}/lib".format(prefix),
            "-lomp",
        ],
    }
    if not check_openmp_support(openmp_flags=flags):
        return False

    ext.extra_compile_args.extend(flags["compiler_flags"])
    ext.extra_link_args.extend(flags["linker_flags"])
    print("Compiling {} with OpenMP via libomp at {}".format(ext.name, prefix))
    return True


# extension arguments
src = [
    "src/erfa.c",
    "src/chealpix.c",
    "src/qp_error.c",
    "src/qp_iers_bulletin_a.c",
    "src/qp_map.c",
    "src/qp_params.c",
    "src/qp_pixel.c",
    "src/qp_pixhash.c",
    "src/qpoint.c",
    "src/quaternion.c",
    "src/sincos.c",
]

extra_args = ["-O3", "-Wall", "-std=c99", "-fPIC"]

ext_qp = Extension("qpoint.libqpoint", src, extra_compile_args=extra_args)

# add openmp support if possible
add_openmp(ext_qp)

# run setup
setup(
    ext_modules=[ext_qp],
    packages=find_packages(),
)
