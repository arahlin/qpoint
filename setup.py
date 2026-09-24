from glob import glob

from setuptools import setup, find_packages, Extension
from extension_helpers import add_openmp_flags_if_available

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
    # exercises the C constructors the ctypes layer never calls; see
    # tests/test_selftest.py
    "src/qp_selftest.c",
]

extra_args = ["-O3", "-Wall", "-std=c99", "-fPIC"]

# Headers are not sources, and setuptools judges an extension stale by
# comparing it against sources + depends only. Without them listed here a
# header-only edit rebuilds nothing at all, says nothing about it, and
# leaves you testing the previous build. The cost is that changing one
# header rebuilds the whole extension: distutils has no per-object
# staleness check, so once build_extension decides to run it recompiles
# every source.
headers = sorted(glob("src/*.h"))

ext_qp = Extension(
    "qpoint.libqpoint", src, depends=headers, extra_compile_args=extra_args
)

# add openmp support if possible
add_openmp_flags_if_available(ext_qp)

# run setup
setup(
    ext_modules=[ext_qp],
    packages=find_packages(),
)
