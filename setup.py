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
]

extra_args = ["-O3", "-Wall", "-std=c99", "-fPIC"]

ext_qp = Extension("qpoint.libqpoint", src, extra_compile_args=extra_args)

# add openmp support if possible
add_openmp_flags_if_available(ext_qp)

# run setup
setup(
    ext_modules=[ext_qp],
    packages=find_packages(),
)
