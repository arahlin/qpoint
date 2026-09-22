from glob import glob

from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext
from extension_helpers import add_openmp_flags_if_available
import pybind11

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

# -ffp-contract=off suppresses FMA contraction, which would otherwise differ between
# the C and C++ builds: they inline differently, so the contraction window differs and
# results diverge in the last bits.
common_args = ["-O3", "-Wall", "-fPIC", "-ffp-contract=off"]

# Headers are not sources, and setuptools judges an extension stale by
# comparing it against sources + depends only. Without them listed here a
# header-only edit rebuilds nothing at all, says nothing about it, and
# leaves you testing the previous build. The cost is that changing one
# header rebuilds the whole extension: distutils has no per-object
# staleness check, so once build_extension decides to run it recompiles
# every source.
headers = sorted(glob("src/*.h"))
headers_cpp = headers + sorted(glob("src/cpp/*.hpp"))

ext_qp = Extension(
    "qpoint.libqpoint",
    src,
    depends=headers,
    extra_compile_args=common_args + ["-std=c99"],
)

# add openmp support if possible
add_openmp_flags_if_available(ext_qp)

# qpoint2 extension: C++ core in src/cpp/, linking only the third-party C
# sources. The qp_*.c files are replaced by the C++ classes; erfa, chealpix
# and sincos are not.
src_cpp2 = ["src/qp2_bindings.cpp"] + sorted(glob("src/cpp/*.cpp"))

ext_qp2 = Extension(
    "qpoint2._libqpoint2",
    src_cpp2 + ["src/erfa.c", "src/chealpix.c", "src/sincos.c"],
    include_dirs=[pybind11.get_include(), "src", "src/cpp"],
    depends=headers_cpp,
    # the -std= flag is applied per file in the custom build_ext below
    extra_compile_args=common_args,
    language="c++",
)

add_openmp_flags_if_available(ext_qp2)

CXX_EXTENSIONS = {"qpoint2._libqpoint2"}


class BuildExt(build_ext):
    """Custom build_ext that applies per-file compiler standards."""

    def build_extension(self, ext):
        if ext.name not in CXX_EXTENSIONS:
            super().build_extension(ext)
            return

        # patch compile() to insert per-file flags
        original_compile = self.compiler.compile

        def patched_compile(sources, **kwargs):
            objects = []
            for src_file in sources:
                if src_file.endswith(".cpp"):
                    std_flag = "-std=c++17"
                else:
                    std_flag = "-std=c99"
                ea = list(kwargs.get("extra_postargs", []))
                ea = [f for f in ea if not f.startswith("-std=")] + [std_flag]
                kw = dict(kwargs, extra_postargs=ea)
                objs = original_compile([src_file], **kw)
                objects.extend(objs)
            return objects

        self.compiler.compile = patched_compile
        try:
            super().build_extension(ext)
        finally:
            self.compiler.compile = original_compile


# run setup
setup(
    ext_modules=[ext_qp, ext_qp2],
    packages=find_packages(),
    cmdclass={"build_ext": BuildExt},
)
