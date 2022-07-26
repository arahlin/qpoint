from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

# build_ext handles intel compiler
# need to pass --compiler=intel (or add to [build_ext] section of setup.cfg)
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext

import os, glob
import sys
import subprocess as sp
import sysconfig

# get version from git tag
version = os.getenv("QPOINT_VERSION")
if version is None:
    if os.path.exists(".git"):
        cmd = 'git describe --abbrev=4 --dirty --always --tags'.split()
        version = sp.check_output(cmd).strip().decode()
    else:
        vline = open('python/_version.py', 'r').read().strip()
        version = vline.split()[-1].strip("'\"")
version_simple = version.split('-')[0]

# write version file
with open('python/_version.py', 'w') as f:
    vtup = tuple(int(x) for x in version_simple.split('.'))
    f.write('__version__ = {!r}\n'.format(version_simple))

print('qpoint version {}'.format(version_simple))
varg = '-DQP_VERSION=\"{}\"'.format(version)

# build dependencies
build_iers = False
try:
    # try to get the latest IERS tables baked in, because why not
    from astropy.utils.iers import IERS_Auto
    build_iers = True
except ImportError:
    if not os.path.exists("src/qp_iers_bulletin_a.c"):
        build_iers = True

if build_iers:
    print('make qp_iers_bulletin_a.c')
    sp.check_call(
        'cd src && {} make_iers_bulletin_a_dot_c.py'.format(sys.executable),
        shell=True,
    )


def read(rel_path):
    # type: (str) -> str
    here = os.path.abspath(os.path.dirname(__file__))
    # intentionally *not* adding an encoding option to open, See:
    #   https://github.com/pypa/virtualenv/issues/201#issuecomment-3145690
    with open(os.path.join(here, rel_path)) as fp:
        return fp.read()


# build classes
class BuildExt(build_ext):
    def build_extensions(self):
        c = self.compiler.compiler[0]
        print("Building extension with {}".format(c))
        plat = sysconfig.get_platform()
        if 'macosx' not in plat:
            if 'gcc' in c:
                for e in self.extensions:
                    e.extra_compile_args.append('-fopenmp')
                    e.libraries.append('gomp')
            elif 'intel' in c:
                for e in self.extensions:
                    e.extra_compile_args.append('-qopenmp')
        build_ext.build_extensions(self)


# setup extension arguments
src = [x for x in glob.glob('erfa/*.c') if not 'test' in x]
src += [x for x in glob.glob('chealpix/*.c') if not 'test' in x]
src += [x for x in glob.glob('src/*.c')]
src = [x for x in src if not x.endswith('iers_bulletin_a.c')]
src += ['src/qp_iers_bulletin_a.c']
incl_dirs = ['src', 'erfa', 'chealpix']
extra_args = ['-O3', '-Wall', '-std=c99', '-fPIC', varg]

# this isn't technically an extension...
# hack to make a shared library to install with the package
ext_qp = Extension(
    'qpoint.libqpoint', src, include_dirs=incl_dirs, extra_compile_args=extra_args
)

# run setup
setup(
    name='qpoint',
    version=version_simple,
    description='A lightweight quaternion-based library for efficient telescope pointing.',
    long_description=read("README.md"),
    long_description_content_type="text/markdown",
    author='Alexandra Rahlin',
    packages=['qpoint'],
    package_dir={'qpoint': 'python'},
    ext_modules=[ext_qp],
    license_files=('LICENSE',),
    cmdclass={
        'build_ext': BuildExt,
    },
    install_requires=['numpy>=1.10.0'],
)
