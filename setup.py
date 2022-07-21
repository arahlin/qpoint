from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

# numpy.distutils handles intel compiler
# need to pass --compiler=intel (or add to [build_ext] section of setup.cfg)
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext

import os, glob
import sys
import subprocess as sp

# get version from git tag
version = os.getenv("QPOINT_VERSION")
if version is None:
    cmd = 'git describe --abbrev=4 --dirty --always --tags'.split()
    version = sp.check_output(cmd).strip().decode()
version_simple = version.split('-')[0]

# write version file
with open('python/_version.py', 'w') as f:
    vtup = tuple(int(x) for x in version_simple.split('.'))
    f.write('__version__ = {:r}\n'.format(version_simple))

print('qpoint version {}'.format(version_simple))
varg = '-DQP_VERSION=\"{}\"'.format(version)

# build dependencies
if not os.path.exists("src/qp_iers_bulletin_a.c"):
    print('make qp_iers_bulletin_a.c')
    sp.check_call(
        'cd src && {} make_iers_bulletin_a_dot_c.py'.format(sys.executable),
        shell=True,
    )

# build classes
class BuildExt(build_ext):
    def build_extensions(self):
        c = self.compiler.compiler[0]
        print("Building extension with {}".format(c))
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
    author='Alexandra Rahlin',
    packages=['qpoint'],
    package_dir={'qpoint': 'python'},
    ext_modules=[ext_qp],
    cmdclass={
        'build_ext': BuildExt,
    },
    install_requires=['numpy>=1.10.0'],
)
