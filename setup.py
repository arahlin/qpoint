from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
# numpy.distutils handles intel compiler
# need to pass --compiler=intel (or add to [build_ext] section of setup.cfg)
from numpy.distutils.core import setup, Extension
import os, glob
from warnings import warn

# get version from git tag
import subprocess as sp
version = sp.check_output(
    'git describe --abbrev=4 --dirty --always --tags'.split()).strip()
vtup = tuple(int(x) for x in version.split('-')[0].split('.'))
print('qpoint version {}'.format(vtup))
with open('python/_version.py', 'w') as f:
    f.write('__version__ = ({:d}, {:d}, {:d})\n'.format(*vtup))
varg = '-DQP_VERSION=\"{}\"'.format(version)

# hack to avoid recompiling sofa every time
libsofa_file = 'sofa/libsofa_c.a'
if not os.path.exists(libsofa_file):
    os.system('make -C sofa')

libchealpix_file = 'chealpix/libchealpix_qp.a'
if not os.path.exists(libchealpix_file):
    os.system('make -C chealpix')

sp.check_call('make -C src qp_iers_bulletin_a.c'.split())
src = [x for x in glob.glob('src/*.c')]
src = [x for x in src if not x.endswith('iers_bulletin_a.c')]
src += ['src/qp_iers_bulletin_a.c']
incl_dirs = ['src','sofa','chealpix']
extra_obj = [libsofa_file, libchealpix_file]
extra_args = ['-O3', '-Wall', '-std=c99', '-fPIC', varg]

# do different stuff if using intel copmilers
if os.getenv("CC", "") == "icc":
    extra_args.append('-qopenmp')
    libs = []
else:
    extra_args.append('-fopenmp')
    libs = ['gomp']

libslarefro_file = 'slarefro/libslarefro.a'
if os.path.exists(libslarefro_file):
    incl_dirs.append('slarefro')
    extra_obj.append(libslarefro_file)
    extra_args.append('-DENABLE_SLAREFRO')

# this isn't technically an extension...
# hack to make a shared library to install with the package
ext_qp = Extension('qpoint.libqpoint', src,
                   include_dirs=incl_dirs,
                   extra_compile_args=extra_args,
                   libraries=libs,
                   extra_objects=extra_obj)

setup(name='qpoint',
      version=version,
      description='Pointing for SPIDER',
      author='ASR',
      packages=['qpoint'],
      package_dir = {'qpoint': 'python'},
      ext_modules=[ext_qp]
      )
