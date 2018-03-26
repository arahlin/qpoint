from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
# numpy.distutils handles intel compiler
# need to pass --compiler=intel (or add to [build_ext] section of setup.cfg)
from numpy.distutils.core import setup, Extension
from numpy.distutils.command.build import build
try:
    import sphinx.setup_command
    sphinx_found = True
except ImportError:
    sphinx_found = False
import os, glob
import sys
import sysconfig

# get version from git tag
import subprocess as sp
version = sp.check_output(
    'git describe --abbrev=4 --dirty --always --tags'.split()).strip().decode()
version_simple = version.split('-')[0]
vtup = tuple(int(x) for x in version_simple.split('.'))
print('qpoint version {}'.format(version_simple))
varg = '-DQP_VERSION=\"{}\"'.format(version)

# build directory names
def build_dir_name(name):
    path = 'build/{name}.{platform}-{version[0]}.{version[1]}'.format(
        name=name, platform=sysconfig.get_platform(), version=sys.version_info)
    return os.path.abspath(path)

# build classes
class BuildLib(build):

    def run(self):
        # build dependencies
        print('make -C sofa')
        sp.check_call('make -C sofa'.split())
        print('make -C chealpix')
        sp.check_call('make -C chealpix'.split())
        print('make -C src')
        sp.check_call('make -C src'.split())

        # write version files
        with open('python/_version.py', 'w') as f:
            f.write('__version__ = ({:d}, {:d}, {:d})\n'.format(*vtup))

        build.run(self)

if sphinx_found:
    class BuildDoc(sphinx.setup_command.BuildDoc):
        def run(self):
            sys.path.insert(0, build_dir_name('lib'))
            sphinx.setup_command.BuildDoc.run(self)
else:
    class BuildDoc(build):
        user_options = [('all', 'a', '')]
        def initialize_options(self):
            raise ImportError('sphinx not found, cannot build docs')

# setup extension arguments
src = [x for x in glob.glob('src/*.c')]
src = [x for x in src if not x.endswith('iers_bulletin_a.c')]
src += ['src/qp_iers_bulletin_a.c']
incl_dirs = ['src','sofa','chealpix']
libsofa_file = 'sofa/libsofa_c.a'
libchealpix_file = 'chealpix/libchealpix_qp.a'
extra_obj = [libsofa_file, libchealpix_file]
extra_args = ['-O3', '-Wall', '-std=c99', '-fPIC', varg]

# do different stuff if using intel compilers
if os.getenv("CC", "") == "icc":
    extra_args.append('-qopenmp')
    libs = []
else:
    extra_args.append('-fopenmp')
    libs = ['gomp']

# this isn't technically an extension...
# hack to make a shared library to install with the package
ext_qp = Extension('qpoint.libqpoint', src,
                   include_dirs=incl_dirs,
                   extra_compile_args=extra_args,
                   libraries=libs,
                   extra_objects=extra_obj)

# run setup
setup(name='qpoint',
      version=version,
      description='Pointing for SPIDER',
      author='Alexandra Rahlin',
      packages=['qpoint'],
      package_dir = {'qpoint': 'python'},
      ext_modules=[ext_qp],
      cmdclass={
          'build': BuildLib,
          'build_sphinx': BuildDoc,
      },
      command_options={
          'build_sphinx': {
              'project': ('setup.py', 'qpoint'),
              'version': ('setup.py', version_simple),
              'release': ('setup.py', version_simple),
              'source_dir': ('setup.py', 'docs'),
              'build_dir': ('setup.py', build_dir_name('docs')),
              }
      }
)
