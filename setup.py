from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

# numpy.distutils handles intel compiler
# need to pass --compiler=intel (or add to [build_ext] section of setup.cfg)
from numpy.distutils.core import setup, Extension
from numpy.distutils.command.build import build
from numpy.distutils.command.build_ext import build_ext

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

version = os.getenv("QPOINT_VERSION")
if version is None:
    cmd = 'git describe --abbrev=4 --dirty --always --tags'.split()
    version = sp.check_output(cmd).strip().decode()
version_simple = version.split('-')[0]
vtup = tuple(int(x) for x in version_simple.split('.'))
print('qpoint version {}'.format(version_simple))
varg = '-DQP_VERSION=\"{}\"'.format(version)

# build directory names
def build_dir_name(name):
    path = 'build/{name}.{platform}-{version[0]}.{version[1]}'.format(
        name=name, platform=sysconfig.get_platform(), version=sys.version_info
    )
    return os.path.abspath(path)


# build classes
class BuildLib(build):
    def run(self):
        # build dependencies
        print('make -C src qp_iers_bulletin_a.c')
        sp.check_call(
            'PYTHON={} make -C src qp_iers_bulletin_a.c'.format(sys.executable),
            shell=True,
        )

        # write version files
        with open('python/_version.py', 'w') as f:
            f.write('__version__ = ({:d}, {:d}, {:d})\n'.format(*vtup))

        build.run(self)


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
    version=version,
    description='A lightweight quaternion-based library for efficient telescope pointing.',
    author='Alexandra Rahlin',
    packages=['qpoint'],
    package_dir={'qpoint': 'python'},
    ext_modules=[ext_qp],
    cmdclass={
        'build': BuildLib,
        'build_ext': BuildExt,
        'build_sphinx': BuildDoc,
    },
    install_requires=['numpy>=1.10.0'],
    command_options={
        'build_sphinx': {
            'project': ('setup.py', 'qpoint'),
            'version': ('setup.py', version_simple),
            'release': ('setup.py', version_simple),
            'source_dir': ('setup.py', 'docs'),
            'build_dir': ('setup.py', build_dir_name('docs')),
        }
    },
)
