from distutils.core import setup, Extension
import os, glob

# hack to avoid recompiling sofa every time
libsofa_file = 'sofa/libsofa_c.a'
if not os.path.exists(libsofa_file):
    os.system('make -C sofa')

libslarefro_file = 'slarefro/libslarefro.a'
if not os.path.exists(libslarefro_file):
    os.system('make -C slarefro')

# this isn't technically an extension...
# hack to make a shared library to install with the package
ext_qp = Extension('qpoint.libqpoint',[x for x in glob.glob('src/*.c')
                                       if 'test_' not in x],
                   include_dirs=['src', 'sofa', 'slarefro'],
                   extra_compile_args=['-O3', '-Wall', '-std=c99', '-fPIC','-fopenmp'],
                   libraries=['gomp','chealpix'],
                   extra_objects=[libsofa_file, libslarefro_file])

setup(name='qpoint',
      version='0.1',
      description='Pointing for SPIDER',
      author='ASR',
      packages=['qpoint'],
      package_dir = {'qpoint': 'python'},
      ext_modules=[ext_qp]
      )
