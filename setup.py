from distutils.core import setup, Extension
import os, glob

# hack to avoid recompiling sofa every time
libsofa_file = 'sofa/libsofa_c.a'
if not os.path.exists(libsofa_file):
    os.system('make -C sofa')

# this isn't technically an extension...
# hack to make a shared library to install with the package
ext_qp = Extension('qpoint.libqpoint',[x for x in glob.glob('src/*.c')
                                       if 'test_' not in x],
                   include_dirs=['src','sofa'],
                   extra_compile_args=['-O3','-Wall','-std=c99','-fPIC'],
                   extra_objects=[libsofa_file])

setup(name='qpoint',
      version='0.1',
      description='Pointing for SPIDER',
      author='ASR',
      packages=['qpoint'],
      ext_modules=[ext_qp]
      )
