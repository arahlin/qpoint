from distutils.core import setup, Extension
import os, glob

# hack to avoid recompiling sofa every time
libsofa_file = 'sofa/libsofa_c.a'
if not os.path.exists(libsofa_file):
    os.system('make -C sofa')

incl_dirs = ['src','sofa']
extra_obj = [libsofa_file]
extra_args = ['-O3', '-Wall', '-std=c99', '-fPIC','-fopenmp']
libs = ['gomp']

libslarefro_file = 'slarefro/libslarefro.a'
if os.path.exists(libslarefro_file):
    incl_dirs.append('slarefro')
    extra_obj.append(libslarefro_file)
    extra_args.append('-DSLAREFRO')

hpx = os.getenv('HEALPIX')
cfits = os.getenv('CFITSIO')
if hpx:
    hpx = hpx.strip()
    extra_obj.append(os.path.join(hpx,'lib/libchealpix.a'))
    if cfits:
        extra_obj.append(os.path.join(cfits, 'lib/libcfitsio.so'))
    else:
        extra_obj.append('/usr/lib64/libcfitsio.so') #location on Feynman if no CFITS variable set
    incl_dirs.append(os.path.join(hpx,'include'))
else:
    libs.append('chealpix')

# this isn't technically an extension...
# hack to make a shared library to install with the package
ext_qp = Extension('qpoint.libqpoint',[x for x in glob.glob('src/*.c')
                                       if 'test_' not in x],
                   include_dirs=incl_dirs,
                   extra_compile_args=extra_args,
                   libraries=libs,
                   extra_objects=extra_obj)

setup(name='qpoint',
      version='0.1',
      description='Pointing for SPIDER',
      author='ASR',
      packages=['qpoint'],
      package_dir = {'qpoint': 'python'},
      ext_modules=[ext_qp]
      )
