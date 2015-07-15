# numpy.distutils handles intel compiler
# need to pass --compiler=intel (or add to [build_ext] section of setup.cfg)
from numpy.distutils.core import setup, Extension
import os, glob
from warnings import warn

# hack to avoid recompiling sofa every time
libsofa_file = 'sofa/libsofa_c.a'
if not os.path.exists(libsofa_file):
    os.system('make -C sofa')

src = [x for x in glob.glob('src/*.c')]
incl_dirs = ['src','sofa']
extra_obj = [libsofa_file]
extra_args = ['-O3', '-Wall', '-std=c99', '-fPIC']

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

hpx = os.getenv('HEALPIX','').strip()
if hpx:
    extra_obj.append(os.path.join(hpx,'lib/libchealpix.a'))
    incl_dirs.append(os.path.join(hpx,'include'))

    # deal with cfitsio library using the standard environment variables
    ext_cfits = os.getenv('EXTERNAL_CFITSIO','').strip()
    if ext_cfits == 'yes':
        ext_pref = os.getenv('CFITSIO_EXT_PREFIX','').strip()
        if ext_pref:
            extra_obj.append(os.path.join(ext_pref,'lib/libcfitsio.a'))
            incl_dirs.append(os.path.join(ext_pref,'include'))
        else:
            ext_lib = os.getenv('CFITSIO_EXT_LIB','').strip()
            if ext_lib:
                extra_obj.append(ext_lib)
            ext_inc = os.getenv('CFITSIO_EXT_INC','').strip()
            if ext_inc:
                incl_dirs.append(ext_inc)
    else:
        # location on Feynman if no CFITS variable set
        warn('Hardcoded path to libcfitsio, consider using the appropriate ' \
             'environment variables (see README)', Warning)
        cflib = '/usr/lib64/libcfitsio.so'
        if os.path.exists(cflib):
            extra_obj.append(cflib)
else:
    warn('Blindly trying to link against chealpix, consider using the ' \
         'appropriate environment variables (see README)', Warning)
    libs.append('chealpix')

# this isn't technically an extension...
# hack to make a shared library to install with the package
ext_qp = Extension('qpoint.libqpoint', src,
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
