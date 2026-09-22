## OpenMP on macOS

The Python extensions build without OpenMP on macOS by default, because
Apple's clang does not ship an OpenMP runtime and linking one causes a hard
conflict: `healpy` bundles its own copy of `libomp.dylib`, and with two LLVM
OpenMP runtimes in a process, whichever loads second segfaults the moment it
spawns threads. There is no safe import order, no diagnostic, and no
environment variable that avoids it.

To get parallel `tod2map`/`map2tod` on macOS anyway, build with GCC, which
uses `libgomp` and does coexist with `healpy`:

```
brew install gcc
CC=gcc-15 CXX=g++-15 pip install -e .
```

This is worth roughly 2.4x on `tod2map`. Note that it changes results in the
last bits or so relative to an Apple clang build, since the two toolchains
differ in codegen and `libm`. Linux is unaffected throughout: OpenMP is
enabled there automatically.

## Building the C library

There is no standalone `make` build. It went away with the move to PEP 517
packaging, which removed the root `Makefile` along with the ones under
`erfa/` and `chealpix/`, and the library is now built only as part of the
Python extension.

The sources are self-contained, though -- `erfa` and `chealpix` are bundled
in `src/`, and nothing outside libc and libm is needed -- so a project that
wants the `C` on its own can compile them directly:

```
cc -std=c99 -O3 -fPIC -Isrc -c src/erfa.c src/chealpix.c src/qp_*.c \
    src/qpoint.c src/quaternion.c src/sincos.c
ar rcs libqpoint.a *.o
```

Add `-fopenmp` for a parallel `tod2map`/`map2tod`, subject to the macOS
caveat above. The headers a caller needs are `src/qpoint.h`,
`src/quaternion.h` and `src/vec3.h`.

## The C tests

Two test programs come with the bundled third-party code and are not built
by anything, so they have to be compiled by hand:

```
cc -std=c99 -O2 -Isrc src/test_erfa.c src/erfa.c -lm -o test_erfa
cc -std=c99 -O2 -Isrc src/test_chealpix.c src/chealpix.c -lm -o test_chealpix
```

Each exits nonzero on failure. They cover the vendored `erfa` and
`chealpix` only; everything else is covered by the Python suite, which
`pytest` runs.
