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

Some users may want to build the `C` library for linking to their own projects,
or have more control over the build process.

To just build and install the `C` library without OpenMP support:

```
make
make install
```

To enable OpenMP support in the `C` library:

```
make ENABLE_OMP=y
```

To build a shared library instead of a static one:

```
make ENABLE_SHARED=y
```

To build a "lite" version of the library without support for OpenMP or any of the healpix mapmaking backend:

```
make ENABLE_LITE=y
```

To install to a user directory (by default this is `$HOME/.local`):

```
make install-user
```

To install to a specific prefix (the directory that will contain the `lib/`
and `include/` subfolders):

```
make PREFIX=/your/install/prefix install
```
