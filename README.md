# qpoint

A lightweight quaternion-based library for telescope pointing.  This library is forked from the `libactpol` pointing library, originally written by M. Nolta.

Written and maintained by A. Rahlin.

Other Contributors:

* Steve Benton
* Anne Gambrel
* Carlo Contaldi

### Requirements

* GNU, Intel or LLVM C compiler
* Python 2.7 (bindings should be compatible with, but have not been tested in, Python 3.)
* `numpy` library, version 1.10.0 or newer (for python bindings)
* `astropy` library, version 1.2 or newer (optional, for baking in IERS-A data)
* [SOFA C library](http://www.iausofa.org/) (issue 2017-04-20 bundled with this package)
* [HEALPix C library](http://healpix.sourceforge.net/) (v. 3.31 bundled with this package)

### Installation

For most users, it should be sufficient to run the following to install the python library:

```
python setup.py install
```

This will install the python bindings and library code compiled with OpenMP support, if found.

To just build and install the `C` library without OpenMP support:

```
make
make install
```

To enable OpenMP support in the `C` library:

```
make ENABLE_OMP=y
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

### Usage

See example c and python code in `examples/`
