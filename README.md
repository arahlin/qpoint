# qpoint

[![](https://badge.fury.io/py/qpoint.svg)](https://badge.fury.io/py/qpoint)
[![](https://github.com/arahlin/qpoint/actions/workflows/documentation.yaml/badge.svg)](https://github.com/arahlin/qpoint/actions/workflows/documentation.yaml)
[![](https://github.com/arahlin/qpoint/actions/workflows/release.yaml/badge.svg)](https://github.com/arahlin/qpoint/actions/workflows/release.yaml)

A lightweight quaternion-based library for telescope pointing.  This library is
forked from the `libactpol` pointing library, originally written by M. Nolta.

Written and maintained by Alexandra Rahlin.

Documentation can be found [here](https://arahlin.github.io/qpoint/).

Other Contributors:

* Steve Benton
* Anne Gambrel
* Carlo Contaldi
* Ivan Padilla
* Matthew Hasselfield

### Requirements

* GNU, Intel or LLVM C compiler
* Python 2.7, Python 3+
* `numpy` library, version 1.10.0 or newer (for python bindings)
* `astropy` library, version 1.2 or newer (optional, for baking in IERS-A data)
* [ERFA C library](https://github.com/liberfa/erfa) (version 2.0.0, based on
  SOFA issue 2021-05-12 bundled with this package)
* [HEALPix C library](http://healpix.sourceforge.net/) (v. 3.31 bundled with
  this package)

### Installation

For most users, it should be sufficient to install the python library from PyPI:

```
pip install qpoint
```

This will install the python bindings and library code compiled with OpenMP
support, if possible (only available with GCC or Intel compilers).

### Usage

To use the pointing library, initialize a `qpoint.QPoint` instance.  When
installed from PyPI, the internal IERS table is left empty.  Use the
`update_iers` argument to update the internal IERS-A table using the IERS
utilities provided by `astropy` (this of course assumes that you have `astropy`
installed on your system):

```
>>> import qpoint as qp
>>> Q = qp.QPoint(update_iers=True)
```

See the [documentation](https://arahlin.github.io/qpoint/) for more details.
See also some example Python code in `examples/`.
