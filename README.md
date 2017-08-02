# qpoint

A lightweight quaternion-based library for telescope pointing.  This library is forked from the `libactpol` pointing library, originally written by M. Nolta.

Written and maintained by A. Rahlin.

Other Contributors:

* Steve Benton
* Anne Gambrel
* Carlo Contaldi

### Requirements

* GNU or Intel C compiler
* Python 2.7 (bindings should be compatible with, but have not been tested in, Python 3.)
* `astropy` library, version 1.2 or newer (optional, for baking in IERS-A data)
* [SOFA C library](http://www.iausofa.org/) (issue 2017-04-20 bundled with this package)
* [HEALPix C library](http://healpix.sourceforge.net/) (v. 3.31 bundled with this package)

### Installation

To install both the `C` library (`libqpoint.a`) and `python` bindings:

```
make
sudo make install
```

For just the `C` library:

```
make qpoint
sudo make install-qpoint
```

For just the `python` bindings (the `C` library is bundled internally):

```
make python
sudo make install-python
```

To install to a user directory (by default this is `$HOME/.local`):

```
make
make install-user # or install-qpoint-user or install-python-user
```

To install to a specific prefix (the directory that will contain the `lib/`
and `include/` subfolders):

```
make
sudo make PREFIX=/your/install/prefix install
```

### Usage

See example c and python code in `examples/`
