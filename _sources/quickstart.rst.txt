Getting Started
===============

Dependencies
------------

The following Python package dependencies should be installed on your system.

* ``numpy`` (required)
* ``healpy`` (optional, for handling healpix maps)
* ``astropy`` (optional, for updating internal IERS-A tables)

Installation
------------

To use the latest version of QPoint, install via ``pip``::

  $ pip install qpoint

Otherwise, for developement work, clone the repository
(`<https://github.com/arahlin/qpoint>`_) and install manually::

  $ python setup.py install

Or, if installing on a shared machine::

  $ python setup.py install --user

Examples
--------

Several example scripts are provided in the `examples/
<https://github.com/arahlin/qpoint/blob/master/examples/>`_ directory.

The script `example_point.py
<https://github.com/arahlin/qpoint/blob/master/examples/example_point.py>`_
shows some typical usage of the :py:class:`~qpoint.qpoint_class.QPoint` class,
e.g. to construct timestreams in celestial coordinates from a set of local
coordinates, and to compute offset detector pointing from a boresight location.

The script `example_map.py
<https://github.com/arahlin/qpoint/blob/master/examples/example_map.py>`_ shows
some typical usage of the :py:class:`~qpoint.qmap_class.QMap` class, e.g. to
simulate detector timestreams with offset pointing from an input HEALpix sky
realization, and to bin individual detector timestreams into a HEALpix map.

Algorithms
----------

`This document <./qpoint.pdf>`_ describes the quaternion-based algorithm used
for computing equatorial right ascension and declination from the local
coordinates of the (potentially moving) observing platform, including the
various approximations that can be made to speed up the calculation.  The
document also covers computing detector pointing offset from the telescope
boresight, using the same quaternion formalism.

`Another document <./mapmaking_hwp.pdf>`_ describes the mapmaking algorithm
(implemented in the :py:class:`~qpoint.qmap_class.QMap` class) for projecting
time-ordered data onto a `HEALPix <https://healpix.sourceforce.io>`_ spherical
coordinate grid to make sky maps of Stokes parameters, typically used by
telescopes designed for observing the cosmic microwave background (e.g. `SPIDER
<https://spider.princeton.edu>`_).
