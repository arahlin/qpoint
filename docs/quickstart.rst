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
<https://github.com/arahlin/qpoint/blob/master/examples/example_point.py>`_ shows some
typical usage of the :py:class:`~qpoint.qpoint_class.QPoint` class, e.g. to construct
timestreams in celestial coordinates from a set of local coordinates, and to
compute offset detector pointing from a boresight location.

The script `example_map.py
<https://github.com/arahlin/qpoint/blob/master/examples/example_map.py>`_ shows some
typical usage of the :py:class:`~qpoint.qmap_class.QMap` class, e.g. to simulate detector
timestreams with offset pointing from an input HEALpix sky realization, and to
bin individual detector timestreams into a HEALpix map.

Memory Initialization
---------------------

Boresight and Detector Pointing
-------------------------------

Mapmaking
---------

Map Simulations
---------------

