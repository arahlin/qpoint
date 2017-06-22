"""qpoint

A lightweight library for efficient pointing.

Based on M. Nolta's libactpol.
Uses the SOFA Software Collection, available from http://www.iausofa.org/
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from ._version import __version__

def version():
    return __version__

from . import tools
from .qpoint_class import *
from .qmap_class import *
