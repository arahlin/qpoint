"""qpoint

A lightweight library for efficient pointing.

Based on M. Nolta's libactpol.
Uses the SOFA Software Collection, available from http://www.iausofa.org/
"""

from _version import __version__

def version():
    return __version__

import tools
from qpoint_class import *
from qmap_class import *
