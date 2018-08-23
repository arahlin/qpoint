"""A lightweight quaternion-based library for efficient telescope pointing.

Based on M. Nolta's libactpol.
Uses the SOFA Software Collection, available from http://www.iausofa.org/
"""

from ._version import __version__

def version():
    """Print version string"""
    return __version__

from . import tools
from .qpoint_class import *
from .qmap_class import *
