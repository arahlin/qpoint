from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import numpy as _np
from ._libqpoint import libqp as _libqp

def refraction(el, temp, press, hum, freq=150.):
    """
    Standalone function for calculating the refraction correction without
    storing any parameters.  Useful for testing, numpy-vectorized.
    
    Arguments:
    
    el           elevation angle, degrees
    temperature  temperature, Celcius
    pressure     pressure, mbar
    humidity     humidity, fraction
    frequency    array frequency, GHz
    
    Output:
    
    delta        refraction correction, in degrees
    """
    
    fvec = _np.vectorize(_libqp.qp_refraction,[_np.double])
    delta = fvec(el, temp, press, hum, freq)
    if delta.shape == ():
        return delta[()]
    return delta
