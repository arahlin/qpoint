from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import numpy as _np
from ._libqpoint import libqp as _libqp

def refraction(el, temp, press, hum, freq=150.):
    """
    Standalone function for calculating the refraction correction without
    storing any parameters.  Useful for testing, numpy-vectorized.

    Arguments
    ---------
    el : array_like
        Observer elevation angle, degrees
    temperature : array_like
        Ambient temperature, Celcius
    pressure : array_like
        Ambient pressure, mbar
    humidity : array_like
        Relative humidity, fraction
    frequency : array_like
        Observing frequency, GHz

    Returns
    -------
    delta : array_like
        Refraction correction, in degrees
    """

    fvec = _np.vectorize(_libqp.qp_refraction,[_np.double])
    delta = fvec(el, temp, press, hum, freq)
    if delta.shape == ():
        return delta[()]
    return delta
