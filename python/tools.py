from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import numpy as np
from ._libqpoint import libqp as qp

__all__ = ['refraction']


def refraction(el, temp, press, hum, freq=150.0):
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

    fvec = np.vectorize(qp.qp_refraction, [np.double])
    delta = fvec(el, temp, press, hum, freq)
    if delta.shape == ():
        return delta[()]
    return delta
