from ._libqpoint2 import refraction as _refraction

__all__ = ["refraction"]


def refraction(el, temp, press, hum, freq=150.0):
    """
    Standalone function for calculating the refraction correction without
    storing any parameters.  Useful for testing, numpy-vectorized.

    Arguments
    ---------
    el : array_like
        Observer elevation angle, degrees
    temp : array_like
        Ambient temperature, Celcius
    press : array_like
        Ambient pressure, mbar
    hum : array_like
        Relative humidity, fraction
    freq : array_like
        Observing frequency, GHz

    Returns
    -------
    delta : array_like
        Refraction correction, in degrees

    Notes
    -----
    Unlike :meth:`QPoint.refraction`, this touches none of the stored
    weather and caches nothing.

    The loop is in the binding rather than in numpy: qpoint reaches the
    same C function through `np.vectorize`, which is a Python-level loop
    over a ctypes call. A scalar argument standing in for a full-length
    one is read at stride 0, so passing an elevation array with scalar
    weather allocates nothing beyond the output.
    """
    return _refraction(el, temp, press, hum, freq)
