import functools
import inspect
from contextlib import contextmanager

import numpy as np

from . import _libqpoint2 as lib

RATE_PARAMS = [
    "daber",
    "lonlat",
    "wobble",
    "dut1",
    "erot",
    "npb",
    "aaber",
    "ref",
    "daber_inv",
    "lonlat_inv",
    "wobble_inv",
    "dut1_inv",
    "erot_inv",
    "npb_inv",
    "aaber_inv",
    "ref_inv",
]

OPTION_PARAMS = [
    "accuracy",
    "mean_aber",
    "fast_aber",
    "fast_math",
    "polconv",
    "pix_order",
    "interp_pix",
    "fast_pix",
    "error_missing",
    "nan_missing",
    "interp_missing",
    "num_threads",
]

WEATHER_PARAMS = ["temperature", "pressure", "humidity", "frequency"]

DOUBLE_PARAMS = ["ref_delta", "dut1"]

# Insertion order matters: get() with no arguments returns the groups in
# this order, and qpoint builds the same dict from qp_funcs.
PARAM_GROUPS = {
    "rates": ["rate_" + r for r in RATE_PARAMS],
    "options": OPTION_PARAMS,
    "weather": WEATHER_PARAMS,
    "params": DOUBLE_PARAMS,
}

ALL_PARAMS = [k for keys in PARAM_GROUPS.values() for k in keys]

__all__ = ["QPoint", "qp_settings"]


# Arguments are passed through to the bindings untouched. Every per-sample
# argument may be a scalar or a full-length float64 array, and the binding
# layer resolves the two against each other with a stride-0 column rather
# than materializing broadcast copies. Arrays are never converted: a wrong
# dtype or layout raises instead of being silently copied.
#
# The one exception is the in-place coordinate rotations below, which write
# their results back into their inputs. Those genuinely do need real,
# writable, full-length arrays, so a scalar has to be materialized.


def _prep_coord_rotation(ra, dec, pa, sin2psi, cos2psi, inplace):
    """
    Validate and produce writable arrays for an in-place coordinate rotation.

    Returns (do_pa, arrays, all_scalar), where do_pa selects pa mode over
    sin2psi/cos2psi mode and all_scalar says every supplied coordinate was
    a true scalar. The rotation needs writable 1-d arrays either way, so
    the caller degrades the result afterwards rather than here.
    """
    do_pa = pa is not None or (sin2psi is None and cos2psi is None)
    if pa is not None and (sin2psi is not None or cos2psi is not None):
        raise KeyError("supply either pa or sin2psi/cos2psi, not both")
    if not do_pa and (sin2psi is None) != (cos2psi is None):
        raise KeyError("both sin2psi and cos2psi are required")

    args = (ra, dec, pa) if do_pa else (ra, dec, sin2psi, cos2psi)
    all_scalar = all(np.ndim(a) == 0 for a in args if a is not None)
    arrays = np.broadcast_arrays(
        *[np.asarray(a if a is not None else 0.0, dtype=float) for a in args]
    )
    if not inplace:
        arrays = [np.array(a) for a in arrays]
    prepped = [np.require(np.atleast_1d(a), float, ["C", "A", "W"]) for a in arrays]
    return do_pa, prepped, all_scalar


def _settings_note(defaults):
    """The Notes text qp_settings appends to a method it wraps."""
    text = (
        "Any parameter accepted by :meth:`QPoint.set` may also be given as a\n"
        "keyword here. It applies for the duration of this call and is put\n"
        "back afterwards."
    )
    if defaults:
        named = ", ".join("``{}={!r}``".format(k, v) for k, v in defaults.items())
        text += "\n\nThis method runs with {} unless the call says otherwise.".format(
            named
        )
    return text


def _append_note(doc, text):
    """
    Append `text` to a docstring as a numpydoc `Note` section.

    `Note`, not `Notes`: napoleon renders the singular as an admonition,
    so the keywords are called out in the built API reference rather than
    trailing a method's own prose, and a method that keeps notes of its
    own keeps them separate.

    The indentation is taken from the docstring itself rather than
    assumed: Python 3.13 strips the common leading whitespace at compile
    time and earlier versions do not, so it is empty on one and a method
    body's worth on the other.
    """
    if not doc:
        return doc
    widths = [len(l) - len(l.lstrip()) for l in doc.splitlines()[1:] if l.strip()]
    pad = " " * min(widths) if widths else ""
    body = "\n".join(pad + l if l else "" for l in text.splitlines())
    return "{d}\n\n{p}Note\n{p}----\n{b}\n{p}".format(d=doc.rstrip(), p=pad, b=body)


def qp_settings(_method=None, **defaults):
    """
    Give a method the per-call parameter keywords.

    Every method of :class:`QPoint` takes the parameters as keywords,
    applied for that call and restored afterwards. This is what does it,
    and it is exported so that a subclass adding its own methods gets the
    same behaviour without repeating the mechanism::

        from qpoint2 import QPoint, qp_settings

        class MyPoint(QPoint):
            @qp_settings
            def my_scan(self, az, el, ctime):
                return self.azel2bore(az, el, None, None, 0.0, 0.0, ctime)

        MyPoint().my_scan(az, el, ctime, accuracy="low", fast_math=True)

    The keywords naming parameters go to :meth:`QPoint.settings`; the
    rest go to the method, which reports an unexpected one the way any
    Python function does. A `Note` section describing them is appended to
    the wrapped method's docstring — napoleon renders it as an admonition
    — so no method writes that out itself.

    Arguments
    ---------
    **defaults
        Parameters to set for this method's calls, under anything the
        caller passes. :meth:`QPoint.rotate_map` is declared
        `@qp_settings(interp_pix=True)`, so it interpolates unless told
        otherwise, whatever the object is set to.

    Notes
    -----
    The split is by name against the parameter list, so **a decorated
    method must not take an argument named after a parameter** -- the
    keyword would be routed to the settings and the argument would never
    be filled. Where the method needs the value, read it back with
    :meth:`QPoint.get_param` inside the call, where the override is in
    force.
    """

    def decorate(method):
        sig = inspect.signature(method)

        @functools.wraps(method)
        def wrapper(self, *args, **kwargs):
            params = dict(defaults)
            params.update({k: kwargs.pop(k) for k in list(kwargs) if k in ALL_PARAMS})
            with self.settings(**params):
                return method(self, *args, **kwargs)

        # the keywords are the decorator's, so the note describing them is
        # its to write -- no method should be repeating it
        wrapper.__doc__ = _append_note(method.__doc__, _settings_note(defaults))
        # the decorator swallowed them, but they are still part of the API
        wrapper.__signature__ = sig.replace(
            parameters=[
                *sig.parameters.values(),
                inspect.Parameter("kwargs", inspect.Parameter.VAR_KEYWORD),
            ]
        )
        return wrapper

    return decorate(_method) if _method is not None else decorate


class QPoint(lib.Pointing):
    """Quaternion-based telescope pointing."""

    def __init__(self, update_iers=False, **kwargs):
        super().__init__()
        if update_iers:
            self.update_bulletin_a()
        self.set(**kwargs)

    # ---- Parameters ----

    def set(self, **kwargs):
        """
        Set any number of parameters.

        Raises KeyError for a name that is not a parameter, so a typo is
        reported rather than silently doing nothing. qpoint ignores it.

        Arguments
        ---------
        rate_daber : {'never', 'once', 'always'}, or float
            Rate at which the diurnal aberration correction is applied in seconds
            (NB: this can only be applied always or never)
        rate_lonlat : {'never', 'once', 'always'}, or float
            Rate at which observer's lon and lat are updated
        rate_wobble : {'never', 'once', 'always'}, or float
            Rate at which the polar motion correction is updated
            (NB: these are not estimated for dates beyond a year from now)
        rate_dut1 : {'never', 'once', 'always'}, or float
            Rate at which the ut1-utc correction is updated
            (NB: this is not estimated for dates beyond a year from now)
        rate_erot : {'never', 'once', 'always'}, or float
            Rate at which the earth's rotation angle is updated
        rate_npb : {'never', 'once', 'always'}, or float
            Rate at which the nutation/precession/frame-bias terms are updated
        rate_aaber : {'never', 'once', 'always'}, or float
            Rate at which the annual aberration correction
            (due to the earth's orbital velocity) is updated
        rate_ref : {'never', 'once', 'always'}, or float
            Rate at which the refaction correction is updated (NB: this
            correction can also be updated manually -- see :meth:`refraction`)
        accuracy : 'low' or 'high'
            If 'low', use a truncated form (2000b) for the NPB correction,
            which is much faster but less accurate. If 'high' (default), use
            the full 2006/2000a form.
        mean_aber : bool
            If True, apply the aberration correction as an average for the
            entire field of view.  This is gives a 1-2 arcsec deviation
            at the edges of the SPIDER field of view.
        fast_aber : bool
            If True, apply the aberration correction using the exact
            angle and quaternion rotation.  Otherwise, use a small-angle
            approximation for the aberration quaternion (default).
        fast_math : bool
            If True, use polynomial approximations for trig functions
        polconv : 'cosmo' or 'iau'
            Specify the 'cosmo' or 'iau' polarization convention
        pix_order : 'nest' or 'ring'
            HEALPix pixel ordering
        interp_pix : bool
            If True, interpolate between pixels in scanning the source map.
        fast_pix : bool
            If True, use `vec2pix` to get pixel number directly from the
            quaternion instead of `ang2pix` from ra/dec.
        error_missing : bool
            If True, raise an error if reading/writing missing pixels.
        nan_missing : bool
            If True, fill samples from missing pixels with NaN.
            Only used if `error_missing` is False.
        interp_missing : bool
            If True and `interp_pix` is True, drop missing neighbors
            and reweight remaining neighbors.  Overrides `nan_missing`.
        num_threads : bool
             Number of openMP threads to use for mapmaking.
        temperature : float
            Ambient temperature, Celcius. For computing refraction corrections.
        pressure : float
            Ambient pressure, mbar. For computing refraction corrections.
        humidity : float
            Relative humidity, fraction. For computing refraction corrections.
        frequency : float
            Observer frequency, GHz. For computing refraction corrections.
        dut1 : float
            UT1 correction
        ref_delta : float
            Refraction correction
        """
        for k, v in kwargs.items():
            self.set_param(k, v)

    def get(self, *args):
        """
        Return parameter values.

        No arguments returns every parameter, grouped: a dict of the four
        group names, each holding that group's parameters. A group name
        ('rates', 'options', 'weather', 'params') returns just that group,
        flat. A single key returns the bare value. Several keys return a
        dict keyed by what was asked for, so a group among them stays
        nested under its own name rather than being merged in.

        Arguments
        ---------
        *args : str
            Parameter names, or the group names 'rates', 'options', 'weather'
            and 'params'. With no arguments, every parameter is returned.

        Returns
        -------
        values : dict or scalar
            The bare value for a single parameter name, otherwise a dict
            keyed by parameter name, or by group name where a group was
            asked for.
        """
        if not args:
            return {
                group: {k: self.get_param(k) for k in keys}
                for group, keys in PARAM_GROUPS.items()
            }
        if len(args) == 1:
            key = args[0]
            if key in PARAM_GROUPS:
                return {k: self.get_param(k) for k in PARAM_GROUPS[key]}
            return self.get_param(key)
        out = {}
        for key in args:
            if key in PARAM_GROUPS:
                out[key] = {k: self.get_param(k) for k in PARAM_GROUPS[key]}
            else:
                out[key] = self.get_param(key)
        return out

    @contextmanager
    def settings(self, **kwargs):
        """
        Temporarily apply parameters, restoring them on exit.

        Arguments
        ---------
        **kwargs
            Parameters to apply, as :meth:`set` accepts them. Each is read
            first and restored on exit; an unrecognized name raises.

        Notes
        -----
        A context manager, so the parameters last as long as the `with`
        block and are put back afterwards even if the body raises. Every
        method taking `**kwargs` uses this internally, which is how a
        per-call parameter override works.
        """
        if not kwargs:
            yield
            return
        # reads every name first, so an unrecognized one raises before
        # anything has been changed
        saved = {k: self.get_param(k) for k in kwargs}
        self.set(**kwargs)
        try:
            yield
        finally:
            self.set(**saved)

    def refraction(self, *args, **kwargs):
        """
        Update or return the refraction correction, in degrees.

        Called with a scalar or ``delta=``, sets the stored correction.
        Called with ``q=``, computes it from that boresight quaternion.
        Weather keywords are applied first.

        Arguments
        ---------
        q : quaternion or array of quaternions
            Observer orientation in horizon coordinates
        temperature : float
            Ambient temperature, Celcius
        pressure : float
            Ambient pressure, mbar
        humidity : float
            Ambient relative humidity, fraction
        frequency : float
            Observing frequency, GHz
        delta : float
            The refraction correction itself, in degrees

        Returns
        -------
        delta : array_like
            Refraction correction computed at each input orientation

        Notes
        -----
        If `q` is given, then the refraction correction in degrees
        is calculated, stored and returned after updating any other given
        parameters. Otherwise, the correction is returned w/out recalculating.

        Alternatively, if a single numerical argument, or the `delta` keyword
        argument is given, then the correction is stored with this value
        instead of being recalculated.

        Numpy-vectorized for the `q` argument.  Note that this is not
        an efficient vectorization, and only the last calculated value is
        stored for use in the coordinate conversion functions.
        """
        if len(args) == 1 and not kwargs and np.isscalar(args[0]):
            self.set_param("ref_delta", args[0])
            return self.get_param("ref_delta")

        arg_names = ["q"] + WEATHER_PARAMS
        for name, val in zip(arg_names, args):
            kwargs[name] = val

        if "delta" in kwargs:
            self.set_param("ref_delta", kwargs.pop("delta"))
            return self.get_param("ref_delta")

        q = kwargs.pop("q", None)
        for k, v in kwargs.items():
            self.set_param(k, v)

        if q is not None:
            return super().update_ref(q)
        return self.get_param("ref_delta")

    # ---- Time ----

    @qp_settings
    def gmst(self, ctime):
        """
        Greenwich mean sidereal time, in hours.

        Arguments
        ---------
        ctime : array_like
            Unix time in seconds UTC

        Returns
        -------
        gmst : array_like
            Greenwich mean sidereal time of the observer
        """
        return super().gmst(ctime)

    @qp_settings
    def lmst(self, ctime, lon):
        """
        Local mean sidereal time, in hours.

        Arguments
        ---------
        ctime : array_like
            Unix time in seconds UTC
        lon : array_like
            Observer longitude (degrees)

        Returns
        -------
        lmst : array_like
            Local mean sidereal time of the observer
        """
        return super().lmst(ctime, lon)

    # ---- Quaternion construction ----

    def det_offset(self, delta_az, delta_el, delta_psi):
        """
        Quaternion for the requested detector centroid offset from boresight.

        Arguments
        ---------
        delta_az : array_like
            Azimuthal centroid offset of the detector in degrees
        delta_el : array_like
            Elevation centroid offset of the detector in degrees
        delta_psi : array_like
            Polarization offset of the detector from vertical in degrees

        Returns
        -------
        q : array_like
            Detector centroid offset quaternion for each detector
        """
        return lib.det_offset(delta_az, delta_el, delta_psi)

    def hwp_quat(self, theta):
        """
        Quaternion for rotation by 2*theta (physical HWP angle).

        Arguments
        ---------
        theta : array_like
            HWP physical angle in degrees

        Returns
        -------
        q : array_like
            Quaternion for each hwp angle
        """
        return lib.hwp_quat(theta)

    @qp_settings
    def radecpa2quat(self, ra, dec, pa):
        """
        Quaternion from ra/dec/position angle, in degrees.

        Arguments
        ---------
        ra : array_like
            Right ascension angle
        dec : array_like
            Declination angle
        pa : array_like
            Position angle

        Returns
        -------
        q : array_like
            Quaternion constructed from the input angles.
        """
        return super().radecpa2quat(ra, dec, pa)

    @qp_settings
    def quat2radecpa(self, quat):
        """
        ra/dec/position angle from a quaternion, in degrees.

        Arguments
        ---------
        quat : quaternion or array of quaternions
            Orientation quaternions, of shape (N, 4).

        Returns
        -------
        ra : array_like
            Right ascension in degrees.
        dec : array_like
            Declination in degrees.
        pa : array_like
            Position angle in degrees.
        """
        return super().quat2radecpa(quat)

    # ---- Boresight ----

    @qp_settings
    def azel2bore(self, az, el, pitch, roll, lon, lat, ctime):
        """
        Boresight quaternion from az/el/pitch/roll/lon/lat/ctime.

        Arguments
        ---------
        az : array_like
            Boresight azimuth in degrees
        el : array_like
            Boresight elevation in degrees
        pitch : array_like
            Boresight pitch in degrees.  If `None`, this term is ignored.
        roll : array_like
            Boresight roll in degrees.  If `None`, this term is ignored.
        lon : array_like
            Observer longitude in degrees
        lat : array_like
            Observer latitude in degrees
        ctime : array_like
            Unix time in seconds UTC

        Returns
        -------
        q : array_like
            Nx4 numpy array of quaternions for each supplied timestamp.
        """
        return super().azel2bore(az, el, None, pitch, roll, lon, lat, ctime)

    @qp_settings
    def azelpsi2bore(self, az, el, psi, pitch, roll, lon, lat, ctime):
        """
        Boresight quaternion, accounting for FPU boresight rotation psi.

        Arguments
        ---------
        az : array_like
            Boresight azimuth in degrees
        el : array_like
            Boresight elevation in degrees
        psi : array_like
            Boresight rotation angle in degrees
        pitch : array_like
            Boresight pitch in degrees.  If `None`, this term is ignored.
        roll : array_like
            Boresight roll in degrees.  If `None`, this term is ignored.
        lon : array_like
            Observer longitude in degrees
        lat : array_like
            Observer latitude in degrees
        ctime : array_like
            Unix time in seconds UTC

        Returns
        -------
        q : array_like
            Nx4 numpy array of quaternions for each supplied timestamp.
        """
        return super().azel2bore(az, el, psi, pitch, roll, lon, lat, ctime)

    # ---- Coordinate conversion ----

    @qp_settings
    def bore2radec(
        self, q_off, ctime, q_bore, q_hwp=None, sindec=False, return_pa=False
    ):
        """
        Sky coordinates for a detector offset and a boresight timestream.

        Returns (ra, dec, sin2psi, cos2psi), or (ra, dec, pa) if return_pa.
        dec is replaced by sin(dec) if sindec. Unlike the C library, the two
        options are independent and may be combined.

        Arguments
        ---------
        q_off : quaternion
            Detector offset quaternion for a single detector, calculated using
            :meth:`det_offset`.
        ctime : array_like
            Unix time in seconds UTC, broadcastable to shape (N,),
            the long dimension of `q_bore`.
        q_bore : quaternion or array of quaternions
            Nx4 array of quaternions encoding the boresight orientation on the
            sky (as output by :meth:`azel2radec` or equivalent)
        q_hwp : quaternion or array of quaternions, optional
            HWP angle quaternions calculated using :meth:`hwp_quat`.  Must be
            broadcastable to the same shape as `q_bore`.
        sindec : bool, optional
            If `True`, return sin(dec) instead of dec in degrees
            (default False).
        return_pa : bool, optional
            If `True`, return pa instead of sin2psi / cos2psi

        Returns
        -------
        ra : array_like
            Detector right ascension in degrees
        dec/sindec : array_like
            Detector declination in degrees or sin(dec) if `sindec` is `True`.
        pa/sin2psi : array_like
            Detector polarization orientation if `return_pa` is `True`, or
            sin(2*pa) if `return_pa` is `False`.
        cos2psi : array_like
            detector polarization orientation cos(2*pa), if `return_pa` is `False`.
        """
        if ctime is None:
            if not self.get_param("mean_aber"):
                raise ValueError("ctime required if mean_aber is False")
            # a stride-0 zero column, so no array is allocated
            ctime = 0.0
        return super().bore2radec(
            q_off, ctime, q_bore, q_hwp, bool(sindec), bool(return_pa)
        )

    @qp_settings
    def bore2azel(self, q_bore, lon, lat, ctime):
        """
        Horizon coordinates from a boresight quaternion timestream.

        Arguments
        ---------
        q_bore : array_like
            Nx4 array of boresight quaternions (as output by :meth:`azel2bore`).
        lon : array_like
            Observer longitude in degrees.
        lat : array_like
            Observer latitude in degrees.
        ctime : array_like
            Unix time in seconds UTC

        Returns
        -------
        az : array_like
            Azimuth in degrees
        el : array_like
            Elevation in degrees
        pa : array_like
            Position angle in horizon coordinates
        """
        return super().bore2azel(q_bore, lon, lat, ctime)

    @qp_settings
    def azel2radec(
        self,
        delta_az,
        delta_el,
        delta_psi,
        az,
        el,
        pitch,
        roll,
        lon,
        lat,
        ctime,
        hwp=None,
        sindec=False,
        return_pa=False,
    ):
        """
        Sky coordinates from az/el boresight and a detector offset.

        Arguments
        ---------
        delta_az : float
            Azimuthal offset of the detector in degrees
        delta_el : float
            Elevation offset of the detector in degrees
        delta_psi : float
            Polarization offset of the detector in degrees
        az : array_like
            Boresight azimuth in degrees
        el : array_like
            Boresight elevation in degrees
        pitch : array_like
            Boresight pitch in degrees.  If None, this term is ignored.
        roll : array_like
            Boresight roll in degrees.  If None, this term is ignored.
        lon : array_like
            Observer longitude in degrees.
        lat : array_like
            Observer latitude in degrees.
        ctime : array_like
            Unix time in seconds UTC
        hwp : array_like, optional
            HWP angles in degrees
        sindec : bool, optional
            If `True`, return sin(dec) instead of dec in degrees (default False)
        return_pa : bool, optional
            If `True`, return pa instead of sin2psi/cos2psi

        Returns
        -------
        ra : array_like
            Detector right ascension in degrees
        dec/sindec : array_like
            Detector declination in degrees
        pa : array_like
            Detector position angle, if `return_pa` is True
        sin2psi : array_like
            Detector polarization orientation, if `return_pa` is False
        cos2psi : array_like
            Detector polarization orientation, if `return_pa` is False
        """
        return super().azel2radec(
            delta_az,
            delta_el,
            delta_psi,
            az,
            el,
            None,
            pitch,
            roll,
            lon,
            lat,
            ctime,
            hwp,
            bool(sindec),
            bool(return_pa),
        )

    @qp_settings
    def azelpsi2radec(
        self,
        delta_az,
        delta_el,
        delta_psi,
        az,
        el,
        psi,
        pitch,
        roll,
        lon,
        lat,
        ctime,
        hwp=None,
        sindec=False,
        return_pa=False,
    ):
        """
        Sky coordinates, accounting for FPU boresight rotation psi.

        Arguments
        ---------
        delta_az : float
            Azimuthal offset of the detector in degrees
        delta_el : float
            Elevation offset of the detector in degrees
        delta_psi : float
            Polarization offset of the detector in degrees
        az : array_like
            Boresight azimuth in degrees
        el : array_like
            Boresight elevation in degrees
        psi : array_like
            Boresight rotation in degrees
        pitch : array_like
            Boresight pitch in degrees.  If None, this term is ignored.
        roll : array_like
            Boresight roll in degrees.  If None, this term is ignored.
        lon : array_like
            Observer longitude in degrees.
        lat : array_like
            Observer latitude in degrees.
        ctime : array_like
            Unix time in seconds UTC
        hwp : array_like, optional
            HWP angles in degrees
        sindec : bool, optional
            If `True`, return sin(dec) instead of dec in degrees (default False)
        return_pa : bool, optional
            If `True`, return pa instead of sin2psi/cos2psi

        Returns
        -------
        ra : array_like
            Detector right ascension in degrees
        dec/sindec : array_like
            Detector declination in degrees
        pa : array_like
            Detector position angle, if `return_pa` is True
        sin2psi : array_like
            Detector polarization orientation, if `return_pa` is False
        cos2psi : array_like
            Detector polarization orientation, if `return_pa` is False
        """
        return super().azel2radec(
            delta_az,
            delta_el,
            delta_psi,
            az,
            el,
            psi,
            pitch,
            roll,
            lon,
            lat,
            ctime,
            hwp,
            bool(sindec),
            bool(return_pa),
        )

    @qp_settings
    def radec2azel(self, ra, dec, pa, lon, lat, ctime):
        """
        Horizon coordinates from equatorial coordinates.

        Arguments
        ---------
        ra : array_like
            Right ascension angle
        dec : array_like
            Declination angle
        pa : array_like
            Position angle in equatorial coordinates
        lon : array_like
            Observer longitude in degrees.
        lat : array_like
            Observer latitude in degrees.
        ctime : array_like
            Unix time in seconds UTC

        Returns
        -------
        az : array_like
            Azimuth in degrees
        el : array_like
            Elevation in degrees
        hpa : array_like
            Position angle in horizon coordinates
        """
        return super().radec2azel(ra, dec, pa, lon, lat, ctime)

    # ---- Pixelization ----

    @qp_settings
    def radec2pix(self, ra, dec, nside=256):
        """
        HEALPix pixel numbers for the given sky coordinates.

        Arguments
        ---------
        ra : array_like
            Right ascension angle
        dec : array_like
            Declination angle
        nside : int
            HEALpix resolution parameter

        Returns
        -------
        pix : array_like
            Pixel number(s) corresponding to the input positions(s).
        """
        return super().radec2pix(ra, dec, nside)

    @qp_settings
    def quat2pix(self, quat, nside=256, pol=True):
        """
        Pixel number and polarization angle for a quaternion.

        Arguments
        ---------
        quat : quaternion or array of quaternions
            Pointing orientation(s)
        nside : int, optional
            HEALpix resolution parameter
        pol : bool, optional
            If True, return sin2psi and cos2psi along with the pixel number(s)

        Returns
        -------
        pix : array_like
            Pixel number(s) for the given input quaternion(s)
        sin2psi : array_like
        cos2psi : array_like
            Polarization coefficients, if `pol` is `True`.
        """
        pix, sin2psi, cos2psi = super().quat2pix(quat, nside, False)
        return (pix, sin2psi, cos2psi) if pol else pix

    @qp_settings
    def quat2pixpa(self, quat, nside=256):
        """
        Pixel number and position angle for a quaternion.

        Arguments
        ---------
        quat : quaternion or array of quaternions
            Orientation quaternions, of shape (N, 4).
        nside : int, optional
            HEALPix resolution of the pixelization.

        Returns
        -------
        pix : array_like
            Pixel number for each quaternion.
        pa : array_like
            Position angle in degrees.
        """
        return super().quat2pix(quat, nside, True)

    @qp_settings
    def bore2pix(
        self,
        q_off,
        ctime,
        q_bore,
        q_hwp=None,
        nside=256,
        pol=True,
        return_pa=False,
    ):
        """
        Pixel and polarization timestreams for a detector offset.

        Arguments
        ---------
        q_off : quaternion
            Detector offset quaternion for a single detector,
            calculated using :meth:`det_offset`.
        ctime : array_like
            Unix times in seconds UTC, broadcastable to shape (N,),
            the long dimenions of `q_bore`.
        q_bore : quaternion or array of quaternions
            Nx4 array of quaternions encoding the boresight orientation on the
            sky (as output by :meth:`azel2radec` or equivalent)
        q_hwp : quaternion or array of quaternions, optional
            HWP angle quaternions calculated using :meth:`hwp_quat`.  Must be
            broadcastable to the same shape as `q_bore`.
        nside : int, optional
            HEALpix map dimension.  Default: 256.
        pol : bool, optional
            If `False`, return only the pixel timestream
        return_pa : bool, optional
            If `True`, return pa instead of sin2psi / cos2psi

        Returns
        -------
        pix : array_like
            Detector pixel number
        pa/sin2psi : array_like
            Detector polarization orientation if `return_pa` is `True`, or
            sin(2*pa) if `return_pa` is `False`.
        cos2psi : array_like
            detector polarization orientation cos(2*pa), if `return_pa` is `False`.
        """
        if ctime is None:
            if not self.get_param("mean_aber"):
                raise ValueError("ctime required if mean_aber is False")
            ctime = 0.0
        out = super().bore2pix(q_off, ctime, q_bore, q_hwp, nside, bool(return_pa))
        if return_pa:
            return out
        pix, sin2psi, cos2psi = out
        return (pix, sin2psi, cos2psi) if pol else pix

    # ---- Galactic rotation ----

    @qp_settings
    def rotate_quat(self, quat, coord=("C", "G"), inplace=True):
        """
        Rotate a quaternion between celestial and galactic coordinates.

        Arguments
        ---------
        quat : array_like
            array of quaternions, of shape (n, 4)
        coord : list, optional
            2-element list of input and output coordinates
        inplace : bool, optional
            If True, apply the rotation in-place on the input quaternion.
            Otherwise, return a copy of the input array.  Default: True.

        Returns
        -------
        quat : array_like
            rotated quaternion array
        """
        if tuple(coord) == ("C", "G"):
            to_gal = True
        elif tuple(coord) == ("G", "C"):
            to_gal = False
        else:
            raise ValueError("Unsupported coord: {}".format(coord))

        arr = np.asarray(quat, dtype=float)
        # A bare (4,) is one quaternion without a sample axis and comes
        # back the same way; an (n, 4) keeps its axis, n of 1 included.
        single = arr.ndim == 1
        quat = np.atleast_2d(arr)
        if not inplace:
            quat = np.array(quat)
        quat = np.require(quat, float, ["C", "A", "W"])
        super().rotate_quat(quat, to_gal)
        return quat[0] if single else quat

    @qp_settings
    def radec2gal(self, ra, dec, pa=None, sin2psi=None, cos2psi=None, inplace=True):
        """
        Rotate equatorial coordinates to galactic.

        Arguments
        ---------
        ra : array_like
            Right ascension in degrees, of shape (N,). Rotated in place unless
            `inplace` is False.
        dec : array_like
            Declination in degrees, of shape (N,). Rotated in place unless
            `inplace` is False.
        pa : array_like, optional
            Position angle in degrees. Supply this or the `sin2psi`/`cos2psi`
            pair, not both.
        sin2psi : array_like, optional
            sin(2*pa), if the polarization angle is carried as a pair.
        cos2psi : array_like, optional
            cos(2*pa), paired with `sin2psi`.
        inplace : bool, optional
            If True, the default, rotate the caller's arrays and return them.
            If False, work on copies and leave the inputs untouched.

        Returns
        -------
        ra : array_like
            Rotated right ascension in degrees.
        dec : array_like
            Rotated declination in degrees.
        pa/sin2psi : array_like
            Rotated position angle, or sin(2*pa) if the pair was supplied.
        cos2psi : array_like
            Rotated cos(2*pa), if the pair was supplied.
        """
        return self._rotate_coord(ra, dec, pa, sin2psi, cos2psi, inplace, True)

    @qp_settings
    def gal2radec(self, ra, dec, pa=None, sin2psi=None, cos2psi=None, inplace=True):
        """
        Rotate galactic coordinates to equatorial.

        Arguments
        ---------
        ra : array_like
            Right ascension in degrees, of shape (N,). Rotated in place unless
            `inplace` is False.
        dec : array_like
            Declination in degrees, of shape (N,). Rotated in place unless
            `inplace` is False.
        pa : array_like, optional
            Position angle in degrees. Supply this or the `sin2psi`/`cos2psi`
            pair, not both.
        sin2psi : array_like, optional
            sin(2*pa), if the polarization angle is carried as a pair.
        cos2psi : array_like, optional
            cos(2*pa), paired with `sin2psi`.
        inplace : bool, optional
            If True, the default, rotate the caller's arrays and return them.
            If False, work on copies and leave the inputs untouched.

        Returns
        -------
        ra : array_like
            Rotated right ascension in degrees.
        dec : array_like
            Rotated declination in degrees.
        pa/sin2psi : array_like
            Rotated position angle, or sin(2*pa) if the pair was supplied.
        cos2psi : array_like
            Rotated cos(2*pa), if the pair was supplied.
        """
        return self._rotate_coord(ra, dec, pa, sin2psi, cos2psi, inplace, False)

    def _rotate_coord(self, ra, dec, pa, sin2psi, cos2psi, inplace, to_gal):
        do_pa, bc, all_scalar = _prep_coord_rotation(
            ra, dec, pa, sin2psi, cos2psi, inplace
        )
        if do_pa:
            super().rotate_coord(bc[0], bc[1], bc[2], None, None, to_gal)
        else:
            super().rotate_coord(bc[0], bc[1], None, bc[2], bc[3], to_gal)
        # Scalars in, scalars out. There was nothing to rotate in place
        # in that case: the arrays above were made here.
        if all_scalar:
            return tuple(a[0] for a in bc)
        return tuple(bc)

    @qp_settings
    def rotate_coord(
        self,
        ra,
        dec,
        pa=None,
        sin2psi=None,
        cos2psi=None,
        coord=("C", "G"),
        inplace=True,
    ):
        """
        Rotate sky coordinates between celestial and galactic.

        Arguments
        ---------
        ra : array_like
            Right ascension in degrees, of shape (N,). Rotated in place unless
            `inplace` is False.
        dec : array_like
            Declination in degrees, of shape (N,). Rotated in place unless
            `inplace` is False.
        pa : array_like, optional
            Position angle in degrees. Supply this or the `sin2psi`/`cos2psi`
            pair, not both.
        sin2psi : array_like, optional
            sin(2*pa), if the polarization angle is carried as a pair.
        cos2psi : array_like, optional
            cos(2*pa), paired with `sin2psi`.
        inplace : bool, optional
            If True, the default, rotate the caller's arrays and return them.
            If False, work on copies and leave the inputs untouched.
        coord : tuple of str, optional
            The frames to rotate from and to, as a pair. ('C', 'G') is
            celestial to galactic and ('G', 'C') the reverse; anything else
            raises ValueError.

        Returns
        -------
        ra : array_like
            Rotated right ascension in degrees.
        dec : array_like
            Rotated declination in degrees.
        pa/sin2psi : array_like
            Rotated position angle, or sin(2*pa) if the pair was supplied.
        cos2psi : array_like
            Rotated cos(2*pa), if the pair was supplied.
        """
        if tuple(coord) == ("C", "G"):
            fn = self.radec2gal
        elif tuple(coord) == ("G", "C"):
            fn = self.gal2radec
        else:
            raise ValueError("Unsupported coord: {}".format(coord))
        return fn(ra, dec, pa=pa, sin2psi=sin2psi, cos2psi=cos2psi, inplace=inplace)

    # ---- IERS Bulletin A ----

    def update_bulletin_a(self, start_year=2000):
        """
        Load IERS Bulletin A from astropy's auto-updating table.

        Arguments
        ---------
        start_year : int, optional
            Oldest year for which data should be stored.

        Returns
        -------
        mjd : array_like
            Modified Julian date
        dut1 : array_like
            UT1-UTC time correction
        x : array_like
        y : array_like
            Polar motion (wobble) corrections
        """
        try:
            from astropy.utils.iers import IERS_Auto
        except ImportError:
            from warnings import warn

            warn("astropy is required to update IERS Bulletin A", ImportWarning)
            return None

        iers = IERS_Auto.open()
        mjds = np.asarray(iers["MJD"], dtype=float)
        dut1 = np.asarray(iers["UT1_UTC"], dtype=float)
        x = np.asarray(iers["PM_x"], dtype=float)
        y = np.asarray(iers["PM_y"], dtype=float)

        year = np.asarray(iers["year"], dtype=int) + 1900
        year[np.concatenate([[False], np.ediff1d(year) < 0]).cumsum() > 0] += 100
        keep = year >= start_year

        mjds, dut1, x, y = mjds[keep], dut1[keep], x[keep], y[keep]
        self.set_bulletin_a(
            int(mjds[0]),
            int(mjds[-1]),
            *[np.require(a, float, ["C", "A"]) for a in (dut1, x, y)],
        )
        return mjds, dut1, x, y

    def load_bulletin_a(self, filename, columns=("mjd", "dut1", "x", "y"), **kwargs):
        """
        Load IERS Bulletin A from a text file.

        Arguments
        ---------
        filename : string
            Name of the text file containing IERS Bulletin A parameters.
        columns : list of strings
            list of columns as they appear in the file.
            A KeyError is raise if the list does not contain
            each of ['mjd', 'dut1', 'x', 'y'].

        Any other keyword arguments are passed to the `numpy.loadtxt` function

        Returns
        -------
        mjd : array_like
            Modified Julian date
        dut1 : array_like
            UT1-UTC time correction
        x : array_like
        y : array_like
            Polar motion corrections
        """
        req = ["mjd", "dut1", "x", "y"]
        missing = [c for c in req if c not in columns]
        if missing:
            raise KeyError("Missing columns {}".format(missing))

        data = np.loadtxt(filename, unpack=True, **kwargs)
        cols = {c: np.asarray(d, dtype=float) for c, d in zip(columns, data)}
        mjd, dut1, x, y = (cols[c] for c in req)
        self.set_bulletin_a(
            int(mjd[0]),
            int(mjd[-1]),
            *[np.require(a, float, ["C", "A"]) for a in (dut1, x, y)],
        )
        return mjd, dut1, x, y
