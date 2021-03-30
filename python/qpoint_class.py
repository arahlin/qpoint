from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
try:
    import itertools.izip as zip
except ImportError:
    pass
import numpy as np
from . import _libqpoint as lib
from ._libqpoint import libqp as qp
from ._libqpoint import check_input, check_inputs, check_output

__all__ = ['QPoint', 'check_input', 'check_inputs', 'check_output']

class QPoint(object):

    def __init__(self, **kwargs):
        """
        Initialize a `QPoint` memory instance for keeping track of pointing
        corrections over time.

        Any keyword arguments are passed to
        :meth:`qpoint.qpoint_class.QPoint.set` to update memory.
        """

        # initialize memory
        self._memory = qp.qp_init_memory()

        # collect all parameter functions
        self._funcs = lib.qp_funcs
        self._all_funcs = dict()
        for k in self._funcs:
            self._all_funcs.update(**self._funcs[k])

        # set any requested parameters
        self.set(**kwargs)

    def print_memory(self):
        """
        Print current memory state in C.
        """
        qp.qp_print_memory(self._memory)

    def __del__(self):
        """
        Free memory before deleting the object
        """
        qp.qp_free_memory(self._memory)

    def _set(self, key, val):
        """
        Set a single parameter to the given value.  See
        :meth:`qpoint.qpoint_class.QPoint.set` for a list of parameter names.
        """
        if key not in self._all_funcs:
            raise KeyError('Unknown parameter {}'.format(key))
        val = self._all_funcs[key]['check_set'](val)
        self._all_funcs[key]['set'](self._memory, val)

    def _get(self, key):
        """
        Get the value for a single parameter.  See
        :meth:`qpoint.qpoint_class.QPoint.set` for a list of parameter names.
        """
        if key not in self._all_funcs:
            raise KeyError('Unknown parameter {}'.format(key))
        val = self._all_funcs[key]['get'](self._memory)
        return self._all_funcs[key]['check_get'](val)

    def set(self, **kwargs):
        """
        Set computation options.

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
            Rate at which the refaction correction is updated
            (NB: this correction can also be updated manually -- see `refraction`)
        accuracy : 'low' or 'high'
            If 'low', use a truncated form (2000b) for the NPB correction,
            which is much faster but less accurate. If 'high' (default), use
            the full 2006/2000a form.
        mean_aber : bool
            If True, apply the aberration correction as an average for the
            entire field of view.  This is gives a 1-2 arcsec deviation
            at the edges of the SPIDER field of view.
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

        for k in kwargs:
            try:
                self._set(k, kwargs[k])
            except KeyError:
                continue

    def get(self, *args):
        """
        Returns a dictionary of the requested state parameters.  If no
        parameters are supplied, then all are returned.  If a single parameter
        is supplied, then just that value is returned.  See
        :meth:`qpoint.qpoint_class.QPoint.set` for a list of parameter names.

        Can also select 'options', 'rates', 'weather', or 'params' to return
        all of that subset of parameters.
        """
        state = dict()
        if not len(args):
            for k in self._funcs:
                state[k] = dict()
                for kk in self._funcs[k]:
                    state[k][kk] = self._get(kk)
            return state

        for arg in args:
            if arg in self._funcs:
                state[arg] = dict()
                for k in self._funcs[arg]:
                    state[arg][k] = self._get(k)
            else:
                state[arg] = self._get(arg)
        if len(args) == 1:
            return state[args[0]]
        return state

    def reset_rates(self):
        """
        Reset update counters for each state.  Useful to force an updated
        correction term at the beginning of each chunk.
        """
        qp.qp_reset_rates(self._memory)

    def reset_inv_rates(self):
        """
        Reset update counters for each inverse state.  Useful to force an updated
        correction term at the beginning of each chunk.
        """
        qp.qp_reset_inv_rates(self._memory)

    def refraction(self, *args, **kwargs):
        """refraction(q, **kwargs)
        refraction(delta)

        Update refraction parameters

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

        if len(args) == 1 and len(kwargs) == 0:
            v = args[0]
            if np.isscalar(v):
                self._set('ref_delta', v)
                return v

        if 'delta' in kwargs:
            v = kwargs.get('delta')
            self._set('ref_delta', v)
            return v

        arg_names = ['q'] + list(self._funcs['weather'])
        for idx, a in enumerate(args):
            kwargs[arg_names[idx]] = a

        for w in self._funcs['weather']:
            if w in kwargs:
                self._set(w, kwargs.get(w))

        q = kwargs.get('q', None)
        if q is not None:
            def func(x0, x1, x2, x3):
                q = np.ascontiguousarray([x0, x1, x2, x3])
                return qp.qp_update_ref(self._memory, q)
            fvec = np.vectorize(func, [np.double])
            if q.size // 4 > 1:
                q = q.transpose()
            delta = fvec(q[0], q[1], q[2], q[3])
            if delta.shape == ():
                return delta[()]
            return delta
        return self._get('ref_delta')

    def gmst(self, ctime, **kwargs):
        """
        Return Greenwich mean sidereal time for given ctimes.
        Vectorized.

        Arguments
        ---------
        ctime : array_like
            Unix time in seconds UTC

        Returns
        -------
        gmst : array_like
            Greenwich mean sidereal time of the observer

        Notes
        -----
        Any keywords accepted by the :meth:`qpoint.qpoint_class.QPoint.set`
        method can also be passed here, and will be processed prior to
        calculation.
        """

        self.set(**kwargs)

        ctime = check_input('ctime', np.atleast_1d(ctime))
        n = ctime.size

        gmst = check_output('gmst', shape=ctime.shape)
        qp.qp_gmstn(self._memory, ctime, gmst, n)

        if n == 1:
            return gmst[0]
        return gmst

    def lmst(self, ctime, lon, **kwargs):
        """
        Return local mean sidereal time for given ctimes and longitudes.
        Vectorized, input arguments must be broadcastable to the same shape.

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

        Notes
        -----
        Any keywords accepted by the :meth:`qpoint.qpoint_class.QPoint.set`
        method can also be passed here, and will be processed prior to
        calculation.
        """

        self.set(**kwargs)

        ctime, lon = check_inputs(ctime, lon)
        n = ctime.size

        lmst = check_output('lmst', shape=ctime.shape)
        qp.qp_lmstn(self._memory, ctime, lon, lmst, n)

        if n == 1:
            return lmst[0]
        return lmst

    def dipole(self, ctime, ra, dec, **kwargs):
        """
        Return dipole amplitude in the given equatorial direction.
        Vectorized, input arguments must be broadcastable to the same shape.

        Arguments
        ---------
        ctime : array_like
            Unix time in seconds UTC
        ra : array_like
            Right ascension on the sky, in degrees.
        dec : array_like
            Declination on the sky, in degrees

        Returns
        -------
        dipole : array_like
            Dipole amplitude in K

        Notes
        -----
        Any keywords accepted by the :meth:`qpoint.qpoint_class.QPoint.set`
        method can also be passed here, and will be processed prior to
        calculation.
        """

        self.set(**kwargs)

        ctime, ra, dec = check_inputs(ctime, ra, dec)
        n = ctime.size

        dipole = check_output('dipole', shape=ctime.shape, **kwargs)
        qp.qp_dipolen(self._memory, ctime, ra, dec, dipole, n)

        if n == 1:
            return dipole[0]
        return dipole

    def bore2dipole(self, q_off, ctime, q_bore, **kwargs):
        """
        Calculate dipole timestream for given offset and boresight pointing.

        Arguments
        ---------
        q_off : quaternion
            Detector offset quaternion for a single detector, calculated using
            `det_offset`
        ctime : array_like
            Array of unix times in seconds UTC
        q_bore : quaternion or array of quaternions
            Array of quaternions encoding the boresight orientation
            on the sky (as output by `azel2radec` or similar).
            Broadcastable to the same length as `ctime`.

        Returns
        -------
        dipole : array_like
            Dipole amplitude in K

        Notes
        -----
        Any keywords accepted by the :meth:`qpoint.qpoint_class.QPoint.set`
        method can also be passed here, and will be processed prior to
        calculation.
        """

        self.set(**kwargs)

        q_off = check_input('q_off', q_off, shape=(4,), quat=True)
        ctime = check_input('ctime', ctime)
        n = ctime.size
        q_bore = check_input('q_bore', q_bore, shape=(n, 4), quat=True)

        dipole = check_output('dipole', shape=ctime.shape, **kwargs)
        qp.qp_bore2dipole(self._memory, q_off, ctime, q_bore, dipole, n)

        if n == 1:
            return dipole[0]
        return dipole

    def det_offset(self, delta_az, delta_el, delta_psi):
        """
        Return quaternion corresponding to the requested detector centroid
        offset from boresight.  Vectorized, input arguments must be
        broadcastable to the same shape.

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

        delta_az, delta_el, delta_psi = \
            check_inputs(delta_az, delta_el, delta_psi)
        ndet = delta_az.size

        quat = check_output('quat', shape=(ndet,4))
        qp.qp_det_offsetn(delta_az, delta_el, delta_psi, quat, ndet)

        if ndet == 1:
            return quat[0]
        return quat

    def bore_offset(self, q_bore, ang1=None, ang2=None, ang3=None,
                    post=False, inplace=False):
        """
        Apply a fixed or variable offset to the boresight quaternion.
        Input arguments must be broadcastable to the same shape.

        Arguments
        ---------
        q_bore : array_like
            boresight pointing quaternion
        ang1 : array_like, optional
            Azimuthal or ra offset in degrees
        ang2 : array_like, optional
            Elevation or dec offset in degrees
        ang3 : array_like, optional
            Position angle offset in degrees
        post : bool, optional
            If False, apply offset as an az/el/pa pre-rotation
            If True, apply offset as an ra/dec/pa post-rotation
        inplace : bool, optional
            If True, apply the rotation in-place in memory.

        Returns
        -------
        q_bore : array_like
            Offset boresight quaternion
        """
        q_bore = check_input('q_bore', np.atleast_2d(q_bore), quat=True,
                             inplace=inplace, output=True)
        n = q_bore.shape[0]
        if all([a is None for a in [ang1, ang2, ang3]]):
            raise ValueError('One of ang1, ang2, ang3 is required')
        ang1, ang2, ang3 = \
            check_inputs(ang1, ang2, ang3, shape=(n,))
        qp.qp_bore_offset(self._memory, q_bore, ang1, ang2, ang3, n, int(post))
        return q_bore

    def hwp_quat(self, theta):
        """
        Return quaternion corresponding to rotation by 2 * theta, where theta is
        the physical waveplate angle. Vectorized.

        Arguments
        ---------
        theta : array_like
            HWP physical angle in degrees

        Returns
        -------
        q : array_like
            Quaternion for each hwp angle
        """
        theta = check_input('theta', np.atleast_1d(theta))
        n = theta.size

        quat = check_output('quat', shape=(n,4))
        qp.qp_hwp_quatn(theta, quat, n)

        if n == 1:
            return quat[0]
        return quat

    def azel2bore(self, az, el, pitch, roll, lon, lat, ctime, q=None,
                  **kwargs):
        """
        Estimate the quaternion for the boresight orientation on the sky given
        the attitude (az/el/pitch/roll), location on the earth (lon/lat) and
        ctime. Input vectors must be numpy-array-like and broadcastable to the
        same shape.

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
        q : array_like, optional
            Output quaternion array initialized by user.  Supply this
            for in-place computation.

        Returns
        -------
        q : array_like
            Nx4 numpy array of quaternions for each supplied timestamp.

        Notes
        -----
        Any keywords accepted by the :meth:`qpoint.qpoint_class.QPoint.set`
        method can also be passed here, and will be processed prior to
        calculation.
        """

        self.set(**kwargs)

        az, el, pitch, roll, lon, lat, ctime = \
            check_inputs(az, el, pitch, roll, lon, lat, ctime)
        n = az.size

        # identity quaternion
        q = check_output('q', q, shape=(n,4), fill=[1,0,0,0])

        qp.qp_azel2bore(self._memory, az, el, pitch, roll, lon, lat,
                        ctime, q, n)

        return q

    def azelpsi2bore(self, az, el, psi, pitch, roll, lon, lat, ctime, q=None,
                  **kwargs):
        """
        Estimate the quaternion for the boresight orientation on the sky given
        the attitude (az/el/psi/pitch/roll), location on the earth (lon/lat) and
        ctime. Input vectors must be numpy-array-like and broadcastable to the
        same shape.

        This is an augmented version of azel2bore to accept rotations of the focal
        plane about the optical axis.

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
        q : array_like, optional
            Output quaternion array initialized by user.  Supply this
            for in-place computation.

        Returns
        -------
        q : array_like
            Nx4 numpy array of quaternions for each supplied timestamp.

        Notes
        -----
        Any keywords accepted by the :meth:`qpoint.qpoint_class.QPoint.set`
        method can also be passed here, and will be processed prior to
        calculation.
        """

        self.set(**kwargs)

        az, el, psi, pitch, roll, lon, lat, ctime = \
            check_inputs(az, el, psi, pitch, roll, lon, lat, ctime)
        n = az.size

        # identity quaternion
        q = check_output('q', q, shape=(n,4), fill=[1,0,0,0])

        qp.qp_azelpsi2bore(self._memory, az, el, psi, pitch, roll, lon, lat,
                        ctime, q, n)

        return q

    def bore2radec(self, q_off, ctime, q_bore, q_hwp=None, sindec=False,
                   return_pa=False, ra=None, dec=None, pa=None,
                   sin2psi=None, cos2psi=None, **kwargs):
        """
        Calculate the orientation on the sky for a detector offset from the
        boresight.  Detector offsets are defined assuming the boresight is
        pointed toward the horizon, and that the boresight polarization axis is
        along the vertical.

        Arguments
        ---------
        q_off : quaternion
            Detector offset quaternion for a single detector, calculated using
            `det_offset`.
        ctime : array_like
            Unix time in seconds UTC, broadcastable to shape (N,),
            the long dimension of `q_bore`.
        q_bore : quaternion or array of quaternions
            Nx4 array of quaternions encoding the boresight orientation
            on the sky (as output by `azel2radec` or equivalent)
        q_hwp : quaternion or array of quaternions, optional
            HWP angle quaternions calculated using `hwp_quat`.
            Must be broadcastable to the same shape as `q_bore`.
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

        Notes
        -----
        Any keywords accepted by the :meth:`qpoint.qpoint_class.QPoint.set`
        method can also be passed here, and will be processed prior to
        calculation.

        Pre-allocated output arguments can also be supplied as input keywords
        for in-place operation.
        """

        self.set(**kwargs)

        q_off  = check_input('q_off', q_off, quat=True)
        q_bore = check_input('q_bore', np.atleast_2d(q_bore), quat=True)
        shape = (q_bore.size // 4,)
        if ctime is None:
            if not self.get('mean_aber'):
                raise ValueError('ctime required if mean_aber is False')
            ctime = np.zeros(shape, dtype=q_bore.dtype)
        ctime  = check_input('ctime', ctime, shape=shape)
        pars = dict(shape=ctime.shape, dtype=np.double)
        ra = check_output('ra', ra, **pars)
        dec = check_output('dec', dec, **pars)
        if return_pa:
            if sindec:
                raise ValueError('Cannot use sindec with return_pa=True')
            pa = check_output('pa', pa, **pars)
        else:
            sin2psi = check_output('sin2psi', sin2psi, **pars)
            cos2psi = check_output('cos2psi', cos2psi, **pars)
        n = ctime.size

        if q_hwp is None:
            if return_pa:
                qp.qp_bore2radecpa(self._memory, q_off, ctime, q_bore,
                                   ra, dec, pa, n)
            elif sindec:
                qp.qp_bore2rasindec(self._memory, q_off, ctime, q_bore,
                                    ra, dec, sin2psi, cos2psi, n)
            else:
                qp.qp_bore2radec(self._memory, q_off, ctime, q_bore,
                                 ra, dec, sin2psi, cos2psi, n)
        else:
            q_hwp = check_input('q_hwp', q_hwp, shape=q_bore.shape)
            if return_pa:
                qp.qp_bore2radecpa_hwp(self._memory, q_off, ctime, q_bore,
                                       q_hwp, ra, dec, pa, n)
            elif sindec:
                qp.qp_bore2rasindec_hwp(self._memory, q_off, ctime, q_bore,
                                        q_hwp, ra, dec, sin2psi, cos2psi, n)
            else:
                qp.qp_bore2radec_hwp(self._memory, q_off, ctime, q_bore,
                                     q_hwp, ra, dec, sin2psi, cos2psi, n)

        if return_pa:
            if n == 1:
                return ra[0], dec[0], pa[0]
            return ra, dec, pa

        if n == 1:
            return ra[0], dec[0], sin2psi[0], cos2psi[0]
        return ra, dec, sin2psi, cos2psi

    def azel2radec(self, delta_az, delta_el, delta_psi,
                   az, el, pitch, roll, lon, lat, ctime,
                   hwp=None, sindec=False, return_pa=False,
                   ra=None, dec=None, pa=None, sin2psi=None,
                   cos2psi=None, **kwargs):
        """
        Estimate the orientation on the sky for a detector offset from
        boresight, given the boresight attitude (az/el/pitch/roll), location on
        the earth (lon/lat) and UTC time.  Input vectors must be
        numpy-array-like and broadcastable to the same shape. Detector offsets
        are defined assuming the boresight is pointed toward the horizon, and
        that the boresight polarization axis is along the horizontal.

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

        Notes
        -----
        Any keywords accepted by the :meth:`qpoint.qpoint_class.QPoint.set`
        method can also be passed here, and will be processed prior to
        calculation.
        """

        self.set(**kwargs)

        az, el, pitch, roll, lon, lat, ctime = \
            check_inputs(az, el, pitch, roll, lon, lat, ctime)

        ra = check_output('ra', ra, shape=az.shape, dtype=np.double)
        dec = check_output('dec', dec, shape=az.shape, dtype=np.double)
        if return_pa:
            pa = check_output('pa', pa, shape=az.shape, dtype=np.double)
        else:
            sin2psi = check_output('sin2psi', sin2psi, shape=az.shape,
                                   dtype=np.double)
            cos2psi = check_output('cos2psi', cos2psi, shape=az.shape,
                                   dtype=np.double)
        n = az.size

        if hwp is None:
            if return_pa:
                qp.qp_azel2radecpa(self._memory, delta_az, delta_el, delta_psi,
                                   az, el, pitch, roll, lon, lat, ctime,
                                   ra, dec, pa, n)
            elif sindec:
                qp.qp_azel2rasindec(self._memory, delta_az, delta_el, delta_psi,
                                    az, el, pitch, roll, lon, lat, ctime,
                                    ra, dec, sin2psi, cos2psi, n)
            else:
                qp.qp_azel2radec(self._memory, delta_az, delta_el, delta_psi,
                                 az, el, pitch, roll, lon, lat, ctime,
                                 ra, dec, sin2psi, cos2psi, n)
        else:
            hwp = check_input('hwp', hwp, shape=az.shape)

            if return_pa:
                qp.qp_azel2radec_hwp(self._memory, delta_az, delta_el, delta_psi,
                                     az, el, pitch, roll, lon, lat, ctime, hwp,
                                     ra, dec, pa, n)
            elif sindec:
                qp.qp_azel2rasindec_hwp(self._memory, delta_az, delta_el, delta_psi,
                                        az, el, pitch, roll, lon, lat, ctime, hwp,
                                        ra, dec, sin2psi, cos2psi, n)
            else:
                qp.qp_azel2radec_hwp(self._memory, delta_az, delta_el, delta_psi,
                                     az, el, pitch, roll, lon, lat, ctime, hwp,
                                     ra, dec, sin2psi, cos2psi, n)

        if return_pa:
            return ra, dec, pa
        return ra, dec, sin2psi, cos2psi

    def azelpsi2radec(self, delta_az, delta_el, delta_psi,
                   az, el, psi, pitch, roll, lon, lat, ctime,
                   hwp=None, sindec=False, return_pa=False,
                   ra=None, dec=None, pa=None, sin2psi=None,
                   cos2psi=None, **kwargs):
        """
        Estimate the orientation on the sky for a detector offset from
        boresight, given the boresight attitude (az/el/pitch/roll), location on
        the earth (lon/lat) and UTC time.  Input vectors must be
        numpy-array-like and broadcastable to the same shape. Detector offsets
        are defined assuming the boresight is pointed toward the horizon, and
        that the boresight polarization axis is along the horizontal.

        This is agumented from azel2radec(), to accept rotations of the focal
        plane about the optical axis.

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

        Notes
        -----
        Any keywords accepted by the :meth:`qpoint.qpoint_class.QPoint.set`
        method can also be passed here, and will be processed prior to
        calculation.
        """

        self.set(**kwargs)

        az, el, psi, pitch, roll, lon, lat, ctime = \
            check_inputs(az, el, psi, pitch, roll, lon, lat, ctime)

        ra = check_output('ra', ra, shape=az.shape, dtype=np.double)
        dec = check_output('dec', dec, shape=az.shape, dtype=np.double)
        if return_pa:
            pa = check_output('pa', pa, shape=az.shape, dtype=np.double)
        else:
            sin2psi = check_output('sin2psi', sin2psi, shape=az.shape,
                                   dtype=np.double)
            cos2psi = check_output('cos2psi', cos2psi, shape=az.shape,
                                   dtype=np.double)
        n = az.size

        if hwp is None:
            if return_pa:
                qp.qp_azelpsi2radecpa(self._memory, delta_az, delta_el, delta_psi,
                                   az, el, psi, pitch, roll, lon, lat, ctime,
                                   ra, dec, pa, n)
            elif sindec:
                qp.qp_azelpsi2rasindec(self._memory, delta_az, delta_el, delta_psi,
                                    az, el, psi, pitch, roll, lon, lat, ctime,
                                    ra, dec, sin2psi, cos2psi, n)
            else:
                qp.qp_azelpsi2radec(self._memory, delta_az, delta_el, delta_psi,
                                 az, el, psi, pitch, roll, lon, lat, ctime,
                                 ra, dec, sin2psi, cos2psi, n)
        else:
            hwp = check_input('hwp', hwp, shape=az.shape)

            if return_pa:
                qp.qp_azelpsi2radec_hwp(self._memory, delta_az, delta_el, delta_psi,
                                     az, el, psi, pitch, roll, lon, lat, ctime, hwp,
                                     ra, dec, pa, n)
            elif sindec:
                qp.qp_azelpsi2rasindec_hwp(self._memory, delta_az, delta_el, delta_psi,
                                        az, el, psi, pitch, roll, lon, lat, ctime, hwp,
                                        ra, dec, sin2psi, cos2psi, n)
            else:
                qp.qp_azelpsi2radec_hwp(self._memory, delta_az, delta_el, delta_psi,
                                     az, el, psi, pitch, roll, lon, lat, ctime, hwp,
                                     ra, dec, sin2psi, cos2psi, n)

        if return_pa:
            return ra, dec, pa
        return ra, dec, sin2psi, cos2psi

    def radec2azel(self, ra, dec, pa, lon, lat, ctime, az=None, el=None, hpa=None,
                   **kwargs):
        """
        Estimate the horizon coordinates for a given set of equatorial coordinates
        (ra/dec/psi), location on the earth (lon/lat) and UTC time.  Input vectors
        must be numpy-array-like and broadcastable to the same shape.

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

        Notes
        -----
        Any keywords accepted by the :meth:`qpoint.qpoint_class.QPoint.set`
        method can also be passed here, and will be processed prior to
        calculation.
        """

        self.set(**kwargs)

        ra, dec, pa, lon, lat, ctime = \
            check_inputs(ra, dec, pa, lon, lat, ctime)

        az = check_output('az', az, shape=ra.shape, dtype=np.double)
        el = check_output('el', el, shape=ra.shape, dtype=np.double)
        hpa = check_output('hpa', hpa, shape=ra.shape, dtype=np.double)
        n = ra.size

        qp.qp_radec2azel(self._memory, ra, dec, pa, lon, lat, ctime,
                         az, el, hpa, n)

        return az, el, hpa

    def radecpa2quat(self, ra, dec, pa, **kwargs):
        """
        Calculate quaternion for input ra/dec/pa. Vectorized, input arguments
        must be broadcastable to the same shape.

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

        Notes
        -----
        Any keywords accepted by the :meth:`qpoint.qpoint_class.QPoint.set`
        method can also be passed here, and will be processed prior to
        calculation.
        """
        self.set(**kwargs)

        ra, dec, pa = check_inputs(ra, dec, pa)
        n = ra.size
        quat = check_output('quat', shape=(n,4), dtype=np.double, **kwargs)
        qp.qp_radecpa2quatn(self._memory, ra, dec, pa, quat, n)

        if n == 1:
            return quat[0]
        return quat

    def quat2radecpa(self, quat, **kwargs):
        """
        Calculate ra/dec/pa for input quaternion(s).

        Arguments
        ---------
        q : array_like
            Pointing quaternion

        Returns
        -------
        ra : array_like
            Right ascension angle
        dec : array_like
            Declination angle
        pa : array_like
            Position angle

        Notes
        -----
        Any keywords accepted by the :meth:`qpoint.qpoint_class.QPoint.set`
        method can also be passed here, and will be processed prior to
        calculation.
        """
        self.set(**kwargs)

        quat = check_input('quat', np.atleast_2d(quat), quat=True)
        n = quat.shape[0]
        ra = check_output('ra', shape=(n,), dtype=np.double, **kwargs)
        dec = check_output('dec', shape=(n,), dtype=np.double, **kwargs)
        pa = check_output('pa', shape=(n,), dtype=np.double, **kwargs)

        qp.qp_quat2radecpan(self._memory, quat, ra, dec, pa, n)
        if n == 1:
            return ra[0], dec[0], pa[0]
        return ra, dec, pa

    def quat2pixpa(self, quat, nside=256, **kwargs):
        """
        Calculate ra/dec/pa for input quaternion(s).

        Arguments
        ---------
        q : array_like
            Pointing quaternion

        Returns
        -------
        pix : array_like
            Pixel number(s) corresponding to the input positions(s).
        pa : array_like
            Position angle

        Notes
        -----
        Any keywords accepted by the :meth:`qpoint.qpoint_class.QPoint.set`
        method can also be passed here, and will be processed prior to
        calculation.
        """
        self.set(**kwargs)

        quat = check_input('quat', np.atleast_2d(quat), quat=True)
        n = quat.shape[0]
        pix = check_output('pix', shape=(n,), dtype=np.int, **kwargs)
        pa = check_output('pa', shape=(n,), dtype=np.double, **kwargs)

        qp.qp_quat2pixpan(self._memory, quat, pix, pa, n)
        if n == 1:
            return pix[0], pa[0]
        return pix, pa

    def radec2pix(self, ra, dec, nside=256, **kwargs):
        """
        Calculate HEALpix pixel number for given ra/dec and nside.  Vectorized,
        input arguments must be broadcastable to the same shape.

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

        Notes
        -----
        Any keywords accepted by the :meth:`qpoint.qpoint_class.QPoint.set`
        method can also be passed here, and will be processed prior to
        calculation.
        """

        self.set(**kwargs)

        ra, dec = check_inputs(ra, dec)
        n = ra.size

        pix = check_output('pix', shape=ra.shape, dtype=np.int, **kwargs)
        qp.qp_radec2pixn(self._memory, ra, dec, nside, pix, n)

        if n == 1:
            return pix[0]
        return pix

    def rotate_quat(self, quat, coord=['C', 'G'], inplace=True, **kwargs):
        """
        Rotate a quaternion from one coordinate system to another.
        Supported coordinates:

        C: celestial (equatorial) coordinates
        G: galactic coordinates

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

        Notes
        -----
        Any keywords accepted by the :meth:`qpoint.qpoint_class.QPoint.set`
        method can also be passed here, and will be processed prior to
        calculation.
        """

        self.set(**kwargs)

        quat = check_input('quat', np.atleast_2d(quat), quat=True,
                           inplace=inplace, output=True)
        n = quat.size // 4

        if coord[0] == 'C' and coord[1] == 'G':
            qp.qp_radec2gal_quatn(self._memory, quat, n)
        elif coord[0] == 'G' and coord[1] == 'C':
            qp.qp_gal2radec_quatn(self._memory, quat, n)
        else:
            raise ValueError('unsupported coord {}'.format(coord))

        return quat.squeeze()

    def rotate_coord(self, ra, dec, pa=None, sin2psi=None, cos2psi=None,
                     coord=['C', 'G'], inplace=True, **kwargs):
        """
        Rotate coordinates from one coordinate system to another. Vectorized,
        input arguments must be broadcastable to the same shape.
        Supported coordinates:

        C: celestial (equatorial) coordinates
        G: galactic coordinates

        Arguments
        ---------
        ra, dec, pa : array_like
          -- or --
        ra, dec, sin2psi, cos2psi : array_like
            arrays of coordinates, of shape (n,)
            If none of pa, sin2psi or cos2psi are supplied,
            pa = 0 is assumed.
        coord : list, optional
            2-element list of input and output coordinates.
            Supported systems: C, G.
        inplace : bool, optional
            If True, apply the rotation in-place on the input coordinates.
            Otherwise, return a copy of the input array.  Default: True.

        Returns
        -------
        ra, dec, pa : array_like
          -- or --
        ra, dec, sin2psi, cos2psi : array_like
            rotated coordinate arrays, same form as input

        Notes
        -----
        Any keywords accepted by the :meth:`qpoint.qpoint_class.QPoint.set`
        method can also be passed here, and will be processed prior to
        calculation.
        """
        self.set(**kwargs)

        do_pa = True
        if pa is None:
            if sin2psi is None and cos2psi is None:
                pass
            elif sin2psi is None or cos2psi is None:
                raise KeyError('both sin2psi and cos2psi arguments required')
            else:
                do_pa = False
        else:
            if sin2psi is not None or cos2psi is not None:
                raise KeyError('ambiguous pol arguments, supply either pa '
                               'or sin2psi/cos2psi only')

        ra, dec, pa, sin2psi, cos2psi = \
            check_inputs(ra, dec, pa, sin2psi, cos2psi, inplace=inplace, output=True)
        n = ra.size

        if coord[0] == 'C' and coord[1] == 'G':
            if do_pa:
                qp.qp_radecpa2galn(self._memory, ra, dec, pa, n)
            else:
                qp.qp_radec2galn(self._memory, ra, dec, sin2psi, cos2psi, n)
        elif coord[0] == 'G' and coord[1] == 'C':
            if do_pa:
                qp.qp_gal2radecpan(self._memory, ra, dec, pa, n)
            else:
                qp.qp_gal2radecn(self._memory, ra, dec, sin2psi, cos2psi, n)

        if n == 1:
            if do_pa:
                return ra[0], dec[0], pa[0]
            return ra[0], dec[0], sin2psi[0], cos2psi[0]
        if do_pa:
            return ra, dec, pa
        return ra, dec, sin2psi, cos2psi

    def radec2gal(self, ra, dec, pa=None, sin2psi=None, cos2psi=None,
                  inplace=True, **kwargs):
        """
        Rotate celestial (equatorial) coordinates to galactic coordinates.
        This is equivalent to calling `rotate_coord(..., coord=['C', 'G'])`.
        Vectorized, input arguments must be broadcastable to the same shape.

        Arguments
        ---------
        ra, dec, pa : array_like
          -- or --
        ra, dec, sin2psi, cos2psi : array_like
            arrays of coordinates, of shape (n,)
            If none of pa, sin2psi or cos2psi are supplied,
            pa = 0 is assumed.
        inplace : bool, optional
            If True, apply the rotation in-place on the input quaternion.
            Otherwise, return a copy of the input array.  Default: True.

        Returns
        -------
        l, b, pa : array_like
          -- or --
        l, b, sin2psi, cos2psi : array_like
            rotated coordinate arrays, same form as input

        Notes
        -----
        Any keywords accepted by the :meth:`qpoint.qpoint_class.QPoint.set`
        method can also be passed here, and will be processed prior to
        calculation.
        """
        return self.rotate_coord(ra, dec, pa, sin2psi, cos2psi,
                                 coord=['C','G'], inplace=inplace, **kwargs)

    def gal2radec(self, ra, dec, pa=None, sin2psi=None, cos2psi=None,
                  inplace=True, **kwargs):
        """
        Rotate galactic coordinates to celestial (equatorial) coordinates.
        This is equivalent to calling `rotate_coord(..., coord=['G', 'C'])`.
        Vectorized, input arguments must be broadcastable to the same shape.

        Arguments
        ---------
        l, b, pa : array_like
          -- or --
        l, b, sin2psi, cos2psi : array_like
            arrays of coordinates, of shape (n,)
            If none of pa, sin2psi or cos2psi are supplied,
            pa = 0 is assumed.
        inplace : bool, optional
            If True, apply the rotation in-place on the input quaternion.
            Otherwise, return a copy of the input array.  Default: True.

        Returns
        -------
        ra, dec, pa : array_like
          -- or --
        ra, dec, sin2psi, cos2psi : array_like
            rotated coordinate arrays, same form as input

        Notes
        -----
        Any keywords accepted by the :meth:`qpoint.qpoint_class.QPoint.set`
        method can also be passed here, and will be processed prior to
        calculation.
        """
        return self.rotate_coord(ra, dec, pa, sin2psi, cos2psi,
                                 coord=['G','C'], inplace=inplace, **kwargs)

    def rotate_map(self, map_in, coord=['C','G'], map_out=None,
                   interp_pix=True, **kwargs):
        """
        Rotate a polarized 3-x-npix map from one coordinate system to another.
        Supported coordinates:

        C = celestial (equatorial J2000)
        G = galactic

        Arguments
        ---------
        map_in : array_like
            Input map, of shape (3, N)
        coord : list, optional
            2-element list of input and output coordinates.
            Supported systems: C, G.
        map_out : array_like, optional
            Rotated output map, for inplace operation.  Same shape as `map_in`.
        interp_pix : bool, optional
            If True, interpolate the rotated map.

        Returns
        -------
        map_out : array_like
            Rotated output map.

        Notes
        -----
        Any keywords accepted by the :meth:`qpoint.qpoint_class.QPoint.set`
        method can also be passed here, and will be processed prior to
        calculation.

        Only full-sky maps are currently supported.
        """

        from warnings import warn
        warn('This code is buggy, use at your own risk', UserWarning)

        interp_orig = self.get('interp_pix')
        self.set(interp_pix=interp_pix, **kwargs)

        from .qmap_class import check_map
        map_in, nside = check_map(map_in)
        map_out = check_output(
            'map_out', map_out, shape=map_in.shape, dtype=map_in.dtype, fill=0)

        try:
            coord_in = coord[0]
            coord_out = coord[1]
        except:
            raise ValueError('unable to parse coord')

        map_in_p = lib.pointer_2d(map_in)
        map_out_p = lib.pointer_2d(map_out)

        qp.qp_rotate_map(self._memory, nside, map_in_p, coord_in,
                         map_out_p, coord_out)

        self.set(interp_pix=interp_orig)
        return map_out

    def quat2pix(self, quat, nside=256, pol=True, **kwargs):
        """
        Calculate HEALpix pixel number and optional polarization angle
        for a given orientation.

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

        Notes
        -----
        Any keywords accepted by the :meth:`qpoint.qpoint_class.QPoint.set`
        method can also be passed here, and will be processed prior to
        calculation.
        """

        self.set(**kwargs)

        quat = check_input('quat', np.atleast_2d(quat), quat=True)

        n = quat.shape[0]
        shape = (n,)
        pix = check_output('pix', shape=shape, dtype=np.int, **kwargs)
        sin2psi = check_output('sin2psi', shape=shape, **kwargs)
        cos2psi = check_output('cos2psi', shape=shape, **kwargs)
        qp.qp_quat2pixn(self._memory, quat, nside, pix, sin2psi, cos2psi, n)

        if n == 1:
            pix, sin2psi, cos2psi = pix[0], sin2psi[0], cos2psi[0]
        if pol:
            return pix, sin2psi, cos2psi
        return pix

    def bore2pix(self, q_off, ctime, q_bore, q_hwp=None, nside=256, pol=True,
                 return_pa=False, **kwargs):
        """
        Calculate the orientation on the sky for a detector offset from the
        boresight.  Detector offsets are defined assuming the boresight is
        pointed toward the horizon, and that the boresight polarization axis is
        along the vertical.

        Arguments
        ---------
        q_off : quaternion
            Detector offset quaternion for a single detector,
            calculated using `det_offset`.
        ctime : array_like
            Unix times in seconds UTC, broadcastable to shape (N,),
            the long dimenions of `q_bore`.
        q_bore : quaternion or array of quaternions
            Nx4 array of quaternions encoding the boresight orientation
            on the sky (as output by `azel2radec` or equivalent)
        q_hwp : quaternion or array of quaternions, optional
            HWP angle quaternions calculated using `hwp_quat`.
            Must be broadcastable to the same shape as `q_bore`.
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

        Notes
        -----
        Any keywords accepted by the :meth:`qpoint.qpoint_class.QPoint.set`
        method can also be passed here, and will be processed prior to
        calculation.
        """

        self.set(**kwargs)

        q_off  = check_input('q_off', q_off, quat=True)
        q_bore = check_input('q_bore', q_bore, quat=True)
        if ctime is None:
            if not self.get('mean_aber'):
                raise ValueError('ctime required if mean_aber is False')
            ctime = np.zeros((q_bore.size // 4,), dtype=q_bore.dtype)
        ctime  = check_input('ctime', ctime)
        pix  = check_output('pix', shape=ctime.shape,
                                dtype=np.int, **kwargs)
        if return_pa:
            pa = check_output('pa', shape=ctime.shape, **kwargs)
        else:
            sin2psi = check_output('sin2psi', shape=ctime.shape,
                                   **kwargs)
            cos2psi = check_output('cos2psi', shape=ctime.shape,
                                   **kwargs)
        n = ctime.size

        if q_hwp is None:
            if return_pa:
                qp.qp_bore2pixpa(self._memory, q_off, ctime, q_bore,
                                 nside, pix, pa, n)
            else:
                qp.qp_bore2pix(self._memory, q_off, ctime, q_bore,
                               nside, pix, sin2psi, cos2psi, n)
        else:
            q_hwp = check_input('q_hwp', q_hwp, shape=q_bore.shape)

            if return_pa:
                qp.qp_bore2pixpa_hwp(self._memory, q_off, ctime, q_bore,
                                     q_hwp, nside, pix, pa, n)
            else:
                qp.qp_bore2pix_hwp(self._memory, q_off, ctime, q_bore,
                                   q_hwp, nside, pix, sin2psi, cos2psi, n)

        if pol is True:
            if return_pa:
                return pix, pa
            return pix, sin2psi, cos2psi
        return pix

    def get_interp_val(self, map_in, ra, dec, nest=False):
        """
        Interpolate map pixels to these coordinates.  Uses a C implementation
        of the bilinear interpolation method `get_interpol()` as implemented
        in the equivalent healpix_cxx / healpy function.

        Arguments
        ---------
        map_in : array_like
            A single healpix map or list of maps which to interpolate from.
        ra, dec: array_like
            Timestreams of coordinates to interpolate to, in degrees,
            broadcastable to a common shape (nsample,)
        nest: bool, optional
            If True, input map is in the nested pixel ordering scheme.
            Otherwise, ring ordering is assumed.
            Default: False.

        Returns
        -------
        values : array_like
            Array of interpolated map values, of shape (nmap, nsample).
        """

        pix_order = self.get('pix_order')
        if nest:
            self.set(pix_order='nest')
        else:
            self.set(pix_order='ring')

        ra, dec = check_inputs(ra, dec)
        n = ra.size

        from .qmap_class import check_map
        map_in, nside = check_map(map_in)

        val = check_output('value', shape=(len(map_in), n))

        for m, v in zip(map_in, val):
            qp.qp_get_interp_valn(self._memory, nside, m, ra, dec, v, n)

        self.set(pix_order=pix_order)

        v = v.squeeze()
        if not v.shape:
            return v[()]
        return v

    def update_bulletin_a(self, start_year=2000):
        """
        Update the IERS Bulletin A database using astropy, and return the stored
        entries.  Issues an ImportWarning if astropy version 1.2 or newer is not
        found.

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
            warn('Compatible Astropy not found.  Install astropy v1.2 or newer '
                 'for accurate polar motion and UT1 corrections', ImportWarning)
            return

        columns = ['MJD', 'UT1_UTC', 'PM_x', 'PM_y', 'year']
        iers_table = IERS_Auto.open()[columns].as_array()

        # check year
        year = iers_table['year'] + 1900
        wraps, = np.where(np.ediff1d(year) < 0)
        for idx in wraps:
            year[idx + 1:] += 100
        iers_table['year'] = year
        iers_table = iers_table[year >= start_year]

        # check MJD
        mjds = iers_table['MJD']
        mjd_min = int(mjds.min())
        mjd_max = int(mjds.max())

        # update table
        dut1 = np.array(iers_table['UT1_UTC'])
        x = np.array(iers_table['PM_x'])
        y = np.array(iers_table['PM_y'])
        qp.qp_set_iers_bulletin_a(self._memory, mjd_min, mjd_max, dut1, x, y)

        return mjds, dut1, x, y

    def load_bulletin_a(self, filename, columns=['mjd','dut1','x','y'], **kwargs):
        """
        Load IERS Bulletin A from file and store in memory.  The file must be
        readable using `numpy.loadtxt` with `unpack=True`, and is assumed to be
        sorted by mjd.

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

        req_columns = ['mjd','dut1','x','y']
        if not set(req_columns) <= set(columns):
            raise KeyError(
                'Missing columns {}'.format(list(set(req_columns)-set(columns))))
        kwargs['unpack'] = True
        data = np.loadtxt(filename, **kwargs)
        mjd, x, y, dut1 = (data[columns.index(x)] for x in req_columns)
        mjd_min, mjd_max = int(mjd[0]), int(mjd[-1])

        try:
            qp.qp_set_iers_bulletin_a(self._memory, mjd_min, mjd_max, dut1, x, y)
        except:
            raise RuntimeError(
                'Error loading Bulletin A data from file {}'.format(filename))

        return mjd, dut1, x, y

    def get_bulletin_a(self, mjd):
        """
        Return dut1/x/y for given mjd. Numpy-vectorized.
        """

        def func(x):
            return lib.qp_get_bulletin_a(self._memory, x)
        fvec = np.vectorize(func, [np.double]*3)

        out = fvec(mjd)
        if out[0].shape == ():
            return tuple(x[()] for x in out)
        return out
