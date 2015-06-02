import numpy as np
import _libqpoint as lib
from _libqpoint import libqp as qp
from _libqpoint import check_input, check_inputs, check_output

__all__ = ['QPoint', 'check_input', 'check_inputs', 'check_output']

class QPoint(object):

    def __init__(self, **kwargs):
        """
        Initialize a QPoint memory instance for keeping track of pointing
        corrections over time.

        Any keyword arguments are passed to QPoint.set to update memory.
        """

        # initialize memory
        self._memory = qp.qp_init_memory()

        # collect all parameter functions
        self._funcs = lib.qp_funcs
        self._all_funcs = dict()
        for k,v in self._funcs.items():
            self._all_funcs.update(**v)

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
        Set a single parameter to the given value.  See QPoint.set for a list
        of parameter names.
        """
        if key not in self._all_funcs:
            raise KeyError,'Unknown parameter {}'.format(key)
        val = self._all_funcs[key]['check_set'](val)
        self._all_funcs[key]['set'](self._memory, val)

    def _get(self, key):
        """
        Get the value for a single parameter.  See QPoint.set for a list of
        parameter names.
        """
        if key not in self._all_funcs:
            raise KeyError,'Unknown parameter {}'.format(key)
        val = self._all_funcs[key]['get'](self._memory)
        return self._all_funcs[key]['check_get'](val)

    def set(self, **kwargs):
        """
        Available keywords are:

        * Correction rates:
          NB: these can be 'never' (-999), 'once' (-1), 'always' (0) or seconds
        rate_daber     Rate at which the diurnal aberration correction is
                       applied (NB: this can only be applied always or never)
        rate_lonlat    Rate at which observer's lon and lat are updated
        rate_wobble    Rate at which the polar motion correction is updated
                       (NB: these are not estimated for dates beyond a year
                       from now)
        rate_dut1      Rate at which the ut1-utc correction is updated
                       (NB: this is not estimated for dates beyond a year from
                       now)
        rate_erot      Rate at which the earth's rotation angle is updated
        rate_npb       Rate at which the nutation/precession/frame-bias terms
                       are updated
        rate_aaber     Rate at which the annual aberration correction
                       (due to the earth's orbital velocity) is updated
        rate_ref       Rate at which the refaction correction is updated
                       (NB: this correction can also be updated manually -- see
                       refraction)

        * Options:
        accuracy       If 'low', use a truncated form (2000b) for the NPB
                       correction, which is much faster but less accurate.
                       If 'high' (default), use the full 2006/2000a form.
        mean_aber      If True, apply the aberration correction as an average
                       for the entire field of view.  This is gives a 1-2
                       arcsec deviation at the edges of the SPIDER field of
                       view.
        fast_math      If True, use polynomial approximations for trig
                       functions
        polconv        Specify the 'cosmo' or 'iau' polarization convention
        pix_order      'nest' or 'ring' for healpix pixel ordering
        interp_pix     If True, interpolate between pixels in scanning the
                       source map.
        fast_pix       If True, use vec2pix to get pixel number directly from
                       the quaternion instead of ang2pix from ra/dec.
        num_threads    Number of threads for openMP bore2map computation

        * Weather:
        height         height above sea level, meters
        temperature    temperature, Celcius
        pressure       pressure, mbar
        humidity       relative humidity, fraction
        frequency      observer frequency, GHz
        lapse_rate     tropospheric lapse rate, K/m

        * Parameters:
        dut1           UT1 correction
        ref_tol        Tolerance on refraction correction, in radians
        ref_delta      Refraction correction
        """

        for k,v in kwargs.items():
            try:
                self._set(k ,v)
            except KeyError:
                continue

    def get(self,*args):
        """
        Returns a dictionary of the requested state parameters.  If no
        parameters are supplied, then all are returned.  If a single parameter
        is supplied, then just that value is returned.  See QPoint.set for a
        list of parameter names.

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

    def refraction(self, *args, **kwargs):
        """
        Update refraction parameters

        Arguments (positional or keyword):

        q            observer orientation in horizon coordinates
        lat          latitude, degrees
        height       height above sea level, meters
        temperature  temperature, Celcius
        pressure     pressure, mbar
        humidity     humidity, fraction
        frequency    array frequency, GHz
        lapse_rate   tropospheric lapse rate, K/m
        tolerance    tolerance on convergence, radians
        delta        the refraction correction itself, in degrees

        If both el and lat are given, then the refraction correction in degrees
        is calculated, stored and returned after updating any other given
        parameters. Otherwise, the correction is returned w/out recalculating.

        Alternatively, if a single numerical argument, or the 'delta' keyword
        argument is given, then the correction is stored with this value
        instead of being recalculated.

        Numpy-vectorized for el and lat arguments.  Note that this is not
        an efficient vectorization, and only the last calculated value is
        stored for use in the coordinate conversion functions.
        """

        if len(args) == 1 and len(kwargs) == 0:
            v = args[0]
            self._set('ref_delta', v)
            return v

        if 'delta' in kwargs:
            v = kwargs.get('delta')
            self._set('ref_delta', v)
            return v

        arg_names = ['q','lat'] + self._funcs['weather'].keys()
        for idx,a in enumerate(args):
            kwargs[arg_names[idx]] = a

        for w in self._funcs['weather']:
            if w in kwargs:
                self._set(w, kwargs.get(w))

        q = kwargs.get('q',None)
        lat = kwargs.get('lat',None)
        if q is not None and lat is not None:
            def func(x0, x1, x2, x3, y):
                q = np.ascontiguousarray([x0,x1,x2,x3])
                return qp.qp_update_ref(self._memory, q, y)
            fvec = np.vectorize(func,[np.double])
            if q.size / 4 > 1:
                q = q.transpose()
            delta = fvec(q[0], q[1], q[2], q[3], lat)
            if delta.shape == ():
                return delta[()]
            return delta
        return self._get('ref_delta')

    def gmst(self, ctime, **kwargs):
        """
        Return Greenwich mean sidereal time for given ctimes and longitudes.
        Vectorized.

        Arguments:

        ctime      unix time in seconds UTC

        Keyword arguments:

        Any keywords accepted by the QPoint.set function can also be passed
        here, and will be processed prior to calculation.

        Outputs:

        gmst       Greenwich mean sidereal time of the observer
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
        Vectorized.

        Arguments:

        ctime      unix time in seconds UTC
        lon        observer longitude (degrees)

        Keyword arguments:

        Any keywords accepted by the QPoint.set function can also be passed
        here, and will be processed prior to calculation.

        Outputs:

        lmst       local mean sidereal time of the observer
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
        Vectorized.

        Arguments:

        ctime      unix time in seconds UTC
        ra         right ascension on the sky
        dec        declination on the sky

        Keyword arguments:

        Any keywords accepted by the QPoint.set function can also be passed
        here, and will be processed prior to calculation.

        Outputs:

        dipole     dipole amplitude in K
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

        Arguments:

        q_off      Detector offset quaternion for a single detector,
                   calculated using det_offset
        ctime      array of unix times in seconds UTC
        q_bore     Nx4 array of quaternions encoding the boresight orientation
                   on the sky (as output by azel2radec)

        Outputs:

        dipole     dipole amplitude in K
        """

        self.set(**kwargs)

        q_off = check_input('q_off', q_off, shape=(4,))
        ctime = check_input('ctime', ctime)
        n = ctime.size
        q_bore = check_input('q_bore', q_bore, shape=(n, 4))

        dipole = check_output('dipole', shape=ctime.shape, **kwargs)
        qp.qp_bore2dipole(self._memory, q_off, ctime, q_bore, dipole, n)

        if n == 1:
            return dipole[0]
        return dipole

    def det_offset(self, delta_az, delta_el, delta_psi):
        """
        Return quaternion corresponding to the requested detector offset.
        Vectorized.

        Arguments:

        delta_az   azimuthal offset of the detector (degrees)
        delta_el   elevation offset of the detector (degrees)
        delta_psi  polarization offset of the detector (degrees)

        Outputs:

        q          detector offset quaternion for each detector
        """

        delta_az, delta_el, delta_psi = \
            check_inputs(delta_az, delta_el, delta_psi)
        ndet = delta_az.size

        quat = check_output('quat', shape=(ndet,4))
        qp.qp_det_offsetn(delta_az, delta_el, delta_psi, quat, ndet)

        if ndet == 1:
            return quat[0]
        return quat

    def hwp_quat(self, theta):
        """
        Return quaternion corresponding to rotation by 2*theta,
        where theta is the physical waveplate angle.
        Vectorized.

        Arguments:

        theta      hwp physical angle (degrees)

        Outputs:

        q          quaternion for each hwp angle
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
        ctime. Input vectors must be numpy-array-like and of the same shape.

        Arguments:

        az         boresight azimuth (degrees)
        el         boresight elevation (degrees)
        pitch      boresight pitch (degrees); can be None
        roll       boresight pitch (degrees); can be None
        lon        observer longitude (degrees)
        lat        observer latitude (degrees)
        ctime      unix time in seconds UTC
        q          output quaternion array initialized by user

        Keywork arguments:

        Any keywords accepted by the QPoint.set function can also be passed
        here, and will be processed prior to calculation.

        Output:

        q          Nx4 numpy array of quaternions for each supplied timestamp.
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

    def bore2radec(self, q_off, ctime, q_bore, q_hwp=None, sindec=False,
                   ra=None, dec=None, sin2psi=None, cos2psi=None, **kwargs):
        """
        Calculate the orientation on the sky for a detector offset from the
        boresight.  Detector offsets are defined assuming the boresight is
        pointed toward the horizon, and that the boresight polarization axis is
        along the vertical.

        Arguments:

        q_off      Detector offset quaternion for a single detector,
                   calculated using det_offset
        ctime      array of unix times in seconds UTC
        q_bore     Nx4 array of quaternions encoding the boresight orientation
                   on the sky (as output by azel2radec)

        Keyword arguments:

        q_hwp      HWP angle quaternions calculated using hwp_quat
                   must be same shape as q_bore
        sindec     If True, return sin(dec) instead of dec in degrees
                   (default False)

        Any keywords accepted by the QPoint.set function can also be passed
        here, and will be processed prior to calculation.

        Outputs:

        ra         detector right ascension (degrees)
        dec/sindec detector declination (degrees) or sin(dec)
        sin2psi    detector polarization orientation
        cos2psi    detector polarization orientation
        """

        self.set(**kwargs)

        q_off  = check_input('q_off', q_off)
        q_bore = check_input('q_bore', q_bore)
        if ctime is None:
            if not self.get('mean_aber'):
                raise ValueError,'ctime required if mean_aber is False'
            ctime = np.zeros((q_bore.size/4,), dtype=q_bore.dtype)
        ctime  = check_input('ctime', ctime)
        ra = check_output('ra', ra, shape=ctime.shape, dtype=np.double)
        dec = check_output('dec', dec, shape=ctime.shape, dtype=np.double)
        sin2psi = check_output('sin2psi', sin2psi, shape=ctime.shape,
                                   dtype=np.double)
        cos2psi = check_output('cos2psi', cos2psi, shape=ctime.shape,
                                   dtype=np.double)
        n = ctime.size

        if q_hwp is None:
            if sindec:
                qp.qp_bore2rasindec(self._memory, q_off, ctime, q_bore,
                                    ra, dec, sin2psi, cos2psi, n)
            else:
                qp.qp_bore2radec(self._memory, q_off, ctime, q_bore,
                                 ra, dec, sin2psi, cos2psi, n)
        else:
            q_hwp = check_input('q_hwp', q_hwp, shape=q_bore.shape)
            if sindec:
                qp.qp_bore2rasindec_hwp(self._memory, q_off, ctime, q_bore,
                                        q_hwp, ra, dec, sin2psi, cos2psi, n)
            else:
                qp.qp_bore2radec_hwp(self._memory, q_off, ctime, q_bore,
                                     q_hwp, ra, dec, sin2psi, cos2psi, n)

        return ra, dec, sin2psi, cos2psi

    def azel2radec(self, delta_az, delta_el, delta_psi,
                   az, el, pitch, roll, lon, lat, ctime,
                   hwp=None, sindec=False, ra=None, dec=None,
                   sin2psi=None, cos2psi=None, **kwargs):
        """
        Estimate the orientation on the sky for a detector offset from
        boresight, given the boresight attitude (az/el/pitch/roll), location on
        the earth (lon/lat) and UTC time.  Input vectors must be
        numpy-array-like and of the same shape. Detector offsets are defined
        assuming the boresight is pointed toward the horizon, and that the
        boresight polarization axis is along the horizontal.

        Arguments:

        delta_az   azimuthal offset of the detector (degrees)
        delta_el   elevation offset of the detector (degrees)
        delta_psi  polarization offset of the detector (degrees)
        az         boresight azimuth (degrees)
        el         boresight elevation (degrees)
        pitch      boresight pitch (degrees); can be None
        roll       boresight roll (degrees); can be None
        lon        observer longitude (degrees)
        lat        observer latitude (degrees)
        ctime      unix time in seconds UTC

        Keyword arguments:

        hwp        HWP angles (degrees)
        sindec     If True, return sin(dec) instead of dec in degrees
                   (default False)

        Any keywords accepted by the QPoint.set function can also be passed
        here, and will be processed prior to calculation.

        Outputs:

        ra         detector right ascension (degrees)
        dec/sindec detector declination (degrees)
        sin2psi    detector polarization orientation
        cos2psi    detector polarization orientation
        """

        self.set(**kwargs)

        az, el, pitch, roll, lon, lat, ctime = \
            check_inputs(az, el, pitch, roll, lon, lat, ctime)

        ra = check_output('ra', ra, shape=az.shape, dtype=np.double)
        dec = check_output('dec', dec, shape=az.shape, dtype=np.double)
        sin2psi = check_output('sin2psi', sin2psi, shape=az.shape,
                                   dtype=np.double)
        cos2psi = check_output('cos2psi', cos2psi, shape=az.shape,
                                   dtype=np.double)
        n = az.size

        if hwp is None:
            if sindec:
                qp.qp_azel2rasindec(self._memory, delta_az, delta_el, delta_psi,
                                    az, el, pitch, roll, lon, lat, ctime,
                                    ra, dec, sin2psi, cos2psi, n)
            else:
                qp.qp_azel2radec(self._memory, delta_az, delta_el, delta_psi,
                                 az, el, pitch, roll, lon, lat, ctime,
                                 ra, dec, sin2psi, cos2psi, n)
        else:
            hwp = check_input('hwp', hwp, shape=az.shape)

            if sindec:
                qp.qp_azel2rasindec_hwp(self._memory, delta_az, delta_el, delta_psi,
                                        az, el, pitch, roll, lon, lat, ctime, hwp,
                                        ra, dec, sin2psi, cos2psi, n)
            else:
                qp.qp_azel2radec_hwp(self._memory, delta_az, delta_el, delta_psi,
                                     az, el, pitch, roll, lon, lat, ctime, hwp,
                                     ra, dec, sin2psi, cos2psi, n)

        return ra, dec, sin2psi, cos2psi

    def radecpa2quat(self, ra, dec, pa, **kwargs):
        """
        Calculate quaternion for input ra/dec/pa.
        """
        self.set(**kwargs)

        ra, dec, pa = check_inputs(ra, dec, pa)
        n = ra.size
        quat = check_output('quat', shape=(n,4), dtype=np.double,
                                **kwargs)
        qp.qp_radecpa2quatn(self._memory, ra, dec, pa, quat, n)

        if n == 1:
            return quat[0]
        return quat

    def quat2radecpa(self, quat, **kwargs):
        """
        Calculate ra/dec/pa for input quaternion(s).
        """
        self.set(**kwargs)

        quat = check_input('quat', np.atleast_2d(quat))
        n = quat.shape[0]
        ra = check_output('ra', shape=(n,), dtype=np.double, **kwargs)
        dec = check_output('dec', shape=(n,), dtype=np.double, **kwargs)
        pa = check_output('pa', shape=(n,), dtype=np.double, **kwargs)

        qp.qp_quat2radecpan(self._memory, quat, ra, dec, pa, n)
        if n == 1:
            return ra[0], dec[0], pa[0]
        return ra, dec, pa

    def radec2pix(self, ra, dec, nside=256, **kwargs):
        """
        Calculate healpix pixel number for given ra/dec and nside
        """

        self.set(**kwargs)

        ra, dec = check_inputs(ra, dec)
        n = ra.size

        pix = check_output('pix', shape=ra.shape, dtype=np.int, **kwargs)
        qp.qp_radec2pixn(self._memory, ra, dec, nside, pix, n)

        if n == 1:
            return pix[0]
        return pix

    def radec2gal(self, ra, dec, sin2psi, cos2psi, inplace=True, **kwargs):
        """
        Rotate celestial coordinates to galactic coordinates.
        """

        self.set(**kwargs)

        ra, dec, sin2psi, cos2psi = \
            check_inputs(ra, dec, sin2psi, cos2psi, inplace=inplace)
        n = ra.size

        qp.qp_radec2galn(self._memory, ra, dec, sin2psi, cos2psi, n)

        if n == 1:
            return ra[0], dec[0], sin2psi[0], cos2psi[0]
        return ra, dec, sin2psi, cos2psi

    def gal2radec(self, ra, dec, sin2psi, cos2psi, inplace=True, **kwargs):
        """
        Rotate celestial coordinates to galactic coordinates.
        """

        self.set(**kwargs)

        ra, dec, sin2psi, cos2psi = \
            check_inputs(ra, dec, sin2psi, cos2psi, inplace=inplace)
        n = ra.size

        qp.qp_gal2radecn(self._memory, ra, dec, sin2psi, cos2psi, n)

        if n == 1:
            return ra[0], dec[0], sin2psi[0], cos2psi[0]
        return ra, dec, sin2psi, cos2psi

    def rotate_map(self, map_in, coord=['C','G'], **kwargs):
        """
        Rotate a polarized npix-x-3 map from one coordinate system to another.
        Supported coordinates:

        C = celestial (J2000)
        G = galactic
        """

        from warnings import warn
        warn('This code is buggy, use at your own risk', UserWarning)

        from qmap_class import check_map
        map_in, nside = check_map(map_in)
        map_out = check_output(
            'map_out', map_out, shape=map_in.shape, fill=0)

        try:
            coord_in = coord[0]
            coord_out = coord[1]
        except:
            raise ValueError,'unable to parse coord'

        qp.qp_rotate_map(self._memory, nside, map_in, coord_in,
                         map_out, coord_out)
        return map_out

    def quat2pix(self, quat, nside=256, pol=True, **kwargs):
        """
        Calculate healpix pixel number and polarization angle given
        quaternion and nside
        """

        self.set(**kwargs)

        quat = check_input('quat', np.atleast_2d(quat))

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
                 **kwargs):
        """
        Calculate the orientation on the sky for a detector offset from the
        boresight.  Detector offsets are defined assuming the boresight is
        pointed toward the horizon, and that the boresight polarization axis is
        along the vertical.

        Arguments:

        q_off      Detector offset quaternion for a single detector,
                   calculated using det_offset
        ctime      array of unix times in seconds UTC
        q_bore     Nx4 array of quaternions encoding the boresight orientation
                   on the sky (as output by azel2radec)

        Keyword arguments:

        q_hwp      HWP angle quaternions calculated using hwp_quat
                   must be same shape as q_bore
        nside      map dimension
        pol        if False, return only the pixel timestream

        Any keywords accepted by the QPoint.set function can also be passed
        here, and will be processed prior to calculation.

        Outputs:

        pix        detector pixel number
        sin2psi    detector polarization orientation (if pol is True)
        cos2psi    detector polarization orientation (if pol is True)
        """

        self.set(**kwargs)

        q_off  = check_input('q_off', q_off)
        q_bore = check_input('q_bore', q_bore)
        if ctime is None:
            if not self.get('mean_aber'):
                raise ValueError,'ctime required if mean_aber is False'
            ctime = np.zeros((q_bore.size/4,), dtype=q_bore.dtype)
        ctime  = check_input('ctime', ctime)
        pix  = check_output('pix', shape=ctime.shape,
                                dtype=np.int, **kwargs)
        sin2psi = check_output('sin2psi', shape=ctime.shape,
                                   **kwargs)
        cos2psi = check_output('cos2psi', shape=ctime.shape,
                                   **kwargs)
        n = ctime.size

        if q_hwp is None:
            qp.qp_bore2pix(self._memory, q_off, ctime, q_bore,
                           nside, pix, sin2psi, cos2psi, n)
        else:
            q_hwp = check_input('q_hwp', q_hwp, shape=q_bore.shape)

            qp.qp_bore2pix_hwp(self._memory, q_off, ctime, q_bore,
                               q_hwp, nside, pix, sin2psi, cos2psi, n)

        if pol is True:
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
            of shape (nsample,)
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

        from qmap_class import check_map
        map_in, nside = check_map(map_in)

        val = check_output('value', shape=(len(map_in), n))

        for m, v in zip(map_in, val):
            qp.qp_get_interp_valn(self._memory, nside, m, ra, dec, v, n)

        self.set(pix_order=pix_order)

        v = v.squeeze()
        if not v.shape:
            return v[()]
        return v

    def load_bulletin_a(self, filename, columns=['mjd','dut1','x','y'], **kwargs):
        """
        Load IERS Bulletin A from file and store in memory.  The file must be
        readable using numpy.loadtxt with unpack=True, and is assumed to be sorted
        by mjd.

        Keyword arguments:

        columns    list of columns as they appear in the file.
                   A KeyError is raise if the list does not contain
                   each of ['mjd', 'dut1', 'x', 'y']

        Any other keyword arguments are passed to the numpy.loadtxt function

        Output:

        mjd, dut1, x, y
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
