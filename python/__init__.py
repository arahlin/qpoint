"""qpoint

A lightweight library for efficient pointing.

Based on M. Nolta's libactpol.
Uses the SOFA Software Collection, available from http://www.iausofa.org/
"""

import numpy as _np
from _libqpoint import libqp as _libqp

class QPoint(object):
    
    def __init__(self, **kwargs):
        """
        Initialize a QPoint memory instance for keeping track of pointing
        corrections over time.
        
        Any keyword arguments are passed to QPoint.set to update memory.
        """
        
        # initialize memory
        self._memory = _libqp.qp_init_memory()
        
        # collect all parameter functions
        from _libqpoint import qp_funcs
        self._funcs = qp_funcs
        self._all_funcs = dict()
        for k,v in qp_funcs.items():
            self._all_funcs.update(**v)
        
        # set any requested parameters
        self.set(**kwargs)
        
    def __del__(self):
        """
        Free memory before deleting the object
        """
        _libqp.qp_free_memory(self._memory)
    
    def _set(self, key, val):
        """
        Set a single parameter to the given value.  See QPoint.set for a list
        of parameter names.
        """
        if key not in self._all_funcs:
            raise KeyError,'Unknown parameter %s' % key
        val = self._all_funcs[key]['check_set'](val)
        self._all_funcs[key]['set'](self._memory, val)
    
    def _get(self, key):
        """
        Get the value for a single parameter.  See QPoint.set for a list of
        parameter names.
        """
        if key not in self._all_funcs:
            raise KeyError,'Unknown parameter %s' % key
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
        pair_dets      If True, A/B detectors are paired in bore2map
                       (factor of 2 speed up in computation, but assumes ideality)
        pix_order      'nest' or 'ring' for healpix pixel ordering
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
            self._set(k ,v)
    
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
        _libqp.qp_reset_rates(self._memory)
    
    def refraction(self, *args, **kwargs):
        """
        Update refraction parameters
        
        Arguments (positional or keyword):
        
        el           elevation angle, degrees
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
            v = kwargs.get('ref_delta')
            self._set('ref_delta', v)
            return v
        
        arg_names = ['el','lat'] + self._funcs['weather'].keys()
        for idx,a in enumerate(args):
            kwargs[arg_names[idx]] = a
        
        for w in self._funcs['weather']:
            if w in kwargs:
                self._set(w, kwargs.get(w))
        
        el = kwargs.get('el',None)
        lat = kwargs.get('lat',None)
        if el is not None and lat is not None:
            def func(x, y): return _libqp.qp_update_ref(self._memory, x, y)
            fvec = _np.vectorize(func,[_np.double])
            delta = fvec(el, lat)
            if delta.shape == ():
                return delta[()]
            return delta
        return self._get('ref_delta')
    
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
        
        ctime = _np.array(_np.atleast_1d(ctime), dtype=_np.double)
        lon   = _np.array(_np.atleast_1d(lon),   dtype=_np.double)
        
        if ctime.shape != lon.shape:
            raise ValueError, "input vectors must have the same shape"
        
        if ctime.size == 1:
            return _libqp.qp_lmst(self._memory, ctime[0], lon[0])
        
        lmst = _np.empty(ctime.shape, dtype=_np.double)
        _libqp.qp_lmstn(self._memory, ctime, lon, lmst, ctime.size)
        return lmst
        
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
        
        delta_az, delta_el, delta_psi \
            = [_np.array(_np.atleast_1d(x), dtype=_np.double) for x in
               _np.broadcast_arrays(delta_az, delta_el, delta_psi)]
        ndet = delta_az.size
        
        for x in (delta_el, delta_psi):
            if x.shape != delta_az.shape:
                raise ValueError, "input offset vectors must have the same shape"
        
        if ndet == 1:
            quat = _np.empty((4,), dtype=_np.double)
            _libqp.qp_det_offset(delta_az[0], delta_el[0], delta_psi[0], quat)
        else:
            quat = _np.empty(delta_az.shape + (4,), dtype=_np.double)
            _libqp.qp_det_offsetn(delta_az, delta_el, delta_psi, quat, ndet)
        
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
        theta = _np.array(_np.atleast_1d(theta), dtype=_np.double)
        if theta.size == 1:
            quat = _np.empty((4,), dtype=_np.double)
            _libqp.qp_hwp_quat(theta[0], quat)
        else:
            quat = _np.empty(theta.shape + (4,), dtype=_np.double)
            _libqp.qp_hwp_quatn(theta, quat, theta.size)
        return quat
    
    def azel2bore(self, az, el, pitch, roll, lon, lat, ctime, **kwargs):
        """
        Estimate the quaternion for the boresight orientation on the sky given
        the attitude (az/el/pitch/roll), location on the earth (lon/lat) and
        ctime. Input vectors must be numpy-array-like and of the same shape.
        
        Arguments:
        
        az         boresight azimuth (degrees)
        el         boresight elevation (degrees)
        pitch      boresight pitch (degrees)
        roll       boresight pitch (degrees)
        lon        observer longitude (degrees)
        lat        observer latitude (degrees)
        ctime      unix time in seconds UTC
        
        Keywork arguments:
            
        Any keywords accepted by the QPoint.set function can also be passed
        here, and will be processed prior to calculation.
        
        Output:
        
        q          Nx4 numpy array of quaternions for each supplied timestamp.
        """
        
        self.set(**kwargs)
        
        az    = _np.asarray(az,    dtype=_np.double)
        el    = _np.asarray(el,    dtype=_np.double)
        pitch = _np.asarray(pitch, dtype=_np.double)
        roll  = _np.asarray(roll,  dtype=_np.double)
        lon   = _np.asarray(lon,   dtype=_np.double)
        lat   = _np.asarray(lat,   dtype=_np.double)
        ctime = _np.asarray(ctime, dtype=_np.double)
        q = _np.empty(az.shape + (4,), dtype=_np.double)
        n = az.size
    
        for x in (el,pitch,roll,lon,lat,ctime):
            if x.shape != az.shape:
                raise ValueError,"input vectors must have the same shape"
        
        _libqp.qp_azel2bore(self._memory, az, el, pitch, roll, lon, lat,
                            ctime, q, n)
        
        return q
    
    def bore2radec(self, q_off, ctime, q_bore, q_hwp=None, sindec=False, **kwargs):
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
        
        q_off  = _np.asarray(q_off,  dtype=_np.double)
        ctime  = _np.asarray(ctime,  dtype=_np.double)
        q_bore = _np.asarray(q_bore, dtype=_np.double)
        ra  = _np.empty(ctime.shape, dtype=_np.double)
        dec = _np.empty(ctime.shape, dtype=_np.double)
        sin2psi = _np.empty(ctime.shape, dtype=_np.double)
        cos2psi = _np.empty(ctime.shape, dtype=_np.double)
        n = ctime.size
        
        if q_bore.shape != ctime.shape + (4,):
            raise ValueError,\
                'ctime and q_bore must have compatible shapes (N,) and (N,4)'
        
        if q_hwp is None:
            if sindec:
                _libqp.qp_bore2rasindec(self._memory, q_off, ctime, q_bore,
                                        ra, dec, sin2psi, cos2psi, n)
            else:
                _libqp.qp_bore2radec(self._memory, q_off, ctime, q_bore,
                                     ra, dec, sin2psi, cos2psi, n)
        else:
            q_hwp = _np.asarray(q_hwp, dtype=_np.double)
            if q_bore.shape != q_hwp.shape:
                raise ValueError,\
                    'q_bore and q_hwp must have the same shape'
            
            if sindec:
                _libqp.qp_bore2rasindec_hwp(self._memory, q_off, ctime, q_bore,
                                            q_hwp, ra, dec, sin2psi, cos2psi, n)
            else:
                _libqp.qp_bore2radec_hwp(self._memory, q_off, ctime, q_bore,
                                         q_hwp, ra, dec, sin2psi, cos2psi, n)
        
        return ra, dec, sin2psi, cos2psi
    
    def azel2radec(self, delta_az, delta_el, delta_psi,
                   az, el, pitch, roll, lon, lat, ctime,
                   hwp=None, sindec=False, **kwargs):
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
        pitch      boresight pitch (degrees)
        roll       boresight roll (degrees)
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
        
        az    = _np.asarray(az,    dtype=_np.double)
        el    = _np.asarray(el,    dtype=_np.double)
        pitch = _np.asarray(pitch, dtype=_np.double)
        roll  = _np.asarray(roll,  dtype=_np.double)
        lon   = _np.asarray(lon,   dtype=_np.double)
        lat   = _np.asarray(lat,   dtype=_np.double)
        ctime = _np.asarray(ctime, dtype=_np.double)
        ra  = _np.empty(az.shape, dtype=_np.double)
        dec = _np.empty(az.shape, dtype=_np.double)
        sin2psi = _np.empty(az.shape, dtype=_np.double)
        cos2psi = _np.empty(az.shape, dtype=_np.double)
        n = az.size
        
        for x in (el, pitch, roll, lon, lat, ctime):
            if x.shape != az.shape:
                raise ValueError,"input vectors must have the same shape"
        
        if hwp is None:
            if sindec:
                _libqp.qp_azel2rasindec(self._memory, delta_az, delta_el, delta_psi,
                                        az, el, pitch, roll, lon, lat, ctime,
                                        ra, dec, sin2psi, cos2psi, n)
            else:
                _libqp.qp_azel2radec(self._memory, delta_az, delta_el, delta_psi,
                                     az, el, pitch, roll, lon, lat, ctime,
                                     ra, dec, sin2psi, cos2psi, n)
        else:
            hwp = _np.asarray(hwp, dtype=_np.double)
            if hwp.shape != az.shape:
                raise ValueError,"input vectors must have the same shape"
            
            if sindec:
                _libqp.qp_azel2rasindec_hwp(self._memory, delta_az, delta_el, delta_psi,
                                            az, el, pitch, roll, lon, lat, ctime, hwp,
                                            ra, dec, sin2psi, cos2psi, n)
            else:
                _libqp.qp_azel2radec_hwp(self._memory, delta_az, delta_el, delta_psi,
                                         az, el, pitch, roll, lon, lat, ctime, hwp,
                                         ra, dec, sin2psi, cos2psi, n)
        
        return ra, dec, sin2psi, cos2psi
        
    def bore2map(self, q_off, ctime, q_bore, nside=256,
                 q_hwp=None, pmap=None, **kwargs):
        """
        Calculate hits map for given detectors and boresight orientations.
        Returns a npix-x-6 array containing (hits, p01, p02, p11, p12, p22).
        Optionally accumulate hits in place by supplying an existing pmap array.
        """
        
        self.set(**kwargs)
        
        q_off = _np.asarray(q_off, dtype=_np.double)
        ndet = q_off.size/4
        
        ctime  = _np.asarray(ctime,  dtype=_np.double)
        q_bore = _np.asarray(q_bore, dtype=_np.double)
        n = ctime.size
        
        if ctime.shape[0] != q_bore.shape[0]:
            raise ValueError,"input ctime and q_bore must have compatible shape"
        
        mshape = (12*nside*nside, 6)
        if pmap is None:
            pmap = _np.zeros(mshape, dtype=_np.double)
        else:
            nside = int(_np.sqrt(pmap.size/6/12))
            mshape = (12*nside*nside, 6)
        if pmap.shape != mshape:
            print pmap.shape
            raise ValueError,"input pmap must have shape %s" % str(mshape)
        if pmap.dtype != _np.double:
            raise TypeError,"input pmap must be of type numpy.double"
        
        if q_hwp is None:
            _libqp.qp_bore2map(self._memory, q_off, ndet, ctime, q_bore, n,
                               pmap, nside)
        else:
            q_hwp = _np.asarray(q_hwp, dtype=_np.double)
            if q_hwp.shape != q_bore.shape:
                raise ValueError, "input q_hwp must have the same shape as q_bore"
            
            _libqp.qp_bore2map_hwp(self._memory, q_off, ndet, ctime, q_bore,
                                   q_hwp, n, pmap, nside)
        
        return pmap
        
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
            raise KeyError,\
                'Missing columns %s' % str(list(set(req_columns)-set(columns)))
        kwargs['unpack'] = True
        data = _np.loadtxt(filename, **kwargs)
        mjd, x, y, dut1 = (data[columns.index(x)] for x in req_columns)
        mjd_min, mjd_max = int(mjd[0]), int(mjd[-1])
        
        try:
            _libqp.set_iers_bulletin_a(self._memory, mjd_min, mjd_max, dut1, x, y)
        except:
            raise RuntimeError, \
                'Error loading Bulletin A data from file %s' % filename
        
        return mjd, dut1, x, y
    
    def get_bulletin_a(self, mjd):
        """
        Return dut1/x/y for given mjd. Numpy-vectorized.
        """
        
        from _libqpoint import get_bulletin_a
        def func(x): return get_bulletin_a(self._memory, x)
        fvec = _np.vectorize(func, [_np.double]*3)
        
        out = fvec(mjd)
        if out[0].shape == ():
            return tuple(x[()] for x in out)
        return out

def refraction(el, lat, height, temp, press, hum,
               freq=150., lapse=0.0065, tol=1e-8):
    """
    Standalone function for calculating the refraction correction without
    storing any parameters.  Useful for testing, numpy-vectorized.
    
    Arguments:
    
    el           elevation angle, degrees
    lat          latitude, degrees
    height       height above sea level, meters
    temperature  temperature, Celcius
    pressure     pressure, mbar
    humidity     humidity, fraction
    frequency    array frequency, GHz
    lapse_rate   tropospheric lapse rate, K/m
    tolerance    tolerance on convergence, radians
    
    Output:
    
    delta        refraction correction, in degrees
    """
    
    fvec = _np.vectorize(_libqp.qp_refraction,[_np.double])
    delta = fvec(el, lat, height, temp, press, hum, freq, lapse, tol)
    if delta.shape == ():
        return delta[()]
    return delta

# for debugging
def _plot_diff(ang1,ang2,asec=True,n=None):
    scale = 3600 if asec else 1
    if n is None:
        n = len(ang1[0])
    
    dra = (ang1[0]-ang2[0])*scale
    ddec = (ang1[1]-ang2[1])*scale
    ds2p = (ang1[2]-ang2[2])
    dc2p = (ang1[3]-ang2[3])
    
    import pylab
    pylab.figure();
    ax = pylab.subplot(411);
    pylab.plot(dra[:n]);
    pylab.subplot(412,sharex=ax);
    pylab.plot(ddec[:n]);
    pylab.subplot(413,sharex=ax);
    pylab.plot(ds2p[:n]);
    pylab.subplot(414,sharex=ax);
    pylab.plot(dc2p[:n]);

