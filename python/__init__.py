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
        
        # initialize memoery
        self._memory = _libqp.qp_init_memory()
        
        # collect all parameter functions
        from _libqpoint import qp_funcs
        self._funcs = qp_funcs
        self._all_funcs = dict()
        for k,v in qp_funcs.items():
            self._all_funcs.update(**v)
        
        # set an reqeusted parameters
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
        
        Note that this is not a vectorized function.
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
            delta = _libqp.qp_update_ref(self._memory, el, lat)
        else:
            delta = self._get('ref_delta')
        return delta
    
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
    
    def bore2radec(self, delta_az, delta_el, delta_psi, ctime, q_bore,
                   **kwargs):
        """
        Calculate the orientation on the sky for a detector offset from the
        boresight.  Detector offsets are defined assuming the boresight is
        pointed toward the horizon, and that the boresight polarization axis is
        along the vertical.
        
        Arguments:
        
        delta_az   azimuthal offset of the detector (degrees)
        delta_el   elevation offset of the detector (degrees)
        delta_psi  polarization angle of the detector (degrees)
        ctime      array of unix times in seconds UTC
        q_bore     Nx4 array of quaternions encoding the boresight orientation 
                   on the sky (as output by azel2radec)
        
        Keyword arguments:
        
        Any keywords accepted by the QPoint.set function can also be passed
        here, and will be processed prior to calculation.    
        
        Outputs:
        
        ra         detector right ascension (degrees)
        dec        detector declination (degrees)
        sin2psi    detector polarization orientation
        cos2psi    detector polarization orientation
        """
        
        self.set(**kwargs)
        
        ctime  = _np.asarray(ctime,  dtype=_np.double)
        q_bore = _np.asarray(q_bore, dtype=_np.double)
        ra  = _np.empty(ctime.shape, dtype=_np.double)
        dec = _np.empty(ctime.shape, dtype=_np.double)
        sin2psi = _np.empty(ctime.shape, dtype=_np.double)
        cos2psi = _np.empty(ctime.shape, dtype=_np.double)
        n = ctime.size
        
        if q_bore.shape != ctime.shape + (4,):
            raise ValueError,\
                'ctime and q must have compatible shapes (N,) and (N,4)'
        
        _libqp.qp_bore2radec(self._memory, delta_az, delta_el, delta_psi,
                             ctime, q_bore, ra, dec, sin2psi, cos2psi, n)
        
        return ra, dec, sin2psi, cos2psi
    
    def bore2rasindec(self, delta_az, delta_el, delta_psi, ctime, q_bore,
                      **kwargs):
        """
        Calculate the orientation on the sky for a detector offset from the
        boresight. Detector offsets are defined assuming the boresight is
        pointed toward the horizon, and that the boresight polarization axis is
        along the vertical.
        
        Arguments:
        
        delta_az   azimuthal offset of the detector (degrees)
        delta_el   elevation offset of the detector (degrees)
        delta_psi  polarization angle of the detector (degrees)
        ctime      array of unix times in seconds UTC
        q_bore     Nx4 array of quaternions encoding the boresight orientation 
                   on the sky (as output by azel2radec)
        
        Keyword arguments:
        
        Any keywords accepted by the QPoint.set function can also be passed
        here, and will be processed prior to calculation.    
        
        Outputs:
        
        ra         detector right ascension (degrees)
        sindec     detector sin(declination)
        sin2psi    detector polarization orientation
        cos2psi    detector polarization orientation
        """
        
        self.set(**kwargs)
        
        ctime  = _np.asarray(ctime,  dtype=_np.double)
        q_bore = _np.asarray(q_bore, dtype=_np.double)
        ra  = _np.empty(ctime.shape, dtype=_np.double)
        sindec = _np.empty(ctime.shape, dtype=_np.double)
        sin2psi = _np.empty(ctime.shape, dtype=_np.double)
        cos2psi = _np.empty(ctime.shape, dtype=_np.double)
        n = ctime.size
        
        if q_bore.shape != ctime.shape + (4,):
            raise ValueError, \
                'ctime and q must have compatible shapes (N,) and (N,4)'
        
        _libqp.qp_bore2rasindec(self._memory, delta_az, delta_el, delta_psi,
                                ctime, q_bore, ra, sindec, sin2psi, cos2psi, n)
        
        return ra, sindec, sin2psi, cos2psi
    
    def azel2radec(self, delta_az, delta_el, delta_psi,
                   az, el, pitch, roll, lon, lat, ctime, **kwargs):
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
        
        Any keywords accepted by the QPoint.set function can also be passed
        here, and will be processed prior to calculation.    
        
        Outputs:
        
        ra         detector right ascension (degrees)
        dec        detector declination (degrees)
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
            
        _libqp.qp_azel2radec(self._memory, delta_az, delta_el, delta_psi,
                             az, el, pitch, roll, lon, lat, ctime,
                             ra, dec, sin2psi, cos2psi, n)
        
        return ra, dec, sin2psi, cos2psi
    
    def azel2rasindec(self, delta_az, delta_el, delta_psi,
                      az, el, pitch, roll, lon, lat, ctime, **kwargs):
        """
        Estimate the orientation on the sky for a detector offset from
        boresight, given the boresight attitude (az/el/pitch/roll), location on
        the earth (lon/lat) and UTC time.  Input vectors must by
        numpy-array-like and of the same shape.  Detector offsets are defined
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
        
        Any keywords accepted by the QPoint.set function can also be passed
        here, and will be processed prior to calculation.
        
        Outputs:
        
        ra         detector right ascension (degrees)
        sindec     detector sin(declination)
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
        sindec = _np.empty(az.shape, dtype=_np.double)
        sin2psi = _np.empty(az.shape, dtype=_np.double)
        cos2psi = _np.empty(az.shape, dtype=_np.double)
        n = az.size
        
        for x in (el, pitch, roll, lon, lat, ctime):
            if x.shape != az.shape:
                raise ValueError,"input vectors must have the same shape"
        
        _libqp.qp_azel2rasindec(self._memory, delta_az, delta_el, delta_psi,
                                az, el, pitch, roll, lon, lat, ctime,
                                ra, sindec, sin2psi, cos2psi, n)
        
        return ra, sindec, sin2psi, cos2psi

    def load_bulletin_a(filename, columns=['mjd','dut1','x','y'], **kwargs):
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
        Return dut1/x/y for given mjd. Vectorized.
        """
        
        from _libqpoint import get_bulletin_a
        def func(x): return get_bulletin_a(self._memory, x)
        fvec = _np.vectorize(func, [_np.double, _np.double, _np.double])
        
        return fvec(mjd)

def refraction(el, lat, height, temp, press, hum,
               freq=150., lapse=0.0065, tol=1e-8):
    """
    Standalone function for calculating the refraction correction without
    storing any parameters.  Useful for testing.  Note that this is not a
    vectorized function.
    
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
    
    return _libqp.qp_refraction(el, lat, height, temp, press, hum,
                                freq, lapse, tol)

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

