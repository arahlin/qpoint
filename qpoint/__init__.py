"""qpoint

A lightweight library for efficient pointing.

Based in M. Nolta's libactpol.
Uses the SOFA Software Collection, available from http://www.iausofa.org/
"""

import ctypes as _ct
import numpy as _np
import os as _os

_qp = _np.ctypeslib.load_library('libqpoint',_os.path.dirname(__file__))
_ndp = _np.ctypeslib.ndpointer

_qp.qp_azel2radec.argtypes = (_ct.c_double, _ct.c_double, _ct.c_double, # offset
                              _ndp(dtype=_np.double), # az
                              _ndp(dtype=_np.double), # el
                              _ndp(dtype=_np.double), # pitch
                              _ndp(dtype=_np.double), # roll
                              _ndp(dtype=_np.double), # lon
                              _ndp(dtype=_np.double), # lat
                              _ndp(dtype=_np.double), # ctime
                              _ndp(dtype=_np.double), # ra
                              _ndp(dtype=_np.double), # dec
                              _ndp(dtype=_np.double), # sin2psi
                              _ndp(dtype=_np.double), # cos2psi
                              _ct.c_int) # n

_qp.qp_azel2rasindec.argtypes = (_ct.c_double, _ct.c_double, _ct.c_double, # offset
                                 _ndp(dtype=_np.double), # az
                                 _ndp(dtype=_np.double), # el
                                 _ndp(dtype=_np.double), # pitch
                                 _ndp(dtype=_np.double), # roll
                                 _ndp(dtype=_np.double), # lon
                                 _ndp(dtype=_np.double), # lat
                                 _ndp(dtype=_np.double), # ctime
                                 _ndp(dtype=_np.double), # ra
                                 _ndp(dtype=_np.double), # sindec
                                 _ndp(dtype=_np.double), # sin2psi
                                 _ndp(dtype=_np.double), # cos2psi
                                 _ct.c_int) # n

_qp.qp_azel2bore.argtypes = (_ndp(dtype=_np.double), # az
                             _ndp(dtype=_np.double), # el
                             _ndp(dtype=_np.double), # pitch
                             _ndp(dtype=_np.double), # roll
                             _ndp(dtype=_np.double), # lon
                             _ndp(dtype=_np.double), # lat
                             _ndp(dtype=_np.double), # ctime
                             _ndp(dtype=_np.double), # q
                             _ct.c_int) # n

_qp.qp_bore2radec.argtypes = (_ct.c_double, _ct.c_double, _ct.c_double, # offset
                              _ndp(dtype=_np.double), # ctime
                              _ndp(dtype=_np.double), # q
                              _ndp(dtype=_np.double), # ra
                              _ndp(dtype=_np.double), # dec
                              _ndp(dtype=_np.double), # sin2psi
                              _ndp(dtype=_np.double), # cos2psi
                              _ct.c_int) # n

_qp.qp_bore2rasindec.argtypes = (_ct.c_double, _ct.c_double, _ct.c_double, # offset
                                 _ndp(dtype=_np.double), # ctime
                                 _ndp(dtype=_np.double), # q
                                 _ndp(dtype=_np.double), # ra
                                 _ndp(dtype=_np.double), # sindec
                                 _ndp(dtype=_np.double), # sin2psi
                                 _ndp(dtype=_np.double), # cos2psi
                                 _ct.c_int) # n

def _set_sfunc(step):
    f = _qp['qp_set_%s'%step]
    f.argtypes = (_ct.c_double,)
    return f
def _reset_sfunc(step):
    f = _qp['qp_reset_%s'%step]
    return f
def _get_sfunc(step):
    f = _qp['qp_get_%s'%step]
    f.restype = _ct.c_double
    return f

_steps = ['lonlat','npb','erot','daber','aaber','wobble','dut1']
_step_funcs = dict()
for _s in _steps:
    _step_funcs[_s] = dict()
    _step_funcs[_s]['set'] = _set_sfunc(_s)
    _step_funcs[_s]['reset'] = _reset_sfunc(_s)
    _step_funcs[_s]['get'] = _get_sfunc(_s)

def _set_pfunc(param):
    f = _qp['qp_set_%s'%param]
    f.argtypes = (_ct.c_int,)
    return f

def _get_pfunc(param):
    f = _qp['qp_get_%s'%param]
    f.restype = _ct.c_int
    return f

def _check_accuracy(acc):
    if acc is None:
        return 0
    if acc in [0,1]:
        return acc
    if isinstance(acc,basestring):
        if acc.lower() in ['high']:
            return 0
        if acc.lower() in ['low']:
            return 1
    return 0

def _check_get_accuracy(ac):
    if ac == 1:
        return 'low'
    return 'high'

def _check_mean_aber(ab):
    if ab is None:
        return 0
    if ab in [False, 0]:
        return 0
    if ab in [True, 1]:
        return 1
    return 0

def _check_get_mean_aber(ab):
    if ab == 1:
        return True
    return False

def _check_fast_math(fast):
    if fast is None:
        return 0
    if fast in [True,1]:
        return 1
    if fast in [False,0]:
        return 0
    return 0

def _check_get_fast_math(fast):
    if fast == 1:
        return True
    return False

def _check_polconv(pol):
    if pol is None:
        return 0
    if isinstance(pol,basestring):
        if pol.lower() in ['cosmo','healpix']:
            return 0
        if pol.lower() in ['iau']:
            return 1
    if pol in [0,1]:
        return pol
    return 0

def _check_get_polconv(pol):
    if pol == 1:
        return 'iau'
    return 'healpix'

_params = ['accuracy','mean_aber','fast_math','polconv']
_param_funcs = dict()
for _p in _params:
    _param_funcs[_p] = dict()
    _param_funcs[_p]['set'] = _set_pfunc(_p)
    _param_funcs[_p]['get'] = _get_pfunc(_p)
    _param_funcs[_p]['check'] = globals()['_check_%s'%_p]
    _param_funcs[_p]['check_get'] = globals()['_check_get_%s'%_p]

def set_params(**kwargs):
    """
    Set the update rates, in seconds for each step of the azel2radec
    coordinate transformation.
    
    For each step the rate can be set to 'always' (0), 'once' (-1),
    'never' (-999), or as a period in seconds.
    
    Available keywords are:
    
    daber     Rate at which the diurnal aberration correction is applied
              (NB: this can only be applied always or never)
    lonlat    Rate at which observer's lon and lat are updated
    wobble    Rate at which the polar motion correction is updated
              (NB: these are not estimated for dates beyond a year from now)
    dut1      Rate at which the ut1-utc correction is updated
              (NB: this is not estimated for dates beyond a year from now)
    erot      Rate at which the earth's rotation angle is updated
    npb       Rate at which the nutation/precession/frame-bias terms
              are updated
    aaber     Rate at which the annual aberration correction
              (due to the earth's orbital velocity) is updated

    accuracy   If 'low', use a truncated form (2000b) for the NPB correction,
               which is much faster but less accurate.
               If 'high' (default), use the full 2006/2000a form.
    mean_aber  If True, apply the aberration correction as an average for the
               entire field of view.  This is gives a 1-2 arcsec deviation
               at the edges of the SPIDER field of view.
    fast_math  If True, use polynomial approximations for trig functions
    polconv    Specify the 'cosmo' or 'iau' polarization convention
    """

    init_params()
    rdict = {'always':0,'once':-1,'never':-999}
    for step in _steps:
        v = kwargs.pop(step,None)
        if v is None: continue
        v = rdict.get(v,v)
        if not _np.isscalar(v):
            raise TypeError,'rate %s must be a scalar value' % step
        _step_funcs[step]['set'](v)
    for param in _params:
        v = kwargs.pop(param,None)
        if v is None: continue
        v = _param_funcs[param]['check'](v)
        if not _np.isscalar(v):
            raise TypeError,'parameter %s must be a scalar value' % param
        _param_funcs[param]['set'](v)

def get_params():
    """
    Returns a dictionary of the current update rates for each step.
    """
    init_params()
    rdict = {0:'always',-1:'once',-999:'never'}
    state = dict()
    for step in _steps:
        v = _step_funcs[step]['get']()
        state[step] = rdict.get(v,v)
    for param in _params:
        v = _param_funcs[param]['get']()
        state[param] = _param_funcs[param]['check_get'](v)
    return state

def reset_params():
    """
    Reset update counters for each step.
    Useful to force an updated correction term at the beginning of each chunk.
    """
    _qp.qp_reset_params()

def init_params():
    """
    Initialize the parameter structure for keeping track of each step.
    If it has already been initialized, then nothing is done.
    """
    _qp.qp_init_params()

def azel2bore(az, el, pitch, roll, lon, lat, ctime, **kwargs):
    """
    Estimate the quaternion for the boresight orientation on the sky
    given the attitude (az/el/pitch/roll), location on the earth (lon/lat)
    and ctime. Input vectors must be numpy-array-like and of the same shape.
    
    Arguments:
    
    az         boresight azimuth (degrees)
    el         boresight elevation (degrees)
    pitch      boresight pitch (degrees)
    roll       boresight pitch (degrees)
    lon        observer longitude (degrees)
    lat        observer latitude (degrees)
    ctime      unix time in seconds UTC
    
    Keywork arguments:
        
    Any keywords accepted by the set_params function can also be passed here,
    and will be processed prior to calculation.
    
    Output:
    
    q          Nx4 numpy array of quaternions for each supplied timestamp.
    """
    
    set_params(**kwargs)
    
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
    
    _qp.qp_azel2bore(az, el, pitch, roll, lon, lat, ctime, q, n)
    
    return q

def bore2radec(delta_az, delta_el, delta_psi, ctime, q_bore, **kwargs):
    """
    Calculate the orientation on the sky for a detector offset from the boresight.
    Detector offsets are defined assuming the boresight is pointed toward the horizon,
    and that the boresight polarization axis is along the vertical.
    
    Arguments:
    
    delta_az   azimuthal offset of the detector (degrees)
    delta_el   elevation offset of the detector (degrees)
    delta_psi  polarization angle of the detector (degrees)
    ctime      array of unix times in seconds UTC
    q_bore     Nx4 array of quaternions encoding the boresight orientation on the sky
               (as output by azel2radec)
    
    Keyword arguments:
    
    Any keywords accepted by the set_params function can also be passed here,
    and will be processed prior to calculation.    
    
    Outputs:
    
    ra         detector right ascension (degrees)
    dec        detector declination (degrees)
    sin2psi    detector polarization orientation
    cos2psi    detector polarization orientation
    """
    
    set_params(**kwargs)
    
    ctime  = _np.asarray(ctime,  dtype=_np.double)
    q_bore = _np.asarray(q_bore, dtype=_np.double)
    ra  = _np.empty(ctime.shape, dtype=_np.double)
    dec = _np.empty(ctime.shape, dtype=_np.double)
    sin2psi = _np.empty(ctime.shape, dtype=_np.double)
    cos2psi = _np.empty(ctime.shape, dtype=_np.double)
    n = ctime.size
    
    if q_bore.shape != ctime.shape + (4,):
        raise ValueError,'ctime and q must have compatible shapes (N,) and (N,4)'
    
    _qp.qp_bore2radec(delta_az, delta_el, delta_psi, ctime, q_bore,
                      ra, dec, sin2psi, cos2psi, n)
    
    return ra, dec, sin2psi, cos2psi
    
def bore2rasindec(delta_az, delta_el, delta_psi, ctime, q_bore, **kwargs):
    """
    Calculate the orientation on the sky for a detector offset from the boresight.
    Detector offsets are defined assuming the boresight is pointed toward the horizon,
    and that the boresight polarization axis is along the vertical.
    
    Arguments:
    
    delta_az   azimuthal offset of the detector (degrees)
    delta_el   elevation offset of the detector (degrees)
    delta_psi  polarization angle of the detector (degrees)
    ctime      array of unix times in seconds UTC
    q_bore     Nx4 array of quaternions encoding the boresight orientation on the sky
               (as output by azel2radec)
    
    Keyword arguments:
    
    Any keywords accepted by the set_params function can also be passed here,
    and will be processed prior to calculation.    
    
    Outputs:
    
    ra         detector right ascension (degrees)
    sindec     detector sin(declination)
    sin2psi    detector polarization orientation
    cos2psi    detector polarization orientation
    """
    
    set_params(**kwargs)
    
    ctime  = _np.asarray(ctime,  dtype=_np.double)
    q_bore = _np.asarray(q_bore, dtype=_np.double)
    ra  = _np.empty(ctime.shape, dtype=_np.double)
    sindec = _np.empty(ctime.shape, dtype=_np.double)
    sin2psi = _np.empty(ctime.shape, dtype=_np.double)
    cos2psi = _np.empty(ctime.shape, dtype=_np.double)
    n = ctime.size
    
    if q_bore.shape != ctime.shape + (4,):
        raise ValueError,'ctime and q must have compatible shapes (N,) and (N,4)'
    
    _qp.qp_bore2rasindec(delta_az, delta_el, delta_psi, ctime, q_bore,
                         ra, sindec, sin2psi, cos2psi, n)
    
    return ra, sindec, sin2psi, cos2psi
    
def azel2radec(delta_az, delta_el, delta_psi,
               az, el, pitch, roll, lon, lat, ctime, **kwargs):
    """
    Estimate the orientation on the sky for a detector offset from boresight,
    given the boresight attitude (az/el/pitch/roll), location on the earth (lon/lat)
    and UTC time.  Input vectors must by numpy-array-like and of the same shape.
    Detector offsets are defined assuming the boresight is pointed toward the horizon,
    and that the boresight polarization axis is along the horizontal.
    
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
    
    Any keywords accepted by the set_params function can also be passed here,
    and will be processed prior to calculation.    

    Outputs:
    
    ra         detector right ascension (degrees)
    dec        detector declination (degrees)
    sin2psi    detector polarization orientation
    cos2psi    detector polarization orientation
    """
    
    set_params(**kwargs)
    
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
        
    _qp.qp_azel2radec(delta_az, delta_el, delta_psi,
                      az, el, pitch, roll, lon, lat, ctime,
                      ra, dec, sin2psi, cos2psi, n)
    
    return ra, dec, sin2psi, cos2psi

def azel2rasindec(delta_az, delta_el, delta_psi,
                  az, el, pitch, roll, lon, lat, ctime, **kwargs):
    """
    Estimate the orientation on the sky for a detector offset from boresight,
    given the boresight attitude (az/el/pitch/roll), location on the earth (lon/lat)
    and UTC time.  Input vectors must by numpy-array-like and of the same shape.
    Detector offsets are defined assuming the boresight is pointed toward the horizon,
    and that the boresight polarization axis is along the horizontal.
    
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
    
    Any keywords accepted by the set_params function can also be passed here,
    and will be processed prior to calculation.    

    Outputs:
    
    ra         detector right ascension (degrees)
    sindec     detector sin(declination)
    sin2psi    detector polarization orientation
    cos2psi    detector polarization orientation
    """
    
    set_params(**kwargs)
    
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
    
    _qp.qp_azel2rasindec(delta_az, delta_el, delta_psi,
                         az, el, pitch, roll, lon, lat, ctime,
                         ra, sindec, sin2psi, cos2psi, n)
    
    return ra, sindec, sin2psi, cos2psi

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
