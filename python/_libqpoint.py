# ctypes wrapper for libqpoint

import ctypes as ct
import numpy as np
import os

# structures

class qp_state_t(ct.Structure):
    _fields_ = [
        ('update_rate', ct.c_double),
        ('ctime_last', ct.c_double)
        ]

class qp_weather_t(ct.Structure):
    _fields_ = [
        ('height', ct.c_double),
        ('temperature', ct.c_double),
        ('pressure', ct.c_double),
        ('humidity', ct.c_double),
        ('frequency', ct.c_double),
        ('lapse_rate', ct.c_double),
        ]

class qp_bulletina_entry_t(ct.Structure):
    _fields_ = [
        ('x', ct.c_float),
        ('y', ct.c_float),
        ('dut1', ct.c_float),
        ]

class qp_bulletina_t(ct.Structure):
    _fields_ = [
        ('entries', ct.POINTER(qp_bulletina_entry_t)),
        ('mjd_min', ct.c_int),
        ('mjd_max', ct.c_int),
        ]

class qp_memory_t(ct.Structure):
    _fields_ = [
        ('initialized', ct.c_int),
        ('state_daber', qp_state_t),
        ('state_lonlat', qp_state_t),
        ('state_wobble', qp_state_t),
        ('state_dut1', qp_state_t),
        ('state_erot', qp_state_t),
        ('state_npb', qp_state_t),
        ('state_aaber', qp_state_t),
        ('state_ref', qp_state_t),
        ('weather', qp_weather_t),
        ('ref_tol', ct.c_double),
        ('ref_delta', ct.c_double),
        ('dut1', ct.c_double),
        ('q_lonlat', ct.c_double * 4),
        ('q_wobble', ct.c_double * 4),
        ('q_npb', ct.c_double * 4),
        ('q_erot', ct.c_double * 4),
        ('beta_earth', ct.c_double * 3),
        ('bulletinA', qp_bulletina_t),
        ('accuracy', ct.c_int),
        ('mean_aber', ct.c_int),
        ('fast_math', ct.c_int),
        ('polconv', ct.c_int),
        ('pair_dets', ct.c_int),
        ('pix_order', ct.c_int),
        ('num_threads', ct.c_int),
        ]

# library functions

libqp = np.ctypeslib.load_library('libqpoint',os.path.dirname(__file__))
NDP = np.ctypeslib.ndpointer
qp_memory_t_p = ct.POINTER(qp_memory_t)

libqp.qp_init_memory.restype = qp_memory_t_p
libqp.qp_free_memory.argtypes = (qp_memory_t_p,)

libqp.qp_azel2radec.argtypes = (qp_memory_t_p, # params
                                ct.c_double, ct.c_double, ct.c_double, # offset
                                NDP(dtype=np.double), # az
                                NDP(dtype=np.double), # el
                                NDP(dtype=np.double), # pitch
                                NDP(dtype=np.double), # roll
                                NDP(dtype=np.double), # lon
                                NDP(dtype=np.double), # lat
                                NDP(dtype=np.double), # ctime
                                NDP(dtype=np.double), # ra
                                NDP(dtype=np.double), # dec
                                NDP(dtype=np.double), # sin2psi
                                NDP(dtype=np.double), # cos2psi
                                ct.c_int) # n

libqp.qp_azel2radec_hwp.argtypes = (qp_memory_t_p, # params
                                    ct.c_double, ct.c_double, ct.c_double, # offset
                                    NDP(dtype=np.double), # az
                                    NDP(dtype=np.double), # el
                                    NDP(dtype=np.double), # pitch
                                    NDP(dtype=np.double), # roll
                                    NDP(dtype=np.double), # lon
                                    NDP(dtype=np.double), # lat
                                    NDP(dtype=np.double), # ctime
                                    NDP(dtype=np.double), # hwp
                                    NDP(dtype=np.double), # ra
                                    NDP(dtype=np.double), # dec
                                    NDP(dtype=np.double), # sin2psi
                                    NDP(dtype=np.double), # cos2psi
                                    ct.c_int) # n

libqp.qp_azel2rasindec.argtypes = (qp_memory_t_p, # params
                                   ct.c_double, ct.c_double, ct.c_double, # offset
                                   NDP(dtype=np.double), # az
                                   NDP(dtype=np.double), # el
                                   NDP(dtype=np.double), # pitch
                                   NDP(dtype=np.double), # roll
                                   NDP(dtype=np.double), # lon
                                   NDP(dtype=np.double), # lat
                                   NDP(dtype=np.double), # ctime
                                   NDP(dtype=np.double), # ra
                                   NDP(dtype=np.double), # sindec
                                   NDP(dtype=np.double), # sin2psi
                                   NDP(dtype=np.double), # cos2psi
                                   ct.c_int) # n

libqp.qp_azel2rasindec_hwp.argtypes = (qp_memory_t_p, # params
                                       ct.c_double, ct.c_double, ct.c_double, # offset
                                       NDP(dtype=np.double), # az
                                       NDP(dtype=np.double), # el
                                       NDP(dtype=np.double), # pitch
                                       NDP(dtype=np.double), # roll
                                       NDP(dtype=np.double), # lon
                                       NDP(dtype=np.double), # lat
                                       NDP(dtype=np.double), # ctime
                                       NDP(dtype=np.double), # hwp
                                       NDP(dtype=np.double), # ra
                                       NDP(dtype=np.double), # sindec
                                       NDP(dtype=np.double), # sin2psi
                                       NDP(dtype=np.double), # cos2psi
                                       ct.c_int) # n

libqp.qp_azel2bore.argtypes = (qp_memory_t_p, # params
                               NDP(dtype=np.double), # az
                               NDP(dtype=np.double), # el
                               NDP(dtype=np.double), # pitch
                               NDP(dtype=np.double), # roll
                               NDP(dtype=np.double), # lon
                               NDP(dtype=np.double), # lat
                               NDP(dtype=np.double), # ctime
                               NDP(dtype=np.double), # q
                               ct.c_int) # n

libqp.qp_det_offset.argtypes = (ct.c_double, ct.c_double, ct.c_double, # offset
                               NDP(dtype=np.double)) # quat

libqp.qp_det_offsetn.argtypes = (NDP(dtype=np.double), # delta_az
                                NDP(dtype=np.double), # delta_el
                                NDP(dtype=np.double), # delta_psi
                                NDP(dtype=np.double), # quat
                                ct.c_int) # n

libqp.qp_hwp_quat.argtypes = (ct.c_double, # ang
                              NDP(dtype=np.double)) # quat

libqp.qp_hwp_quatn.argtypes = (NDP(dtype=np.double), # ang
                               NDP(dtype=np.double), # quat
                               ct.c_int) # n

libqp.qp_gmst.argtypes = (qp_memory_t_p, # params
                          ct.c_double) # ctime
libqp.qp_gmst.restype = ct.c_double # gmst

libqp.qp_gmstn.argtypes = (qp_memory_t_p, # params
                           NDP(dtype=np.double), # ctime
                           NDP(dtype=np.double), # gmst
                           ct.c_int) # n

libqp.qp_lmst.argtypes = (qp_memory_t_p, # params
                          ct.c_double, # ctime
                          ct.c_double) # lon
libqp.qp_lmst.restype = ct.c_double # lmst

libqp.qp_lmstn.argtypes = (qp_memory_t_p, # params
                           NDP(dtype=np.double), # ctime
                           NDP(dtype=np.double), # lon
                           NDP(dtype=np.double), # lmst
                           ct.c_int) # n

libqp.qp_dipole.argtypes = (qp_memory_t_p, # params
                            ct.c_double, # ctime
                            ct.c_double, # ra
                            ct.c_double) # dec
libqp.qp_dipole.restype = ct.c_double # dipole

libqp.qp_dipolen.argtypes = (qp_memory_t_p, # params
                             NDP(dtype=np.double), # ctime
                             NDP(dtype=np.double), # ra
                             NDP(dtype=np.double), # dec
                             NDP(dtype=np.double), # dipole
                             ct.c_int) # n

libqp.qp_bore2radec.argtypes = (qp_memory_t_p, # params
                                NDP(dtype=np.double), # offset
                                NDP(dtype=np.double), # ctime
                                NDP(dtype=np.double), # q
                                NDP(dtype=np.double), # ra
                                NDP(dtype=np.double), # dec
                                NDP(dtype=np.double), # sin2psi
                                NDP(dtype=np.double), # cos2psi
                                ct.c_int) # n

libqp.qp_bore2radec_hwp.argtypes = (qp_memory_t_p, # params
                                    NDP(dtype=np.double), # offset
                                    NDP(dtype=np.double), # ctime
                                    NDP(dtype=np.double), # q_bore
                                    NDP(dtype=np.double), # q_hwp
                                    NDP(dtype=np.double), # ra
                                    NDP(dtype=np.double), # dec
                                    NDP(dtype=np.double), # sin2psi
                                    NDP(dtype=np.double), # cos2psi
                                    ct.c_int) # n

libqp.qp_bore2rasindec.argtypes = (qp_memory_t_p, # params
                                   NDP(dtype=np.double), # offset
                                   NDP(dtype=np.double), # ctime
                                   NDP(dtype=np.double), # q
                                   NDP(dtype=np.double), # ra
                                   NDP(dtype=np.double), # sindec
                                   NDP(dtype=np.double), # sin2psi
                                   NDP(dtype=np.double), # cos2psi
                                   ct.c_int) # n

libqp.qp_bore2rasindec_hwp.argtypes = (qp_memory_t_p, # params
                                       NDP(dtype=np.double), # offset
                                       NDP(dtype=np.double), # ctime
                                       NDP(dtype=np.double), # q_bore
                                       NDP(dtype=np.double), # q_hwp
                                       NDP(dtype=np.double), # ra
                                       NDP(dtype=np.double), # sindec
                                       NDP(dtype=np.double), # sin2psi
                                       NDP(dtype=np.double), # cos2psi
                                       ct.c_int) # n

libqp.qp_radec2pix.argtypes = (qp_memory_t_p, # params
                               ct.c_int, # nside
                               ct.c_double, # ra
                               ct.c_double) # dec
libqp.qp_radec2pix.restype = ct.c_int

libqp.qp_radec2pixn.argtypes = (qp_memory_t_p, # params
                                ct.c_int, # nside
                                NDP(dtype=np.double), # ra
                                NDP(dtype=np.double), # dec
                                NDP(dtype=np.int), # pix
                                ct.c_int) # n

libqp.qp_bore2map.argtypes = (qp_memory_t_p, # params
                              NDP(dtype=np.double), # offsets
                              ct.c_int, # ndet
                              NDP(dtype=np.double), # ctime
                              NDP(dtype=np.double), # q_bore
                              ct.c_int, # n
                              NDP(dtype=np.double), # pmap
                              ct.c_int) # nside

libqp.qp_bore2map_hwp.argtypes = (qp_memory_t_p, # params
                                  NDP(dtype=np.double), # offsets
                                  ct.c_int, # ndet
                                  NDP(dtype=np.double), # ctime
                                  NDP(dtype=np.double), # q_bore
                                  NDP(dtype=np.double), # q_hwp
                                  ct.c_int, # n
                                  NDP(dtype=np.double), # pmap
                                  ct.c_int) # nside

libqp.set_iers_bulletin_a.argtypes = (qp_memory_t_p,
                                      ct.c_int, ct.c_int, # mjd_min, mjd_max
                                      NDP(dtype=np.double), # dut1
                                      NDP(dtype=np.double), # x
                                      NDP(dtype=np.double)) # y
libqp.set_iers_bulletin_a.restype = ct.c_int

libqp.get_iers_bulletin_a.argtypes = (qp_memory_t_p,
                                      ct.c_double, # mjd
                                      ct.POINTER(ct.c_double), # dut1
                                      ct.POINTER(ct.c_double), # x
                                      ct.POINTER(ct.c_double)) # y
libqp.get_iers_bulletin_a.restype = ct.c_int

def get_bulletin_a(mem, mjd):
    dut1 = ct.c_double()
    x = ct.c_double()
    y = ct.c_double()
    
    libqp.get_iers_bulletin_a(mem, mjd, ct.byref(dut1),
                              ct.byref(x), ct.byref(y))
    return dut1.value, x.value, y.value

libqp.qp_refraction.argtypes = (ct.c_double, # elevation angle
                                ct.c_double, # height
                                ct.c_double, # temperature
                                ct.c_double, # pressure
                                ct.c_double, # humidity
                                ct.c_double, # frequency
                                ct.c_double, # latitude
                                ct.c_double, # lapse_rate
                                 ct.c_double) # tolerance
libqp.qp_refraction.restype = ct.c_double

libqp.qp_update_ref.argtypes = (qp_memory_t_p, # params
                                  ct.c_double, # el
                                  ct.c_double) # lat
libqp.qp_update_ref.restype = ct.c_double

# parameters and options

def check_set_float(val):
    if not np.isscalar(val):
        raise TypeError,'val must be a scalar value'
    return float(val)

def check_set_int(val):
    if not np.isscalar(val):
        raise TypeError,'val must be a scalar value'
    return int(val)

def check_pass(val):
    return val

def set_rfunc(state):
    f = libqp['qp_set_rate_%s'%state]
    f.argtypes = (qp_memory_t_p,ct.c_double)
    return f
def reset_rfunc(state):
    f = libqp['qp_reset_rate_%s'%state]
    f.argtypes = (qp_memory_t_p,)
    return f
def get_rfunc(state):
    f = libqp['qp_get_rate_%s'%state]
    f.argtypes = (qp_memory_t_p,)
    f.restype = ct.c_double
    return f

def check_set_state(state):
    rdict = {'always':0,'once':-1,'never':-999}
    state = rdict.get(state,state)
    return check_set_float(state)

def check_get_state(state):
    rdict = {0:'always',-1:'once',-999:'never'}
    state = rdict.get(state,state)
    return state

states = ['lonlat','npb','erot','daber','aaber','wobble','dut1','ref']
state_funcs = dict()
for s in states:
    k = 'rate_%s' % s
    state_funcs[k] = dict()
    state_funcs[k]['set'] = set_rfunc(s)
    state_funcs[k]['reset'] = reset_rfunc(s)
    state_funcs[k]['get'] = get_rfunc(s)
    state_funcs[k]['check_set'] = check_set_state
    state_funcs[k]['check_get'] = check_get_state

def set_wfunc(par):
    f = libqp['qp_set_weather_%s'%par]
    f.argtypes = (qp_memory_t_p,ct.c_double)
    return f

def get_wfunc(par):
    f = libqp['qp_get_weather_%s'%par]
    f.argtypes = (qp_memory_t_p,)
    f.restype = ct.c_double
    return f

weather_params = ['height','temperature','pressure','humidity',
                  'frequency','lapse_rate']
weather_funcs = dict()
for w in weather_params:
    weather_funcs[w] = dict()
    weather_funcs[w]['set'] = set_wfunc(w)
    weather_funcs[w]['get'] = get_wfunc(w)
    weather_funcs[w]['check_set'] = check_set_float
    weather_funcs[w]['check_get'] = check_pass

def set_ofunc(option):
    f = libqp['qp_set_opt_%s'%option]
    f.argtypes = (qp_memory_t_p,ct.c_int)
    return f

def get_ofunc(option):
    f = libqp['qp_get_opt_%s'%option]
    f.argtypes = (qp_memory_t_p,)
    f.restype = ct.c_int
    return f

def check_set_accuracy(acc):
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

def check_get_accuracy(ac):
    if ac == 1:
        return 'low'
    return 'high'

def check_set_mean_aber(ab):
    if ab is None:
        return 0
    if ab in [False, 0]:
        return 0
    if ab in [True, 1]:
        return 1
    return 0

def check_get_mean_aber(ab):
    if ab == 1:
        return True
    return False

def check_set_fast_math(fast):
    if fast is None:
        return 0
    if fast in [True,1]:
        return 1
    if fast in [False,0]:
        return 0
    return 0

def check_get_fast_math(fast):
    if fast == 1:
        return True
    return False

def check_set_polconv(pol):
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

def check_get_polconv(pol):
    if pol == 1:
        return 'iau'
    return 'healpix'

def check_set_pair_dets(pair):
    if pair is None:
        return 0
    if pair in [1,True]:
        return 1
    if pair in [0,False]:
        return 0
    return 0

def check_get_pair_dets(pair):
    if pair == 1:
        return True
    return False

def check_set_pix_order(order):
    if order is None:
        return 0
    if isinstance(order,basestring):
        if order.lower() in ['ring']:
            return 0
        if order.lower() in ['nest','nested']:
            return 1
    if order in [0,1]:
        return order
    return 0

def check_get_pix_order(order):
    if order == 1:
        return 'nest'
    return 'ring'

def check_set_num_threads(nt):
    if nt is None:
        return 0
    return int(nt)

def check_get_num_threads(nt):
    return nt

options = ['accuracy','mean_aber','fast_math','polconv','pair_dets',
           'pix_order','num_threads']
option_funcs = dict()
for p in options:
    option_funcs[p] = dict()
    option_funcs[p]['set'] = set_ofunc(p)
    option_funcs[p]['get'] = get_ofunc(p)
    option_funcs[p]['check_set'] = globals()['check_set_%s'%p]
    option_funcs[p]['check_get'] = globals()['check_get_%s'%p]

def set_pfunc(par):
    f = libqp['qp_set_%s'%par]
    f.argtypes = (qp_memory_t_p,ct.c_double)
    return f

def get_pfunc(par):
    f = libqp['qp_get_%s'%par]
    f.argtypes = (qp_memory_t_p,)
    f.restype = ct.c_double
    return f

double_params = ['ref_tol','ref_delta','dut1']
double_funcs = dict()
for p in double_params:
    double_funcs[p] = dict()
    double_funcs[p]['set'] = set_pfunc(p)
    double_funcs[p]['get'] = get_pfunc(p)
    double_funcs[p]['check_set'] = check_set_float
    double_funcs[p]['check_get'] = check_pass

qp_funcs = dict()
qp_funcs['rates'] = state_funcs
qp_funcs['options'] = option_funcs
qp_funcs['weather'] = weather_funcs
qp_funcs['params'] = double_funcs
