# ctypes wrapper for libqpoint

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import ctypes as ct
import numpy as np
import os

libqp = np.ctypeslib.load_library('libqpoint',os.path.dirname(__file__))

# **********************************************************************
# Types
# **********************************************************************

NDP = np.ctypeslib.ndpointer

quat_t = NDP(np.double, ndim=1, shape=(4,), flags=['A','C'])
quat_t_p = NDP(np.double, ndim=2, flags=['A','C'])
mueller_t = NDP(np.double, ndim=1, shape=(4,), flags=['A','C'])
mueller_t_p = NDP(np.double, ndim=2, flags=['A', 'C'])
vec3_t = NDP(np.double, ndim=1, shape=(3,), flags=['A','C'])
vec3_t_p = NDP(np.double, ndim=2, flags=['A','C'])

wquat_t = NDP(np.double, ndim=1, shape=(4,), flags=['A','C','W'])
wquat_t_p = NDP(np.double, ndim=2, flags=['A','C','W'])
wvec3_t = NDP(np.double, ndim=1, shape=(3,), flags=['A','C','W'])
wvec3_t_p = NDP(np.double, ndim=2, flags=['A','C','W'])

arr = NDP(np.double, ndim=1, flags=['A','C'])
arrf = NDP(np.uint8, ndim=1, flags=['A','C'])
warr = NDP(np.double, ndim=1, flags=['A','C','W'])
warri = NDP(np.int, ndim=1, flags=['A','C','W'])

arr2 = NDP(np.uintp, ndim=1, flags=['A','C'])
larr = NDP(np.long, ndim=1, flags=['A','C'])

QP_DO_ALWAYS = ct.c_int.in_dll(libqp, "QP_DO_ALWAYS").value
QP_DO_ONCE = ct.c_int.in_dll(libqp, "QP_DO_ONCE").value
QP_DO_NEVER = ct.c_int.in_dll(libqp, "QP_DO_NEVER").value

class qp_state_t(ct.Structure):
    _fields_ = [
        ('update_rate', ct.c_double),
        ('ctime_last', ct.c_double)
        ]

class qp_weather_t(ct.Structure):
    _fields_ = [
        ('temperature', ct.c_double),
        ('pressure', ct.c_double),
        ('humidity', ct.c_double),
        ('frequency', ct.c_double),
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
        ('init', ct.c_int),

        ('state_daber', qp_state_t),
        ('state_lonlat', qp_state_t),
        ('state_wobble', qp_state_t),
        ('state_dut1', qp_state_t),
        ('state_erot', qp_state_t),
        ('state_npb', qp_state_t),
        ('state_aaber', qp_state_t),
        ('state_ref', qp_state_t),

        ('state_daber_inv', qp_state_t),
        ('state_lonlat_inv', qp_state_t),
        ('state_wobble_inv', qp_state_t),
        ('state_dut1_inv', qp_state_t),
        ('state_erot_inv', qp_state_t),
        ('state_npb_inv', qp_state_t),
        ('state_aaber_inv', qp_state_t),
        ('state_ref_inv', qp_state_t),

        ('weather', qp_weather_t),
        ('ref_delta', ct.c_double),
        ('q_ref', ct.c_double * 4),
        ('q_ref_inv', ct.c_double * 4),
        ('dut1', ct.c_double),
        ('q_lonlat', ct.c_double * 4),
        ('q_lonlat_inv', ct.c_double * 4),
        ('q_wobble', ct.c_double * 4),
        ('q_wobble_inv', ct.c_double * 4),
        ('q_npb', ct.c_double * 4),
        ('q_npb_inv', ct.c_double * 4),
        ('q_erot', ct.c_double * 4),
        ('q_erot_inv', ct.c_double * 4),
        ('q_gal', ct.c_double * 4),
        ('q_gal_inv', ct.c_double * 4),
        ('gal_init', ct.c_int),
        ('v_dipole', ct.c_double * 3),
        ('dipole_init', ct.c_int),
        ('beta_earth', ct.c_double * 3),
        ('beta_rot', ct.c_double * 3),
        ('bulletinA', qp_bulletina_t),

        ('accuracy', ct.c_int),
        ('mean_aber', ct.c_int),
        ('fast_math', ct.c_int),
        ('polconv', ct.c_int),
        ('pix_order', ct.c_int),
        ('interp_pix', ct.c_int),
        ('fast_pix', ct.c_int),
        ('error_missing', ct.c_int),
        ('nan_missing', ct.c_int),
        ('interp_missing', ct.c_int),
        ('num_threads', ct.c_int),
        ('thread_num', ct.c_int),
        ]

qp_memory_t_p = ct.POINTER(qp_memory_t)

QP_STRUCT_INIT = 1
QP_STRUCT_MALLOC = 2
QP_ARR_INIT_PTR = 4
QP_ARR_MALLOC_1D = 8
QP_ARR_MALLOC_2D = 16

class qp_det_t(ct.Structure):
    _fields_ = [
        ('init', ct.c_int),
        ('q_off', ct.c_double * 4),
        ('weight', ct.c_double),
        ('gain', ct.c_double),
        ('mueller', ct.c_double * 4),
        ('n', ct.c_size_t),
        ('tod_init', ct.c_int),
        ('tod', ct.POINTER(ct.c_double)),
        ('flag_init', ct.c_int),
        ('flag', ct.POINTER(ct.c_uint8)),
        ('weights_init', ct.c_int),
        ('weights', ct.POINTER(ct.c_double)),
        ]
qp_det_t_p = ct.POINTER(qp_det_t)

class qp_detarr_t(ct.Structure):
    _fields_ = [
        ('init', ct.c_int),
        ('n', ct.c_size_t),
        ('arr_init', ct.c_int),
        ('diff', ct.c_size_t),
        ('arr', qp_det_t_p)
        ]
qp_detarr_t_p = ct.POINTER(qp_detarr_t)

class qp_point_t(ct.Structure):
    _fields_ = [
        ('init', ct.c_int),
        ('n', ct.c_size_t),
        ('q_bore_init', ct.c_int),
        ('q_bore', ct.POINTER(ct.c_double * 4)),
        ('ctime_init', ct.c_int),
        ('ctime', ct.POINTER(ct.c_double)),
        ('q_hwp_init', ct.c_int),
        ('q_hwp', ct.POINTER(ct.c_double * 4))
        ]
qp_point_t_p = ct.POINTER(qp_point_t)

qp_vec_mode = ct.c_uint
QP_VEC_NONE = 0
QP_VEC_TEMP = 1
QP_VEC_POL = 2
QP_VEC_VPOL = 3
QP_VEC_D1 = 4
QP_VEC_D1_POL = 5
QP_VEC_D2 = 6
QP_VEC_D2_POL = 7
vec_modes = {1  : QP_VEC_TEMP,
             3  : {True  : QP_VEC_POL,
                   False : QP_VEC_D1},
             4  : QP_VEC_VPOL,
             6  : QP_VEC_D2,
             9  : QP_VEC_D1_POL,
             18 : QP_VEC_D2_POL}

qp_proj_mode = ct.c_uint
QP_PROJ_NONE = 0
QP_PROJ_TEMP = 1
QP_PROJ_POL = 2
QP_PROJ_VPOL = 3
proj_modes = {1  : QP_PROJ_TEMP,
              6  : QP_PROJ_POL,
              10 : QP_PROJ_VPOL}

class qp_map_t(ct.Structure):
    _fields_ = [
        ('init', ct.c_int),
        ('partial', ct.c_int),
        ('nside', ct.c_size_t),
        ('npix', ct.c_size_t),
        ('pixinfo_init', ct.c_int),
        ('pixinfo', ct.c_void_p),
        ('pixhash_init', ct.c_int),
        ('pixhash', ct.c_void_p),
        ('num_vec', ct.c_size_t),
        ('vec_mode', qp_vec_mode),
        ('vec1d_init', ct.c_int),
        ('vec1d', ct.POINTER(ct.c_double)),
        ('vec_init', ct.c_int),
        ('vec', ct.POINTER(ct.POINTER(ct.c_double))),
        ('num_proj', ct.c_size_t),
        ('proj_mode', qp_proj_mode),
        ('proj1d_init', ct.c_int),
        ('proj1d', ct.POINTER(ct.c_double)),
        ('proj_init', ct.c_int),
        ('proj', ct.POINTER(ct.POINTER(ct.c_double)))
        ]
qp_map_t_p = ct.POINTER(qp_map_t)

def pointer_2d(d):
    return (d.__array_interface__['data'][0] +
            np.arange(d.shape[0]) * d.strides[0]).astype(np.uintp)

def as_ctypes(d):
    return np.ctypeslib.as_ctypes(d)

def setargs(fname, arg=None, res=None):
    func = getattr(libqp, fname)
    if arg is not None and not isinstance(arg, (tuple, list)):
        arg = (arg,)
    func.argtypes = arg
    func.restype = res

# **********************************************************************
# Functions
# **********************************************************************

setargs('qp_init_memory', res=qp_memory_t_p)
setargs('qp_free_memory', arg=qp_memory_t_p)
setargs('qp_print_memory', arg=qp_memory_t_p)

setargs('qp_reset_rates', arg=qp_memory_t_p)
setargs('qp_reset_inv_rates', arg=qp_memory_t_p)

setargs('qp_get_error_code', arg=qp_memory_t_p, res=ct.c_int)
setargs('qp_get_error_string', arg=qp_memory_t_p, res=ct.c_char_p)

setargs('qp_azel2radec',
        arg=(qp_memory_t_p, # params
             ct.c_double, ct.c_double, ct.c_double, # offset
             arr, arr, arr, arr, arr, arr, arr, # a/e/p/r/l/l/t
             warr, warr, warr, warr, ct.c_int))
setargs('qp_azelpsi2radec',
        arg=(qp_memory_t_p, # params
             ct.c_double, ct.c_double, ct.c_double, # offset
             arr, arr, arr, arr, arr, arr, arr, arr, # a/e/p/p/r/l/l/t
             warr, warr, warr, warr, ct.c_int))
setargs('qp_azel2radecpa',
        arg=(qp_memory_t_p, # params
             ct.c_double, ct.c_double, ct.c_double, # offset
             arr, arr, arr, arr, arr, arr, arr, # a/e/p/r/l/l/t
             warr, warr, warr, ct.c_int))
setargs('qp_azelpsi2radecpa',
        arg=(qp_memory_t_p, # params
             ct.c_double, ct.c_double, ct.c_double, # offset
             arr, arr, arr, arr, arr, arr, arr, arr, # a/e/p/p/r/l/l/t
             warr, warr, warr, ct.c_int))
setargs('qp_radec2azel',
        arg=(qp_memory_t_p, # params
             arr, arr, arr, arr, arr, arr, arr, arr, arr, #r/d/p/l/l/t/a/e/p
             ct.c_int))
setargs('qp_azel2radec_hwp',
        arg=(qp_memory_t_p, # params
             ct.c_double, ct.c_double, ct.c_double, # offset
             arr, arr, arr, arr, arr, arr, arr, arr, # a/e/p/r/l/l/t/hwp
             warr, warr, warr, warr, ct.c_int))
setargs('qp_azelpsi2radec_hwp',
        arg=(qp_memory_t_p, # params
             ct.c_double, ct.c_double, ct.c_double, # offset
             arr, arr, arr, arr, arr, arr, arr, arr, arr, # a/e/p/p/r/l/l/t/hwp
             warr, warr, warr, warr, ct.c_int))
setargs('qp_azel2radecpa_hwp',
        arg=(qp_memory_t_p, # params
             ct.c_double, ct.c_double, ct.c_double, # offset
             arr, arr, arr, arr, arr, arr, arr, arr, # a/e/p/r/l/l/t/hwp
             warr, warr, warr, ct.c_int))
setargs('qp_azelpsi2radecpa_hwp',
        arg=(qp_memory_t_p, # params
             ct.c_double, ct.c_double, ct.c_double, # offset
             arr, arr, arr, arr, arr, arr, arr, arr, arr, # a/e/p/p/r/l/l/t/hwp
             warr, warr, warr, ct.c_int))
setargs('qp_azel2rasindec',
        arg=(qp_memory_t_p, # params
             ct.c_double, ct.c_double, ct.c_double, # offset
             arr, arr, arr, arr, arr, arr, arr, # a/e/p/r/l/l/t
             warr, warr, warr, warr, ct.c_int))
setargs('qp_azelpsi2rasindec',
        arg=(qp_memory_t_p, # params
             ct.c_double, ct.c_double, ct.c_double, # offset
             arr, arr, arr, arr, arr, arr, arr, arr, # a/e/p/p/r/l/l/t
             warr, warr, warr, warr, ct.c_int))
setargs('qp_azel2rasindec_hwp',
        arg=(qp_memory_t_p, # params
             ct.c_double, ct.c_double, ct.c_double, # offset
             arr, arr, arr, arr, arr, arr, arr, arr, # a/e/p/r/l/l/t/hwp
             warr, warr, warr, warr, ct.c_int))
setargs('qp_azelpsi2rasindec_hwp',
        arg=(qp_memory_t_p, # params
             ct.c_double, ct.c_double, ct.c_double, # offset
             arr, arr, arr, arr, arr, arr, arr, arr, arr, # a/e/p/p/r/l/l/t/hwp
             warr, warr, warr, warr, ct.c_int))
setargs('qp_azel2bore',
        arg=(qp_memory_t_p, # params
             arr, arr, arr, arr, arr, arr, arr, # a/e/p/r/l/l/t
             wquat_t_p, ct.c_int))
setargs('qp_azelpsi2bore',
        arg=(qp_memory_t_p, # params
             arr, arr, arr, arr, arr, arr, arr, arr, # a/e/p/p/r/l/l/t
             wquat_t_p, ct.c_int))

setargs('qp_det_offsetn', arg=(arr, arr, arr, wquat_t_p, ct.c_int))
setargs('qp_bore_offset', arg=(qp_memory_t_p, wquat_t_p, arr, arr, arr,
                               ct.c_int, ct.c_int))
setargs('qp_hwp_quatn', arg=(arr, wquat_t_p, ct.c_int))
setargs('qp_gmstn', arg=(qp_memory_t_p, arr, warr, ct.c_int))
setargs('qp_lmstn', arg=(qp_memory_t_p, arr, arr, warr, ct.c_int))
setargs('qp_dipolen', arg=(qp_memory_t_p, arr, arr, arr, warr, ct.c_int))
setargs('qp_bore2dipole',
        arg=(qp_memory_t_p, quat_t, arr, quat_t_p, warr, ct.c_int))

setargs('qp_bore2radec',
        arg=(qp_memory_t_p, quat_t, arr, quat_t_p,
             warr, warr, warr, warr, ct.c_int))
setargs('qp_bore2radec_hwp',
        arg=(qp_memory_t_p, quat_t, arr, quat_t_p, quat_t_p,
             warr, warr, warr, warr, ct.c_int))
setargs('qp_bore2rasindec',
        arg=(qp_memory_t_p, quat_t, arr, quat_t_p,
             warr, warr, warr, warr, ct.c_int))
setargs('qp_bore2rasindec_hwp',
        arg=(qp_memory_t_p, quat_t, arr, quat_t_p, quat_t_p,
             warr, warr, warr, warr, ct.c_int))
setargs('qp_bore2radecpa',
        arg=(qp_memory_t_p, quat_t, arr, quat_t_p,
             warr, warr, warr, ct.c_int))
setargs('qp_bore2radecpa_hwp',
        arg=(qp_memory_t_p, quat_t, arr, quat_t_p, quat_t_p,
             warr, warr, warr, ct.c_int))

setargs('qp_radecpa2quatn',
        arg=(qp_memory_t_p, arr, arr, arr, wquat_t_p, ct.c_int))
setargs('qp_quat2radecpan',
        arg=(qp_memory_t_p, quat_t_p, warr, warr, warr, ct.c_int))

setargs('qp_radec2pixn',
        arg=(qp_memory_t_p, arr, arr, ct.c_int, warri, ct.c_int))

setargs('qp_radec2gal_quatn', arg=(qp_memory_t_p, wquat_t_p, ct.c_int))
setargs('qp_gal2radec_quatn', arg=(qp_memory_t_p, wquat_t_p, ct.c_int))

setargs('qp_radec2galn',
        arg=(qp_memory_t_p, warr, warr, warr, warr, ct.c_int))
setargs('qp_gal2radecn',
        arg=(qp_memory_t_p, warr, warr, warr, warr, ct.c_int))
setargs('qp_radecpa2galn',
        arg=(qp_memory_t_p, warr, warr, warr, ct.c_int))
setargs('qp_gal2radecpan',
        arg=(qp_memory_t_p, warr, warr, warr, ct.c_int))
setargs('qp_rotate_map',
        arg=(qp_memory_t_p, ct.c_int, arr2, ct.c_char,
             arr2, ct.c_char))

setargs('qp_quat2pixn',
        arg=(qp_memory_t_p, quat_t_p, ct.c_int, warri, warr, warr, ct.c_int))
setargs('qp_quat2pixpan',
        arg=(qp_memory_t_p, quat_t_p, ct.c_int, warri, warr, ct.c_int))
setargs('qp_bore2pix',
        arg=(qp_memory_t_p, quat_t, arr, quat_t_p, ct.c_int,
             warri, warr, warr, ct.c_int))
setargs('qp_bore2pix_hwp',
        arg=(qp_memory_t_p, quat_t, arr, quat_t_p, quat_t_p, ct.c_int,
             warri, warr, warr, ct.c_int))
setargs('qp_bore2pixpa',
        arg=(qp_memory_t_p, quat_t, arr, quat_t_p, ct.c_int,
             warri, warr, ct.c_int))
setargs('qp_bore2pixpa_hwp',
        arg=(qp_memory_t_p, quat_t, arr, quat_t_p, quat_t_p, ct.c_int,
             warri, warr, ct.c_int))

setargs('qp_get_interp_valn',
        arg=(qp_memory_t_p, ct.c_int, arr, arr, arr, warr, ct.c_int))

setargs('qp_set_iers_bulletin_a',
        arg=(qp_memory_t_p, ct.c_int, ct.c_int, arr, arr, arr),
        res=ct.c_int)
setargs('qp_get_iers_bulletin_a',
        arg=(qp_memory_t_p, ct.c_double, ct.POINTER(ct.c_double),
             ct.POINTER(ct.c_double), ct.POINTER(ct.c_double)),
        res=ct.c_int)

def qp_get_bulletin_a(mem, mjd):
    dut1 = ct.c_double()
    x = ct.c_double()
    y = ct.c_double()

    libqp.qp_get_iers_bulletin_a(mem, mjd, ct.byref(dut1),
                              ct.byref(x), ct.byref(y))
    return dut1.value, x.value, y.value

setargs('qp_refraction', arg=(ct.c_double,) * 5, res=ct.c_double)
setargs('qp_update_ref', arg=(qp_memory_t_p, quat_t),
        res=ct.c_double)

def get_vec_mode(map_in=None, pol=True, vpol=False):
    if pol is None:
        pol = True
    if map_in is False:
        return QP_VEC_NONE
    if map_in is None:
        return QP_VEC_VPOL if vpol else QP_VEC_POL if pol else QP_VEC_TEMP
    n = len(map_in)
    if n not in vec_modes:
        raise ValueError('Unrecognized map')
    mode = vec_modes[n]
    if isinstance(mode, dict):
        mode = mode[bool(pol)]
    return mode

def get_proj_mode(proj_in=None, pol=True, vpol=False):
    if pol is None:
        pol = True
    if proj_in is False:
        return QP_PROJ_NONE
    if proj_in is None:
        return QP_PROJ_VPOL if vpol else QP_PROJ_POL if pol else QP_PROJ_TEMP
    n = len(proj_in)
    if n not in proj_modes:
        raise ValueError('Unrecognized proj')
    mode = proj_modes[n]
    if isinstance(mode, dict):
        mode = mode[bool(pol)]
    return mode

# **********************************************************************
# Mapping functions
# **********************************************************************

# initialize detectors
setargs('qp_init_det', arg=(quat_t, ct.c_double, ct.c_double, mueller_t),
        res=qp_det_t_p)
setargs('qp_default_det', res=qp_det_t_p)
setargs('qp_init_det_tod', arg=(qp_det_t_p, ct.c_size_t))
setargs('qp_init_det_tod_from_array',
        arg=(qp_det_t_p, arr, ct.c_size_t, ct.c_int))
setargs('qp_init_det_flag', arg=(qp_det_t_p, ct.c_size_t))
setargs('qp_init_det_flag_from_array',
        arg=(qp_det_t_p, arrf, ct.c_size_t, ct.c_int))
setargs('qp_free_det', arg=qp_det_t_p)
setargs('qp_init_detarr', arg=(quat_t_p, arr, arr, mueller_t_p, ct.c_size_t),
        res=qp_detarr_t_p)
setargs('qp_init_detarr_tod', arg=(qp_detarr_t_p, ct.c_size_t));
setargs('qp_init_detarr_tod_from_array_1d',
        arg=(qp_detarr_t_p, arr, ct.c_size_t, ct.c_int));
setargs('qp_init_detarr_flag', arg=(qp_detarr_t_p, ct.c_size_t));
setargs('qp_init_detarr_flag_from_array_1d',
        arg=(qp_detarr_t_p, arr, ct.c_size_t, ct.c_int));
setargs('qp_free_detarr', arg=qp_detarr_t_p);

# initialize pointing
setargs('qp_init_point', arg=(ct.c_size_t, ct.c_int, ct.c_int),
        res=qp_point_t_p)
setargs('qp_init_point_from_arrays',
        arg=(quat_t_p, arr, quat_t_p, ct.c_size_t, ct.c_int),
        res=qp_point_t_p)
setargs('qp_free_point', arg=qp_point_t_p)

# initialize maps
setargs('qp_init_map', arg=(ct.c_size_t, ct.c_size_t, qp_vec_mode, qp_proj_mode),
        res=qp_map_t_p)
setargs('qp_init_map_from_arrays_1d',
        arg=(arr, arr, ct.c_size_t, ct.c_size_t,
             qp_vec_mode, qp_proj_mode, ct.c_int),
        res=qp_map_t_p)
setargs('qp_init_map_from_map',
        arg=(qp_map_t_p, ct.c_int, ct.c_int), res=qp_map_t_p)
setargs('qp_free_map', arg=qp_map_t_p)
setargs('qp_reshape_map', arg=qp_map_t_p, res=ct.c_int)
setargs('qp_init_map_pixhash',
        arg=(qp_map_t_p, larr, ct.c_size_t), res=ct.c_int)

# tod -> map
setargs('qp_add_map', arg=(qp_map_t_p, qp_map_t_p), res=ct.c_int)
setargs('qp_tod2map1',
        arg=(qp_memory_t_p, qp_det_t_p, qp_point_t_p, qp_map_t_p),
        res=ct.c_int)
setargs('qp_tod2map',
        arg=(qp_memory_t_p, qp_detarr_t_p, qp_point_t_p, qp_map_t_p),
        res=ct.c_int)

# map -> tod
setargs('qp_map2tod1',
        arg=(qp_memory_t_p, qp_det_t_p, qp_point_t_p, qp_map_t_p),
        res=ct.c_int)
setargs('qp_map2tod',
        arg=(qp_memory_t_p, qp_detarr_t_p, qp_point_t_p, qp_map_t_p),
        res=ct.c_int)

# **********************************************************************
# Parameters
# **********************************************************************

def check_set_float(val):
    if not np.isscalar(val):
        raise TypeError('val must be a scalar value')
    return float(val)

def check_set_int(val):
    if not np.isscalar(val):
        raise TypeError('val must be a scalar value')
    return int(val)

def check_pass(val):
    return val

def set_rfunc(state):
    f = libqp['qp_set_rate_%s'%state]
    f.argtypes = (qp_memory_t_p,ct.c_double)
    f.restype = None
    return f
def reset_rfunc(state):
    f = libqp['qp_reset_rate_%s'%state]
    f.argtypes = (qp_memory_t_p,)
    f.restype = None
    return f
def get_rfunc(state):
    f = libqp['qp_get_rate_%s'%state]
    f.argtypes = (qp_memory_t_p,)
    f.restype = ct.c_double
    return f

def check_set_state(state):
    rdict = {'always':QP_DO_ALWAYS, 'once':QP_DO_ONCE, 'never':QP_DO_NEVER}
    state = rdict.get(state,state)
    return check_set_float(state)

def check_get_state(state):
    rdict = {QP_DO_ALWAYS:'always', QP_DO_ONCE:'once', QP_DO_NEVER:'never'}
    state = rdict.get(state,state)
    return state

states = ['lonlat', 'npb', 'erot', 'daber', 'aaber', 'wobble', 'dut1', 'ref']
inv_states = [k + '_inv' for k in states]
state_funcs = dict()
for s in states + inv_states:
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
    f.restype = None
    return f

def get_wfunc(par):
    f = libqp['qp_get_weather_%s'%par]
    f.argtypes = (qp_memory_t_p,)
    f.restype = ct.c_double
    return f

weather_params = ['temperature', 'pressure', 'humidity', 'frequency']
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
    f.restype = None
    return f

def get_ofunc(option):
    f = libqp['qp_get_opt_%s'%option]
    f.argtypes = (qp_memory_t_p,)
    f.restype = ct.c_int
    return f

opts = {'accuracy': {0: 'high',
                     1: 'low'},
        'polconv': {0: ['healpix', 'cosmo'],
                    1: 'iau'},
        'pix_order': {0: 'ring',
                      1: ['nest', 'nested']}}
defaults = {'accuracy': 0,
            'polconv': 0,
            'pix_order': 0}

def check_set_dict(opt):
    def func(val):
        if val is None:
            return defaults[opt]
        if isinstance(val, str):
            val = val.lower()
        for k in opts[opt]:
            v = opts[opt][k]
            if isinstance(v, list):
                if val in v:
                    return k
            else:
                if val == v:
                    return k
            if val == k:
                return k
        return defaults[opt]
    return func

def check_get_dict(opt):
    def func(val):
        for k in opts[opt]:
            v = opts[opt][k]
            if val == k:
                if isinstance(v, list):
                    return v[0]
                return v
        raise ValueError(
            'unrecognized value {} for option {}'.format(val, opt))
    return func

check_set_accuracy = check_set_dict('accuracy')
check_get_accuracy = check_get_dict('accuracy')

check_set_polconv = check_set_dict('polconv')
check_get_polconv = check_get_dict('polconv')

check_set_pix_order = check_set_dict('pix_order')
check_get_pix_order = check_get_dict('pix_order')

def check_get_bool(val):
    return bool(val)

def check_set_bool(val):
    return int(bool(val))

check_set_mean_aber = check_set_bool
check_get_mean_aber = check_get_bool

check_set_fast_math = check_set_bool
check_get_fast_math = check_get_bool

check_set_fast_pix = check_set_bool
check_get_fast_pix = check_get_bool

check_set_error_missing = check_set_bool
check_get_error_missing = check_get_bool

check_set_nan_missing = check_set_bool
check_get_nan_missing = check_get_bool

check_set_interp_missing = check_set_bool
check_get_interp_missing = check_get_bool

check_set_interp_pix = check_set_bool
check_get_interp_pix = check_get_bool

def check_set_num_threads(nt):
    if nt is None:
        return 0
    return int(nt)

def check_get_num_threads(nt):
    return nt

def check_set_thread_num(tn):
    if tn is None:
        return 0
    return int(nt)

def check_get_thread_num(tn):
    return tn

options = ['accuracy', 'mean_aber', 'fast_math', 'polconv', 'pix_order',
           'interp_pix', 'fast_pix', 'error_missing', 'nan_missing',
           'interp_missing', 'num_threads', 'thread_num']
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
    f.restype = None
    return f

def get_pfunc(par):
    f = libqp['qp_get_%s'%par]
    f.argtypes = (qp_memory_t_p,)
    f.restype = ct.c_double
    return f

double_params = ['ref_delta','dut1']
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

# **********************************************************************
# Argument checking
# **********************************************************************

def check_flags(arg):
    return (arg.flags['C_CONTIGUOUS'] |
            arg.flags['OWNDATA'] << 1 |
            arg.flags['WRITEABLE'] << 2 |
            arg.flags['ALIGNED'] << 3)

def check_input(name, arg, shape=None, quat=False, dtype=np.double,
                inplace=True, fill=0, allow_transpose=True,
                allow_tuple=True, output=False):
    """
    Ensure input argument is an aligned array of the right type and shape.

    Arguments
    ---------
    name : string
        Name of the argument
    arg : array_like
        The argument itself
    shape : tuple, optional
        If supplied, ensure `arg` has this shape.
    quat : bool, optional
        If True, ensure that the last dimension of the array has length 4.
    dtype : numpy.dtype, optional
        Ensure `arg` is of this dtype.  Default: numpy.double.
    inplace : bool, optional
        If False make sure that a copy of the input `arg` is made prior to
        returning.  Otherwise, the output `arg` may share memory with the
        input.
    fill : scalar, optional
        If the input is missing or empty, fill with this value.
        If None, leave the array empty.
    allow_transpose : bool, optional
        If True, transpose the input array if the transposed shape matches
        `shape`.
    allow_tuple : bool, optional
        If True, `numpy.vstack` the input `arg` if it is a tuple.
    output : bool, optional
        If True, ensure that the output `arg` is a writeable array.
        Otherwise, the array is only ensured to be aligned and
        C-contiguous.

    Returns
    -------
    arg :
        Aligned, contiguous and properly typed and shaped array for
        passing to the C library.
    """
    if arg is None:
        if shape is None:
            raise ValueError('need shape to initialize input!')
        if fill is None:
            arg = np.empty(shape, dtype=dtype)
        else:
            arg = fill * np.ones(shape, dtype=dtype)
    if isinstance(arg, tuple) and allow_tuple:
        arg = np.vstack(arg)
    if np.isscalar(arg):
        arg = np.array(arg)
    if not isinstance(arg, np.ndarray):
        raise TypeError('input {} must be of type numpy.ndarray'.format(name))
    if quat and arg.shape[-1] != 4:
        raise ValueError('input {} is not a valid quaternion array')
    if shape is not None:
        if arg.shape != shape:
            if arg.T.shape == shape and allow_transpose:
                arg = arg.T
            else:
                try:
                    _, arg = np.broadcast_arrays(np.empty(shape), arg)
                    # ensure writeable flag is set properly if necessary
                    # broadcast arrays will not be writable in the future
                    if inplace or output:
                        arg = arg.copy()
                    else:
                        arg.flags['WRITEABLE'] = False
                except ValueError:
                    s = 'input {} of shape {} cannot be broadcast to shape {}'
                    raise ValueError(s.format(name, arg.shape, shape))
    istat = check_flags(arg)
    arg = np.require(arg, dtype, list('AC' + 'W'*output))
    ostat = check_flags(arg)
    if istat == ostat and inplace is False:
        return arg.copy()
    return arg

def check_inputs(*args, **kwargs):
    """
    Ensure that a group of input arguments have the same shape by broadcasting.

    Arguments
    ---------
    args :
        A list of broadcastable array_like arguments.
    kwargs :
        Dictionary of arguments to pass to `check_input`.

    Returns
    -------
    args :
        A list of broadcast and properly aligned/shaped/typed arrays for
        passing to the C library as a set of timestreams.
    """
    args = [arg if arg is not None else 0 for arg in args]
    if 'shape' not in kwargs:
        kwargs['shape'] = np.broadcast(*args).shape
    return [
        check_input('argument {}'.format(ix), np.atleast_1d(x), **kwargs)
        for ix, x in enumerate(args)
    ]

def check_output(name, arg=None, shape=None, quat=False, dtype=np.double,
                 inplace=True, fill=None, allow_transpose=True,
                 allow_tuple=True, **kwargs):
    """
    Ensure that the output argument is properly aligned/shaped/typed.
    Check input kwargs to see if a pointer to the output array
    has been passed in.

    Arguments
    ---------
    name : string
        Name of the argument
    arg : array_like, optional
        The argument itself.  If not supplied, kwargs is checked.
    shape : tuple, optional
        If supplied, ensure `arg` has this shape.  If `arg` is None,
        this argument is required.
    dtype : numpy.dtype, optional
        Ensure `arg` is of this dtype.  Default: numpy.double.
    inplace : bool, optional
        If False make sure that a copy of the input `arg` is made prior to
        returning.  Otherwise, the output `arg` may share memory with the
        input.
    fill : scalar, optional
        If the input is missing or empty, fill with this value.
        If None, leave the array empty.
    allow_transpose : bool, optional
        If True, transpose the input array if the transposed shape matches
        `shape`.
    allow_tuple : bool, optional
        If True, `numpy.vstack` the input `arg` if it is a tuple.
    kwargs :
        Any remaining input arguments.  If `arg` is None,
        `kwargs` is searched for the `name` key.  If not found, a
        an empty array is created of the appropriate shape.

    Returns
    -------
    arg :
        Aligned, contiguous, writeable and properly typed and shaped array
        for passing to the C library.
    """
    if arg is None:
        arg = kwargs.pop(name, None)
    return check_input(name, arg, shape, quat, dtype, inplace, fill,
                       allow_transpose, allow_tuple, True)
