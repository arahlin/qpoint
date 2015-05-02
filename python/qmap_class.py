import numpy as np
from qpoint_class import QPoint
import healpy as hp
import ctypes as ct

import _libqpoint as lib
from _libqpoint import libqp as qp

__all__ = ['QMap']

def check_map(map_in):
    """
    Return a properly transposed and memory-aligned map and its nside.
    """
    map_in = np.atleast_2d(map_in)
    if np.argmax(map_in.shape) == 0:
        map_in = map_in.T
    nside = hp.get_nside(map_in)
    return lib.check_input('map', map_in), nside

class QMap(QPoint):
    """
    Quaternion-based mapmaker that generates per-channel pointing
    on-the-fly.

    Input and output maps are initialized once
    """

    def __init__(self, nside=256, pol=True,
                 source_map=None, source_pol=True,
                 q_bore=None, ctime=None, q_hwp=None, **kwargs):
        """
        Initialize the internal structures and data depo for
        mapmaking and/or timestream generation.

        Arguments
        ---------
        nside : int, optional
            Resolution of the output maps. Default: 256
        pol : bool, optional
            If True, output maps are polarized.
        source_map : array_like, optional
            If supplied, passed to `init_source()` to initialize
            the source map structure.
        source_pol : bool, optional
            If True, source_map is polarized.  See `init_source()`
            for details.
        q_bore, ctime, q_hwp : array_like, optional
            Boresight pointing data.  See `init_point()` for details.
            If not supplied, the pointing structure is left
            uninitialized.

        Remaining arguments are passed to `QPoint.set()`.

        Notes
        -----
        An internal `depo` dictionary attribute stores the source and output
        maps, timestreams, and pointing data for retrieval by the user.
        Only pointers to these arrays in memory are passed to the C
        library.  To ensure that extraneous copies of data are not made,
        supply these methods with C-contiguous arrays of the correct shape.
        """
        super(QMap, self).__init__(**kwargs)

        self.reset()

        self.init_dest(nside=nside, pol=pol)

        if source_map is not None:
            self.init_source(source_map, pol=source_pol)

        if q_bore is not None:
            self.init_point(q_bore, ctime=ctime, q_hwp=q_hwp)

    def init_source(self, source_map, pol=True, reset=False):
        """
        Initialize the source map structure.  Timestreams are
        produced by scanning this map.

        Arguments
        ---------
        source_map : array_like
            Input map.  Must be of shape (N, npix), where N can be
            1, 3, 6, 9, or 18.
        pol : bool, optional
            If True, and the map shape is (3, npix), then input is a
            polarized map (and not T + 1st derivatives).
        reset : bool, optional
            If True, and if the structure has already been initialized,
            it is reset and re-initialized with the new map.  If False,
            a RuntimeError is raised if the structure has already been
            initialized.

        Notes
        -----
        This method will automatically determine the type of map
        given its shape.  Note that for N=3, the pol keyword argument
        should be used to disambiguate the two map types.  By default,
        a polarized map with (T,Q,U) components is assumed.

        N     map_in contents
        1  :  T
        3  :  (T, Q, U) --or-- (T, dTdt, dTdp)
        6  :  (T, dTdt, dTdp) + (dT2dt2, dT2dpdt, dT2dp2)
        9  :  (T, Q, U) + (dTdt, dQdt, dUdt, dTdp, dQdp, dUdp)
        18 :  (N=9) + (dT2dt2, dQ2dt2, dU2dt2, dT2dpdt, dQ2dpdt, dU2dpdt,
                       dT2dp2, dQ2dp2, dU2dp2)
        """

        if self._source.contents.init:
            if reset:
                self.reset_source()
            else:
                raise RuntimeError, 'source already initialized'

        # check map shape and create pointer
        smap, nside = check_map(source_map)

        # store map
        self.depo['source_map'] = smap

        # initialize
        source = self._source.contents
        source.nside = nside
        source.npix = hp.nside2npix(nside)
        source.num_vec = len(source_map)
        source.vec_mode = lib.get_vec_mode(smap, pol)
        source.vec1d = lib.as_ctypes(smap.ravel())
        source.vec1d_init = lib.QP_ARR_INIT_PTR
        source.vec = None
        source.vec_init = 0
        source.num_proj = 0
        source.proj_mode = 0
        source.proj = None
        source.proj_init = 0
        source.proj1d = None
        source.proj1d_init = 0
        source.init = lib.QP_STRUCT_INIT

        if qp.qp_reshape_map(self._source):
            raise RuntimeError, qp.qp_get_error_string(self._memory)

    def reset_source(self):
        """
        Reset the source map structure.  Must be reinitialized to
        produce more timestreams.
        """
        if hasattr(self, '_source'):
            qp.qp_free_map(self._source)
        self.depo.pop('source_map', None)
        self._source = ct.pointer(lib.qp_map_t())

    def init_dest(self, nside=256, pol=True, vec=None, proj=None, copy=False,
                  reset=False):
        """
        Initialize the destination map structure.  Timestreams are binned
        and projection matrices accumulated into this structure.

        Arguments
        ---------
        nside : int, optional
            map dimension
        pol : bool, optional
            If True, a polarized map will be created.
        vec : array_like or bool, optional
            If supplied, nside and pol are determined from this map, and
            the vector (binned signal) map is initialized from this.
            If False, mapping from timestreams is disabled.
        proj : array_like or bool, optional
            If supplied, the projection matrix map is initialized from this.
            If False, accumulation of the projection matrix is disabled.
        copy : bool, optional
            If True and vec/proj are supplied, make copies of these inputs
            to avoid in-place operations.
        reset : bool, optional
            If True, and if the structure has already been initialized,
            it is reset and re-initialized with the new map.  If False,
            a RuntimeError is raised if the structure has already been
            initialized.
        """

        if vec is False and proj is False:
            raise ValueError('one of vec or proj must not be False')

        if self._dest.contents.init:
            if reset:
                self.reset_dest()
            else:
                raise RuntimeError,'dest already initialized'

        npix = hp.nside2npix(nside)

        if vec is None:
            if pol:
                vec = np.zeros((3, npix), dtype=np.double)
            else:
                vec = np.zeros((1, npix), dtype=np.double)
        elif vec is not False:
            veci = vec
            vec, nside = check_map(veci)
            npix = hp.nside2npix(nside)
            if len(vec) == 1:
                pol = False
            elif len(vec) == 3:
                pol = True
            else:
                raise ValueError('vec has incompatible shape')
            if copy and np.may_share_memory(veci, vec):
                vec = vec.copy()

        if proj is None:
            if pol:
                proj = np.zeros((6, npix), dtype=np.double)
            else:
                proj = np.zeros((1, npix), dtype=np.double)
        elif proj is not False:
            proji = proj
            proj, pnside = check_map(proji)

            if vec is not False:
                if len(proj) != [1,6][pol]:
                    raise ValueError('proj has incompatible shape')
                if pnside != nside:
                    raise ValueError('proj has incompatible nside')
            else:
                if len(proj) == 1:
                    pol = False
                elif len(proj) == 6:
                    pol = True
                else:
                    raise ValueError('proj has incompatible shape')
                nside = pnside
                npix = hp.nside2npix(nside)

            if copy and np.may_share_memory(proji, proj):
                proj = proj.copy()

        # store arrays for later retrieval
        self.depo['vec'] = vec
        self.depo['proj'] = proj

        # initialize
        dest = self._dest.contents
        dest.nside = nside
        dest.npix = npix
        if vec is not False:
            dest.num_vec = len(vec)
            dest.vec_mode = lib.get_vec_mode(vec, pol)
            dest.vec1d = lib.as_ctypes(vec.ravel())
            dest.vec1d_init = lib.QP_ARR_INIT_PTR
        if proj is not False:
            dest.num_proj = len(proj)
            dest.proj_mode = lib.get_proj_mode(proj, pol)
            dest.proj1d = lib.as_ctypes(proj.ravel())
            dest.proj1d_init = lib.QP_ARR_INIT_PTR
        dest.vec = None
        dest.vec_init = 0
        dest.proj = None
        dest.proj_init = 0
        dest.init = lib.QP_STRUCT_INIT

        if qp.qp_reshape_map(self._dest):
            raise RuntimeError, qp.qp_get_error_string(self._memory)

    def reset_dest(self):
        """
        Reset the destination map structure.
        Must be reinitialized to continue mapmaking.
        """
        if hasattr(self, '_dest'):
            qp.qp_free_map(self._dest)
        self.depo.pop('vec', None)
        self.depo.pop('proj', None)
        self._dest = ct.pointer(lib.qp_map_t())

    def init_point(self, q_bore=None, ctime=None, q_hwp=None):
        """
        Initialize or update the boresight pointing data structure.

        Arguments
        ---------
        q_bore : array_like, optional
             Boresight pointing quaternion, of shape (nsamp, 4).
             If supplied, the pointing structure is reset if already
             initialized.
        ctime : array_like, optional
             time since the UTC epoch.  If not None, the time array
             is updated to this. Shape must be (nsamp,)
        q_hwp : array_like, optional
             Waveplate quaternion.  If not None, the quaternion is
             updated to this. Shape must be (nsamp, 4)
        """

        point = self._point.contents

        if q_bore is not None:
            if point.init:
                self.reset_point()
            point = self._point.contents
            q_bore = lib.check_input('q_bore', q_bore)
            n = q_bore.size / 4
            point.n = n
            self.depo['q_bore'] = q_bore
            point.q_bore = lib.as_ctypes(q_bore)
            point.q_bore_init = lib.QP_ARR_INIT_PTR
            point.init = lib.QP_STRUCT_INIT

        if not point.init:
            raise RuntimeError, 'point not initialized'

        n = point.n

        if ctime is False:
            point.ctime_init = 0
            point.ctime = None
        elif ctime is not None:
            if not point.init:
                raise RuntimeError, 'point not initialized'
            ctime = lib.check_input('ctime', ctime, shape=(n,))
            self.depo['ctime'] = ctime
            point.ctime_init = lib.QP_ARR_INIT_PTR
            point.ctime = lib.as_ctypes(ctime)

        if q_hwp is False:
            point.q_hwp_init = 0
            point.q_hwp = None
        elif q_hwp is not None:
            if not point.init:
                raise RuntimeError, 'point not initialized'
            q_hwp = lib.check_input('q_hwp', q_hwp, shape=(n, 4))
            self.depo['q_hwp'] = q_hwp
            point.q_hwp_init = lib.QP_ARR_INIT_PTR
            point.q_hwp = lib.as_ctypes(q_hwp)

    def reset_point(self):
        """
        Reset the pointing data structure.
        """
        if hasattr(self, '_point'):
            qp.qp_free_point(self._point)
        self.depo.pop('q_bore', None)
        self.depo.pop('ctime', None)
        self.depo.pop('q_hwp', None)
        self._point = ct.pointer(lib.qp_point_t())

    def init_detarr(self, q_off, weight=None, pol_eff=None, tod=None, flag=None,
                    write=False):
        """
        Initialize the detector listing structure.  Detector properties and
        timestreams are passed to and from the mapmaker through this structure.

        Arguments
        ---------
        q_off : array_like
            Array of offset quaternions, of shape (ndet, 4).
        weight : array_like, optional
            Per-channel mapmaking weights, of shape (ndet,) or a constant.
            Default : 1.
        pol_eff : array_like, optional
            Per-channel polarization efficiencies, of shape(ndet,) or a constant.
            Default: 1.
        tod : array_like, optional
            Timestream array, of shape (ndet, nsamp).  nsamp must match that of
            the pointing structure.  If not supplied and `write` is True, then
            a zero-filled timestream array is initialized.
        flag : array_like, optional
            Flag array, of shape (ndet, nsamp), for excluding data from
            mapmaking.  nsamp must match that of the pointing structure.
            If not supplied, a zero-filled array is initialized (i.e. no 
            flagged samples).
        write : bool, optional
             If True, the timestreams are ensured writable and created if
             necessary.
        """

        self.reset_detarr()

        # check inputs
        q_off = lib.check_input('q_off', q_off)
        n = q_off.size / 4
        weight = lib.check_input('weight', weight, shape=(n,), fill=1)
        pol_eff = lib.check_input('pol_eff', pol_eff, shape=(n,), fill=1)

        ns = self._point.contents.n
        shape = (n, ns)

        if write:
            tod = lib.check_output('tod', tod, shape=shape, fill=0)
            self.depo['tod'] = tod
        elif tod is not None:
            tod = lib.check_input('tod', tod, shape=shape)
            self.depo['tod'] = tod
        if flag is not None:
            flag = lib.check_input('flag', flag, dtype=np.uint8, shape=shape)
            self.depo['flag'] = flag

        # populate array
        dets = (lib.qp_det_t * n)()
        for idx, (q, w, p) in enumerate(zip(q_off, weight, pol_eff)):
            dets[idx].init = lib.QP_STRUCT_INIT
            dets[idx].q_off = lib.as_ctypes(q)
            dets[idx].weight = w
            dets[idx].pol_eff = p
            if tod is not None:
                dets[idx].n = ns
                dets[idx].tod_init = lib.QP_ARR_INIT_PTR
                dets[idx].tod = lib.as_ctypes(tod[idx])
            else:
                dets[idx].tod_init = 0
            if flag is not None:
                dets[idx].n = ns
                dets[idx].flag_init = lib.QP_ARR_INIT_PTR
                dets[idx].flag = lib.as_ctypes(flag[idx])
            else:
                dets[idx].flag_init = 0

        detarr = lib.qp_detarr_t()
        detarr.n = n
        detarr.init = lib.QP_STRUCT_INIT
        detarr.arr_init = lib.QP_ARR_INIT_PTR
        detarr.arr = dets

        self._detarr = ct.pointer(detarr)

    def reset_detarr(self):
        """
        Reset the detector array structure.
        """
        if hasattr(self, '_detarr') and self._detarr is not None:
            qp.qp_free_detarr(self._detarr)
        self._detarr = None
        self.depo.pop('tod', None)
        self.depo.pop('flag', None)

    def reset(self):
        """
        Reset the internal data structures, and clear the data depo.
        """
        if not hasattr(self, 'depo'):
            self.depo = dict()
        self.reset_source()
        self.reset_dest()
        self.reset_point()
        self.reset_detarr()
        for k in self.depo:
            self.depo.pop(k)

    def from_tod(self, q_off, tod=None, count_hits=True, weight=None,
                 pol_eff=None, flag=None, **kwargs):
        """
        Calculate signal and hits maps for given detectors.

        Arguments
        ---------
        q_off : array_like
          quaternion offset array, of shape (ndet, 4)
        tod : array_like, optional
          output array for timestreams, of shape (ndet, nsamp)
          if not supplied, only the projection map is populated.
        count_hits : bool, optional
          if True (default), populate projection map.
        weight : array_like, optional
          array of channel weights, of shape (ndet,).  Defaults to 1 if not
          supplied.
        pol_eff : array_like, optional
          array of polarization efficiencies, of shape (ndet,).  Defaults to 1
          if not supplied.
        flag : array_like, optional
          array of flag timestreams for each channel, of shape (ndet, nsamp).

        The remaining keyword arguments are passed to the QPoint.set method.

        Returns
        -------
        vec : array_like, optional
            binned signal map, if tod is supplied
        proj : array_like, optional
            binned projection matrix map, if count_hits is True
        """

        self.set(**kwargs)

        # initialize detectors
        self.init_detarr(q_off, weight=weight, pol_eff=pol_eff,
                         tod=tod, flag=flag)

        # check and cache modes
        dest = self._dest.contents
        return_map = True
        vec_mode = dest.vec_mode
        if tod is None or tod is False:
            return_map = False
            dest.vec_mode = 0
        return_proj = True
        proj_mode = dest.proj_mode
        if not count_hits:
            if tod is None:
                raise ValueError, 'must supply either tod or count_hits'
            return_proj = False
            dest.proj_mode = 0

        # run
        if qp.qp_tod2map(self._memory, self._detarr, self._point, self._dest):
            raise RuntimeError, qp.qp_get_error_string(self._memory)

        # reset modes
        dest.vec_mode = vec_mode
        dest.proj_mode = proj_mode

        # clean up
        self.reset_detarr()

        # return
        ret = ((self.depo['vec'],) * return_map +
               (self.depo['proj'],) * return_proj)
        if len(ret) == 1:
            return ret[0]
        return ret

    def to_tod(self, q_off, pol_eff=None, tod=None, **kwargs):
        """
        Calculate signal TOD from source map for multiple channels.

        Arguments
        ---------
        q_off : array_like
          quaternion offset array, of shape (ndet, 4)
        pol_eff : array_like, optional
          array of polarization efficiencies, of shape (ndet,).  Defaults to 1
          if not supplied.
        tod : array_like, optional
          output array for timestreams, of shape (ndet, nsamp)
          use this keyword argument for in-place computation.

        The remaining keyword arguments are passed to the QPoint.set method.

        Returns
        -------
        tod : array_like
          A timestream sampled from the input map for each requested detector.
          The output array shape is (ndet, nsamp).
        """

        self.set(**kwargs)

        # initialize detectors
        self.init_detarr(q_off, pol_eff=pol_eff, tod=tod, write=True)

        # run
        if qp.qp_map2tod(self._memory, self._detarr, self._point, self._source):
            raise RuntimeError, qp.qp_get_error_string(self._memory)
        tod = self.depo.pop('tod')

        # clean up
        self.reset_detarr()

        # return
        return tod
