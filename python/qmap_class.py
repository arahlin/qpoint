import numpy as _np
from qpoint_class import QPoint

class QMap(QPoint):

    def __init__(self, *args, **kwargs):
        super(QMap, self).__init__(*args, **kwargs)

    def _check_map(self, m):
        as_tuple = False
        if isinstance(m, tuple):
            as_tuple = True
            m = _np.vstack(m)
        if not isinstance(m, _np.ndarray):
            raise TypeError, 'map must be of type numpy.ndarray'
        shape_in = m.shape
        if m.ndim == 2 and _np.argmin(shape_in) == 0:
            m = m.T
        return m, shape_in, as_tuple

    def _return_map(self, m, shape=None, as_tuple=True):
        if shape is not None:
            if m.T.shape == shape:
                m = m.T
            if m.shape != shape:
                raise ValueError,'map shape has changed!'
        if as_tuple:
            if m.ndim == 2 and _np.argmax(m.shape) == 0:
                m = m.T
        return m

    def rotate_map(self, map_in, coord=['C','G'], map_out=None,
                   tuple_out=None, **kwargs):
        """
        Rotate a polarized npix-x-3 map from one coordinate system to another.
        Supported coordinates:

        C = celestial (J2000)
        G = galactic
        """

        from warnings import warn
        warn('This code is buggy, use at your own risk', UserWarning)

        map_in, shape_in, tuple_in = self._check_map(map_in)
        if tuple_out is None:
            tuple_out = tuple_in

        nside = int(_np.sqrt(max(map_in.shape) / 12))
        shape = (12*nside*nside, 3)
        map_in = self._check_input('map_in', map_in, shape=shape)
        map_out = self._check_output('map_out', map_out, shape=shape, fill=0)

        try:
            coord_in = coord[0]
            coord_out = coord[1]
        except:
            raise ValueError,'unable to parse coord'

        _libqp.qp_rotate_map(self._memory, nside, map_in, coord_in, map_out, coord_out)
        return self._return_map(map_out, shape_in, tuple_out)

    def tod2map(self, q_off, ctime, q_bore, nside=256, pol=True,
                q_hwp=None, tod=None, smap=None, pmap=None,
                tuple_out=None, **kwargs):
        """
        Calculate signal and hits maps for given detectors and boresight orientations.

        Arguments
        ---------
        q_off : array_like
          quaternion offset array, of shape ndex-x-4
        ctime : array_like
          UTC time in seconds
        q_bore : array_like
          boresight quaternion timestream, of shape nsamp-x-4

        Keyword Arguments
        -----------------
        nside : int, optional
          healpix map size.  Determined from output map, if given.  Default: 256.
        pol : bool, optional
          write to polarized map and/or pointing matrix if True.  Determined
          from output maps if given.  Default: True.
        q_hwp : array_like, optional
          HWP quaternion timestream, same shape as q_bore
        tod : array_like, optional
          output array for timestreams, of shape ndet-x-nsamp
          if not supplied, smap is not computed.
        smap : array_like or bool, optional
          output signal map, of shape npix-x-3 (if pol) or npix (if not pol).
          use this argument for in-place computation. if False, do not compute.
        pmap : array_like or bool, optional
          output pointing matrix map, of shape npix-x-6 (if pol) or npix (if not pol).
          use this argument for in-place computation. if False, do not compute.
        tuple_out : bool, optional
          if True, return the map as a tuple if npix-sized arrays.
          set this to False if these maps will be added to my more processes.

        The remaining keyword arguments are passed to the QPoint.set method.

        Returns
        -------
        smap : array_like
           accumulated signal map, if requested
        pmap: array_like
           accumulated pointing matrix map, if requested
        """

        self.set(**kwargs)

        q_off = self._check_input('q_off', q_off)
        ndet = q_off.size/4

        q_bore = self._check_input('q_bore', q_bore)
        if ctime is None:
            if not self.get('mean_aber'):
                raise ValueError,'ctime required if mean_aber is False'
            ctime = _np.zeros((q_bore.size/4,), dtype=q_bore.dtype)
        ctime  = self._check_input('ctime', ctime)
        n = ctime.size

        do_pnt = not (pmap is False)
        do_sig = not (tod is None or tod is False or smap is False)

        if not (do_pnt or do_sig):
            raise KeyError, 'Either smap or pmap must not be False'

        npix = 12 * nside * nside
        mshapep = (npix, 6) if pol else (npix,)

        if do_pnt:
            if pmap is None:
                pmap = _np.zeros(mshapep, dtype=_np.double)
            pmap, pshape_in, ptuple_in = self._check_map(pmap)
            npix = max(pmap.shape)
            pol = True if min(pmap.shape) == 6 else False
            nside = int(_np.sqrt(npix/12))
            mshapep = (npix, 6) if pol else (npix,)
            pmap = self._check_output('pmap', pmap, shape=mshapep)

        mshapes = (npix, 3) if pol else (npix,)

        if do_sig:
            tod = self._check_input('tod', tod, shape=(ndet, n))
            todp = (tod.__array_interface__['data'][0] +
                    _np.arange(tod.shape[0]) * tod.strides[0]).astype(_np.uintp)

            if smap is None:
                smap = _np.zeros(mshapes, dtype=_np.double)
            smap, sshape_in, stuple_in = self._check_map(smap)
            npix = max(smap.shape)
            pol = True if min(smap.shape) == 3 else False
            nside = int(_np.sqrt(npix/12))
            mshapes = (npix, 3) if pol else (npix,)
            smap = self._check_output('smap', smap, shape=mshapes)

            if do_pnt:
                if mshapes[0] != mshapep[0]:
                    raise ValueError,'smap and pmap have incompatible shapes'
                if len(mshapes) == 2 and len(mshapep) == 2:
                    if (mshapes[1] > 1) != (mshapep[1] > 1):
                        raise ValueError,'smap and pmap have incompatible shapes'

        if tuple_out is None:
            if do_pnt:
                ptuple_out = ptuple_in
            if do_sig:
                stuple_out = stuple_in
        else:
            if do_pnt:
                ptuple_out = tuple_out
            if do_sig:
                stuple_out = tuple_out

        if pol is False:
            if do_pnt and do_sig:
                _libqp.qp_tod2map_sigpnt_nopol(self._memory, q_off, ndet, ctime, q_bore,
                                               todp, n, smap, pmap, nside)
            elif do_sig:
                _libqp.qp_tod2map_sig_nopol(self._memory, q_off, ndet, ctime, q_bore,
                                            todp, n, smap, nside)
            elif do_pnt:
                _libqp.qp_tod2map_pnt_nopol(self,_memory, q_off, ndet, ctime, q_bore,
                                            n, pmap, nside)
        elif q_hwp is None:
            if do_pnt and do_sig:
                _libqp.qp_tod2map_sigpnt(self._memory, q_off, ndet, ctime, q_bore,
                                         todp, n, smap, pmap, nside)
            elif do_sig:
                _libqp.qp_tod2map_sig(self._memory, q_off, ndet, ctime, q_bore,
                                      todp, n, smap, nside)
            elif do_pnt:
                _libqp.qp_tod2map_pnt(self._memory, q_off, ndet, ctime, q_bore,
                                      n, pmap, nside)
        else:
            q_hwp = self._check_input('q_hwp', q_hwp, shape=q_bore.shape)

            if do_pnt and do_sig:
                _libqp.qp_tod2map_sigpnt_hwp(self._memory, q_off, ndet, ctime, q_bore,
                                             q_hwp, todp, n, smap, pmap, nside)
            elif do_sig:
                _libqp.qp_tod2map_sig_hwp(self._memory, q_off, ndet, ctime, q_bore,
                                          q_hwp, todp, n, smap, nside)
            elif do_pnt:
                _libqp.qp_tod2map_pnt_hwp(self._memory, q_off, ndet, ctime, q_bore,
                                          q_hwp, n, pmap, nside)

        if do_pnt and do_sig:
            return (self._return_map(smap, sshape_in, stuple_out),
                    self._return_map(pmap, pshape_in, ptuple_out))
        elif do_sig:
            return self._return_map(smap, sshape_in, stuple_out)
        elif do_pnt:
            return self._return_map(pmap, pshape_in, ptuple_out)

    def bore2map(self, *args, **kwargs):
        """
        *** This method is deprecated, use QPoint.tod2map ***
        """
        from warnings import warn
        warn("Use of QPoint.bore2map is deprecated, use QPoint.tod2map instead",
             DeprecationWarning)
        return self.tod2map(*args, **kwargs)
    bore2map.__doc__ += tod2map.__doc__

    def map2tod(self, q_off, ctime, q_bore, map_in, q_hwp=None, tod=None,
                pol=True, **kwargs):
        """
        Calculate signal TOD from input map given detector offsets,
        boresight orientation and hwp orientation (optional).

        Arguments
        ---------
        q_off : array_like
          quaternion offset array, of shape ndex-x-4
        ctime : array_like
          UTC time in seconds
        q_bore : array_like
          boresight quaternion timestream, of shape nsamp-x-4
        map_in : array_like
          input map of shape npix-x-N.  See notes for allowed values of N

        Keyword Arguments
        -----------------
        pol : bool, optional
          if True, input map is polarized.  this argument is only used if
          the shape of map_in is npix-x-3.  Default: True
        q_hwp : array_like, optional
          HWP quaternion timestream, same shape as q_bore
        tod : array_like, optional
          output array for timestreams, of shape ndet-x-nsamp
          use this keyword argument for in-place computation.

        The remaining keyword arguments are passed to the QPoint.set method.

        Returns
        -------
        tod : array_like
          A timestream sampled from the input map for each requested detector.
          The output array shape is ndet-x-nsamp.

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

        self.set(**kwargs)

        q_off = self._check_input('q_off', q_off)
        ndet = q_off.size/4

        q_bore = self._check_input('q_bore', q_bore)
        if ctime is None:
            if not self.get('mean_aber'):
                raise ValueError,'ctime required if mean_aber is False'
            ctime = _np.zeros((q_bore.size/4,), dtype=q_bore.dtype)
        ctime  = self._check_input('ctime', ctime)
        n = ctime.size

        map_in, _, _ = self._check_map(map_in)
        npix = max(map_in.shape)
        nside = _np.sqrt(npix / 12)
        if _np.floor(_np.log2(nside)) != _np.log2(nside):
            raise ValueError,'invalid nside'
        nside = int(nside)
        ncol = map_in.size / npix
        mshape = (npix, ncol) if ncol > 1 else (npix,)
        map_in = self._check_input('map_in', map_in, shape=mshape)

        if ncol not in [1,3,6,9,18]:
            raise ValueError,'map must have 1, 3, 6, 9 or 18 columns'

        tod = self._check_output('tod', tod, shape=(ndet, n))
        # array of pointers
        todp = (tod.__array_interface__['data'][0] +
                _np.arange(tod.shape[0]) * tod.strides[0]).astype(_np.uintp)

        if q_hwp is not None and ncol in [3,9,18]:
            q_hwp = self._check_input('q_hwp', q_hwp, shape=q_bore.shape)

        if ncol == 1:
            _libqp.qp_map2tod_nopol(self._memory, q_off, ndet, ctime, q_bore,
                                    map_in, nside, todp, n)
        elif ncol == 3:
            if not pol:
                _libqp.qp_map2tod_der1_nopol(self._memory, q_off, ndet, ctime, q_bore,
                                             map_in, nside, todp, n)
            elif q_hwp is None:
                _libqp.qp_map2tod(self._memory, q_off, ndet, ctime, q_bore,
                                  map_in, nside, todp, n)
            else:
                _libqp.qp_map2tod_hwp(self._memory, q_off, ndet, ctime, q_bore, q_hwp,
                                      map_in, nside, todp, n)
        elif ncol == 6:
            _libqp.qp_map2tod_der2_nopol(self._memory, q_off, ndet, ctime, q_bore,
                                         map_in, nside, todp, n)
        elif ncol == 9:
            if q_hwp is None:
                _libqp.qp_map2tod_der1(self._memory, q_off, ndet, ctime, q_bore,
                                       map_in, nside, todp, n)
            else:
                _libqp.qp_map2tod_der1_hwp(self._memory, q_off, ndet, ctime, q_bore,
                                           q_hwp, map_in, nside, todp, n)
        elif ncol == 18:
            if q_hwp is None:
                _libqp.qp_map2tod_der2(self._memory, q_off, ndet, ctime, q_bore,
                                       map_in, nside, todp, n)
            else:
                _libqp.qp_map2tod_der2_hwp(self._memory, q_off, ndet, ctime, q_bore,
                                           q_hwp, map_in, nside, todp, n)

        return tod
