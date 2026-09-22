import numpy as np
from .qpoint_class import QPoint, qp_settings
from . import _libqpoint2 as lib

__all__ = ["QMap", "check_map", "check_proj"]


# ---------------------------------------------------------------------------
# HEALPix helpers
# ---------------------------------------------------------------------------


def nside2npix(nside):
    """
    Number of pixels in a full-sky map of this resolution.

    Arguments
    ---------
    nside : int
        HEALPix resolution parameter.

    Returns
    -------
    npix : int
        Number of pixels in a full-sky map at that resolution.
    """
    return 12 * nside * nside


def npix2nside(npix):
    """
    Resolution parameter for a full-sky map of this length.

    Arguments
    ---------
    npix : int
        Number of pixels in a full-sky map.

    Returns
    -------
    nside : int
        The resolution parameter it corresponds to. Raises ValueError if
        `npix` is not 12 * nside ** 2 for an integer power of two.
    """
    nside = np.sqrt(npix / 12.0)
    if nside != np.floor(nside):
        raise ValueError("Invalid npix, must be 12 * nside**2")
    return int(nside)


def check_map(map_in, copy=False, partial=False, dtype=np.double):
    """
    Return a properly shaped, C-contiguous map and its nside (or npix if partial).

    A map is read as (nrow, npix) as given and never transposed: guessing
    read a (3, 2) partial map of two pixels as three bogus pixels, and
    quietly accepted (npix, nrow) input that nothing else in the API
    produces. Nothing here tests the aspect ratio in its place -- a
    transposed map fails on npix, or on the row counts the modes allow, or
    against the length of the pixel list for a partial map.

    Arguments
    ---------
    map_in : map or list of maps
        Input map(s)
    copy : bool, optional
        If True, ensure that output map does not share memory with
        the input map.   Use if you do not want in-place operations to
        modify the map contents.
    partial : bool, optional
        If True, the map is not checked to ensure a proper healpix nside,
        and the number of pixels is returned instead.
    dtype : dtype, optional
        Data type required of the map. Default: np.double.

    Returns
    -------
    map_out : numpy.ndarray
        Properly shaped and memory-aligned map, copied from the input
        if necessary.
    nside or npix : int
        If partial is False, the map nside. Otherwise, the number of pixels.
    """
    map_out = np.atleast_2d(map_in)
    dim2 = map_out.shape[1] if partial else npix2nside(map_out.shape[1])
    map_out = np.require(map_out, dtype, ["A", "C"])
    if copy and np.may_share_memory(map_in, map_out):
        map_out = map_out.copy()
    return map_out, dim2


def check_proj(proj_in, copy=False, partial=False):
    """
    Return a properly shaped projection map, its nside, and the solution nmap.

    Arguments
    ---------
    proj_in : map or list of maps
        Input projection matrix map
    copy : bool, optional
        If True, ensure that output map does not share memory with
        the input map.   Use if you do not want in-place operations to
        modify the map contents.
    partial : bool, optional
        If True, the map is not checked to ensure a proper healpix nside,
        and the number of pixels is returned instead.

    Returns
    -------
    map_out : numpy.ndarray
        Properly shaped and memory-aligned map, copied from the input
        if necessary.
    nside or npix : int
        If partial is False, the map nside. Otherwise, the number of pixels.
    nmap : int
        The map size this projection matrix is intended to invert, i.e.
        the solution to `len(proj) = nmap * (nmap + 1) / 2`.  Raises an
        error if an integer solution is not found.
    """
    proj_out, dim2 = check_map(proj_in, copy=copy, partial=partial)
    nmap = int((np.sqrt(8 * len(proj_out) + 1) - 1)) // 2
    if nmap * (nmap + 1) // 2 != len(proj_out):
        raise ValueError("proj has incompatible shape")
    return proj_out, dim2, nmap


# ---------------------------------------------------------------------------
# QMap
# ---------------------------------------------------------------------------


class QMap(QPoint):
    """Quaternion-based mapmaker over the C++ core."""

    def __init__(
        self,
        nside=None,
        pol=None,
        vpol=None,
        source_map=None,
        source_pol=None,
        source_vpol=None,
        q_bore=None,
        ctime=None,
        q_hwp=None,
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.reset()

        if nside is not None:
            self.init_dest(nside=nside, pol=pol, vpol=vpol)

        if source_map is not None:
            self.init_source(source_map, pol=source_pol, vpol=source_vpol)

        if q_bore is not None:
            self.init_point(q_bore, ctime=ctime, q_hwp=q_hwp)

    # ---- Internal structures ----

    def reset(self):
        """
        Reset all internal structures.
        """
        self._source = lib.QpMap()
        self._dest = lib.QpMap()
        self._point = lib.QpPoint()
        self._detarr = None

    # ---- Source map ----

    def source_is_init(self):
        """
        Whether a source map has been initialized.

        Returns
        -------
        bool
            True once :meth:`init_source` has been called.
        """
        return self._source.is_init()

    def init_source(
        self,
        source_map,
        pol=None,
        pixels=None,
        nside=None,
        vpol=None,
        reset=False,
        update=False,
    ):
        """
        Initialize the source map structure for map2tod.

        Arguments
        ---------
        source_map : array_like
            Input map.  Must be of shape `(N, npix)`, where `N` can be
            1, 3, 6, 9, or 18.
        pol : bool, optional
            If `True`, and the map shape is `(3, npix)`, then input is a
            polarized map (and not T + first derivatives).  Only consulted
            where the row count leaves the mode open, which for a source
            map is 3 rows; defaults to `True`.
        pixels : 1D array_like, optional
            Array of pixel numbers for each map index, if `source_map` is
            a partial map.
        nside : int, optional
            map dimension.  If `pixels` is supplied, this argument is required.
            Otherwise, the nside is determined from the input map.
        vpol : bool, optional
            If `True`, and the input map shape is `(4, npix)`, then input is
            a polarized map that includes V polarization.  Defaults to
            `False`, so a 4-row map is rejected unless this is set.
        reset : bool, optional
            If `True`, and if the structure has already been initialized,
            it is reset and re-initialized with the new map.  If `False`,
            a `RuntimeError` is raised if the structure has already been
            initialized.
        update : bool, optional
            If `True`, and if the structure has already been initialized,
            the supplied `source_map` is replaced in the existing source
            structure rather than reinitializing from scratch.

        Notes
        -----
        This method will automatically determine the type of map
        given its shape.  Note that for `N=3`, the `pol` keyword argument
        should be used to disambiguate the two map types.  By default,
        a polarized map with `(T,Q,U)` components is assumed.

        * A map of shape `(1, npix)` or `(npix,)` contains only a `T` map.
        * A map of shape `(3, npix)` contains `(T, Q, U)` if `pol` is `True`,
          or `(T, dTdt, dTdp)` if `pol` is `False`.
        * A map of shape `(4, npix)` contains `(T, Q, U, V)`.
        * A map of shape `(6, npix)` contains `(T, dTdt, dTdp, dT2dt2,
          dT2dpdt, dT2dp2)`.
        * A map of shape `(9, npix)` contains `(T, Q, U, dTdt, dQdt, dUdt,
          dTdp, dQdp, dUdp)`.
        * A map of shape `(18, npix)` contains all the columns of the
          9-column map, followed by `(dT2dt2, dQ2dt2, dU2dt2, dT2dpdt,
          dQ2dpdt, dU2dpdt, dT2dp2, dQ2dp2, dU2dp2)`.
        """
        # The binding layer takes arrays as they come and never reshapes
        # them, so shaping the map is this layer's job -- and doing it once
        # here is what keeps it from happening twice.
        smap, _ = check_map(source_map, partial=pixels is not None)

        if self.source_is_init():
            if reset:
                self._source = lib.QpMap()
            elif update:
                # The mode is fixed at setup; update only swaps buffers.
                self._source.update_vec(smap)
                return
            else:
                raise RuntimeError("source already initialized")

        # False, not None: a source map has no projection component at all,
        # where None would have setup allocate an empty one. setup infers
        # nside and the mode from the map.
        self._source.setup(smap, False, pixels, nside, pol, vpol)

    def reset_source(self):
        """
        Release the source map.
        """
        self._source = lib.QpMap()

    def source_is_pol(self):
        """
        Whether the source map is polarized.

        Returns
        -------
        bool
            True if the source map carries polarization.
        """
        if not self.source_is_init():
            raise RuntimeError("source map not initialized")
        return self._source.is_pol()

    def source_is_vpol(self):
        """
        Whether the source map carries a V component.

        Returns
        -------
        bool
            True if the source map carries the V component as well.
        """
        if not self.source_is_init():
            raise RuntimeError("source map not initialized")
        return self._source.is_vpol()

    # ---- Dest map ----

    def dest_is_init(self):
        """
        Whether a destination map has been initialized.

        Returns
        -------
        bool
            True once :meth:`init_dest` has been called.
        """
        return self._dest.is_init()

    def init_dest(
        self,
        nside=None,
        pol=None,
        vec=None,
        proj=None,
        pixels=None,
        vpol=None,
        copy=False,
        reset=False,
        update=False,
    ):
        """
        Initialize the destination map structure for tod2map.

        Arguments
        ---------
        nside : int, optional
            map dimension.  If `pixels` is supplied, this argument is required.
            Otherwise, the default is 256.
        pol : bool, optional
            If True, a polarized map will be created.  Consulted only when
            neither `vec` nor `proj` is supplied, since either of those
            settles the mode by its shape.  Defaults to True in that case.
        vec : array_like or bool, optional, shape (N, npix)
            If supplied, nside and pol are determined from this map, and
            the vector (binned signal) map is initialized from this.
            If False, accumulation of this map from timestreams is disabled.
        proj : array_like or bool, optional, shape (N*(N+1)/2, npix)
            Array of upper-triangular elements of the projection matrix
            for each pixel.  If not supplied, a blank map of the appropriate
            shape is created. If False, accumulation of the projection matrix
            is disabled.
        pixels : 1D array_like, optional
            Array of pixel numbers for each map index, if `vec` and `proj` are
            partial maps.
        vpol : bool, optional
            If True, a polarized map including V polarization will be created.
            Read the same way as `pol`, and defaults to False.
        copy : bool, optional
            If True and vec/proj are supplied, make copies of these inputs
            to avoid in-place operations.
        reset : bool, optional
            If True, and if the structure has already been initialized,
            it is reset and re-initialized with the new map.  If False,
            a RuntimeError is raised if the structure has already been
            initialized.
        update : bool, optional
            If True, and if the structure has already been initialized,
            the supplied vec and proj are replaced in the existing dest
            structure rather than reinitializing from scratch.
        """
        if vec is False and proj is False:
            raise ValueError("one of vec or proj must not be False")

        # Shape each supplied map exactly once: the binding layer takes
        # arrays as they come and never reshapes them, and normalizing in
        # both places is what used to transpose a small partial map back.
        partial = pixels is not None
        if vec is not None and vec is not False:
            vec, _ = check_map(vec, copy=copy, partial=partial)
        if proj is not None and proj is not False:
            proj, _ = check_map(proj, copy=copy, partial=partial)

        if self.dest_is_init():
            if reset:
                self._dest = lib.QpMap()
            elif update:
                # Same three forms as below, against the installed maps: a
                # map replaces one, None replaces it with zeros, False
                # disables it. None here means "not supplied", so a
                # component that is switched off stays off -- the binding
                # would read it as a request for a fresh one.
                if vec is not None or self._dest.has_vec():
                    self._dest.update_vec(vec)
                if proj is not None or self._dest.has_proj():
                    self._dest.update_proj(proj)
                return
            else:
                raise RuntimeError("dest already initialized")

        # pol and vpol decide the mode only where nothing else does. A
        # supplied vec settles it here, and setup settles it from a
        # supplied proj, sizing the defaulted vec to match and reading the
        # mode back off that. Both flags default to None rather than to
        # T,Q,U so that the signature says which case is which: the
        # default is what you get when neither component was supplied, not
        # something imposed over a map whose shape already answered.
        if vec is not None and vec is not False:
            pol = len(vec) >= 3
            vpol = len(vec) == 4
        elif proj is None or proj is False:
            pol = True if pol is None else pol
            vpol = False if vpol is None else vpol

        # setup takes it from here: it allocates zeros for a component
        # passed as None, drops one passed as False, and infers nside and
        # npix from whichever of the maps, pixels or nside it is given.
        self._dest.setup(vec, proj, pixels, nside, pol, vpol)

    def reset_dest(self):
        """
        Release the destination map and its projection matrix.
        """
        self._dest = lib.QpMap()

    def dest_is_pol(self):
        """
        Whether the destination map is polarized.

        Returns
        -------
        bool
            True if the destination map carries polarization.
        """
        if not self.dest_is_init():
            raise RuntimeError("dest map not initialized")
        return self._dest.is_pol()

    def dest_is_vpol(self):
        """
        Whether the destination map carries a V component.

        Returns
        -------
        bool
            True if the destination map carries the V component as well.
        """
        if not self.dest_is_init():
            raise RuntimeError("dest map not initialized")
        return self._dest.is_vpol()

    # ---- Pointing ----

    def point_is_init(self):
        """
        Whether the boresight pointing has been initialized.

        Returns
        -------
        bool
            True once :meth:`init_point` has been called.
        """
        return self._point.is_init()

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
        if q_bore is not None:
            self._point.set_bore(q_bore)

        if not self.point_is_init():
            raise RuntimeError("point not initialized")

        if ctime is False:
            self._point.clear_ctime()
        elif ctime is not None:
            self._point.set_ctime(ctime)

        if q_hwp is False:
            self._point.clear_hwp()
        elif q_hwp is not None:
            self._point.set_hwp(q_hwp)

    def reset_point(self):
        """
        Release the boresight pointing.
        """
        self._point = lib.QpPoint()

    # ---- Detarr ----

    def init_detarr(
        self,
        q_off,
        weight=None,
        gain=None,
        mueller=None,
        tod=None,
        flag=None,
        weights=None,
        do_diff=False,
        write=False,
    ):
        """
        Initialize the detector array structure.

        Arguments
        ---------
        q_off : array_like
            Array of offset quaternions, of shape (ndet, 4).
        weight : array_like, optional
            Per-channel mapmaking weights, of shape (ndet,) or a constant.
            Default : 1.
        gain : array_like, optional
            Per-channel gains, of shape (ndet,) or a constant.
            Default : 1.
        mueller : array_like, optional
            Per-channel polarization efficiencies, of shape(ndet, 4).
            Default : [1., 1., 0., 1.] per det.
        tod : array_like, optional
            Timestream array, of shape (ndet, nsamp).  nsamp must match that of
            the pointing structure.  If not supplied and `write` is True, then
            a zero-filled timestream array is initialized.
        flag : array_like, optional
            Flag array, of shape (ndet, nsamp), for excluding data from
            mapmaking.  nsamp must match that of the pointing structure.
            If not supplied, a zero-filled array is initialized (i.e. no
            flagged samples).
        weights : array_like, optional
            Weight array, of shape (ndet, nsamp), for weighting each sample of
            data.  nsamp must match that of the pointing structure.  If not
            supplied, this option is not used.
        do_diff : bool, optional
            If True, initialize pairs of arrays for pair-differenced mapmaking.
        write : bool, optional
            If True, the timestreams are ensured writable and created if
            necessary.
        """
        self.reset_detarr()
        ns = self._point.n_samples()
        self._detarr = lib.QpDetArr()
        self._detarr.setup(
            q_off, weight, gain, mueller, tod, flag, weights, ns, do_diff, write
        )

    def reset_detarr(self):
        """
        Release the detector array and the timestreams it holds.
        """
        self._detarr = None

    # ---- Mapmaking ----

    @qp_settings
    def from_tod(
        self,
        q_off,
        tod=None,
        count_hits=True,
        weight=None,
        gain=None,
        mueller=None,
        flag=None,
        weights=None,
        do_diff=False,
    ):
        """
        Accumulate signal and projection maps from TOD.

        Returns (vec,), (proj,), or (vec, proj) depending on tod and count_hits.

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
        gain : array_like, optional
            Per-channel gains, of shape (ndet,) or a constant.
            Default : 1.
        mueller : array_like, optional
            array of Mueller matrix A/B/C elements, of shape (ndet,3).  Defaults to
            [1, 1, 0] per channel if not supplied.
        flag : array_like, optional
            array of flag timestreams for each channel, of shape (ndet, nsamp).
        weights : array_like, optional
            array of weight timestreams for each channel, of shape (ndet, nsamp).
        do_diff: do timestream differencing. Assumes first half of tods are
            one pair and the second half are the other.
        do_diff : bool, optional
            If True, difference detector pairs: detector i is paired with
            detector i + ndet/2, their sum accumulated into temperature and
            their difference into polarization. A sample flagged on either
            detector of a pair is dropped.

        Returns
        -------
        vec : array_like, optional
            binned signal map, if tod is supplied
        proj : array_like, optional
            binned projection matrix map, if count_hits is True
        """
        self.init_detarr(
            q_off,
            weight=weight,
            gain=gain,
            mueller=mueller,
            tod=tod,
            flag=flag,
            weights=weights,
            do_diff=do_diff,
        )

        return_vec = tod is not None and tod is not False and self._dest.has_vec()
        return_proj = count_hits and self._dest.has_proj()
        if not return_vec and not return_proj:
            raise RuntimeError("Nothing to do")

        try:
            lib.tod2map(self, self._detarr, self._point, self._dest)
        finally:
            self.reset_detarr()

        out = ()
        if return_vec:
            out += (self._dest.get_vec().squeeze(),)
        if return_proj:
            out += (self._dest.get_proj().squeeze(),)
        return out[0] if len(out) == 1 else out

    @qp_settings
    def to_tod(self, q_off, gain=None, mueller=None, tod=None, flag=None):
        """
        Sample timestreams from the source map.

        Returns tod array of shape (ndet, nsamp).

        Arguments
        ---------
        q_off : array_like
            quaternion offset array, of shape (ndet, 4)
        gain : array_like, optional
            Per-channel gains, of shape (ndet,) or a constant.
            Default : 1.
        mueller : array_like, optional
            array of Mueller matrix A/B/C/D elements, of shape (ndet, 4).  Defaults t
            [1, 1, 0, 1] per channel if not supplied.
        tod : array_like, optional
            output array for timestreams, of shape (ndet, nsamp)
            use this keyword argument for in-place computation.
        flag : array_like, optional
            Flag array of shape (ndet, nsamp). A flagged sample is left as it
            was found rather than written.

        Returns
        -------
        tod : array_like
            A timestream sampled from the input map for each requested detector.
            The output array shape is (ndet, nsamp).
        """
        self.init_detarr(
            q_off, gain=gain, mueller=mueller, tod=tod, flag=flag, write=True
        )

        try:
            lib.map2tod(self, self._detarr, self._point, self._source)
            tod = self._detarr.get_tod()
        finally:
            self.reset_detarr()

        return tod

    # ---- Projection matrix utilities (pure numpy, identical to qpoint) ----

    def proj_cond(self, proj=None, mode=None, partial=False):
        """
        Hits-normalized condition number per pixel.

        mode selects how it is computed: any order accepted by
        numpy.linalg.cond uses that norm, via an SVD. 'eigh' uses the
        symmetric eigenvalues instead, which is faster but not a drop-in:
        for a nearly singular pixel the smallest eigenvalue is rounding
        noise, so the ratio can differ from the SVD by tens of percent.

        Arguments
        ---------
        proj : array_like
            Projection matrix, of shape (N*(N+1)/2, npix).
            If None, the projection matrix installed by
            :meth:`init_dest`.
        mode : {None, 1, -1, 2, -2, inf, -inf, 'fro'}, optional
            condition number order.  See `numpy.linalg.cond`.
            Default: None (2-norm from SVD)
        partial : bool, optional
            If True, the map is not checked to ensure a proper healpix nside.

        Returns
        -------
        cond : array_like
            Condition number of each pixel.
        """
        if proj is None:
            proj = (
                self._dest.get_proj()
                if self._dest.is_init() and self._dest.has_proj()
                else None
            )
        if proj is None or proj is False:
            raise ValueError("missing proj")
        # Not copied: nothing below writes to it. The normalization is
        # applied to the hit columns pulled out of it, which is a fresh
        # small array -- at nside 512 a full-sky proj is 150 MB and a real
        # scan hits a fraction of a percent of it.
        proj, _, nmap = check_proj(proj, copy=False, partial=partial)
        nproj = len(proj)

        m = proj[0].astype(bool)
        sel = np.flatnonzero(m)
        cond = np.full(len(m), np.inf)

        if nmap == 1:
            # Hits-normalized, so a 1x1 matrix is its own hit count over
            # itself. Written as the division rather than as 1.0 so a
            # degenerate hit count gives what it always gave.
            cond[sel] = proj[0, sel] / proj[0, sel]
            return cond

        idx = np.zeros((nmap, nmap), dtype=int)
        rtri, ctri = np.triu_indices(nmap)
        idx[rtri, ctri] = idx[ctri, rtri] = np.arange(nproj)

        # Only the hit pixels: an unhit pixel's condition number is
        # infinite by definition, and a full-sky map of a real scan is
        # almost all unhit. Each pixel's matrix is independent, so the
        # ones that are kept are unchanged.
        mats = (proj[:, sel] / proj[0, sel])[idx].transpose(2, 0, 1)
        if mode == "eigh":
            # A projection matrix is symmetric and positive semi-definite,
            # so its singular values are the moduli of its eigenvalues.
            w = np.abs(np.linalg.eigvalsh(mats))
            with np.errstate(divide="ignore", invalid="ignore"):
                cond[sel] = w.max(axis=-1) / w.min(axis=-1)
        else:
            cond[sel] = np.linalg.cond(mats, p=mode)
        cond[cond > 1.0 / np.finfo(float).eps] = np.inf
        return cond

    def solve_map(
        self,
        vec=None,
        proj=None,
        mask=None,
        copy=True,
        return_proj=False,
        return_mask=False,
        partial=None,
        fill=0,
        cond=None,
        cond_thresh=1e6,
        cond_mode=None,
        method="exact",
    ):
        """
        Solve for a map from the binned signal and projection matrix.

        Arguments
        ---------
        vec : array_like, optional
            A map or list of N maps.  If None, the signal map installed
            by :meth:`init_dest`.
        proj : array_like, optional
            An array of upper-triangular projection matrices for each pixel, of
            shape (N*(N+1)/2, npix).  If None, the projection matrix
            installed by :meth:`init_dest`.
        mask : array_like, optional
            A mask of shape (npix,), evaluates to True where pixels are valid.
            The input mask in converted to a boolean array if supplied.
        copy : bool, optional
            if False, do the computation in-place so that the input maps are
            modified.  Otherwise, a copy is created prior to solving.
            Default: False.
        return_proj : bool, optional
            if True, return the Cholesky-decomposed projection matrix.
            if False, and inplace is True, the input projection matrix
            is not modified.
        return_mask : bool, optional
            if True, return the mask array, updated with any pixels
            that could not be solved.
        partial : bool, optional
            If True, the map is not checked to ensure a proper healpix nside.
        fill : scalar, optional
            Fill the solved map where proj == 0 with this value.  Default: 0.
        cond : array_like, optional
            A map of condition number per pixel.  If not supplied, this will be
            calculated using :meth:`proj_cond`
        cond_thresh : scalar, optional
            A threshold to place on the condition number to exclude pixels
            prior to solving.  Reduce this to avoid `LinAlgError` due to
            singular matrices.
        method : string, optional
            Map inversion method.  If "exact", invert the pointing matrix directly
            If "cho", use Cholesky decomposition to solve.  Default: "exact".
        cond_mode : str, optional
            Passed to :meth:`proj_cond` when the condition number has to be
            computed here. Default: None, the singular values.

        Returns
        -------
        map : array_like
            A solved map or set of maps, in shape (N, npix).
        proj_out : array_like
            The upper triangular elements of the decomposed projection matrix,
            (if method is 'cho') or of the matrix inverse (if method is 'exact'),
            if requested, in shape (N*(N+1)/2, npix).
        mask : array_like
            1-d array, True for valid pixels, if `return_mask` is True
        """
        if partial is None:
            partial = self._dest.is_init() and self._dest.is_partial()

        if vec is None:
            vec = (
                self._dest.get_vec()
                if self._dest.is_init() and self._dest.has_vec()
                else None
            )
        # False is how this API spells "no such component", so a caller who
        # passes one through gets the same message qpoint gives rather than
        # a complaint about npix from further down
        if vec is None or vec is False:
            raise ValueError("missing vec")
        vec, nside = check_map(vec, copy=copy, partial=partial)

        if proj is None:
            proj = (
                self._dest.get_proj()
                if self._dest.is_init() and self._dest.has_proj()
                else None
            )
        if proj is None or proj is False:
            raise ValueError("missing proj")
        # proj is written only where it is returned, so where it is not
        # there is nothing to protect the caller's array from.
        pcopy = copy if return_proj else False
        proj, pnside, nmap = check_proj(proj, copy=pcopy, partial=partial)

        if pnside != nside or nmap != len(vec):
            raise ValueError("vec and proj have incompatible shapes")
        nproj = len(proj)

        if mask is None:
            # the hit pixels are the whole mask; no need for a ones array
            # and a second pass to and it away
            mask = proj[0] != 0
        else:
            mask, mnside = check_map(mask, copy=copy, partial=partial)
            if mnside != nside:
                raise ValueError("mask has incompatible shape")
            mask = mask.squeeze().astype(bool)
        mask &= proj[0].astype(bool)

        if len(vec) == 1:
            vec = vec.squeeze()
            proj = proj.squeeze()
            vec[mask] /= proj[mask]
            vec[~mask] = fill
            ret = (vec,) + (proj,) * return_proj + (mask,) * return_mask
            return ret[0] if len(ret) == 1 else ret

        idx = np.zeros((nmap, nmap), dtype=int)
        rtri, ctri = np.triu_indices(nmap)
        idx[rtri, ctri] = idx[ctri, rtri] = np.arange(nproj)

        # An ill-conditioned pixel is excluded whichever solver runs, so
        # this sits above the choice of one: cho_factor succeeds on a
        # rank-deficient matrix rather than raising, and would otherwise
        # return whatever it made of a one- or two-hit pixel.
        if cond is None and cond_thresh is not None:
            cond = self.proj_cond(proj=proj, mode=cond_mode, partial=partial)
        if cond is not None:
            mask &= cond < cond_thresh

        if method == "exact":
            # Solve where there is something to solve. Everything else
            # used to be handed an identity matrix to invert, which for a
            # full-sky map of a real scan is almost every pixel.
            sel = np.flatnonzero(mask)
            vec[:, sel] = np.linalg.solve(
                proj[:, sel][idx].transpose(2, 0, 1),
                vec[:, sel].transpose()[..., np.newaxis],
            )[..., 0].transpose()
            vec[:, ~mask] = fill
            # Only worth doing to a proj the caller gets back; otherwise it
            # is a throwaway copy and this is a pass over the whole map.
            if return_proj:
                proj[:, ~mask] = 0
            ret = (vec,) + (proj,) * return_proj + (mask,) * return_mask
            return ret[0] if len(ret) == 1 else ret

        if method != "cho":
            raise ValueError(f"Unrecognized method {method!r}")

        from scipy.linalg import LinAlgError, cho_factor, cho_solve

        # Only the pixels worth solving, as above. This used to expand
        # proj[idx] for the whole map -- an (nmap, nmap, npix) temporary --
        # and then walk all of it in Python to skip the masked ones.
        sel = np.flatnonzero(mask)
        vec[:, ~mask] = fill
        if return_proj:
            proj[:, ~mask] = 0
        mats = proj[:, sel][idx].transpose(2, 0, 1)

        for jj, ii in enumerate(sel):
            A = mats[jj]
            try:
                # cho_factor overwrites A, which is why proj takes the
                # decomposition from it afterwards rather than from proj
                vec[:, ii] = cho_solve(cho_factor(A, False, True), vec[:, ii], True)
            # LinAlgError is the factorization giving up on a matrix that is
            # not positive definite; ValueError is cho_factor refusing one
            # that holds an inf or a nan. Both mean this pixel cannot be
            # solved; anything else is a bug and should not be swallowed.
            except (LinAlgError, ValueError):
                mask[ii] = False
                vec[:, ii] = fill
                if return_proj:
                    proj[:, ii] = 0
            else:
                if return_proj:
                    proj[:, ii] = A[rtri, ctri]

        ret = (vec,) + (proj,) * return_proj + (mask,) * return_mask
        return ret[0] if len(ret) == 1 else ret

    def solve_map_cho(self, *args, **kwargs):
        """
        Solve using Cholesky decomposition.

        Arguments
        ---------
        vec : array_like, optional
            A map or list of N maps.  If None, the signal map installed
            by :meth:`init_dest`.
        proj : array_like, optional
            An array of upper-triangular projection matrices for each pixel, of
            shape (N*(N+1)/2, npix).  If None, the projection matrix
            installed by :meth:`init_dest`.
        mask : array_like, optional
            A mask of shape (npix,), evaluates to True where pixels are valid.
            The input mask in converted to a boolean array if supplied.
        copy : bool, optional
            if False, do the computation in-place so that the input maps are
            modified.  Otherwise, a copy is created prior to solving.
            Default: False.
        return_proj : bool, optional
            if True, return the Cholesky-decomposed projection matrix.
            if False, and inplace is True, the input projection matrix
            is not modified.
        return_mask : bool, optional
            if True, return the mask array, updated with any pixels
            that could not be solved.
        partial : bool, optional
            If True, the map is not checked to ensure a proper healpix nside.
        fill : scalar, optional
            Fill the solved map where proj == 0 with this value.  Default: 0.
        cond : array_like, optional
            A map of condition number per pixel.  If not supplied, this will be
            calculated using :meth:`proj_cond`
        cond_thresh : scalar, optional
            A threshold to place on the condition number to exclude pixels
            prior to solving.  Reduce this to avoid `LinAlgError` due to
            singular matrices.

        Returns
        -------
        map : array_like
            A solved map or set of maps, in shape (N, npix).
        proj_out : array_like
            The upper triangular elements of the decomposed projection matrix,
            if requested, in shape (N*(N+1)/2, npix).
        mask : array_like
            1-d array, True for valid pixels, if `return_mask` is True
        """
        kwargs["method"] = "cho"
        return self.solve_map(*args, **kwargs)

    def unsolve_map(
        self,
        map_in,
        proj=None,
        mask=None,
        copy=True,
        return_proj=False,
        return_mask=False,
        partial=None,
        fill=0,
    ):
        """
        Invert the solved map to recover the binned vec array.

        Arguments
        ---------
        map_in : array_like
            A map or list of N maps.
        proj : array_like, optional
            An array of upper-triangular projection matrices for each pixel, of
            shape (N*(N+1)/2, npix).  If None, the projection matrix
            installed by :meth:`init_dest`.
        mask : array_like, optional
            A mask of shape (npix,), evaluates to True where pixels are valid.
            The input mask in converted to a boolean array if supplied.
        copy : bool, optional
            if False, do the computation in-place so that the input maps are
            modified.  Otherwise, a copy is created prior to solving.
            Default: False.
        return_proj : bool, optional
            if True, return the Cholesky-decomposed projection matrix.
            if False, and inplace is True, the input projection matrix
            is not modified.
        return_mask : bool, optional
            if True, return the mask array, updated with any pixels
            that could not be solved.
        partial : bool, optional
            If True, the map is not checked to ensure a proper healpix nside.
        fill : scalar, optional
            Fill the solved map where proj == 0 with this value.  Default: 0.

        Returns
        -------
        vec : array_like
            A binned map or set of maps, in shape (N, npix).
        proj_out : array_like
            The upper triangular elements of the projection matrix,
            if requested, in shape (N*(N+1)/2, npix).
        mask : array_like
            1-d array, True for valid pixels, if `return_mask` is True
        """
        if partial is None:
            partial = self._dest.is_init() and self._dest.is_partial()

        map_in, nside = check_map(map_in, copy=copy, partial=partial)

        if proj is None:
            proj = (
                self._dest.get_proj()
                if self._dest.is_init() and self._dest.has_proj()
                else None
            )
        if proj is None or proj is False:
            raise ValueError("missing proj")
        # proj is written only where it is returned, so where it is not
        # there is nothing to protect the caller's array from.
        pcopy = copy if return_proj else False
        proj, pnside, nmap = check_proj(proj, copy=pcopy, partial=partial)

        if pnside != nside or nmap != len(map_in):
            raise ValueError("map_in and proj have incompatible shapes")
        nproj = len(proj)

        if mask is None:
            mask = proj[0] != 0
        else:
            mask, mnside = check_map(mask, copy=copy, partial=partial)
            if mnside != nside:
                raise ValueError("mask has incompatible shape")
            mask = mask.squeeze().astype(bool)
        mask &= proj[0].astype(bool)

        if len(map_in) == 1:
            map_in = map_in.squeeze()
            proj = proj.squeeze()
            map_in[mask] *= proj[mask]
            map_in[~mask] = fill
        else:
            idx = np.zeros((nmap, nmap), dtype=int)
            rtri, ctri = np.triu_indices(nmap)
            idx[rtri, ctri] = idx[ctri, rtri] = np.arange(nproj)
            # Only the masked-in pixels, as solve_map does: proj[idx] over
            # the whole map is an (nmap, nmap, npix) temporary, and every
            # pixel it computes outside the mask is then overwritten.
            sel = np.flatnonzero(mask)
            map_in[:, sel] = np.einsum(
                "ij...,j...->i...", proj[:, sel][idx], map_in[:, sel]
            )
            map_in[:, ~mask] = fill

        ret = (map_in,) + (proj,) * return_proj + (mask,) * return_mask
        return ret[0] if len(ret) == 1 else ret
