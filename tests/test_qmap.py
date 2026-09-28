"""Tests for qpoint.QMap mapmaking class and helper functions."""

import os
import sys

import numpy as np
import pytest
import qpoint
from qpoint._libqpoint import (
    QP_VEC_D1,
    QP_VEC_D1_POL,
    QP_VEC_D2,
    QP_VEC_POL,
    QP_VEC_TEMP,
)
from qpoint.qmap_class import nside2npix, npix2nside, check_map, check_proj

# Small nside for fast tests
NSIDE = 8
NPIX = 12 * NSIDE * NSIDE
N = 200  # enough samples to hit several pixels

# Reference observing parameters
CTIME = 1418662800.0
LON = 165.7
LAT = -77.6


def make_bore_and_ctime(qm, n=N):
    """Return q_bore and ctime covering a wide az range at fixed el."""
    az = np.linspace(0, 360, n, endpoint=False)
    el = 45.0 * np.ones(n)
    ctime = CTIME + np.arange(n, dtype=float)
    q_bore = qm.azel2bore(az, el, None, None, LON, LAT, ctime)
    return q_bore, ctime


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------


class TestNsideNpix:
    def test_nside2npix(self):
        assert nside2npix(8) == 768
        assert nside2npix(16) == 3072
        assert nside2npix(64) == 49152

    def test_npix2nside(self):
        assert npix2nside(768) == 8
        assert npix2nside(3072) == 16
        assert npix2nside(49152) == 64

    def test_roundtrip(self):
        for nside in (1, 2, 4, 8, 16, 32, 64, 128):
            assert npix2nside(nside2npix(nside)) == nside

    def test_npix2nside_invalid(self):
        with pytest.raises(ValueError):
            npix2nside(100)


class TestCheckMap:
    def test_1d_map(self):
        m = np.ones(NPIX)
        m_out, nside = check_map(m)
        assert m_out.shape == (1, NPIX)
        assert nside == NSIDE

    def test_2d_map(self):
        m = np.ones((3, NPIX))
        m_out, nside = check_map(m)
        assert m_out.shape == (3, NPIX)
        assert nside == NSIDE

    def test_transposed_map(self):
        m = np.ones((NPIX, 3))  # should be auto-transposed
        m_out, nside = check_map(m)
        assert m_out.shape == (3, NPIX)

    def test_copy(self):
        m = np.ones(NPIX)
        m_out, _ = check_map(m, copy=True)
        assert not np.may_share_memory(m, m_out)

    def test_partial(self):
        npix = 100
        m = np.ones(npix)
        m_out, n = check_map(m, partial=True)
        assert n == npix

    def test_invalid_nside(self):
        m = np.ones(100)
        with pytest.raises(ValueError):
            check_map(m)


class TestCheckProj:
    def test_temp_proj(self):
        proj = np.ones((1, NPIX))
        proj_out, nside, nmap = check_proj(proj)
        assert nmap == 1
        assert nside == NSIDE

    def test_pol_proj(self):
        proj = np.ones((6, NPIX))
        proj_out, nside, nmap = check_proj(proj)
        assert nmap == 3  # 3*(3+1)/2 = 6

    def test_vpol_proj(self):
        proj = np.ones((10, NPIX))
        proj_out, nside, nmap = check_proj(proj)
        assert nmap == 4  # 4*(4+1)/2 = 10

    def test_invalid_proj(self):
        proj = np.ones((5, NPIX))  # 5 is not triangular
        with pytest.raises(ValueError):
            check_proj(proj)


# ---------------------------------------------------------------------------
# QMap initialization
# ---------------------------------------------------------------------------


class TestQMapInit:
    def test_default_init(self):
        qm = qpoint.QMap(mean_aber=True)
        assert not qm.dest_is_init()
        assert not qm.source_is_init()
        assert not qm.point_is_init()

    def test_init_with_nside(self):
        qm = qpoint.QMap(nside=NSIDE, pol=False, mean_aber=True)
        assert qm.dest_is_init()

    def test_init_with_source(self):
        source = np.ones((1, NPIX))
        qm = qpoint.QMap(source_map=source, source_pol=False, mean_aber=True)
        assert qm.source_is_init()

    def test_init_with_bore(self):
        qm = qpoint.QMap(mean_aber=True)
        q_bore, ctime = make_bore_and_ctime(qm)
        qm2 = qpoint.QMap(q_bore=q_bore, ctime=ctime, mean_aber=True)
        assert qm2.point_is_init()


# ---------------------------------------------------------------------------
# init_dest
# ---------------------------------------------------------------------------


class TestInitDest:
    @pytest.mark.parametrize(
        "kwargs, vec_shape, proj_shape",
        [
            ({"pol": False}, (1, NPIX), (1, NPIX)),
            ({"pol": True}, (3, NPIX), (6, NPIX)),
            ({"vpol": True}, (4, NPIX), (10, NPIX)),
        ],
        ids=["temperature", "polarized", "vpol"],
    )
    def test_installed_shapes(self, kwargs, vec_shape, proj_shape):
        """
        The shapes init_dest installs, read out of the depo: it returns
        nothing, so the depo is how they are reached.
        """
        qm = qpoint.QMap(mean_aber=True)
        assert qm.init_dest(nside=NSIDE, **kwargs) is None
        assert np.shape(qm.depo["vec"]) == vec_shape
        assert np.shape(qm.depo["proj"]) == proj_shape

    def test_init_twice_raises(self):
        qm = qpoint.QMap(mean_aber=True)
        qm.init_dest(nside=NSIDE)
        with pytest.raises(RuntimeError):
            qm.init_dest(nside=NSIDE)

    def test_reset_and_reinit(self):
        qm = qpoint.QMap(mean_aber=True)
        qm.init_dest(nside=NSIDE, pol=False)
        qm.reset_dest()
        assert not qm.dest_is_init()
        qm.init_dest(nside=NSIDE, pol=True)
        assert qm.dest_is_init()

    def test_vec_false_proj_false_raises(self):
        qm = qpoint.QMap(mean_aber=True)
        with pytest.raises(ValueError):
            qm.init_dest(nside=NSIDE, vec=False, proj=False)


# ---------------------------------------------------------------------------
# init_source
# ---------------------------------------------------------------------------


class TestInitSource:
    def test_temperature_map(self):
        qm = qpoint.QMap(mean_aber=True)
        source = np.ones((1, NPIX))
        qm.init_source(source, pol=False)
        assert qm.source_is_init()
        assert not qm.source_is_pol()

    def test_polarized_map(self):
        qm = qpoint.QMap(mean_aber=True)
        source = np.ones((3, NPIX))
        qm.init_source(source, pol=True)
        assert qm.source_is_init()
        assert qm.source_is_pol()

    def test_vpol_map(self):
        qm = qpoint.QMap(mean_aber=True)
        source = np.ones((4, NPIX))
        qm.init_source(source, vpol=True)
        assert qm.source_is_init()
        assert qm.source_is_vpol()

    def test_init_twice_raises(self):
        qm = qpoint.QMap(mean_aber=True)
        source = np.ones((1, NPIX))
        qm.init_source(source, pol=False)
        with pytest.raises(RuntimeError):
            qm.init_source(source, pol=False)

    def test_reset_and_reinit(self):
        qm = qpoint.QMap(mean_aber=True)
        source = np.ones((1, NPIX))
        qm.init_source(source, pol=False)
        qm.reset_source()
        assert not qm.source_is_init()
        qm.init_source(source, pol=False)
        assert qm.source_is_init()


class TestVpolIsChecked:
    """
    A row count settles V polarization on its own -- four rows and nothing
    else -- so init_source never consulted vpol to read a map. It is not
    redundant but upstream: a caller uses it to pick which fields to read off
    disk, so it is verified against the map that arrived.
    """

    @pytest.mark.parametrize("nrow", [1, 3, 4])
    @pytest.mark.parametrize("vpol", [None, True, False])
    def test_vpol_against_the_row_count(self, nrow, vpol):
        qm = qpoint.QMap(mean_aber=True)
        smap = np.zeros((nrow, NPIX))
        if vpol is not None and bool(vpol) != (nrow == 4):
            with pytest.raises(ValueError, match="vpol="):
                qm.init_source(smap, vpol=vpol)
            return
        qm.init_source(smap, vpol=vpol)
        assert qm.source_is_vpol() == (nrow == 4)

    def test_pol_still_decides_where_the_rows_do_not(self):
        """Three rows is the only count it is asked about."""
        for pol, mode in ((None, QP_VEC_POL), (True, QP_VEC_POL), (False, QP_VEC_D1)):
            qm = qpoint.QMap(mean_aber=True)
            qm.init_source(np.zeros((3, NPIX)), pol=pol)
            assert qm._source.contents.vec_mode == mode


class TestPolAndVpolDefaultToNone:
    """
    Neither was ever the last word -- a supplied vec or proj settles the
    mode in init_dest, and the map settles it in init_source -- so a
    concrete default named something the shapes discard. None says "ask
    the shapes" and falls back to T,Q,U where nothing else answers.
    """

    @pytest.mark.parametrize(
        "kwargs, rows",
        [
            ({}, 3),
            ({"pol": None}, 3),
            ({"pol": True}, 3),
            ({"pol": False}, 1),
            ({"vpol": None}, 3),
            ({"vpol": True}, 4),
        ],
    )
    def test_init_dest_defaults(self, kwargs, rows):
        qm = qpoint.QMap(mean_aber=True)
        qm.init_dest(nside=NSIDE, **kwargs)
        assert qm.depo["vec"].shape == (rows, NPIX)

    def test_a_supplied_map_settles_it(self):
        """Which is the whole reason the defaults are None."""
        qm = qpoint.QMap(mean_aber=True)
        qm.init_dest(nside=NSIDE, vec=np.zeros((3, NPIX)))
        assert qm.depo["vec"].shape == (3, NPIX)
        assert qm._dest.contents.vec_mode == QP_VEC_POL

    # a dest has three modes, and one shape per component for each
    ROWS = {
        "vec": {"temp": 1, "pol": 3, "vpol": 4},
        "proj": {"temp": 1, "pol": 6, "vpol": 10},
    }

    @pytest.mark.parametrize("component", ["vec", "proj"])
    @pytest.mark.parametrize(
        "kwargs, mode",
        [
            ({"pol": False}, "pol"),
            ({"pol": True}, "temp"),
            ({"vpol": True}, "pol"),
            ({"vpol": False}, "vpol"),
        ],
    )
    def test_contradicting_it_is_an_error(self, component, kwargs, mode):
        """
        Either component settles the mode on its own, so a flag that
        disagrees is a contradiction rather than an override -- which is
        what it used to be, silently. None leaves it to the shapes.
        """
        rows = self.ROWS[component][mode]
        qm = qpoint.QMap(mean_aber=True)
        with pytest.raises(ValueError, match="contradicts the supplied map"):
            qm.init_dest(nside=NSIDE, **{component: np.zeros((rows, NPIX))}, **kwargs)

    @pytest.mark.parametrize("mode", ["temp", "pol", "vpol"])
    def test_a_supplied_proj_settles_the_mode(self, mode):
        """
        Read the same way as a vec, and enough on its own: a proj-only call
        sizes the vec from it, with or without a flag that agrees.
        """
        rows, vrows = self.ROWS["proj"][mode], self.ROWS["vec"][mode]
        flags = {"temp": {"pol": False}, "pol": {"pol": True}, "vpol": {"vpol": True}}[
            mode
        ]
        for kwargs in ({}, flags):
            qm = qpoint.QMap(mean_aber=True)
            qm.init_dest(nside=NSIDE, proj=np.zeros((rows, NPIX)), **kwargs)
            assert qm.depo["vec"].shape == (vrows, NPIX)


# ---------------------------------------------------------------------------
# init_dest / init_source with update=True
# ---------------------------------------------------------------------------


def update_mapper(**kwargs):
    """A pointed, polarized QMap and a one-detector offset array."""
    qm = qpoint.QMap(nside=NSIDE, pol=True, mean_aber=True, **kwargs)
    q_bore, ctime = make_bore_and_ctime(qm)
    qm.init_point(q_bore, ctime=ctime)
    return qm, np.atleast_2d(qm.det_offset(1.0, 2.0, 0.0))


class TestUpdateDest:
    """
    init_dest(update=True) replaces the maps in a dest that is already set up,
    which is what a scan wants between chunks: the pointing, the pixel hash,
    nside and the modes stay, and only the accumulators are swapped or zeroed.

    Several of its behaviours are surprising, and are pinned as they stand.
    """

    def test_both_components_are_replaced(self):
        """The depo is how the installed maps are reached; nothing is
        returned."""
        qm, _ = update_mapper()
        vec = np.full((3, NPIX), 2.0)
        proj = np.full((6, NPIX), 3.0)
        assert qm.init_dest(vec=vec, proj=proj, update=True) is None
        assert qm.depo["vec"] is vec
        assert qm.depo["proj"] is proj

    def test_the_mapmaker_accumulates_into_the_replacement(self):
        """The ctypes struct is rewritten, not just the depo entry."""
        qm, q_off = update_mapper()
        seed = np.zeros((3, NPIX))
        seed[0, :] = 3.0
        qm.init_dest(vec=seed, update=True)
        vec, _ = qm.from_tod(q_off, tod=np.ones((1, N)))
        assert vec[0].sum() > 3.0 * NPIX
        # Nothing was copied, so the accumulation landed in the caller's array.
        assert np.shares_memory(vec, seed)

    def test_defaults_zero_the_accumulators(self):
        """
        Supplying neither component is the per-chunk reset: each is replaced
        by zeros of its own shape, so the next chunk starts clean without
        re-deriving nside, the modes or the pixel hash.
        """
        qm, q_off = update_mapper()
        _, proj = qm.from_tod(q_off, tod=np.ones((1, N)))
        hits = proj[0].sum()
        assert hits > 0

        qm.init_dest(update=True)
        assert not qm.depo["vec"].any()
        assert not qm.depo["proj"].any()

        _, proj = qm.from_tod(q_off, tod=np.ones((1, N)))
        assert proj[0].sum() == hits

    def test_copy_keeps_the_callers_array_out_of_it(self):
        qm, _ = update_mapper()
        caller = np.ones((3, NPIX))
        qm.init_dest(vec=caller, update=True, copy=True)
        assert not np.shares_memory(qm.depo["vec"], caller)

        qm, _ = update_mapper()
        caller = np.ones((3, NPIX))
        qm.init_dest(vec=caller, update=True)
        assert np.shares_memory(qm.depo["vec"], caller)

    def test_wrong_npix_raises(self):
        qm, _ = update_mapper()
        with pytest.raises(ValueError, match="vec shape does not match"):
            qm.init_dest(vec=np.ones((3, NPIX // 2)), update=True)
        with pytest.raises(ValueError, match="proj shape does not match"):
            qm.init_dest(proj=np.ones((6, NPIX // 2)), update=True)

    def test_a_partial_dest_updates_too(self):
        """The shape check is against the pixel list, not against nside."""
        qm = qpoint.QMap(mean_aber=True)
        pixels = np.arange(50, dtype=np.int64)
        qm.init_dest(nside=NSIDE, pol=True, pixels=pixels)
        qm.init_dest(vec=np.full((3, 50), 4.0), proj=np.full((6, 50), 5.0), update=True)
        assert qm.depo["vec"].shape == (3, 50)
        assert qm.depo["proj"].shape == (6, 50)
        with pytest.raises(ValueError, match="vec shape does not match"):
            qm.init_dest(vec=np.ones((3, NPIX)), update=True)

    def test_update_on_a_fresh_qmap_initializes_instead(self):
        """
        update is only consulted once dest_is_init(), so on an uninitialized
        QMap the flag is ignored and this is an ordinary init_dest.
        """
        qm = qpoint.QMap(mean_aber=True)
        assert not qm.dest_is_init()
        qm.init_dest(nside=NSIDE, pol=True, update=True)
        assert qm.dest_is_init()
        assert qm.depo["vec"].shape == (3, NPIX)
        assert qm.depo["proj"].shape == (6, NPIX)

    def test_a_switched_off_component_stays_off(self):
        """
        And the map supplied for it is dropped without a word. Each branch
        is guarded by what the dest already has, so there is nowhere for a
        vec to go on a proj-only dest. Nothing signals it, now that there
        is no return value whose arity used to hint at it.
        """
        qm, _ = update_mapper()
        qm.reset_dest()
        qm.init_dest(nside=NSIDE, pol=True, vec=False)
        assert qm.depo["vec"] is False

        mine = np.full((3, NPIX), 99.0)
        qm.init_dest(vec=mine, update=True)
        assert qm.depo["vec"] is False
        assert (mine == 99.0).all()

    def test_the_row_count_must_match(self):
        """
        The whole shape is compared, not just the pixel axis, so the row
        count cannot change under update. It used to: num_vec followed the
        replacement while proj kept its own rows, leaving a dest that
        described a T map and a polarized one at once, which from_tod
        accumulated into and solve_map then refused.
        """
        qm, _ = update_mapper()
        with pytest.raises(ValueError, match="vec shape does not match"):
            qm.init_dest(vec=np.ones((1, NPIX)), update=True)
        assert qm._dest.contents.num_vec == 3
        assert qm._dest.contents.vec_mode == QP_VEC_POL

    def test_pol_is_ignored_on_update(self):
        """
        update swaps arrays; it does not reinterpret them. pol used to be
        consulted here, which let it write QP_VEC_D1 onto a *destination*
        map -- a mode nothing in tod2map reads, and one init_dest will not
        build, since a three-row vec forces pol True there.
        """
        qm, _ = update_mapper()
        assert qm._dest.contents.vec_mode == QP_VEC_POL
        qm.init_dest(vec=np.ones((3, NPIX)), pol=False, update=True)
        assert qm._dest.contents.vec_mode == QP_VEC_POL
        assert qm._dest.contents.vec_mode != QP_VEC_D1

    def test_false_is_not_how_a_component_is_switched_off(self):
        """
        It used to reach the shape check as if it were a map and fail on
        the attribute. Switching one off mid-scan still has no spelling
        here, but the error now names the argument and what to use
        instead.
        """
        qm, _ = update_mapper()
        with pytest.raises(ValueError, match="cannot switch vec off"):
            qm.init_dest(vec=False, update=True)


class TestUpdateSource:
    """
    The same for init_source, which is simpler: one map rather than two.

    The mode carries more weight here. qp_map2tod1 reads it with an exact
    switch where qp_tod2map1 only asks whether it is at least POL, so a
    source's mode decides which rows are read and whether the offset from the
    pixel centre is computed at all. And a source map may change how many
    derivative terms it carries, where a dest's vec may not.
    """

    def test_the_source_map_is_replaced(self):
        qm, q_off = update_mapper()
        qm.init_source(np.zeros((3, NPIX)), pol=True)
        assert qm.to_tod(q_off).sum() == 0.0

        source = np.zeros((3, NPIX))
        source[0] = 5.0
        assert qm.init_source(source, update=True) is None
        assert qm.depo["source_map"] is source
        assert qm.to_tod(q_off).sum() == pytest.approx(5.0 * N)

    def test_a_partial_source_updates_too(self):
        qm = qpoint.QMap(mean_aber=True)
        pixels = np.arange(50, dtype=np.int64)
        qm.init_source(np.ones((3, 50)), pol=True, pixels=pixels, nside=NSIDE)
        qm.init_source(np.full((3, 50), 2.0), update=True)
        assert qm.depo["source_map"][0, 0] == 2.0

    def test_update_on_a_fresh_qmap_initializes_instead(self):
        qm = qpoint.QMap(mean_aber=True)
        assert not qm.source_is_init()
        qm.init_source(np.ones((3, NPIX)), pol=True, update=True)
        assert qm.source_is_init()
        assert qm.source_is_pol()

    @pytest.mark.parametrize(
        "start, pol, replacement, mode",
        [
            (1, False, 3, QP_VEC_D1),
            (1, False, 6, QP_VEC_D2),
            (3, True, 9, QP_VEC_D1_POL),
            (3, False, 6, QP_VEC_D2),
            (9, True, 3, QP_VEC_POL),
            (18, True, 1, QP_VEC_TEMP),
        ],
    )
    def test_the_number_of_derivative_terms_may_change(
        self, start, pol, replacement, mode
    ):
        """
        Unlike a dest, whose vec and proj row counts have to correspond, a
        source map may carry a different number of derivative terms than
        the one it replaces -- so the row count is free and the mode
        follows it.

        This is also what reaches qp_reshape_map with the row count grown,
        which overran the row table until the fix two commits below: it was
        sized for the map being replaced and never resized.
        """
        qm = qpoint.QMap(mean_aber=True)
        qm.init_source(np.ones((start, NPIX)), pol=pol)
        qm.init_source(np.full((replacement, NPIX), 2.0), update=True)
        assert qm._source.contents.num_vec == replacement
        assert qm._source.contents.vec_mode == mode
        assert qm.depo["source_map"].shape == (replacement, NPIX)
        # The row pointers have to follow the new buffer, which is the
        # part that used to be written past the end of its table.
        assert (qm.depo["source_map"] == 2.0).all()

    def test_npix_must_still_match(self):
        """The pixelization is what update exists to hold on to."""
        qm = qpoint.QMap(mean_aber=True)
        qm.init_source(np.ones((3, NPIX)), pol=True)
        with pytest.raises(ValueError, match="source_map npix does not match"):
            qm.init_source(np.ones((3, NPIX // 2)), update=True)
        assert qm._source.contents.num_vec == 3

    @pytest.mark.parametrize("pol, mode", [(False, QP_VEC_D1), (True, QP_VEC_POL)])
    def test_three_rows_keep_their_reading_across_a_same_size_update(self, pol, mode):
        """
        Three rows are T,Q,U or T with first derivatives, and pol is what
        picks at init. An update of the same size leaves the mode alone, so
        the reading survives without the caller repeating pol -- which is
        the case that used to flip.
        """
        qm = qpoint.QMap(mean_aber=True)
        qm.init_source(np.ones((3, NPIX)), pol=pol)
        assert qm._source.contents.vec_mode == mode
        qm.init_source(np.full((3, NPIX), 2.0), update=True)
        assert qm._source.contents.vec_mode == mode
        # Even asking for the other reading does not get it.
        qm.init_source(np.ones((3, NPIX)), pol=not pol, update=True)
        assert qm._source.contents.vec_mode == mode
        # reset is how it changes.
        qm.init_source(np.ones((3, NPIX)), pol=not pol, reset=True)
        assert qm._source.contents.vec_mode != mode

    def test_a_round_trip_through_nine_rows_reads_back_polarized(self):
        """
        The polarization comes from the installed mode, through
        source_is_pol, rather than from a remembered copy of the original
        argument. So a D1 source taken to 9 rows -- D1_POL, which carries Q
        and U -- and back to 3 returns as T,Q,U rather than as D1. The
        installed mode is what the structure actually holds; a private copy
        able to disagree with it would be worse. reset is the way back.
        """
        qm = qpoint.QMap(mean_aber=True)
        qm.init_source(np.ones((3, NPIX)), pol=False)
        assert qm._source.contents.vec_mode == QP_VEC_D1

        qm.init_source(np.ones((9, NPIX)), update=True)
        assert qm._source.contents.vec_mode == QP_VEC_D1_POL
        assert qm.source_is_pol()

        qm.init_source(np.ones((3, NPIX)), update=True)
        assert qm._source.contents.vec_mode == QP_VEC_POL

    def test_the_map_mode_survives_an_update(self):
        """
        The consequential one. D1 is T with its two first derivatives, and
        map2tod evaluates the map away from the pixel centre with them; POL
        is T,Q,U and does not. Three rows are either, and pol is what picks
        at init -- so recomputing the mode from pol here, where it defaults
        to True, turned a derivative source into a polarized one and
        changed the timestream with no diagnostic at all.
        """
        smap = np.zeros((3, NPIX))
        smap[0], smap[1], smap[2] = 1.0, 0.5, -0.25

        qm, q_off = update_mapper()
        qm.init_source(smap, pol=False)
        assert qm._source.contents.vec_mode == QP_VEC_D1
        before = qm.to_tod(q_off).copy()

        qm.init_source(smap, update=True)
        assert qm._source.contents.vec_mode == QP_VEC_D1
        assert np.array_equal(qm.to_tod(q_off), before)

        # And the polarized reading of the same map really is different,
        # so the assertion above is not vacuous.
        qm2, q_off2 = update_mapper()
        qm2.init_source(smap, pol=True)
        assert not np.allclose(qm2.to_tod(q_off2), before)

    def test_pol_is_ignored_on_update(self):
        """Asking for the other reading does not get it; reset does."""
        qm, _ = update_mapper()
        qm.init_source(np.zeros((3, NPIX)), pol=False)
        qm.init_source(np.zeros((3, NPIX)), pol=True, update=True)
        assert qm._source.contents.vec_mode == QP_VEC_D1
        assert not qm.source_is_pol()

        qm.init_source(np.zeros((3, NPIX)), pol=True, reset=True)
        assert qm._source.contents.vec_mode == QP_VEC_POL
        assert qm.source_is_pol()


# ---------------------------------------------------------------------------
# init_point
# ---------------------------------------------------------------------------


class TestInitPoint:
    def test_init_with_qbore(self):
        qm = qpoint.QMap(mean_aber=True)
        q_bore, ctime = make_bore_and_ctime(qm)
        qm.init_point(q_bore, ctime=ctime)
        assert qm.point_is_init()

    def test_init_without_point_raises_when_mapping(self):
        qm = qpoint.QMap(nside=NSIDE, pol=False, mean_aber=True)
        q_off = np.atleast_2d(qm.det_offset(0.0, 0.0, 0.0))
        tod = np.ones((1, N))
        with pytest.raises(Exception):
            qm.from_tod(q_off, tod=tod)

    def test_init_with_hwp(self):
        qm = qpoint.QMap(mean_aber=True)
        q_bore, ctime = make_bore_and_ctime(qm)
        q_hwp = qm.hwp_quat(np.zeros(N))
        qm.init_point(q_bore, ctime=ctime, q_hwp=q_hwp)
        assert qm.point_is_init()


# ---------------------------------------------------------------------------
# to_tod: constant map -> all-ones TOD
# ---------------------------------------------------------------------------


class TestToTod:
    def test_temperature_constant_map(self):
        """Scanning a constant T=1 map should produce a TOD of all ones."""
        qm = qpoint.QMap(mean_aber=True)
        q_bore, ctime = make_bore_and_ctime(qm)
        source = np.ones((1, NPIX))
        qm.init_source(source, pol=False)
        qm.init_point(q_bore, ctime=ctime)
        q_off = np.atleast_2d(qm.det_offset(0.0, 0.0, 0.0))
        tod = qm.to_tod(q_off)
        assert tod.shape == (1, N)
        assert np.allclose(tod, 1.0, atol=1e-5)

    def test_temperature_zero_map(self):
        """Scanning a zero map should produce a TOD of all zeros."""
        qm = qpoint.QMap(mean_aber=True)
        q_bore, ctime = make_bore_and_ctime(qm)
        source = np.zeros((1, NPIX))
        qm.init_source(source, pol=False)
        qm.init_point(q_bore, ctime=ctime)
        q_off = np.atleast_2d(qm.det_offset(0.0, 0.0, 0.0))
        tod = qm.to_tod(q_off)
        assert tod.shape == (1, N)
        assert np.allclose(tod, 0.0, atol=1e-10)

    def test_tod_shape_multi_det(self):
        ndet = 4
        qm = qpoint.QMap(mean_aber=True)
        q_bore, ctime = make_bore_and_ctime(qm)
        source = np.ones((1, NPIX))
        qm.init_source(source, pol=False)
        qm.init_point(q_bore, ctime=ctime)
        q_off = np.atleast_2d(
            qm.det_offset([0, 0.5, -0.5, 0], [0, 0, 0, 0.5], [0] * ndet)
        )
        tod = qm.to_tod(q_off)
        assert tod.shape == (ndet, N)

    def test_polarized_constant_T_Q_U(self):
        """Scanning a pol map (T=1, Q=0, U=0) with a single det gives TOD of 1."""
        qm = qpoint.QMap(mean_aber=True)
        q_bore, ctime = make_bore_and_ctime(qm)
        source = np.zeros((3, NPIX))
        source[0] = 1.0  # T = 1 everywhere
        qm.init_source(source, pol=True)
        qm.init_point(q_bore, ctime=ctime)
        q_off = np.atleast_2d(qm.det_offset(0.0, 0.0, 0.0))
        tod = qm.to_tod(q_off)
        assert np.allclose(tod, 1.0, atol=1e-5)


# ---------------------------------------------------------------------------
# from_tod: binning TOD into map
# ---------------------------------------------------------------------------


class TestFromTod:
    def test_projection_map_nonzero_after_from_tod(self):
        """from_tod with any TOD should fill the hits map."""
        qm = qpoint.QMap(nside=NSIDE, pol=False, mean_aber=True)
        q_bore, ctime = make_bore_and_ctime(qm)
        qm.init_point(q_bore, ctime=ctime)
        q_off = np.atleast_2d(qm.det_offset(0.0, 0.0, 0.0))
        tod = np.ones((1, N))
        vec, proj = qm.from_tod(q_off, tod=tod)
        assert np.any(proj > 0)

    def test_proj_equals_hits(self):
        """With unit TOD weights, proj[pix] = number of hits on that pixel."""
        qm = qpoint.QMap(nside=NSIDE, pol=False, mean_aber=True)
        q_bore, ctime = make_bore_and_ctime(qm)
        qm.init_point(q_bore, ctime=ctime)
        q_off = np.atleast_2d(qm.det_offset(0.0, 0.0, 0.0))
        tod = np.ones((1, N))
        vec, proj = qm.from_tod(q_off, tod=tod)
        # For T-only: proj = hits, vec = sum(tod) per pixel
        # vec / proj should equal 1 everywhere observed
        mask = proj.squeeze() > 0
        assert np.allclose(vec.squeeze()[mask] / proj.squeeze()[mask], 1.0, atol=1e-10)

    def test_from_tod_pol_shape(self):
        qm = qpoint.QMap(nside=NSIDE, pol=True, mean_aber=True)
        q_bore, ctime = make_bore_and_ctime(qm)
        qm.init_point(q_bore, ctime=ctime)
        q_off = np.atleast_2d(qm.det_offset(0.0, 0.0, 0.0))
        tod = np.ones((1, N))
        vec, proj = qm.from_tod(q_off, tod=tod)
        assert vec.shape == (3, NPIX)
        assert proj.shape == (6, NPIX)

    def test_count_hits_false(self):
        """count_hits=False should not accumulate the projection map."""
        qm = qpoint.QMap(nside=NSIDE, pol=False, mean_aber=True)
        q_bore, ctime = make_bore_and_ctime(qm)
        qm.init_point(q_bore, ctime=ctime)
        q_off = np.atleast_2d(qm.det_offset(0.0, 0.0, 0.0))
        tod = np.ones((1, N))
        # First call accumulates vec and proj
        vec, proj = qm.from_tod(q_off, tod=tod)
        proj_copy = proj.squeeze().copy()
        # Second call with count_hits=False should not update proj
        tod2 = 2.0 * np.ones((1, N))
        vec = qm.from_tod(q_off, tod=tod2, count_hits=False)
        assert np.allclose(qm.depo["proj"].squeeze(), proj_copy, atol=1e-10)


# ---------------------------------------------------------------------------
# to_tod -> from_tod -> solve_map round-trip
# ---------------------------------------------------------------------------


class TestTodMapRoundtrip:
    def test_constant_temp_roundtrip(self):
        """Scan a constant T map, bin TOD back to map, solve -> recover constant."""
        qm = qpoint.QMap(nside=NSIDE, pol=False, mean_aber=True)
        q_bore, ctime = make_bore_and_ctime(qm)
        source = np.ones((1, NPIX))
        qm.init_source(source, pol=False)
        qm.init_point(q_bore, ctime=ctime)
        q_off = np.atleast_2d(qm.det_offset(0.0, 0.0, 0.0))

        # Produce TOD
        tod = qm.to_tod(q_off)

        # Bin TOD
        vec, proj = qm.from_tod(q_off, tod=tod)

        # Solve
        solved = qm.solve_map()

        # Every observed pixel should be 1
        mask = qm.depo["proj"].squeeze() > 0
        assert np.any(mask), "No pixels were observed"
        assert np.allclose(solved[mask], 1.0, atol=1e-5)

    def test_scaled_temp_roundtrip(self):
        """Scan a constant T=42 map and verify solved map is also 42."""
        scale = 42.0
        qm = qpoint.QMap(nside=NSIDE, pol=False, mean_aber=True)
        q_bore, ctime = make_bore_and_ctime(qm)
        source = scale * np.ones((1, NPIX))
        qm.init_source(source, pol=False)
        qm.init_point(q_bore, ctime=ctime)
        q_off = np.atleast_2d(qm.det_offset(0.0, 0.0, 0.0))

        tod = qm.to_tod(q_off)
        qm.from_tod(q_off, tod=tod)
        solved = qm.solve_map()

        mask = qm.depo["proj"].squeeze() > 0
        assert np.allclose(solved[mask], scale, atol=1e-5)


# ---------------------------------------------------------------------------
# solve_map
# ---------------------------------------------------------------------------


class TestSolveMap:
    def test_temperature_solve(self):
        # vec and proj have equal values so solved = vec/proj = 1 everywhere observed
        vec = np.array([[1.0, 2.0, 0.0, 3.0]])
        proj = np.array([[1.0, 2.0, 0.0, 3.0]])
        qm = qpoint.QMap(mean_aber=True)
        solved = qm.solve_map(vec=vec.copy(), proj=proj.copy(), partial=True)
        assert np.allclose(solved[[0, 1, 3]], 1.0, atol=1e-10)
        assert solved[2] == 0.0  # zero proj -> fill=0

    def test_returns_correct_shape_temperature(self):
        qm = qpoint.QMap(nside=NSIDE, pol=False, mean_aber=True)
        q_bore, ctime = make_bore_and_ctime(qm)
        qm.init_point(q_bore, ctime=ctime)
        q_off = np.atleast_2d(qm.det_offset(0.0, 0.0, 0.0))
        tod = np.ones((1, N))
        qm.from_tod(q_off, tod=tod)
        solved = qm.solve_map()
        assert solved.shape == (NPIX,)

    def test_returns_correct_shape_polarized(self):
        qm = qpoint.QMap(nside=NSIDE, pol=True, mean_aber=True)
        q_bore, ctime = make_bore_and_ctime(qm)
        qm.init_point(q_bore, ctime=ctime)
        q_off = np.atleast_2d(qm.det_offset(0.0, 0.0, 0.0))
        tod = np.ones((1, N))
        qm.from_tod(q_off, tod=tod)
        solved = qm.solve_map()
        assert solved.shape == (3, NPIX)

    def test_return_mask(self):
        qm = qpoint.QMap(nside=NSIDE, pol=False, mean_aber=True)
        q_bore, ctime = make_bore_and_ctime(qm)
        qm.init_point(q_bore, ctime=ctime)
        q_off = np.atleast_2d(qm.det_offset(0.0, 0.0, 0.0))
        tod = np.ones((1, N))
        qm.from_tod(q_off, tod=tod)
        solved, mask = qm.solve_map(return_mask=True)
        assert mask.shape == (NPIX,)
        assert mask.dtype == bool

    def test_fill_value(self):
        qm = qpoint.QMap(nside=NSIDE, pol=False, mean_aber=True)
        q_bore, ctime = make_bore_and_ctime(qm)
        qm.init_point(q_bore, ctime=ctime)
        q_off = np.atleast_2d(qm.det_offset(0.0, 0.0, 0.0))
        tod = np.ones((1, N))
        qm.from_tod(q_off, tod=tod)
        solved = qm.solve_map(fill=np.nan)
        mask = qm.depo["proj"].squeeze() > 0
        assert np.all(np.isnan(solved[~mask]))

    def test_cho_method(self):
        qm = qpoint.QMap(nside=NSIDE, pol=False, mean_aber=True)
        q_bore, ctime = make_bore_and_ctime(qm)
        qm.init_point(q_bore, ctime=ctime)
        q_off = np.atleast_2d(qm.det_offset(0.0, 0.0, 0.0))
        tod = np.ones((1, N))
        qm.from_tod(q_off, tod=tod)
        solved = qm.solve_map(method="cho")
        assert solved.shape == (NPIX,)

    def test_invalid_method_raises(self):
        # Method validation only applies to polarized maps (T-only returns early)
        qm = qpoint.QMap(nside=NSIDE, pol=True, mean_aber=True)
        q_bore, ctime = make_bore_and_ctime(qm)
        qm.init_point(q_bore, ctime=ctime)
        q_off = np.atleast_2d(qm.det_offset(0.0, 0.0, 0.0))
        tod = np.ones((1, N))
        qm.from_tod(q_off, tod=tod)
        with pytest.raises(ValueError):
            qm.solve_map(method="unknown")


# ---------------------------------------------------------------------------
# unsolve_map
# ---------------------------------------------------------------------------


class TestUnsolveMap:
    def test_unsolve_inverts_solve(self):
        """unsolve_map(solve_map(vec, proj), proj) should recover the original vec."""
        vec_orig = np.array([[1.0, 2.0, 0.0, 3.0]])
        proj = np.array([[1.0, 2.0, 0.0, 3.0]])
        qm = qpoint.QMap(mean_aber=True)
        solved = qm.solve_map(vec=vec_orig.copy(), proj=proj.copy(), partial=True)
        unsolved = qm.unsolve_map(solved.reshape(1, -1), proj=proj.copy(), partial=True)
        mask = proj.squeeze() > 0
        assert np.allclose(
            unsolved.squeeze()[mask], vec_orig.squeeze()[mask], atol=1e-10
        )


# ---------------------------------------------------------------------------
# proj_cond
# ---------------------------------------------------------------------------


class TestProjCond:
    def test_temperature_uniform_proj(self):
        """For T-only with uniform hits, condition number should be 1 everywhere."""
        nside = NSIDE
        npix = nside2npix(nside)
        proj = np.ones((1, npix))
        qm = qpoint.QMap(mean_aber=True)
        cond = qm.proj_cond(proj=proj, partial=True)
        assert cond.shape == (npix,)
        assert np.allclose(cond, 1.0, atol=1e-10)

    def test_empty_pixels_have_inf(self):
        """Pixels with zero hits should have infinite condition number."""
        nside = NSIDE
        npix = nside2npix(nside)
        proj = np.ones((1, npix))
        proj[0, : npix // 2] = 0.0  # half the pixels have no hits
        qm = qpoint.QMap(mean_aber=True)
        cond = qm.proj_cond(proj=proj, partial=True)
        assert np.all(np.isinf(cond[: npix // 2]))
        assert np.all(np.isfinite(cond[npix // 2 :]))

    def test_from_depo(self):
        """proj_cond() without args uses depo['proj']."""
        qm = qpoint.QMap(nside=NSIDE, pol=False, mean_aber=True)
        q_bore, ctime = make_bore_and_ctime(qm)
        qm.init_point(q_bore, ctime=ctime)
        q_off = np.atleast_2d(qm.det_offset(0.0, 0.0, 0.0))
        tod = np.ones((1, N))
        qm.from_tod(q_off, tod=tod)
        cond = qm.proj_cond()
        assert cond.shape == (NPIX,)


# ---------------------------------------------------------------------------
# reset methods
# ---------------------------------------------------------------------------


class TestResets:
    def test_reset(self):
        qm = qpoint.QMap(nside=NSIDE, pol=False, mean_aber=True)
        q_bore, ctime = make_bore_and_ctime(qm)
        source = np.ones((1, NPIX))
        qm.init_source(source, pol=False)
        qm.init_point(q_bore, ctime=ctime)
        qm.reset()
        assert not qm.dest_is_init()
        assert not qm.source_is_init()
        assert not qm.point_is_init()
        assert len(qm.depo) == 0

    def test_reset_dest(self):
        qm = qpoint.QMap(nside=NSIDE, pol=False, mean_aber=True)
        assert qm.dest_is_init()
        qm.reset_dest()
        assert not qm.dest_is_init()

    def test_reset_source(self):
        source = np.ones((1, NPIX))
        qm = qpoint.QMap(source_map=source, source_pol=False, mean_aber=True)
        assert qm.source_is_init()
        qm.reset_source()
        assert not qm.source_is_init()

    def test_reset_point(self):
        qm = qpoint.QMap(mean_aber=True)
        q_bore, ctime = make_bore_and_ctime(qm)
        qm.init_point(q_bore, ctime=ctime)
        assert qm.point_is_init()
        qm.reset_point()
        assert not qm.point_is_init()


# ---------------------------------------------------------------------------
# dest_is_pol / source_is_pol
# ---------------------------------------------------------------------------


class TestPolFlags:
    def test_dest_is_pol_false_for_T(self):
        qm = qpoint.QMap(nside=NSIDE, pol=False, mean_aber=True)
        assert not qm.dest_is_pol()

    def test_dest_is_pol_true_for_TQU(self):
        qm = qpoint.QMap(nside=NSIDE, pol=True, mean_aber=True)
        assert qm.dest_is_pol()

    def test_dest_is_vpol(self):
        qm = qpoint.QMap(nside=NSIDE, vpol=True, mean_aber=True)
        assert qm.dest_is_vpol()

    def test_source_is_pol(self):
        source = np.ones((3, NPIX))
        qm = qpoint.QMap(source_map=source, source_pol=True, mean_aber=True)
        assert qm.source_is_pol()

    def test_source_is_not_pol(self):
        source = np.ones((1, NPIX))
        qm = qpoint.QMap(source_map=source, source_pol=False, mean_aber=True)
        assert not qm.source_is_pol()


# ---------------------------------------------------------------------------
# Detector properties: flags, weights, gains
# ---------------------------------------------------------------------------


def make_mapper(ndet=1, nside=NSIDE, pol=True, **kwargs):
    """A pointed QMap and a (ndet, 4) offset array, ready for from_tod."""
    qm = qpoint.QMap(nside=nside, pol=pol, mean_aber=True, **kwargs)
    q_bore, ctime = make_bore_and_ctime(qm)
    qm.init_point(q_bore, ctime=ctime)
    delta = np.arange(ndet, dtype=float)
    q_off = np.atleast_2d(qm.det_offset(1.0 + delta, 2.0 + delta, 30.0 * delta))
    return qm, q_off


class TestFlags:
    """Flagged samples are dropped, which nothing else here exercises."""

    def test_flagged_samples_do_not_accumulate(self):
        qm, q_off = make_mapper()
        _, proj = qm.from_tod(q_off, tod=np.ones((1, N)))
        unflagged = proj[0].sum()

        qm, q_off = make_mapper()
        flag = np.zeros((1, N), dtype=np.uint8)
        flag[0, :10] = 1
        _, proj = qm.from_tod(q_off, tod=np.ones((1, N)), flag=flag)
        assert proj[0].sum() == unflagged - 10

    def test_flagging_everything_leaves_the_map_empty(self):
        qm, q_off = make_mapper()
        flag = np.ones((1, N), dtype=np.uint8)
        vec, proj = qm.from_tod(q_off, tod=np.ones((1, N)), flag=flag)
        assert not np.any(proj)
        assert not np.any(vec)


class TestWeightsAndGain:
    def test_per_channel_weight_scales_the_map(self):
        qm, q_off = make_mapper()
        vec, proj = qm.from_tod(q_off, tod=np.ones((1, N)))
        v1, p1 = vec.copy(), proj.copy()

        qm, q_off = make_mapper()
        vec, proj = qm.from_tod(q_off, tod=np.ones((1, N)), weight=2.5)
        assert np.allclose(vec, 2.5 * v1)
        assert np.allclose(proj, 2.5 * p1)

    def test_per_sample_weights_scale_the_map(self):
        qm, q_off = make_mapper()
        vec, proj = qm.from_tod(q_off, tod=np.ones((1, N)))
        v1, p1 = vec.copy(), proj.copy()

        qm, q_off = make_mapper()
        weights = np.full((1, N), 3.0)
        vec, proj = qm.from_tod(q_off, tod=np.ones((1, N)), weights=weights)
        assert np.allclose(vec, 3.0 * v1)
        assert np.allclose(proj, 3.0 * p1)

    def test_gain_scales_the_signal_but_not_the_hits(self):
        qm, q_off = make_mapper()
        vec, proj = qm.from_tod(q_off, tod=np.ones((1, N)))
        v1, p1 = vec.copy(), proj.copy()

        qm, q_off = make_mapper()
        vec, proj = qm.from_tod(q_off, tod=np.ones((1, N)), gain=4.0)
        assert np.allclose(vec, 4.0 * v1)
        assert np.allclose(proj, p1)


# ---------------------------------------------------------------------------
# Pair differencing
# ---------------------------------------------------------------------------


def make_pair():
    """A pointed QMap and two offsets 90 degrees apart in polarization."""
    qm = qpoint.QMap(nside=NSIDE, pol=True, mean_aber=True)
    q_bore, ctime = make_bore_and_ctime(qm)
    qm.init_point(q_bore, ctime=ctime)
    q_off = qm.det_offset([1.0, 1.0], [2.0, 2.0], [0.0, 90.0])
    return qm, q_off


class TestPairDifference:
    """
    do_diff pairs the first half of the detectors with the second half and
    accumulates their difference into the polarization rows and their sum
    into temperature. Nothing else in this suite reaches that kernel.
    """

    def test_common_mode_cancels_in_polarization(self):
        qm, q_off = make_pair()
        vec, _ = qm.from_tod(q_off, tod=np.ones((2, N)), do_diff=True)
        assert not np.any(vec[1])
        assert not np.any(vec[2])

    def test_common_mode_survives_in_temperature(self):
        qm, q_off = make_pair()
        vec, _ = qm.from_tod(q_off, tod=np.ones((2, N)), do_diff=True)
        assert np.isclose(vec[0].sum(), N)

    def test_a_differential_signal_appears_in_polarization(self):
        qm, q_off = make_pair()
        tod = np.vstack([np.ones(N), -np.ones(N)])
        vec, _ = qm.from_tod(q_off, tod=tod, do_diff=True)
        assert not np.any(vec[0])
        assert np.any(vec[1])
        assert np.any(vec[2])

    def test_one_hit_per_sample_not_per_detector(self):
        qm, q_off = make_pair()
        _, proj = qm.from_tod(q_off, tod=np.ones((2, N)), do_diff=True)
        assert proj[0].sum() == N

    def test_a_flag_on_either_detector_skips_the_sample(self):
        """
        The pair is dropped whichever half carries the flag. The kernel
        used to test flag_init on either detector and then read both
        flag arrays, so this also pins that each is read behind its own.
        """
        for det in (0, 1):
            qm, q_off = make_pair()
            flag = np.zeros((2, N), dtype=np.uint8)
            flag[det, :10] = 1
            _, proj = qm.from_tod(q_off, tod=np.ones((2, N)), do_diff=True, flag=flag)
            assert proj[0].sum() == N - 10


# ---------------------------------------------------------------------------
# Partial maps
# ---------------------------------------------------------------------------


def hit_pixels(count=5):
    """The first `count` pixels a single detector actually lands on."""
    qm, q_off = make_mapper()
    _, proj = qm.from_tod(q_off, tod=np.ones((1, N)))
    return np.nonzero(proj[0])[0][:count].astype(np.int64), proj[0].copy()


class TestPartialDest:
    """
    A destination map can cover a list of pixels rather than the sphere,
    which routes every lookup through the pixel hash.
    """

    def test_shapes_follow_the_pixel_list(self):
        pixels, _ = hit_pixels()
        qm = qpoint.QMap(mean_aber=True, error_missing=False)
        qm.init_dest(nside=NSIDE, pol=True, pixels=pixels)
        q_bore, ctime = make_bore_and_ctime(qm)
        qm.init_point(q_bore, ctime=ctime)
        vec, proj = qm.from_tod(
            np.atleast_2d(qm.det_offset(1.0, 2.0, 0.0)), tod=np.ones((1, N))
        )
        assert vec.shape == (3, len(pixels))
        assert proj.shape == (6, len(pixels))

    def test_hits_match_the_full_sky_map_on_those_pixels(self):
        pixels, full_hits = hit_pixels()
        qm = qpoint.QMap(mean_aber=True, error_missing=False)
        qm.init_dest(nside=NSIDE, pol=True, pixels=pixels)
        q_bore, ctime = make_bore_and_ctime(qm)
        qm.init_point(q_bore, ctime=ctime)
        _, proj = qm.from_tod(
            np.atleast_2d(qm.det_offset(1.0, 2.0, 0.0)), tod=np.ones((1, N))
        )
        assert np.allclose(proj[0], full_hits[pixels])

    def test_a_sample_off_the_map_raises_by_default(self):
        pixels, _ = hit_pixels()
        qm = qpoint.QMap(mean_aber=True)
        qm.init_dest(nside=NSIDE, pol=True, pixels=pixels)
        q_bore, ctime = make_bore_and_ctime(qm)
        qm.init_point(q_bore, ctime=ctime)
        with pytest.raises(RuntimeError, match="out of bounds"):
            qm.from_tod(
                np.atleast_2d(qm.det_offset(1.0, 2.0, 0.0)), tod=np.ones((1, N))
            )

    def test_error_missing_false_drops_it_instead(self):
        pixels, full_hits = hit_pixels()
        qm = qpoint.QMap(mean_aber=True, error_missing=False)
        qm.init_dest(nside=NSIDE, pol=True, pixels=pixels)
        q_bore, ctime = make_bore_and_ctime(qm)
        qm.init_point(q_bore, ctime=ctime)
        _, proj = qm.from_tod(
            np.atleast_2d(qm.det_offset(1.0, 2.0, 0.0)), tod=np.ones((1, N))
        )
        assert proj[0].sum() == full_hits[pixels].sum()
        assert proj[0].sum() < N


# ---------------------------------------------------------------------------
# Map modes
# ---------------------------------------------------------------------------

# (row count, pol) -> the vec mode the row count selects
MODE_ROWS = [
    (1, False, "TEMP"),
    (3, True, "POL"),
    (3, False, "D1"),
    (4, True, "VPOL"),
    (6, False, "D2"),
    (9, True, "D1_POL"),
    (18, True, "D2_POL"),
]


class TestMapModes:
    """
    The row count of a source map selects which reader map2tod uses, and
    the readers index rows directly -- a mode that claims fewer rows than
    its reader touches walks off the end of the map.
    """

    @pytest.mark.parametrize("nrow, pol, name", MODE_ROWS)
    def test_row_count_is_accepted(self, nrow, pol, name):
        qm = qpoint.QMap(mean_aber=True)
        q_bore, ctime = make_bore_and_ctime(qm)
        qm.init_point(q_bore, ctime=ctime)
        qm.init_source(np.zeros((nrow, NPIX)), pol=pol, vpol=(nrow == 4))
        assert qm._source.contents.num_vec == nrow

    @pytest.mark.parametrize("nrow, pol, name", MODE_ROWS)
    def test_every_row_reaches_the_tod(self, nrow, pol, name):
        """
        Perturb each row in turn and require the timestream to change. A
        row that never matters means the mode was inferred too small, and
        the rows past it are read from beyond the map.
        """
        rng = np.random.default_rng(0)
        base = rng.normal(size=(nrow, NPIX))

        def tod_for(source):
            qm = qpoint.QMap(mean_aber=True)
            q_bore, ctime = make_bore_and_ctime(qm)
            qm.init_point(q_bore, ctime=ctime)
            qm.init_source(source, pol=pol, vpol=(nrow == 4))
            q_off = np.atleast_2d(qm.det_offset(1.0, 2.0, 0.0))
            return np.asarray(qm.to_tod(q_off)).copy()

        reference = tod_for(base.copy())
        for row in range(nrow):
            bumped = base.copy()
            bumped[row] += 5.0
            assert not np.array_equal(
                tod_for(bumped), reference
            ), f"row {row} of {name}"


class TestMissingPixelsInToTod:
    def test_nan_missing_marks_samples_off_the_map(self):
        pixels, _ = hit_pixels()
        qm = qpoint.QMap(mean_aber=True, error_missing=False, nan_missing=True)
        q_bore, ctime = make_bore_and_ctime(qm)
        qm.init_point(q_bore, ctime=ctime)
        qm.init_source(np.ones((3, len(pixels))), pol=True, nside=NSIDE, pixels=pixels)
        tod = np.asarray(qm.to_tod(np.atleast_2d(qm.det_offset(1.0, 2.0, 0.0))))
        assert np.isnan(tod).any()
        assert np.isfinite(tod).any()

    def test_without_nan_missing_they_are_left_at_zero(self):
        pixels, _ = hit_pixels()
        qm = qpoint.QMap(mean_aber=True, error_missing=False, nan_missing=False)
        q_bore, ctime = make_bore_and_ctime(qm)
        qm.init_point(q_bore, ctime=ctime)
        qm.init_source(np.ones((3, len(pixels))), pol=True, nside=NSIDE, pixels=pixels)
        tod = np.asarray(qm.to_tod(np.atleast_2d(qm.det_offset(1.0, 2.0, 0.0))))
        assert not np.isnan(tod).any()
        assert np.count_nonzero(tod)


class TestNumThreads:
    """
    Threading is over detectors, each accumulating into its own map before the
    merge, which happens in whatever order the threads finish -- so only
    quantities that sum exactly come back identical. Where OpenMP is off there
    is one thread and everything is trivially identical.
    """

    def run(self, num_threads):
        rng = np.random.default_rng(1)
        qm, q_off = make_mapper(ndet=4, num_threads=num_threads)
        vec, proj = qm.from_tod(q_off, tod=rng.normal(size=(4, N)))
        return np.asarray(vec).copy(), np.asarray(proj).copy()

    def test_the_hit_count_per_pixel_is_exact(self):
        """Whole hits, so the sum is exact whatever order they arrive in."""
        (_, one), (_, four) = self.run(1), self.run(4)
        assert np.array_equal(one[0], four[0])

    def test_the_total_weight_is_exact(self):
        (_, one), (_, four) = self.run(1), self.run(4)
        assert one[0].sum() == four[0].sum() == 4 * N

    def test_the_rest_agrees_to_rounding(self):
        """
        The polarization cross terms and the signal are sums of products,
        which reassociate: identical to about 1e-15 here, checked well
        inside that so the test is not a rounding tripwire.
        """
        (vec_one, one), (vec_four, four) = self.run(1), self.run(4)
        assert np.allclose(one, four, rtol=1e-10, atol=1e-12)
        assert np.allclose(vec_one, vec_four, rtol=1e-10, atol=1e-12)

    def test_the_build_is_threaded_where_it_should_be(self, omp_runtime):
        """
        Everything above passes on a serial build, because one thread
        agrees with itself, so none of it would notice parallelism going
        missing. On the platform that is supposed to have it, say so.

        Linux only, and only under CI: Apple clang ships no OpenMP
        runtime, and linking LLVM's collides with the one healpy bundles,
        so macOS is deliberately serial. A developer's Linux box without
        libgomp is their business; the runners are not.
        """
        if not (os.environ.get("CI") and sys.platform.startswith("linux")):
            pytest.skip("only pinned for CI on Linux")
        assert omp_runtime, "expected an OpenMP build; the threading tests are vacuous"


# ---------------------------------------------------------------------------
# Polarization convention, and the solver that scipy backs
# ---------------------------------------------------------------------------


class TestPolConv:
    """
    polconv picks the sign convention for U. Elsewhere it is only set and
    read back; this is where it shows up, which is in the map rather than
    in bore2radec -- the kernel negates beta, the pointing does not.
    """

    def maps(self, polconv):
        qm, q_off = make_mapper(polconv=polconv)
        vec, _ = qm.from_tod(q_off, tod=np.ones((1, N)))
        return np.asarray(vec).copy()

    def test_u_flips_between_the_conventions(self):
        assert np.allclose(self.maps("cosmo")[2], -self.maps("iau")[2])

    def test_t_and_q_do_not(self):
        cosmo, iau = self.maps("cosmo"), self.maps("iau")
        assert np.allclose(cosmo[0], iau[0])
        assert np.allclose(cosmo[1], iau[1])


class TestSolveMapCho:
    """
    solve_map_cho is the Cholesky path and had no test of its own. It
    excludes the same pixels as the default solver -- cho_factor does not
    fail on a rank-deficient matrix, it returns nonsense, so the condition
    number is what keeps a one- or two-hit pixel out of both of them.
    """

    def solved(self, pol):
        qm, q_off = make_mapper(ndet=4, pol=pol)
        rng = np.random.default_rng(2)
        vec, proj = qm.from_tod(q_off, tod=rng.normal(size=(4, N)))
        direct, mask = qm.solve_map(vec=vec.copy(), proj=proj.copy(), return_mask=True)
        cho = qm.solve_map_cho(vec=vec.copy(), proj=proj.copy())
        return np.asarray(direct), np.asarray(cho), np.asarray(mask, dtype=bool)

    def test_agrees_with_the_default_solver_temperature(self):
        direct, cho, mask = self.solved(pol=False)
        assert mask.any()
        assert np.allclose(direct, cho)

    def test_agrees_with_the_default_solver_polarized(self):
        pytest.importorskip("scipy")
        direct, cho, mask = self.solved(pol=True)
        assert mask.any()
        assert np.allclose(direct, cho)

    def test_both_zero_the_pixels_they_cannot_determine(self):
        pytest.importorskip("scipy")
        direct, cho, mask = self.solved(pol=True)
        assert (~mask).any()
        assert not np.any(direct[:, ~mask])
        assert not np.any(cho[:, ~mask])

    @pytest.mark.parametrize("fault", ["singular", "not finite"])
    def test_a_pixel_the_factorization_refuses_is_masked(self, fault):
        """
        The two faults cho_factor reports, which it does not report alike: a
        matrix that is not positive definite raises LinAlgError, one holding an
        inf or a nan raises ValueError. Both are caught narrowly, so a third kind
        stays a bug rather than a quietly dropped pixel.

        cond is supplied rather than computed, or the conditioning test would
        remove these pixels before the solver saw them.
        """
        pytest.importorskip("scipy")
        qm = qpoint.QMap(nside=NSIDE, pol=True)
        proj = np.zeros((6, NPIX))
        proj[0], proj[3], proj[5] = 4.0, 2.0, 2.0
        proj[1], proj[2], proj[4] = 0.1, 0.1, 0.1
        vec = np.random.default_rng(0).normal(size=(3, NPIX))
        bad = 7
        if fault == "singular":
            proj[:, bad] = [4.0, 4.0, 0.0, 4.0, 0.0, 0.0]
        else:
            proj[1, bad] = np.nan

        out, mask = qm.solve_map(
            vec=vec.copy(),
            proj=proj.copy(),
            method="cho",
            cond=np.zeros(NPIX),
            return_mask=True,
        )
        out, mask = np.asarray(out), np.asarray(mask, dtype=bool)
        assert not mask[bad], "the unsolvable pixel should come back masked"
        assert not np.any(out[:, bad])
        assert mask.sum() == NPIX - 1, "and it should not take the others with it"


class TestSolversLeaveTheInputAlone:
    """
    None of the solvers writes to a proj it does not hand back.

    That is what lets them skip copying it: a full-sky proj is 150 MB at
    nside 512, and the solvers read a fraction of a percent of it. A
    write introduced on one of those paths would corrupt the caller's
    array instead of a private copy, and do it silently, so it is pinned
    here rather than left to the reader.
    """

    def inputs(self, npix=NPIX, nmap=3):
        rng = np.random.default_rng(0)
        nproj = nmap * (nmap + 1) // 2
        proj = np.abs(rng.normal(size=(nproj, npix))) + 1.0
        # a realistic footprint: most of the map never hit
        proj[:, rng.permutation(npix)[: npix // 2]] = 0.0
        return rng.normal(size=(nmap, npix)), proj

    @pytest.mark.parametrize("kwargs", [{}, {"method": "cho"}])
    def test_solve_map(self, kwargs):
        pytest.importorskip("scipy") if kwargs else None
        vec, proj = self.inputs()
        before = proj.copy()
        qpoint.QMap(nside=NSIDE, pol=True).solve_map(vec=vec, proj=proj, **kwargs)
        assert np.array_equal(proj, before)

    def test_proj_cond(self):
        _, proj = self.inputs()
        before = proj.copy()
        qpoint.QMap(nside=NSIDE, pol=True).proj_cond(proj=proj)
        assert np.array_equal(proj, before)

    def test_unsolve_map(self):
        map_in, proj = self.inputs()
        before = proj.copy()
        qpoint.QMap(nside=NSIDE, pol=True).unsolve_map(map_in=map_in, proj=proj)
        assert np.array_equal(proj, before)

    def test_returned_proj_is_still_written(self):
        """
        The other half: asking for it back does give the solved form. The
        Cholesky path is the one that rewrites proj, replacing each hit
        pixel's matrix with its decomposition.
        """
        pytest.importorskip("scipy")
        vec, proj = self.inputs()
        before = proj.copy()
        _, out = qpoint.QMap(nside=NSIDE, pol=True).solve_map(
            vec=vec, proj=proj, return_proj=True, method="cho"
        )
        assert np.array_equal(proj, before), "the input must still be intact"
        assert not np.array_equal(np.asarray(out), before)


class TestCtimeIsRequiredWithoutMeanAber:
    """
    With mean_aber off, aberration is applied per detector and needs the time
    of each sample, so the kernels refuse to run without it. Each checks
    separately -- the binned and differenced paths through tod2map, and
    map2tod -- and the message says which one refused.
    """

    def mapper(self, ndet=1, with_ctime=False):
        """
        A QMap with mean_aber off, pointed with or without ctime.
        make_mapper cannot serve here: it fixes mean_aber=True.
        """
        qm = qpoint.QMap(nside=NSIDE, pol=True, mean_aber=False)
        q_bore, ctime = make_bore_and_ctime(qm)
        qm.init_point(q_bore, ctime=ctime if with_ctime else None)
        delta = np.arange(ndet, dtype=float)
        q_off = np.atleast_2d(qm.det_offset(1.0 + delta, 2.0 + delta, 30.0 * delta))
        return qm, q_off

    def test_from_tod(self):
        qm, q_off = self.mapper()
        with pytest.raises(RuntimeError, match="ctime required"):
            qm.from_tod(q_off, tod=np.ones((1, N)))

    def test_from_tod_differenced(self):
        """The differencing kernel has its own copy of the check."""
        qm, q_off = self.mapper(ndet=2)
        with pytest.raises(RuntimeError, match="ctime required"):
            qm.from_tod(q_off, tod=np.ones((2, N)), do_diff=True)

    def test_to_tod(self):
        qm, q_off = self.mapper()
        qm.init_source(np.zeros((3, NPIX)), pol=True)
        with pytest.raises(RuntimeError, match="ctime required"):
            qm.to_tod(q_off)

    def test_ctime_makes_all_three_work(self):
        """
        The control: the same calls succeed once ctime is supplied, so
        the tests above are pinning the missing time and not some other
        unfinished setup.
        """
        qm, q_off = self.mapper(ndet=2, with_ctime=True)
        qm.from_tod(q_off, tod=np.ones((2, N)))
        qm, q_off = self.mapper(ndet=2, with_ctime=True)
        qm.from_tod(q_off, tod=np.ones((2, N)), do_diff=True)
        qm, q_off = self.mapper(with_ctime=True)
        qm.init_source(np.zeros((3, NPIX)), pol=True)
        assert np.asarray(qm.to_tod(q_off)).shape == (1, N)


class TestDetarrLifetime:
    """
    from_tod builds the detector array and tears it down again, so the
    structure is only alive between init_detarr and reset_detarr.
    """

    def test_from_tod_leaves_no_detector_array_behind(self):
        qm, q_off = make_mapper()
        qm.from_tod(q_off, tod=np.ones((1, N)))
        assert qm._detarr is None

    def test_init_detarr_then_reset(self):
        qm, q_off = make_mapper()
        qm.init_detarr(q_off, tod=np.ones((1, N)))
        assert qm._detarr is not None
        assert "tod" in qm.depo
        qm.reset_detarr()
        assert qm._detarr is None
        assert "tod" not in qm.depo

    def test_reset_is_safe_when_there_is_nothing_to_reset(self):
        qm, _ = make_mapper()
        qm.reset_detarr()
        qm.reset_detarr()

    def test_mapping_accumulates_across_calls(self):
        qm, q_off = make_mapper()
        _, proj = qm.from_tod(q_off, tod=np.ones((1, N)))
        first = proj[0].sum()
        _, proj = qm.from_tod(q_off, tod=np.ones((1, N)))
        assert proj[0].sum() == 2 * first


# ---------------------------------------------------------------------------
# The depo
# ---------------------------------------------------------------------------


DEPO_KEYS = {
    "vec",
    "proj",
    "dest_nside",
    "dest_pixels",
    "source_map",
    "source_nside",
    "source_pixels",
    "q_bore",
    "ctime",
    "q_hwp",
    "tod",
    "flag",
    "weights",
}


class TestDepo:
    """
    Every array QMap hands to C goes into a depo dict, which is also what
    keeps those buffers alive: the ctypes layer takes a bare pointer and
    holds no reference. Thirteen keys, each appearing and disappearing with
    the structure it belongs to. spider_tools subclasses QMap and reads
    them, so the set and its presence rules are interface.
    """

    def test_the_value_is_the_callers_array(self):
        """Not a copy, in either package."""
        vec = np.zeros((3, NPIX))
        qm = qpoint.QMap(mean_aber=True)
        qm.init_dest(nside=NSIDE, vec=vec, pol=True)
        assert qm.depo["vec"] is vec

    def test_a_switched_off_component_is_False_not_absent(self):
        """
        qpoint stores the False sentinel rather than dropping the key, and
        its from_tod tests for it, so the key has to be there either way.
        """
        qm = qpoint.QMap(mean_aber=True)
        qm.init_dest(nside=NSIDE, vec=False, proj=np.zeros((1, NPIX)))
        assert "vec" in qm.depo
        assert qm.depo["vec"] is False

    def test_a_held_depo_follows_a_reinstall(self):
        """
        The liveness property: reset_dest replaces what the depo holds, so
        a caller that cached the dict still sees the new arrays.
        """
        depo = qm_depo = qpoint.QMap(mean_aber=True)
        qm_depo.init_dest(nside=NSIDE, pol=True)
        depo = qm_depo.depo
        other = np.ones((3, NPIX))
        qm_depo.init_dest(nside=NSIDE, vec=other, pol=True, reset=True)
        assert depo["vec"] is other

    @pytest.mark.parametrize(
        "reset, gone",
        [
            ("reset_dest", ("vec", "proj", "dest_nside")),
            ("reset_source", ("source_map", "source_nside")),
            ("reset_point", ("q_bore", "ctime")),
        ],
    )
    def test_each_reset_drops_its_own_keys(self, reset, gone):
        qm, _ = make_mapper()
        qm.init_source(np.ones((3, NPIX)), pol=True)
        for key in gone:
            assert key in qm.depo, key
        getattr(qm, reset)()
        for key in gone:
            assert key not in qm.depo, key

    def test_empty_after_a_full_reset(self):
        qm, _ = make_mapper()
        qm.init_source(np.ones((3, NPIX)), pol=True)
        assert len(qm.depo)
        qm.reset()
        assert len(qm.depo) == 0

    def test_a_partial_map_carries_its_pixel_list(self):
        pixels = np.arange(NPIX // 2, dtype=np.int64)
        qm = qpoint.QMap(mean_aber=True)
        qm.init_dest(nside=NSIDE, pixels=pixels, pol=True)
        assert qm.depo["dest_pixels"] is pixels
        qm = qpoint.QMap(mean_aber=True)
        qm.init_dest(nside=NSIDE, pol=True)
        assert "dest_pixels" not in qm.depo

    def test_the_detector_arrays_come_and_go_with_the_call(self):
        qm, q_off = make_mapper()
        for key in ("tod", "flag", "weights"):
            assert key not in qm.depo
        qm.init_detarr(
            q_off,
            tod=np.ones((1, N)),
            flag=np.zeros((1, N), dtype=np.uint8),
            weights=np.ones((1, N)),
        )
        for key in ("tod", "flag", "weights"):
            assert key in qm.depo, key
        qm.reset_detarr()
        for key in ("tod", "flag", "weights"):
            assert key not in qm.depo, key

    def test_only_a_supplied_flag_appears(self):
        qm, q_off = make_mapper()
        qm.init_detarr(q_off, tod=np.ones((1, N)))
        assert "tod" in qm.depo
        assert "flag" not in qm.depo
        assert "weights" not in qm.depo

    def test_assigning_proj_reaches_the_solver(self):
        """
        The solvers read the depo, so assigning into it changes what they
        use -- which is what makes the dict part of the interface rather
        than a record of it.
        """
        qm, q_off = make_mapper(pol=False)
        qm.from_tod(q_off, tod=np.ones((1, N)))
        hits = np.asarray(qm.depo["proj"]).copy()
        qm.depo["proj"] = 4.0 * hits
        assert np.array_equal(np.asarray(qm.depo["proj"]), 4.0 * hits)
        solved = qm.solve_map()
        assert np.isfinite(solved).any()

    def test_the_key_set_cannot_grow(self):
        qm, q_off = make_mapper()
        qm.init_source(np.ones((3, NPIX)), pol=True)
        qm.init_detarr(q_off, tod=np.ones((1, N)))
        assert set(qm.depo) <= DEPO_KEYS

    def test_arbitrary_keys_are_allowed(self):
        """
        Callers keep their own state here -- spider_tools' UnifileMap has
        24 keys of its own alongside qpoint's 13 -- so the depo has to take
        a key it has never heard of, and pop it again.
        """
        qm, _ = make_mapper()
        qm.depo["dest_coord"] = "C"
        assert qm.depo["dest_coord"] == "C"
        assert "dest_coord" in qm.depo
        assert qm.depo.get("nope") is None
        cache = qm.depo["source_cache"] = {}
        cache["a"] = 1
        assert qm.depo["source_cache"] is cache
        assert qm.depo.pop("dest_coord") == "C"
        assert qm.depo.pop("dest_coord", "dflt") == "dflt"

    def test_an_extra_survives_a_partial_reset_and_not_a_full_one(self):
        qm, _ = make_mapper()
        qm.depo["dest_coord"] = "C"
        qm.reset_dest()
        assert qm.depo["dest_coord"] == "C"
        qm.reset()
        assert "dest_coord" not in qm.depo


class TestDepoDownstreamPatterns:
    """
    The operations spider_tools' UnifileMap and MPIUnifileMap actually
    perform on the depo, replayed here so the contract they depend on is
    pinned in this repo rather than only in theirs.
    """

    def test_the_reduce_dest_cycle(self):
        """
        MPIUnifileMap.reduce_dest reads vec, proj, dest_nside and
        dest_pixels out of the depo, reduces them across ranks and re-inits
        the dest from the result. Without MPI that is the same sequence
        with the reduction as identity.
        """
        # More pixels than proj rows: qpoint still flips a map taller than
        # it is wide, and a 5-pixel TQU proj is (6, 5), which it would read
        # as 5 bogus rows. That divergence has nothing to do with the depo.
        pixels, _ = hit_pixels(count=20)
        qm, q_off = make_mapper(pol=True, error_missing=False)
        qm.init_dest(nside=NSIDE, pol=True, pixels=pixels, reset=True)
        qm.from_tod(q_off, tod=np.ones((1, N)))

        nside = qm.depo["dest_nside"]
        pix = qm.depo.get("dest_pixels", None)
        assert qm.depo["vec"] is not False
        assert qm.depo["proj"] is not False
        vec = np.array(qm.depo["vec"], copy=True)
        proj = np.array(qm.depo["proj"], copy=True)
        zeroed = np.zeros_like(qm.depo["vec"])

        qm.init_dest(nside=nside, vec=vec, proj=proj, pixels=pix, pol=True, reset=True)
        assert np.array_equal(np.asarray(qm.depo["vec"]), vec)
        assert np.shape(zeroed) == np.shape(vec)
        assert qm.depo["dest_nside"] == nside

    def test_the_namespace_it_keeps_alongside(self):
        """
        24 of the 37 keys UnifileMap uses are its own, reached with get,
        pop and `in` rather than plain subscripting.
        """
        extras = {
            "dest_coord": "C",
            "source_coord": "G",
            "point_coord": "C",
            "dest_label": "map01",
            "source_label": "sim",
            "dest_reduced": True,
            "poly_order": 3,
            "poly_scans": [1, 2],
            "turn_idx": np.arange(4),
            "data_range": (0, 100),
        }
        qm, _ = make_mapper()
        for key, value in extras.items():
            qm.depo[key] = value
        for key, value in extras.items():
            assert key in qm.depo
            assert qm.depo.get(key) is value or qm.depo[key] == value
        assert qm.depo.get("point_shift") is None
        for key in extras:
            qm.depo.pop(key, None)
            assert key not in qm.depo

    def test_a_nested_cache_is_mutated_in_place(self):
        """
        MPIUnifileMap builds depo['source_cache'] once and then pops
        entries out of it, so the same object has to come back each time.
        """
        qm, _ = make_mapper()
        if "source_cache" not in qm.depo:
            qm.depo["source_cache"] = dict()
        qm.depo["source_cache"]["label"] = {"local": True}
        assert qm.depo["source_cache"]["label"]["local"] is True
        qm.depo["source_cache"].pop("label", None)
        assert qm.depo["source_cache"] == {}

    def test_assigning_proj_after_solving(self):
        """unimap_exec stores the proj it solved with, on the root rank."""
        qm, q_off = make_mapper(pol=False)
        qm.from_tod(q_off, tod=np.ones((1, N)))
        proj = np.atleast_2d(np.asarray(qm.depo["proj"]))
        qm.solve_map()
        qm.depo["proj"] = np.atleast_2d(proj)
        assert np.array_equal(np.asarray(qm.depo["proj"]), proj)

    def test_the_source_map_is_copied_out(self):
        """unimap_exec seeds its PCG from depo['source_map'].copy()."""
        qm, _ = make_mapper()
        source = np.ones((3, NPIX))
        qm.init_source(source, pol=True)
        x = qm.depo["source_map"].copy()
        assert np.array_equal(x, source)
        assert x is not source


# ---------------------------------------------------------------------------
# fast_pix through the mapmaking kernels
# ---------------------------------------------------------------------------


class TestFastPixInMapmaking:
    """
    fast_pix reaches the mapmaking kernels too, and leaves them
    bit-identical to the two-step path it replaces.
    """

    def source_map(self):
        return np.random.default_rng(0).normal(size=(3, NPIX))

    def test_to_tod(self):
        """Through map2tod's kernel."""
        out = []
        for fast in (False, True):
            qm, q_off = make_mapper(fast_pix=fast)
            qm.init_source(self.source_map(), pol=True)
            out.append(np.asarray(qm.to_tod(q_off)))
        assert np.array_equal(out[0], out[1]), "to_tod fast vs two-step"

    def test_from_tod(self):
        """And through tod2map's."""
        tod = np.random.default_rng(1).normal(size=(1, N))
        out = []
        for fast in (False, True):
            qm, q_off = make_mapper(fast_pix=fast)
            out.append([np.asarray(x) for x in qm.from_tod(q_off, tod=tod.copy())])
        for a, b, name in zip(out[0], out[1], ("vec", "proj")):
            assert np.array_equal(a, b), "from_tod {} fast vs two-step".format(name)


# Resolutions for the derivative round trip below: a low-resolution map
# carrying derivatives, reobserved at the pixel centres of a finer one.
DERIV_NSIDE_LO = 64
DERIV_NSIDE_HI = 256
DERIV_LMAX = 128


@pytest.fixture(scope="module")
def deriv_sky():
    """A smooth band-limited sky at both resolutions, with derivatives."""
    hp = pytest.importorskip("healpy")

    ell = np.arange(DERIV_LMAX + 1)
    cl = np.zeros(DERIV_LMAX + 1)
    cl[2:] = 1.0 / ell[2:] ** 2
    # smoothed to a couple of low-resolution pixels, so the field is actually
    # resolved there and the expansion has something to do
    bl = hp.gauss_beam(2.0 * hp.nside2resol(DERIV_NSIDE_LO), lmax=DERIV_LMAX)
    np.random.seed(0)
    alm = hp.almxfl(hp.synalm(cl, lmax=DERIV_LMAX, new=True), bl)

    m, dt, dp = hp.alm2map_der1(alm, nside=DERIV_NSIDE_LO, lmax=DERIV_LMAX)
    _, dtdt, dtdp = hp.alm2map_der1(
        hp.map2alm(dt, lmax=DERIV_LMAX), nside=DERIV_NSIDE_LO, lmax=DERIV_LMAX
    )
    _, _, dpdp = hp.alm2map_der1(
        hp.map2alm(dp, lmax=DERIV_LMAX), nside=DERIV_NSIDE_LO, lmax=DERIV_LMAX
    )
    source = np.ascontiguousarray([m, dt, dp, dtdt, dtdp, dpdp])

    truth = hp.alm2map(alm, DERIV_NSIDE_HI, lmax=DERIV_LMAX)
    theta, phi = hp.pix2ang(DERIV_NSIDE_HI, np.arange(len(truth)))
    # every fourth fine pixel is plenty, and keeps this quick
    sub = np.arange(0, len(truth), 4)
    return source, truth[sub], theta[sub], phi[sub]


def reobserve_deriv(sky, rows, pol=False, **kwargs):
    """TOD from the first `rows` columns of `sky`, at the fine pixel centres."""
    source, _, theta, phi = sky
    qm = qpoint.QMap(nside=DERIV_NSIDE_LO, pol=False, **kwargs)
    q_bore = np.ascontiguousarray(
        np.atleast_2d(
            qm.radecpa2quat(
                np.rad2deg(phi), 90.0 - np.rad2deg(theta), np.zeros(len(theta))
            )
        )
    )
    qm.init_point(q_bore)
    qm.init_source(np.ascontiguousarray(source[:rows]), pol=pol)
    return np.asarray(qm.to_tod(qm.det_offset(0.0, 0.0, 0.0)))[0]


def rms(x):
    return float(np.sqrt(np.mean(np.asarray(x) ** 2)))


def deriv_series(sky):
    """The Taylor series map2tod should be evaluating, written out here."""
    hp = pytest.importorskip("healpy")
    source, _, theta, phi = sky
    pix = hp.ang2pix(DERIV_NSIDE_LO, theta, phi)
    th_pix, ph_pix = hp.pix2ang(DERIV_NSIDE_LO, pix)
    dth = theta - th_pix
    dph = (phi - ph_pix + np.pi) % (2.0 * np.pi) - np.pi
    arc = dph * np.sin(th_pix)
    return (
        source[0][pix]
        + dth * source[1][pix]
        + arc * source[2][pix]
        + 0.5 * dth**2 * source[3][pix]
        + dth * arc * source[4][pix]
        + 0.5 * arc**2 * source[5][pix]
    )


class TestDerivativeMapsRecoverTheSky:
    """
    map2tod evaluates a derivative source map as a Taylor series about the
    pixel centre, so a low-resolution map carrying its own derivatives should
    reproduce the same sky sampled at higher resolution.

    The derivative columns are gradient components, as `alm2map_der1`
    returns them: the phi derivatives are divided by ``sin(theta)``, so the
    offset that multiplies them is an arc length, and the second derivatives
    are derivatives rather than series coefficients, so the pure second-order
    terms carry a half.  Getting either wrong leaves an error the size of the
    term it belongs to, which shows up here as D2 failing to improve on D1.

    The interpolated modes are measured on the same sky, as the alternative
    way of evaluating a map away from a pixel centre.
    """

    def test_each_order_improves_on_the_last(self, deriv_sky):
        """
        T, then first derivatives, then second: each order should cut the
        residual by a factor of a few.  A D2 map that does no better than D1
        is the signature of a mis-weighted second-order term.
        """
        _, truth, _, _ = deriv_sky
        err = {n: rms(reobserve_deriv(deriv_sky, n) - truth) for n in (1, 3, 6)}
        sky_rms = rms(truth)

        assert err[1] < 0.2 * sky_rms, "pixel-centre sampling is implausibly good"
        assert err[3] < err[1] / 2.0, "first derivatives did not help"
        assert err[6] < err[3] / 2.0, "second derivatives did not help"
        assert err[6] < 0.02 * sky_rms, "D2 residual larger than expected"

    def test_the_expansion_is_the_documented_one(self, deriv_sky):
        """
        Against the series written out here, rather than against the sky, so
        a wrong coefficient cannot hide behind the band limit.
        """
        got = reobserve_deriv(deriv_sky, 6)
        np.testing.assert_allclose(got, deriv_series(deriv_sky), rtol=0, atol=1e-12)

    def test_the_phi_offset_is_an_arc_length(self, deriv_sky):
        """
        Without the sin(theta) the phi terms are wrong by 1/sin(theta), so the
        residual would grow away from the equator.  Here it must not.
        """
        _, truth, theta, _ = deriv_sky
        resid = reobserve_deriv(deriv_sky, 6) - truth
        lat = 90.0 - np.rad2deg(theta)
        near, far = np.abs(lat) < 15.0, np.abs(lat) > 60.0
        assert rms(resid[far]) < 3.0 * rms(resid[near])

    def test_the_polarized_columns_expand_the_same_way(self, deriv_sky):
        """
        An 18-column map with only its temperature columns filled drives the
        D2_POL branch, which carries its own copy of the expansion.
        """
        source = deriv_sky[0]
        wide = np.zeros((18, source.shape[1]))
        wide[::3] = source
        got = reobserve_deriv((wide,) + deriv_sky[1:], 18, pol=True)
        np.testing.assert_allclose(got, deriv_series(deriv_sky), rtol=0, atol=1e-12)
        np.testing.assert_allclose(
            got, reobserve_deriv(deriv_sky, 6), rtol=0, atol=1e-12
        )

    @pytest.mark.parametrize("rows", [1, 3, 4])
    def test_interpolation_recovers_the_sky(self, deriv_sky, rows):
        """
        The interpolated branches, which take four neighbouring pixels rather
        than a series about one centre.  Q, U and V are left at zero, so the
        temperature sky is still the answer, and IDATUM, IPOLDATUM and
        IVPOLDATUM are each driven in turn.
        """
        source, truth, _, _ = deriv_sky
        src = np.zeros((rows, source.shape[1]))
        src[0] = source[0]
        sky = (np.ascontiguousarray(src),) + deriv_sky[1:]

        err = rms(reobserve_deriv(sky, rows, pol=rows > 1, interp_pix=True) - truth)
        centre = rms(reobserve_deriv(sky, rows, pol=rows > 1) - truth)
        assert err < 0.6 * centre, "interpolation did not improve on the centre"
        assert err < 0.1 * rms(truth)

    def test_the_series_beats_interpolation(self, deriv_sky):
        """
        Second derivatives are third order in the offset where interpolation
        is second, so a D2 map should be the better of the two on a sky this
        smooth.  This is the reason the derivative modes exist.
        """
        _, truth, _, _ = deriv_sky
        interp = rms(reobserve_deriv(deriv_sky, 1, interp_pix=True) - truth)
        series = rms(reobserve_deriv(deriv_sky, 6) - truth)
        assert series < 0.5 * interp

    @pytest.mark.parametrize("interp_pix", [False, True])
    def test_interpolation_does_not_reach_these_modes(self, deriv_sky, interp_pix):
        """
        map2tod interpolates for T, POL and VPOL only, so a derivative map
        takes the series either way.  That exclusivity is what lets the
        kernel skip the pixel offset when it is interpolating.
        """
        got = reobserve_deriv(deriv_sky, 6, interp_pix=interp_pix)
        np.testing.assert_allclose(got, deriv_series(deriv_sky), rtol=0, atol=1e-12)

    @pytest.mark.parametrize("fast_math", [False, True])
    def test_fast_math_takes_the_same_path(self, deriv_sky, fast_math):
        """
        The arc length uses poly_sin under fast_math, so both settings have to
        reach the same accuracy.
        """
        _, truth, _, _ = deriv_sky
        err = rms(reobserve_deriv(deriv_sky, 6, fast_math=fast_math) - truth)
        assert err < 0.02 * rms(truth)
