"""Tests for qpoint.QMap mapmaking class and helper functions."""

import os
import sys

import numpy as np
import pytest
import qpoint
from qpoint.qmap_class import nside2npix, npix2nside, check_map, check_proj

qpoint2 = pytest.importorskip("qpoint2")

# Both packages wherever the test is about behaviour. Some reach into the
# wrapper each keeps its maps in, and those differ by design -- qpoint
# holds ctypes structs and a depo dict, qpoint2 holds C++ objects -- so
# those stay with qpoint, and tests/test_parity.py covers the other side.
IMPLS = [qpoint, qpoint2]


@pytest.fixture(params=IMPLS, ids=lambda m: m.__name__)
def mod(request):
    return request.param


# Small nside for fast tests
NSIDE = 8
NPIX = 12 * NSIDE * NSIDE
N = 200  # enough samples to hit several pixels

# Reference observing parameters
CTIME = 1418662800.0
LON = 165.7
LAT = -77.6


def make_bore_and_ctime(mod, qm, n=N):
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
    def test_default_init(self, mod):
        qm = mod.QMap(mean_aber=True)
        assert not qm.dest_is_init()
        assert not qm.source_is_init()
        assert not qm.point_is_init()

    def test_init_with_nside(self, mod):
        qm = mod.QMap(nside=NSIDE, pol=False, mean_aber=True)
        assert qm.dest_is_init()

    def test_init_with_source(self, mod):
        source = np.ones((1, NPIX))
        qm = mod.QMap(source_map=source, source_pol=False, mean_aber=True)
        assert qm.source_is_init()

    def test_init_with_bore(self, mod):
        qm = mod.QMap(mean_aber=True)
        q_bore, ctime = make_bore_and_ctime(mod, qm)
        qm2 = mod.QMap(q_bore=q_bore, ctime=ctime, mean_aber=True)
        assert qm2.point_is_init()


# ---------------------------------------------------------------------------
# init_dest
# ---------------------------------------------------------------------------


class TestInitDest:
    def test_temperature_shapes(self):
        # qpoint only: as above
        qm = qpoint.QMap(mean_aber=True)
        result = qm.init_dest(nside=NSIDE, pol=False)
        if isinstance(result, tuple):
            vec, proj = result
        else:
            vec = result
            proj = qm.depo["proj"].squeeze()
        assert vec.shape == (NPIX,) or vec.shape == (1, NPIX)
        assert proj.shape == (NPIX,) or proj.shape == (1, NPIX)

    def test_polarized_shapes(self):
        # qpoint only: init_dest returns the maps here; qpoint2 returns None and is read with get_vec
        qm = qpoint.QMap(mean_aber=True)
        vec, proj = qm.init_dest(nside=NSIDE, pol=True)
        assert vec.shape == (3, NPIX)
        assert proj.shape == (6, NPIX)

    def test_vpol_shapes(self):
        # qpoint only: as above
        qm = qpoint.QMap(mean_aber=True)
        vec, proj = qm.init_dest(nside=NSIDE, vpol=True)
        assert vec.shape == (4, NPIX)
        assert proj.shape == (10, NPIX)

    def test_init_twice_raises(self, mod):
        qm = mod.QMap(mean_aber=True)
        qm.init_dest(nside=NSIDE)
        with pytest.raises(RuntimeError):
            qm.init_dest(nside=NSIDE)

    def test_reset_and_reinit(self, mod):
        qm = mod.QMap(mean_aber=True)
        qm.init_dest(nside=NSIDE, pol=False)
        qm.reset_dest()
        assert not qm.dest_is_init()
        qm.init_dest(nside=NSIDE, pol=True)
        assert qm.dest_is_init()

    def test_vec_false_proj_false_raises(self, mod):
        qm = mod.QMap(mean_aber=True)
        with pytest.raises(ValueError):
            qm.init_dest(nside=NSIDE, vec=False, proj=False)


# ---------------------------------------------------------------------------
# init_source
# ---------------------------------------------------------------------------


class TestInitSource:
    def test_temperature_map(self, mod):
        qm = mod.QMap(mean_aber=True)
        source = np.ones((1, NPIX))
        qm.init_source(source, pol=False)
        assert qm.source_is_init()
        assert not qm.source_is_pol()

    def test_polarized_map(self, mod):
        qm = mod.QMap(mean_aber=True)
        source = np.ones((3, NPIX))
        qm.init_source(source, pol=True)
        assert qm.source_is_init()
        assert qm.source_is_pol()

    def test_vpol_map(self, mod):
        qm = mod.QMap(mean_aber=True)
        source = np.ones((4, NPIX))
        qm.init_source(source, vpol=True)
        assert qm.source_is_init()
        assert qm.source_is_vpol()

    def test_init_twice_raises(self, mod):
        qm = mod.QMap(mean_aber=True)
        source = np.ones((1, NPIX))
        qm.init_source(source, pol=False)
        with pytest.raises(RuntimeError):
            qm.init_source(source, pol=False)

    def test_reset_and_reinit(self, mod):
        qm = mod.QMap(mean_aber=True)
        source = np.ones((1, NPIX))
        qm.init_source(source, pol=False)
        qm.reset_source()
        assert not qm.source_is_init()
        qm.init_source(source, pol=False)
        assert qm.source_is_init()


# ---------------------------------------------------------------------------
# init_point
# ---------------------------------------------------------------------------


class TestInitPoint:
    def test_init_with_qbore(self, mod):
        qm = mod.QMap(mean_aber=True)
        q_bore, ctime = make_bore_and_ctime(mod, qm)
        qm.init_point(q_bore, ctime=ctime)
        assert qm.point_is_init()

    def test_init_without_point_raises_when_mapping(self, mod):
        qm = mod.QMap(nside=NSIDE, pol=False, mean_aber=True)
        q_off = np.atleast_2d(qm.det_offset(0.0, 0.0, 0.0))
        tod = np.ones((1, N))
        with pytest.raises(Exception):
            qm.from_tod(q_off, tod=tod)

    def test_reset_point(self, mod):
        qm = mod.QMap(mean_aber=True)
        q_bore, ctime = make_bore_and_ctime(mod, qm)
        qm.init_point(q_bore, ctime=ctime)
        assert qm.point_is_init()
        qm.reset_point()
        assert not qm.point_is_init()

    def test_init_with_hwp(self, mod):
        qm = mod.QMap(mean_aber=True)
        q_bore, ctime = make_bore_and_ctime(mod, qm)
        q_hwp = qm.hwp_quat(np.zeros(N))
        qm.init_point(q_bore, ctime=ctime, q_hwp=q_hwp)
        assert qm.point_is_init()


# ---------------------------------------------------------------------------
# to_tod: constant map -> all-ones TOD
# ---------------------------------------------------------------------------


class TestToTod:
    def test_temperature_constant_map(self, mod):
        """Scanning a constant T=1 map should produce a TOD of all ones."""
        qm = mod.QMap(mean_aber=True)
        q_bore, ctime = make_bore_and_ctime(mod, qm)
        source = np.ones((1, NPIX))
        qm.init_source(source, pol=False)
        qm.init_point(q_bore, ctime=ctime)
        q_off = np.atleast_2d(qm.det_offset(0.0, 0.0, 0.0))
        tod = qm.to_tod(q_off)
        assert tod.shape == (1, N)
        assert np.allclose(tod, 1.0, atol=1e-5)

    def test_temperature_zero_map(self, mod):
        """Scanning a zero map should produce a TOD of all zeros."""
        qm = mod.QMap(mean_aber=True)
        q_bore, ctime = make_bore_and_ctime(mod, qm)
        source = np.zeros((1, NPIX))
        qm.init_source(source, pol=False)
        qm.init_point(q_bore, ctime=ctime)
        q_off = np.atleast_2d(qm.det_offset(0.0, 0.0, 0.0))
        tod = qm.to_tod(q_off)
        assert tod.shape == (1, N)
        assert np.allclose(tod, 0.0, atol=1e-10)

    def test_tod_shape_multi_det(self, mod):
        ndet = 4
        qm = mod.QMap(mean_aber=True)
        q_bore, ctime = make_bore_and_ctime(mod, qm)
        source = np.ones((1, NPIX))
        qm.init_source(source, pol=False)
        qm.init_point(q_bore, ctime=ctime)
        q_off = np.atleast_2d(
            qm.det_offset([0, 0.5, -0.5, 0], [0, 0, 0, 0.5], [0] * ndet)
        )
        tod = qm.to_tod(q_off)
        assert tod.shape == (ndet, N)

    def test_polarized_constant_T_Q_U(self, mod):
        """Scanning a pol map (T=1, Q=0, U=0) with a single det gives TOD of 1."""
        qm = mod.QMap(mean_aber=True)
        q_bore, ctime = make_bore_and_ctime(mod, qm)
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
    def test_projection_map_nonzero_after_from_tod(self, mod):
        """from_tod with any TOD should fill the hits map."""
        qm = mod.QMap(nside=NSIDE, pol=False, mean_aber=True)
        q_bore, ctime = make_bore_and_ctime(mod, qm)
        qm.init_point(q_bore, ctime=ctime)
        q_off = np.atleast_2d(qm.det_offset(0.0, 0.0, 0.0))
        tod = np.ones((1, N))
        vec, proj = qm.from_tod(q_off, tod=tod)
        assert np.any(proj > 0)

    def test_proj_equals_hits(self, mod):
        """With unit TOD weights, proj[pix] = number of hits on that pixel."""
        qm = mod.QMap(nside=NSIDE, pol=False, mean_aber=True)
        q_bore, ctime = make_bore_and_ctime(mod, qm)
        qm.init_point(q_bore, ctime=ctime)
        q_off = np.atleast_2d(qm.det_offset(0.0, 0.0, 0.0))
        tod = np.ones((1, N))
        vec, proj = qm.from_tod(q_off, tod=tod)
        # For T-only: proj = hits, vec = sum(tod) per pixel
        # vec / proj should equal 1 everywhere observed
        mask = proj.squeeze() > 0
        assert np.allclose(vec.squeeze()[mask] / proj.squeeze()[mask], 1.0, atol=1e-10)

    def test_from_tod_pol_shape(self, mod):
        qm = mod.QMap(nside=NSIDE, pol=True, mean_aber=True)
        q_bore, ctime = make_bore_and_ctime(mod, qm)
        qm.init_point(q_bore, ctime=ctime)
        q_off = np.atleast_2d(qm.det_offset(0.0, 0.0, 0.0))
        tod = np.ones((1, N))
        vec, proj = qm.from_tod(q_off, tod=tod)
        assert vec.shape == (3, NPIX)
        assert proj.shape == (6, NPIX)

    def test_count_hits_false(self):
        # qpoint only: depo is qpoint's own bookkeeping
        """count_hits=False should not accumulate the projection map."""
        qm = qpoint.QMap(nside=NSIDE, pol=False, mean_aber=True)
        q_bore, ctime = make_bore_and_ctime(qpoint, qm)
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
        # qpoint only: depo is qpoint's own bookkeeping
        """Scan a constant T map, bin TOD back to map, solve -> recover constant."""
        qm = qpoint.QMap(nside=NSIDE, pol=False, mean_aber=True)
        q_bore, ctime = make_bore_and_ctime(qpoint, qm)
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
        # qpoint only: depo is qpoint's own bookkeeping
        """Scan a constant T=42 map and verify solved map is also 42."""
        scale = 42.0
        qm = qpoint.QMap(nside=NSIDE, pol=False, mean_aber=True)
        q_bore, ctime = make_bore_and_ctime(qpoint, qm)
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
    def test_temperature_solve(self, mod):
        # vec and proj have equal values so solved = vec/proj = 1 everywhere observed
        vec = np.array([[1.0, 2.0, 0.0, 3.0]])
        proj = np.array([[1.0, 2.0, 0.0, 3.0]])
        qm = mod.QMap(mean_aber=True)
        solved = qm.solve_map(vec=vec.copy(), proj=proj.copy(), partial=True)
        assert np.allclose(solved[[0, 1, 3]], 1.0, atol=1e-10)
        assert solved[2] == 0.0  # zero proj -> fill=0

    def test_returns_correct_shape_temperature(self, mod):
        qm = mod.QMap(nside=NSIDE, pol=False, mean_aber=True)
        q_bore, ctime = make_bore_and_ctime(mod, qm)
        qm.init_point(q_bore, ctime=ctime)
        q_off = np.atleast_2d(qm.det_offset(0.0, 0.0, 0.0))
        tod = np.ones((1, N))
        qm.from_tod(q_off, tod=tod)
        solved = qm.solve_map()
        assert solved.shape == (NPIX,)

    def test_returns_correct_shape_polarized(self, mod):
        qm = mod.QMap(nside=NSIDE, pol=True, mean_aber=True)
        q_bore, ctime = make_bore_and_ctime(mod, qm)
        qm.init_point(q_bore, ctime=ctime)
        q_off = np.atleast_2d(qm.det_offset(0.0, 0.0, 0.0))
        tod = np.ones((1, N))
        qm.from_tod(q_off, tod=tod)
        solved = qm.solve_map()
        assert solved.shape == (3, NPIX)

    def test_return_mask(self, mod):
        qm = mod.QMap(nside=NSIDE, pol=False, mean_aber=True)
        q_bore, ctime = make_bore_and_ctime(mod, qm)
        qm.init_point(q_bore, ctime=ctime)
        q_off = np.atleast_2d(qm.det_offset(0.0, 0.0, 0.0))
        tod = np.ones((1, N))
        qm.from_tod(q_off, tod=tod)
        solved, mask = qm.solve_map(return_mask=True)
        assert mask.shape == (NPIX,)
        assert mask.dtype == bool

    def test_fill_value(self):
        # qpoint only: depo is qpoint's own bookkeeping
        qm = qpoint.QMap(nside=NSIDE, pol=False, mean_aber=True)
        q_bore, ctime = make_bore_and_ctime(qpoint, qm)
        qm.init_point(q_bore, ctime=ctime)
        q_off = np.atleast_2d(qm.det_offset(0.0, 0.0, 0.0))
        tod = np.ones((1, N))
        qm.from_tod(q_off, tod=tod)
        solved = qm.solve_map(fill=np.nan)
        mask = qm.depo["proj"].squeeze() > 0
        assert np.all(np.isnan(solved[~mask]))

    def test_cho_method(self, mod):
        qm = mod.QMap(nside=NSIDE, pol=False, mean_aber=True)
        q_bore, ctime = make_bore_and_ctime(mod, qm)
        qm.init_point(q_bore, ctime=ctime)
        q_off = np.atleast_2d(qm.det_offset(0.0, 0.0, 0.0))
        tod = np.ones((1, N))
        qm.from_tod(q_off, tod=tod)
        solved = qm.solve_map(method="cho")
        assert solved.shape == (NPIX,)

    def test_invalid_method_raises(self, mod):
        # Method validation only applies to polarized maps (T-only returns early)
        qm = mod.QMap(nside=NSIDE, pol=True, mean_aber=True)
        q_bore, ctime = make_bore_and_ctime(mod, qm)
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
    def test_unsolve_inverts_solve(self, mod):
        """unsolve_map(solve_map(vec, proj), proj) should recover the original vec."""
        vec_orig = np.array([[1.0, 2.0, 0.0, 3.0]])
        proj = np.array([[1.0, 2.0, 0.0, 3.0]])
        qm = mod.QMap(mean_aber=True)
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
    def test_temperature_uniform_proj(self, mod):
        """For T-only with uniform hits, condition number should be 1 everywhere."""
        nside = NSIDE
        npix = nside2npix(nside)
        proj = np.ones((1, npix))
        qm = mod.QMap(mean_aber=True)
        cond = qm.proj_cond(proj=proj, partial=True)
        assert cond.shape == (npix,)
        assert np.allclose(cond, 1.0, atol=1e-10)

    def test_empty_pixels_have_inf(self, mod):
        """Pixels with zero hits should have infinite condition number."""
        nside = NSIDE
        npix = nside2npix(nside)
        proj = np.ones((1, npix))
        proj[0, : npix // 2] = 0.0  # half the pixels have no hits
        qm = mod.QMap(mean_aber=True)
        cond = qm.proj_cond(proj=proj, partial=True)
        assert np.all(np.isinf(cond[: npix // 2]))
        assert np.all(np.isfinite(cond[npix // 2 :]))

    def test_from_depo(self, mod):
        """proj_cond() without args uses depo['proj']."""
        qm = mod.QMap(nside=NSIDE, pol=False, mean_aber=True)
        q_bore, ctime = make_bore_and_ctime(mod, qm)
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
        # qpoint only: depo is qpoint's own bookkeeping
        qm = qpoint.QMap(nside=NSIDE, pol=False, mean_aber=True)
        q_bore, ctime = make_bore_and_ctime(qpoint, qm)
        source = np.ones((1, NPIX))
        qm.init_source(source, pol=False)
        qm.init_point(q_bore, ctime=ctime)
        qm.reset()
        assert not qm.dest_is_init()
        assert not qm.source_is_init()
        assert not qm.point_is_init()
        assert len(qm.depo) == 0

    def test_reset_dest(self, mod):
        qm = mod.QMap(nside=NSIDE, pol=False, mean_aber=True)
        assert qm.dest_is_init()
        qm.reset_dest()
        assert not qm.dest_is_init()

    def test_reset_source(self, mod):
        source = np.ones((1, NPIX))
        qm = mod.QMap(source_map=source, source_pol=False, mean_aber=True)
        assert qm.source_is_init()
        qm.reset_source()
        assert not qm.source_is_init()

    def test_reset_point(self, mod):
        qm = mod.QMap(mean_aber=True)
        q_bore, ctime = make_bore_and_ctime(mod, qm)
        qm.init_point(q_bore, ctime=ctime)
        assert qm.point_is_init()
        qm.reset_point()
        assert not qm.point_is_init()


# ---------------------------------------------------------------------------
# dest_is_pol / source_is_pol
# ---------------------------------------------------------------------------


class TestPolFlags:
    def test_dest_is_pol_false_for_T(self, mod):
        qm = mod.QMap(nside=NSIDE, pol=False, mean_aber=True)
        assert not qm.dest_is_pol()

    def test_dest_is_pol_true_for_TQU(self, mod):
        qm = mod.QMap(nside=NSIDE, pol=True, mean_aber=True)
        assert qm.dest_is_pol()

    def test_dest_is_vpol(self, mod):
        qm = mod.QMap(nside=NSIDE, vpol=True, mean_aber=True)
        assert qm.dest_is_vpol()

    def test_source_is_pol(self, mod):
        source = np.ones((3, NPIX))
        qm = mod.QMap(source_map=source, source_pol=True, mean_aber=True)
        assert qm.source_is_pol()

    def test_source_is_not_pol(self, mod):
        source = np.ones((1, NPIX))
        qm = mod.QMap(source_map=source, source_pol=False, mean_aber=True)
        assert not qm.source_is_pol()


# ---------------------------------------------------------------------------
# Detector properties: flags, weights, gains
# ---------------------------------------------------------------------------


def make_mapper(mod, ndet=1, nside=NSIDE, pol=True, **kwargs):
    """A pointed QMap and a (ndet, 4) offset array, ready for from_tod."""
    qm = mod.QMap(nside=nside, pol=pol, mean_aber=True, **kwargs)
    q_bore, ctime = make_bore_and_ctime(mod, qm)
    qm.init_point(q_bore, ctime=ctime)
    delta = np.arange(ndet, dtype=float)
    q_off = np.atleast_2d(qm.det_offset(1.0 + delta, 2.0 + delta, 30.0 * delta))
    return qm, q_off


class TestFlags:
    """Flagged samples are dropped, which nothing else here exercises."""

    def test_flagged_samples_do_not_accumulate(self, mod):
        qm, q_off = make_mapper(mod)
        _, proj = qm.from_tod(q_off, tod=np.ones((1, N)))
        unflagged = proj[0].sum()

        qm, q_off = make_mapper(mod)
        flag = np.zeros((1, N), dtype=np.uint8)
        flag[0, :10] = 1
        _, proj = qm.from_tod(q_off, tod=np.ones((1, N)), flag=flag)
        assert proj[0].sum() == unflagged - 10

    def test_flagging_everything_leaves_the_map_empty(self, mod):
        qm, q_off = make_mapper(mod)
        flag = np.ones((1, N), dtype=np.uint8)
        vec, proj = qm.from_tod(q_off, tod=np.ones((1, N)), flag=flag)
        assert not np.any(proj)
        assert not np.any(vec)


class TestWeightsAndGain:
    def test_per_channel_weight_scales_the_map(self, mod):
        qm, q_off = make_mapper(mod)
        vec, proj = qm.from_tod(q_off, tod=np.ones((1, N)))
        v1, p1 = vec.copy(), proj.copy()

        qm, q_off = make_mapper(mod)
        vec, proj = qm.from_tod(q_off, tod=np.ones((1, N)), weight=2.5)
        assert np.allclose(vec, 2.5 * v1)
        assert np.allclose(proj, 2.5 * p1)

    def test_per_sample_weights_scale_the_map(self, mod):
        qm, q_off = make_mapper(mod)
        vec, proj = qm.from_tod(q_off, tod=np.ones((1, N)))
        v1, p1 = vec.copy(), proj.copy()

        qm, q_off = make_mapper(mod)
        weights = np.full((1, N), 3.0)
        vec, proj = qm.from_tod(q_off, tod=np.ones((1, N)), weights=weights)
        assert np.allclose(vec, 3.0 * v1)
        assert np.allclose(proj, 3.0 * p1)

    def test_gain_scales_the_signal_but_not_the_hits(self, mod):
        qm, q_off = make_mapper(mod)
        vec, proj = qm.from_tod(q_off, tod=np.ones((1, N)))
        v1, p1 = vec.copy(), proj.copy()

        qm, q_off = make_mapper(mod)
        vec, proj = qm.from_tod(q_off, tod=np.ones((1, N)), gain=4.0)
        assert np.allclose(vec, 4.0 * v1)
        assert np.allclose(proj, p1)


# ---------------------------------------------------------------------------
# Pair differencing
# ---------------------------------------------------------------------------


def make_pair(mod):
    """A pointed QMap and two offsets 90 degrees apart in polarization."""
    qm = mod.QMap(nside=NSIDE, pol=True, mean_aber=True)
    q_bore, ctime = make_bore_and_ctime(mod, qm)
    qm.init_point(q_bore, ctime=ctime)
    q_off = qm.det_offset([1.0, 1.0], [2.0, 2.0], [0.0, 90.0])
    return qm, q_off


class TestPairDifference:
    """
    do_diff pairs the first half of the detectors with the second half and
    accumulates their difference into the polarization rows and their sum
    into temperature. Nothing else in this suite reaches that kernel.
    """

    def test_common_mode_cancels_in_polarization(self, mod):
        qm, q_off = make_pair(mod)
        vec, _ = qm.from_tod(q_off, tod=np.ones((2, N)), do_diff=True)
        assert not np.any(vec[1])
        assert not np.any(vec[2])

    def test_common_mode_survives_in_temperature(self, mod):
        qm, q_off = make_pair(mod)
        vec, _ = qm.from_tod(q_off, tod=np.ones((2, N)), do_diff=True)
        assert np.isclose(vec[0].sum(), N)

    def test_a_differential_signal_appears_in_polarization(self, mod):
        qm, q_off = make_pair(mod)
        tod = np.vstack([np.ones(N), -np.ones(N)])
        vec, _ = qm.from_tod(q_off, tod=tod, do_diff=True)
        assert not np.any(vec[0])
        assert np.any(vec[1])
        assert np.any(vec[2])

    def test_one_hit_per_sample_not_per_detector(self, mod):
        qm, q_off = make_pair(mod)
        _, proj = qm.from_tod(q_off, tod=np.ones((2, N)), do_diff=True)
        assert proj[0].sum() == N

    def test_a_flag_on_either_detector_skips_the_sample(self, mod):
        """
        The pair is dropped whichever half carries the flag. The kernel
        used to test flag_init on either detector and then read both
        flag arrays, so this also pins that each is read behind its own.
        """
        for det in (0, 1):
            qm, q_off = make_pair(mod)
            flag = np.zeros((2, N), dtype=np.uint8)
            flag[det, :10] = 1
            _, proj = qm.from_tod(q_off, tod=np.ones((2, N)), do_diff=True, flag=flag)
            assert proj[0].sum() == N - 10


# ---------------------------------------------------------------------------
# Partial maps
# ---------------------------------------------------------------------------


def hit_pixels(mod, count=5):
    """The first `count` pixels a single detector actually lands on."""
    qm, q_off = make_mapper(mod)
    _, proj = qm.from_tod(q_off, tod=np.ones((1, N)))
    return np.nonzero(proj[0])[0][:count].astype(np.int64), proj[0].copy()


class TestPartialDest:
    """
    A destination map can cover a list of pixels rather than the sphere,
    which routes every lookup through the pixel hash.
    """

    def test_shapes_follow_the_pixel_list(self, mod):
        pixels, _ = hit_pixels(mod)
        qm = mod.QMap(mean_aber=True, error_missing=False)
        qm.init_dest(nside=NSIDE, pol=True, pixels=pixels)
        q_bore, ctime = make_bore_and_ctime(mod, qm)
        qm.init_point(q_bore, ctime=ctime)
        vec, proj = qm.from_tod(
            np.atleast_2d(qm.det_offset(1.0, 2.0, 0.0)), tod=np.ones((1, N))
        )
        assert vec.shape == (3, len(pixels))
        assert proj.shape == (6, len(pixels))

    def test_hits_match_the_full_sky_map_on_those_pixels(self, mod):
        pixels, full_hits = hit_pixels(mod)
        qm = mod.QMap(mean_aber=True, error_missing=False)
        qm.init_dest(nside=NSIDE, pol=True, pixels=pixels)
        q_bore, ctime = make_bore_and_ctime(mod, qm)
        qm.init_point(q_bore, ctime=ctime)
        _, proj = qm.from_tod(
            np.atleast_2d(qm.det_offset(1.0, 2.0, 0.0)), tod=np.ones((1, N))
        )
        assert np.allclose(proj[0], full_hits[pixels])

    def test_a_sample_off_the_map_raises_by_default(self, mod):
        pixels, _ = hit_pixels(mod)
        qm = mod.QMap(mean_aber=True)
        qm.init_dest(nside=NSIDE, pol=True, pixels=pixels)
        q_bore, ctime = make_bore_and_ctime(mod, qm)
        qm.init_point(q_bore, ctime=ctime)
        with pytest.raises(RuntimeError, match="out of bounds"):
            qm.from_tod(
                np.atleast_2d(qm.det_offset(1.0, 2.0, 0.0)), tod=np.ones((1, N))
            )

    def test_error_missing_false_drops_it_instead(self, mod):
        pixels, full_hits = hit_pixels(mod)
        qm = mod.QMap(mean_aber=True, error_missing=False)
        qm.init_dest(nside=NSIDE, pol=True, pixels=pixels)
        q_bore, ctime = make_bore_and_ctime(mod, qm)
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
        # qpoint only: reads the ctypes struct; test_parity pins the same inference for qpoint2
        qm = qpoint.QMap(mean_aber=True)
        q_bore, ctime = make_bore_and_ctime(qpoint, qm)
        qm.init_point(q_bore, ctime=ctime)
        qm.init_source(np.zeros((nrow, NPIX)), pol=pol, vpol=(nrow == 4))
        assert qm._source.contents.num_vec == nrow

    @pytest.mark.parametrize("nrow, pol, name", MODE_ROWS)
    def test_every_row_reaches_the_tod(self, mod, nrow, pol, name):
        """
        Perturb each row in turn and require the timestream to change. A
        row that never matters means the mode was inferred too small, and
        the rows past it are read from beyond the map.
        """
        rng = np.random.default_rng(0)
        base = rng.normal(size=(nrow, NPIX))

        def tod_for(source):
            qm = mod.QMap(mean_aber=True)
            q_bore, ctime = make_bore_and_ctime(mod, qm)
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
    def test_nan_missing_marks_samples_off_the_map(self, mod):
        pixels, _ = hit_pixels(mod)
        qm = mod.QMap(mean_aber=True, error_missing=False, nan_missing=True)
        q_bore, ctime = make_bore_and_ctime(mod, qm)
        qm.init_point(q_bore, ctime=ctime)
        qm.init_source(np.ones((3, len(pixels))), pol=True, nside=NSIDE, pixels=pixels)
        tod = np.asarray(qm.to_tod(np.atleast_2d(qm.det_offset(1.0, 2.0, 0.0))))
        assert np.isnan(tod).any()
        assert np.isfinite(tod).any()

    def test_without_nan_missing_they_are_left_at_zero(self, mod):
        pixels, _ = hit_pixels(mod)
        qm = mod.QMap(mean_aber=True, error_missing=False, nan_missing=False)
        q_bore, ctime = make_bore_and_ctime(mod, qm)
        qm.init_point(q_bore, ctime=ctime)
        qm.init_source(np.ones((3, len(pixels))), pol=True, nside=NSIDE, pixels=pixels)
        tod = np.asarray(qm.to_tod(np.atleast_2d(qm.det_offset(1.0, 2.0, 0.0))))
        assert not np.isnan(tod).any()
        assert np.count_nonzero(tod)


class TestNumThreads:
    """
    Threading is over detectors, and each thread accumulates into its own
    map before they are merged. The merge happens in whatever order the
    threads finish, so only the quantities that sum exactly come back
    identical.

    Where OpenMP is off -- Apple clang, by default -- there is one thread
    and everything is trivially identical. These assertions are written
    for the builds where it is on, which is Linux and any gcc build.
    """

    def run(self, mod, num_threads):
        rng = np.random.default_rng(1)
        qm, q_off = make_mapper(mod, ndet=4, num_threads=num_threads)
        vec, proj = qm.from_tod(q_off, tod=rng.normal(size=(4, N)))
        return np.asarray(vec).copy(), np.asarray(proj).copy()

    def test_the_hit_count_per_pixel_is_exact(self, mod):
        """Whole hits, so the sum is exact whatever order they arrive in."""
        (_, one), (_, four) = self.run(mod, 1), self.run(mod, 4)
        assert np.array_equal(one[0], four[0])

    def test_the_total_weight_is_exact(self, mod):
        (_, one), (_, four) = self.run(mod, 1), self.run(mod, 4)
        assert one[0].sum() == four[0].sum() == 4 * N

    def test_the_rest_agrees_to_rounding(self, mod):
        """
        The polarization cross terms and the signal are sums of products,
        which reassociate: identical to about 1e-15 here, checked well
        inside that so the test is not a rounding tripwire.
        """
        (vec_one, one), (vec_four, four) = self.run(mod, 1), self.run(mod, 4)
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

    def maps(self, mod, polconv):
        qm, q_off = make_mapper(mod, polconv=polconv)
        vec, _ = qm.from_tod(q_off, tod=np.ones((1, N)))
        return np.asarray(vec).copy()

    def test_u_flips_between_the_conventions(self, mod):
        assert np.allclose(self.maps(mod, "cosmo")[2], -self.maps(mod, "iau")[2])

    def test_t_and_q_do_not(self, mod):
        cosmo, iau = self.maps(mod, "cosmo"), self.maps(mod, "iau")
        assert np.allclose(cosmo[0], iau[0])
        assert np.allclose(cosmo[1], iau[1])


class TestSolveMapCho:
    """
    solve_map_cho is the Cholesky path and had no test of its own. It
    excludes the same pixels as the default solver -- cho_factor does not
    fail on a rank-deficient matrix, it returns nonsense, so the condition
    number is what keeps a one- or two-hit pixel out of both of them.
    """

    def solved(self, mod, pol):
        qm, q_off = make_mapper(mod, ndet=4, pol=pol)
        rng = np.random.default_rng(2)
        vec, proj = qm.from_tod(q_off, tod=rng.normal(size=(4, N)))
        direct, mask = qm.solve_map(vec=vec.copy(), proj=proj.copy(), return_mask=True)
        cho = qm.solve_map_cho(vec=vec.copy(), proj=proj.copy())
        return np.asarray(direct), np.asarray(cho), np.asarray(mask, dtype=bool)

    def test_agrees_with_the_default_solver_temperature(self, mod):
        direct, cho, mask = self.solved(mod, pol=False)
        assert mask.any()
        assert np.allclose(direct, cho)

    def test_agrees_with_the_default_solver_polarized(self, mod):
        pytest.importorskip("scipy")
        direct, cho, mask = self.solved(mod, pol=True)
        assert mask.any()
        assert np.allclose(direct, cho)

    def test_both_zero_the_pixels_they_cannot_determine(self, mod):
        pytest.importorskip("scipy")
        direct, cho, mask = self.solved(mod, pol=True)
        assert (~mask).any()
        assert not np.any(direct[:, ~mask])
        assert not np.any(cho[:, ~mask])

    @pytest.mark.parametrize("fault", ["singular", "not finite"])
    def test_a_pixel_the_factorization_refuses_is_masked(self, fault):
        """
        The two faults cho_factor reports, and it does not report them
        the same way: a matrix that is not positive definite raises
        LinAlgError, one holding an inf or a nan raises ValueError. Both
        mean the pixel cannot be solved and both are caught, narrowly, so
        that a third kind of error stays a bug rather than becoming a
        quietly dropped pixel.

        cond is supplied rather than computed, because the conditioning
        test would otherwise remove these pixels before the solver ever
        saw them -- an infinite condition number for the singular one and
        a nan for the other, both of which fail `cond < cond_thresh`.
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
    def test_solve_map(self, mod, kwargs):
        pytest.importorskip("scipy") if kwargs else None
        vec, proj = self.inputs()
        before = proj.copy()
        mod.QMap(nside=NSIDE, pol=True).solve_map(vec=vec, proj=proj, **kwargs)
        assert np.array_equal(proj, before)

    def test_proj_cond(self, mod):
        _, proj = self.inputs()
        before = proj.copy()
        mod.QMap(nside=NSIDE, pol=True).proj_cond(proj=proj)
        assert np.array_equal(proj, before)

    def test_unsolve_map(self, mod):
        map_in, proj = self.inputs()
        before = proj.copy()
        mod.QMap(nside=NSIDE, pol=True).unsolve_map(map_in=map_in, proj=proj)
        assert np.array_equal(proj, before)

    def test_returned_proj_is_still_written(self, mod):
        """
        The other half: asking for it back does give the solved form. The
        Cholesky path is the one that rewrites proj, replacing each hit
        pixel's matrix with its decomposition.
        """
        pytest.importorskip("scipy")
        vec, proj = self.inputs()
        before = proj.copy()
        _, out = mod.QMap(nside=NSIDE, pol=True).solve_map(
            vec=vec, proj=proj, return_proj=True, method="cho"
        )
        assert np.array_equal(proj, before), "the input must still be intact"
        assert not np.array_equal(np.asarray(out), before)


class TestCtimeIsRequiredWithoutMeanAber:
    """
    With mean_aber off, aberration is applied per detector, which needs
    the time of each sample -- so the kernels refuse to run without it.

    Every kernel checks this separately and every check was unexercised:
    the binned and differenced paths through tod2map, and map2tod. The
    rest of the guards in those functions test structures the Python
    layer always fills in, so this is the one a caller can actually
    reach, and the message says which kernel refused.
    """

    def mapper(self, mod, ndet=1, with_ctime=False):
        """
        A QMap with mean_aber off, pointed with or without ctime.
        make_mapper cannot serve here: it fixes mean_aber=True.
        """
        qm = mod.QMap(nside=NSIDE, pol=True, mean_aber=False)
        q_bore, ctime = make_bore_and_ctime(mod, qm)
        qm.init_point(q_bore, ctime=ctime if with_ctime else None)
        delta = np.arange(ndet, dtype=float)
        q_off = np.atleast_2d(qm.det_offset(1.0 + delta, 2.0 + delta, 30.0 * delta))
        return qm, q_off

    def test_from_tod(self, mod):
        qm, q_off = self.mapper(mod)
        with pytest.raises(RuntimeError, match="ctime required"):
            qm.from_tod(q_off, tod=np.ones((1, N)))

    def test_from_tod_differenced(self, mod):
        """The differencing kernel has its own copy of the check."""
        qm, q_off = self.mapper(mod, ndet=2)
        with pytest.raises(RuntimeError, match="ctime required"):
            qm.from_tod(q_off, tod=np.ones((2, N)), do_diff=True)

    def test_to_tod(self, mod):
        qm, q_off = self.mapper(mod)
        qm.init_source(np.zeros((3, NPIX)), pol=True)
        with pytest.raises(RuntimeError, match="ctime required"):
            qm.to_tod(q_off)

    def test_ctime_makes_all_three_work(self, mod):
        """
        The control: the same calls succeed once ctime is supplied, so
        the tests above are pinning the missing time and not some other
        unfinished setup.
        """
        qm, q_off = self.mapper(mod, ndet=2, with_ctime=True)
        qm.from_tod(q_off, tod=np.ones((2, N)))
        qm, q_off = self.mapper(mod, ndet=2, with_ctime=True)
        qm.from_tod(q_off, tod=np.ones((2, N)), do_diff=True)
        qm, q_off = self.mapper(mod, with_ctime=True)
        qm.init_source(np.zeros((3, NPIX)), pol=True)
        assert np.asarray(qm.to_tod(q_off)).shape == (1, N)


class TestDetarrLifetime:
    """
    from_tod builds the detector array and tears it down again, so the
    structure is only alive between init_detarr and reset_detarr.
    """

    def test_from_tod_leaves_no_detector_array_behind(self, mod):
        qm, q_off = make_mapper(mod)
        qm.from_tod(q_off, tod=np.ones((1, N)))
        assert qm._detarr is None

    def test_init_detarr_then_reset(self):
        # qpoint only: depo is qpoint's own bookkeeping
        qm, q_off = make_mapper(qpoint)
        qm.init_detarr(q_off, tod=np.ones((1, N)))
        assert qm._detarr is not None
        assert "tod" in qm.depo
        qm.reset_detarr()
        assert qm._detarr is None
        assert "tod" not in qm.depo

    def test_reset_is_safe_when_there_is_nothing_to_reset(self, mod):
        qm, _ = make_mapper(mod)
        qm.reset_detarr()
        qm.reset_detarr()

    def test_mapping_accumulates_across_calls(self, mod):
        qm, q_off = make_mapper(mod)
        _, proj = qm.from_tod(q_off, tod=np.ones((1, N)))
        first = proj[0].sum()
        _, proj = qm.from_tod(q_off, tod=np.ones((1, N)))
        assert proj[0].sum() == 2 * first


# ---------------------------------------------------------------------------
# qpoint2 only
# ---------------------------------------------------------------------------


class TestFastPixInMapmaking:
    """
    fast_pix reaches the mapmaking kernels too, and has to leave them
    bit-identical to the two-step path it replaces.

    qpoint2 only, and named rather than taken from the mod fixture for
    the same reason as TestFastPix in test_qpoint.py: the C's fast path
    is the less accurate of its two, so it is not held to this.
    """

    def source_map(self):
        return np.random.default_rng(0).normal(size=(3, NPIX))

    def test_to_tod(self):
        """Through map2tod's kernel."""
        out = []
        for fast in (False, True):
            qm, q_off = make_mapper(qpoint2, fast_pix=fast)
            qm.init_source(self.source_map(), pol=True)
            out.append(np.asarray(qm.to_tod(q_off)))
        assert np.array_equal(out[0], out[1]), "to_tod fast vs two-step"

    def test_from_tod(self):
        """And through tod2map's."""
        tod = np.random.default_rng(1).normal(size=(1, N))
        out = []
        for fast in (False, True):
            qm, q_off = make_mapper(qpoint2, fast_pix=fast)
            out.append([np.asarray(x) for x in qm.from_tod(q_off, tod=tod.copy())])
        for a, b, name in zip(out[0], out[1], ("vec", "proj")):
            assert np.array_equal(a, b), "from_tod {} fast vs two-step".format(name)
