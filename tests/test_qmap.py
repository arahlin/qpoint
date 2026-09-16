"""Tests for qpoint.QMap mapmaking class and helper functions."""

import numpy as np
import pytest
import qpoint
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
    def test_temperature_shapes(self):
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
        qm = qpoint.QMap(mean_aber=True)
        vec, proj = qm.init_dest(nside=NSIDE, pol=True)
        assert vec.shape == (3, NPIX)
        assert proj.shape == (6, NPIX)

    def test_vpol_shapes(self):
        qm = qpoint.QMap(mean_aber=True)
        vec, proj = qm.init_dest(nside=NSIDE, vpol=True)
        assert vec.shape == (4, NPIX)
        assert proj.shape == (10, NPIX)

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

    def test_reset_point(self):
        qm = qpoint.QMap(mean_aber=True)
        q_bore, ctime = make_bore_and_ctime(qm)
        qm.init_point(q_bore, ctime=ctime)
        assert qm.point_is_init()
        qm.reset_point()
        assert not qm.point_is_init()

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
