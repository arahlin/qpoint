"""
Cross-implementation parity tests.

qpoint (ctypes over the C library) and qpoint2 (pybind11 over the C++
rewrite) must agree bit for bit: same ERFA, same algorithms, same order of
operations. Everything here asserts exact
equality rather than a tolerance -- a tolerance would hide exactly the kind
of drift this suite exists to catch.
"""

import numpy as np
import pytest

import qpoint

qpoint2 = pytest.importorskip("qpoint2")

# Each implementation is compared against qpoint, the reference.
IMPLS = [pytest.param(qpoint2, id="qpoint2")]

# Option sets that select materially different code paths.
OPTIONS = [
    pytest.param({}, id="defaults"),
    pytest.param({"accuracy": "low"}, id="accuracy-low"),
    pytest.param({"fast_math": True}, id="fast-math"),
    pytest.param({"fast_aber": False}, id="slow-aber"),
    pytest.param({"mean_aber": False}, id="per-det-aber"),
    pytest.param({"polconv": "iau"}, id="polconv-iau"),
]

N = 50
CTIME = 1418662800.0 + np.arange(N, dtype=float)
AZ = np.linspace(0.0, 360.0, N)
EL = np.linspace(30.0, 70.0, N)
PSI = np.linspace(-30.0, 30.0, N)
PITCH = np.linspace(-1.0, 1.0, N)
ROLL = np.linspace(-2.0, 2.0, N)
LON = np.full(N, 165.7)
LAT = np.full(N, -77.6)
RA = np.linspace(0.0, 350.0, N)
DEC = np.linspace(-80.0, 80.0, N)
PA = np.linspace(-170.0, 170.0, N)
HWP = np.linspace(0.0, 180.0, N)


def identical(a, b):
    """Exact equality, treating NaNs in the same slots as equal."""
    a, b = np.asarray(a), np.asarray(b)
    if a.shape != b.shape:
        return False
    return np.array_equal(a, b, equal_nan=np.issubdtype(a.dtype, np.floating))


def assert_identical(expected, actual, label):
    if isinstance(expected, tuple):
        assert len(expected) == len(actual), "{}: arity differs".format(label)
        for i, (e, a) in enumerate(zip(expected, actual)):
            assert identical(e, a), "{}[{}] differs".format(label, i)
    else:
        assert identical(expected, actual), "{} differs".format(label)


def qp(mod, **options):
    """A freshly initialized QPoint, so no rate-cache state carries over."""
    return mod.QPoint(**options)


def bore(mod, **options):
    q = qp(mod, **options)
    return q, np.asarray(q.azel2bore(AZ, EL, PITCH, ROLL, LON, LAT, CTIME))


@pytest.mark.parametrize("mod", IMPLS)
@pytest.mark.parametrize("options", OPTIONS)
class TestPointing:
    def test_azel2bore(self, mod, options):
        _, ref = bore(qpoint, **options)
        _, got = bore(mod, **options)
        assert_identical(ref, got, "azel2bore")

    def test_azelpsi2bore(self, mod, options):
        args = (AZ, EL, PSI, PITCH, ROLL, LON, LAT, CTIME)
        ref = np.asarray(qp(qpoint, **options).azelpsi2bore(*args))
        got = np.asarray(qp(mod, **options).azelpsi2bore(*args))
        assert_identical(ref, got, "azelpsi2bore")

    @pytest.mark.parametrize(
        "kwargs",
        [
            pytest.param({}, id="sincos"),
            pytest.param({"return_pa": True}, id="pa"),
            pytest.param({"sindec": True}, id="sindec"),
        ],
    )
    def test_bore2radec(self, mod, options, kwargs):
        out = []
        for m in (qpoint, mod):
            q, qb = bore(m, **options)
            q_off = q.det_offset(1.0, 2.0, 3.0)
            out.append(tuple(q.bore2radec(q_off, CTIME, qb, **kwargs)))
        assert_identical(out[0], out[1], "bore2radec")

    @pytest.mark.parametrize(
        "kwargs",
        [
            pytest.param({}, id="sincos"),
            pytest.param({"return_pa": True}, id="pa"),
            pytest.param({"sindec": True}, id="sindec"),
        ],
    )
    def test_bore2radec_hwp(self, mod, options, kwargs):
        """
        A HWP crossed with each output mode. The C has a separate entry
        point per combination, so the modes have to be crossed rather than
        tried one at a time: qp_bore2rasindec_hwp is reached only by
        sindec together with a HWP.
        """
        out = []
        for m in (qpoint, mod):
            q, qb = bore(m, **options)
            q_off = q.det_offset(1.0, 2.0, 3.0)
            q_hwp = np.asarray(q.hwp_quat(HWP))
            out.append(tuple(q.bore2radec(q_off, CTIME, qb, q_hwp=q_hwp, **kwargs)))
        assert_identical(out[0], out[1], "bore2radec_hwp")

    def test_bore2azel(self, mod, options):
        out = []
        for m in (qpoint, mod):
            q, qb = bore(m, **options)
            out.append(tuple(q.bore2azel(qb, LON, LAT, CTIME)))
        assert_identical(out[0], out[1], "bore2azel")

    def test_radec2azel(self, mod, options):
        args = (RA, DEC, PA, LON, LAT, CTIME)
        ref = tuple(qp(qpoint, **options).radec2azel(*args))
        got = tuple(qp(mod, **options).radec2azel(*args))
        assert_identical(ref, got, "radec2azel")

    def test_radecpa2quat(self, mod, options):
        ref = np.asarray(qp(qpoint, **options).radecpa2quat(RA, DEC, PA))
        got = np.asarray(qp(mod, **options).radecpa2quat(RA, DEC, PA))
        assert_identical(ref, got, "radecpa2quat")

    def test_quat2radecpa(self, mod, options):
        out = []
        for m in (qpoint, mod):
            q = qp(m, **options)
            out.append(tuple(q.quat2radecpa(q.radecpa2quat(RA, DEC, PA))))
        assert_identical(out[0], out[1], "quat2radecpa")

    @pytest.mark.parametrize(
        "kwargs",
        [
            pytest.param({}, id="sincos"),
            pytest.param({"hwp": HWP}, id="hwp"),
            pytest.param({"return_pa": True}, id="pa"),
            pytest.param({"sindec": True}, id="sindec"),
            pytest.param({"hwp": HWP, "sindec": True}, id="hwp-sindec"),
        ],
    )
    def test_azel2radec(self, mod, options, kwargs):
        args = (1.0, 2.0, 3.0, AZ, EL, PITCH, ROLL, LON, LAT, CTIME)
        ref = tuple(qp(qpoint, **options).azel2radec(*args, **kwargs))
        got = tuple(qp(mod, **options).azel2radec(*args, **kwargs))
        assert_identical(ref, got, "azel2radec")

    @pytest.mark.parametrize(
        "kwargs",
        [
            pytest.param({}, id="sincos"),
            pytest.param({"hwp": HWP}, id="hwp"),
            pytest.param({"return_pa": True}, id="pa"),
            pytest.param({"sindec": True}, id="sindec"),
            pytest.param({"hwp": HWP, "sindec": True}, id="hwp-sindec"),
        ],
    )
    def test_azelpsi2radec(self, mod, options, kwargs):
        args = (1.0, 2.0, 3.0, AZ, EL, PSI, PITCH, ROLL, LON, LAT, CTIME)
        ref = tuple(qp(qpoint, **options).azelpsi2radec(*args, **kwargs))
        got = tuple(qp(mod, **options).azelpsi2radec(*args, **kwargs))
        assert_identical(ref, got, "azelpsi2radec")

    def test_azel2radec_restores_mean_aber(self, mod, options):
        """The call forces mean_aber on; it must put it back."""
        q = qp(mod, **options)
        before = q.get("mean_aber")
        q.azel2radec(1.0, 2.0, 3.0, AZ, EL, PITCH, ROLL, LON, LAT, CTIME)
        assert q.get("mean_aber") == before

    def test_gmst(self, mod, options):
        ref = qp(qpoint, **options).gmst(CTIME)
        got = qp(mod, **options).gmst(CTIME)
        assert_identical(ref, got, "gmst")

    def test_lmst(self, mod, options):
        ref = qp(qpoint, **options).lmst(CTIME, LON)
        got = qp(mod, **options).lmst(CTIME, LON)
        assert_identical(ref, got, "lmst")


OMEGA = np.random.default_rng(0).normal(size=(3, N)) * 0.01


@pytest.mark.parametrize("mod", IMPLS)
@pytest.mark.parametrize("options", OPTIONS)
class TestDipole:
    def test_dipole(self, mod, options):
        ref = qp(qpoint, **options).dipole(CTIME, RA, DEC)
        got = qp(mod, **options).dipole(CTIME, RA, DEC)
        assert_identical(ref, got, "dipole")

    def test_bore2dipole(self, mod, options):
        out = []
        for m in (qpoint, mod):
            q, qb = bore(m, **options)
            out.append(q.bore2dipole(q.det_offset(1.0, 2.0, 3.0), CTIME, qb))
        assert_identical(out[0], out[1], "bore2dipole")

    def test_scalar_matches_array(self, mod, options):
        """The dipole is smooth, so a scalar direction must match element 0."""
        q = qp(mod, **options)
        assert identical(
            q.dipole(CTIME[0], RA[0], DEC[0]),
            np.asarray(q.dipole(CTIME, RA, DEC))[0],
        )


@pytest.mark.parametrize("mod", IMPLS)
class TestOmega2AzElPsi:
    @pytest.mark.parametrize("options", OPTIONS)
    def test_matches_reference(self, mod, options):
        args = (10.0, 45.0, 0.0, OMEGA[0], OMEGA[1], OMEGA[2], 0.01)
        ref = tuple(qp(qpoint, **options).omega2azelpsi(*args))
        got = tuple(qp(mod, **options).omega2azelpsi(*args))
        assert_identical(ref, got, "omega2azelpsi")

    def test_zero_rates_hold_position(self, mod):
        zero = np.zeros(N)
        az, el, psi = qp(mod).omega2azelpsi(10.0, 45.0, 0.0, zero, zero, zero, 0.01)
        assert np.allclose(az, 10.0)
        assert np.allclose(el, 45.0)
        assert np.allclose(psi, 0.0)


@pytest.mark.parametrize("mod", IMPLS)
class TestBoreOffset:
    @pytest.mark.parametrize("post", [False, True])
    def test_matches_reference(self, mod, post):
        out = []
        for m in (qpoint, mod):
            _, qb = bore(m)
            out.append(qp(m).bore_offset(qb.copy(), 1.0, 2.0, 3.0, post=post))
        assert_identical(out[0], out[1], "bore_offset")

    @pytest.mark.parametrize("post", [False, True])
    def test_per_sample_angles(self, mod, post):
        out = []
        for m in (qpoint, mod):
            _, qb = bore(m)
            out.append(qp(m).bore_offset(qb.copy(), RA / 100, DEC / 100, PA, post=post))
        assert_identical(out[0], out[1], "bore_offset per-sample")

    def test_no_angle_raises(self, mod):
        _, qb = bore(mod)
        with pytest.raises(ValueError, match="ang1"):
            qp(mod).bore_offset(qb)

    def test_inplace_defaults_off(self, mod):
        """bore_offset is the one rotation that does not mutate by default."""
        _, qb = bore(mod)
        before = qb.copy()
        qp(mod).bore_offset(qb, 1.0, 2.0, 3.0)
        assert identical(qb, before)

    def test_inplace_true_mutates(self, mod):
        _, qb = bore(mod)
        before = qb.copy()
        out = qp(mod).bore_offset(qb, 1.0, 2.0, 3.0, inplace=True)
        assert not np.array_equal(qb, before)
        assert np.shares_memory(out, qb)


# Option sets that matter for pixelization specifically.
PIX_OPTIONS = [
    pytest.param({}, id="defaults"),
    pytest.param({"fast_pix": True}, id="fast-pix"),
    pytest.param({"pix_order": "nest"}, id="nest"),
    pytest.param({"fast_pix": True, "pix_order": "nest"}, id="fast-pix-nest"),
    pytest.param({"fast_math": True}, id="fast-math"),
    pytest.param({"polconv": "iau"}, id="polconv-iau"),
]

NSIDE = 128


@pytest.mark.parametrize("mod", IMPLS)
@pytest.mark.parametrize("options", PIX_OPTIONS)
class TestPixelization:
    def test_radec2pix(self, mod, options):
        ref = qp(qpoint, **options).radec2pix(RA, DEC, nside=NSIDE)
        got = qp(mod, **options).radec2pix(RA, DEC, nside=NSIDE)
        assert_identical(ref, got, "radec2pix")

    @pytest.mark.parametrize(
        "kwargs",
        [pytest.param({}, id="pol"), pytest.param({"pol": False}, id="no-pol")],
    )
    def test_quat2pix(self, mod, options, kwargs):
        out = []
        for m in (qpoint, mod):
            q = qp(m, **options)
            res = q.quat2pix(q.radecpa2quat(RA, DEC, PA), nside=NSIDE, **kwargs)
            out.append(res if isinstance(res, tuple) else (res,))
        assert_identical(out[0], out[1], "quat2pix")

    def test_quat2pixpa(self, mod, options):
        out = []
        for m in (qpoint, mod):
            q = qp(m, **options)
            out.append(tuple(q.quat2pixpa(q.radecpa2quat(RA, DEC, PA), nside=NSIDE)))
        assert_identical(out[0], out[1], "quat2pixpa")

    @pytest.mark.parametrize(
        "kwargs",
        [
            pytest.param({}, id="pol"),
            pytest.param({"pol": False}, id="no-pol"),
            pytest.param({"return_pa": True}, id="pa"),
        ],
    )
    def test_bore2pix(self, mod, options, kwargs):
        out = []
        for m in (qpoint, mod):
            q, qb = bore(m, **options)
            q_off = q.det_offset(1.0, 2.0, 3.0)
            res = q.bore2pix(q_off, CTIME, qb, nside=NSIDE, **kwargs)
            out.append(res if isinstance(res, tuple) else (res,))
        assert_identical(out[0], out[1], "bore2pix")

    @pytest.mark.parametrize(
        "kwargs",
        [
            pytest.param({}, id="pol"),
            pytest.param({"return_pa": True}, id="pa"),
        ],
    )
    def test_bore2pix_hwp(self, mod, options, kwargs):
        """
        Pixelization with a HWP, which nothing else exercises. Crossed
        with return_pa because that combination is its own C entry point,
        qp_bore2pixpa_hwp, and is not reached by either on its own.
        """
        out = []
        for m in (qpoint, mod):
            q, qb = bore(m, **options)
            q_off = q.det_offset(1.0, 2.0, 3.0)
            q_hwp = np.asarray(q.hwp_quat(HWP))
            res = q.bore2pix(q_off, CTIME, qb, q_hwp=q_hwp, nside=NSIDE, **kwargs)
            out.append(res if isinstance(res, tuple) else (res,))
        assert_identical(out[0], out[1], "bore2pix_hwp")


MAP_RNG = np.random.default_rng(0)
INTERP_DEC = np.linspace(-85.0, 85.0, N)


@pytest.mark.parametrize("mod", IMPLS)
@pytest.mark.parametrize("nest", [False, True])
class TestInterpVal:
    @pytest.mark.parametrize("nside", [16, 64])
    def test_single_map(self, mod, nest, nside):
        m = MAP_RNG.normal(size=12 * nside * nside)
        ref = qp(qpoint).get_interp_val(m, RA, INTERP_DEC, nest=nest)
        got = qp(mod).get_interp_val(m, RA, INTERP_DEC, nest=nest)
        assert_identical(ref, got, "get_interp_val")

    def test_multi_map(self, mod, nest):
        m = MAP_RNG.normal(size=(3, 12 * 16 * 16))
        ref = qp(qpoint).get_interp_val(m, RA, INTERP_DEC, nest=nest)
        got = qp(mod).get_interp_val(m, RA, INTERP_DEC, nest=nest)
        assert np.asarray(got).shape == (3, N)
        assert_identical(ref, got, "get_interp_val multi")

    def test_constant_map_interpolates_to_constant(self, mod, nest):
        m = np.full(12 * 16 * 16, 2.5)
        got = qp(mod).get_interp_val(m, RA, INTERP_DEC, nest=nest)
        assert np.allclose(got, 2.5)


@pytest.mark.parametrize("mod", IMPLS)
@pytest.mark.parametrize("options", PIX_OPTIONS)
class TestGalacticRotation:
    def test_radec2gal_pa(self, mod, options):
        out = []
        for m in (qpoint, mod):
            out.append(
                tuple(qp(m, **options).radec2gal(RA.copy(), DEC.copy(), PA.copy()))
            )
        assert_identical(out[0], out[1], "radec2gal")

    def test_gal2radec_pa(self, mod, options):
        out = []
        for m in (qpoint, mod):
            out.append(
                tuple(qp(m, **options).gal2radec(RA.copy(), DEC.copy(), PA.copy()))
            )
        assert_identical(out[0], out[1], "gal2radec")

    def test_radec2gal_sincos(self, mod, options):
        out = []
        for m in (qpoint, mod):
            out.append(
                tuple(
                    qp(m, **options).radec2gal(
                        RA.copy(),
                        DEC.copy(),
                        sin2psi=np.sin(PA).copy(),
                        cos2psi=np.cos(PA).copy(),
                    )
                )
            )
        assert_identical(out[0], out[1], "radec2gal sin/cos")

    def test_rotate_quat(self, mod, options):
        out = []
        for m in (qpoint, mod):
            q = qp(m, **options)
            out.append(np.asarray(q.rotate_quat(q.radecpa2quat(RA, DEC, PA))))
        assert_identical(out[0], out[1], "rotate_quat")

    def test_rotate_coord_gc(self, mod, options):
        out = []
        for m in (qpoint, mod):
            out.append(
                tuple(
                    qp(m, **options).rotate_coord(
                        RA.copy(), DEC.copy(), PA.copy(), coord=("G", "C")
                    )
                )
            )
        assert_identical(out[0], out[1], "rotate_coord G->C")


@pytest.mark.parametrize("mod", IMPLS)
class TestRotateMap:

    @staticmethod
    def smooth_map(nside=32):
        """Band-limited, so pixel resampling error stays small."""
        hp = pytest.importorskip("healpy")
        rng = np.random.default_rng(0)
        npix = 12 * nside * nside
        alms = tuple(hp.map2alm(rng.normal(size=npix), lmax=24) for _ in range(3))
        return np.array(hp.alm2map(alms, nside))

    @pytest.mark.parametrize("coord", [("C", "G"), ("G", "C")])
    @pytest.mark.parametrize("interp", [True, False])
    def test_matches_reference(self, mod, coord, interp):
        m = self.smooth_map()
        ref = np.asarray(qp(qpoint).rotate_map(m, coord=coord, interp_pix=interp))
        got = np.asarray(qp(mod).rotate_map(m, coord=coord, interp_pix=interp))
        assert_identical(ref, got, "rotate_map")

    def test_nest_ordering(self, mod):
        hp = pytest.importorskip("healpy")
        m = np.array([hp.reorder(x, r2n=True) for x in self.smooth_map()])
        ref = np.asarray(qp(qpoint, pix_order="nest").rotate_map(m))
        got = np.asarray(qp(mod, pix_order="nest").rotate_map(m))
        assert_identical(ref, got, "rotate_map nest")

    @pytest.mark.parametrize("coord", [("C", "G"), ("G", "C")])
    def test_matches_healpy(self, mod, coord):
        """Interpolated rotation should agree with healpy's own."""
        hp = pytest.importorskip("healpy")
        m = self.smooth_map()
        got = np.asarray(qp(mod).rotate_map(m, coord=coord, interp_pix=True))
        ref = np.asarray(hp.Rotator(coord=list(coord)).rotate_map_pixel(m))
        assert np.allclose(got, ref, rtol=0, atol=1e-5)

    def test_round_trip_is_lossy_like_healpy(self, mod):
        """
        Pixel-space rotation does not invert. That is inherent, not a bug:
        healpy loses the same amount, so pin them against each other rather
        than against the input.
        """
        hp = pytest.importorskip("healpy")
        m = self.smooth_map()
        q = qp(mod)
        back = np.asarray(
            q.rotate_map(
                np.asarray(q.rotate_map(m, coord=("C", "G"))), coord=("G", "C")
            )
        )
        hg = hp.Rotator(coord=["C", "G"]).rotate_map_pixel(m)
        hback = np.asarray(hp.Rotator(coord=["G", "C"]).rotate_map_pixel(hg))
        assert not np.allclose(back, m, atol=1e-3)  # genuinely lossy
        assert np.allclose(back, hback, rtol=0, atol=1e-5)  # but lossy the same way

    @pytest.mark.parametrize("nrow", [1, 2, 4])
    def test_rejects_wrong_row_count(self, mod, nrow):
        """Fewer than three rows used to read off the end and segfault."""
        npix = 12 * 16 * 16
        with pytest.raises(ValueError, match="3 rows"):
            qp(mod).rotate_map(np.ones((nrow, npix)))

    @pytest.mark.parametrize("coord", [("C", "E"), ("X", "G"), (1, 2)])
    def test_unrecognized_coord_raises(self, mod, coord):
        """Used to return an all-zero map, silently destroying the input."""
        with pytest.raises(ValueError):
            qp(mod).rotate_map(self.smooth_map(), coord=coord)

    @pytest.mark.parametrize("coord", [("c", "g"), ("C", "g"), ("g", "C")])
    def test_coord_is_case_insensitive(self, mod, coord):
        m = self.smooth_map()
        upper = tuple(c.upper() for c in coord)
        ref = np.asarray(qp(mod).rotate_map(m, coord=upper))
        got = np.asarray(qp(mod).rotate_map(m, coord=coord))
        assert_identical(ref, got, "rotate_map case")

    def test_no_warning(self, mod):
        """The blanket "this code is buggy" warning is gone."""
        import warnings

        with warnings.catch_warnings():
            warnings.simplefilter("error", UserWarning)
            qp(mod).rotate_map(self.smooth_map(), coord=("C", "G"))

    @pytest.mark.parametrize("coord", [("C", "C"), ("G", "G"), ("c", "c")])
    def test_same_coord_returns_the_input(self, mod, coord):
        """Also used to return zeros. Passes the map straight through now."""
        m = self.smooth_map()
        out = qp(mod).rotate_map(m, coord=coord)
        assert identical(out, m)
        assert np.shares_memory(out, m)


@pytest.mark.parametrize("mod", IMPLS)
class TestRotationIsInPlace:
    """
    The rotations write through to the caller's arrays. For qpoint2 this
    also pins the zero-copy contract: the binding writes into the very
    buffer it was handed.
    """

    def test_inplace_mutates_input(self, mod):
        ra, dec, pa = RA.copy(), DEC.copy(), PA.copy()
        before = ra.copy()
        out = qp(mod).radec2gal(ra, dec, pa, inplace=True)
        assert not np.array_equal(ra, before)
        assert np.shares_memory(out[0], ra)

    def test_not_inplace_leaves_input_alone(self, mod):
        ra, dec, pa = RA.copy(), DEC.copy(), PA.copy()
        before = ra.copy()
        out = qp(mod).radec2gal(ra, dec, pa, inplace=False)
        assert identical(ra, before)
        assert not np.shares_memory(out[0], ra)

    def test_unsupported_coord_raises(self, mod):
        with pytest.raises(ValueError, match="[Uu]nsupported coord"):
            qp(mod).rotate_coord(RA.copy(), DEC.copy(), PA.copy(), coord=("C", "E"))


class TestQpSettingsIsPublic:
    """
    The per-call parameter decorator is exported, so a subclass adding a
    method gets the same keyword handling as the built-in ones.
    """

    def subclass(self):
        class MyPoint(qpoint2.QPoint):
            @qpoint2.qp_settings
            def scan(self, ctime, offset=0.0):
                return self.gmst(ctime) + offset

        return MyPoint()

    def test_importable_from_the_package(self):
        assert qpoint2.qp_settings is not None
        assert "qp_settings" in qpoint2.__all__

    def test_parameters_are_applied_and_restored(self):
        q = self.subclass()
        before = q.get("accuracy")
        assert q.scan(CTIME[0], accuracy="low") != q.scan(CTIME[0])
        assert q.get("accuracy") == before

    def test_the_methods_own_arguments_still_reach_it(self):
        q = self.subclass()
        assert q.scan(CTIME[0], offset=1.0) == q.scan(CTIME[0]) + 1.0
        assert q.scan(CTIME[0], offset=1.0, accuracy="low") == (
            q.scan(CTIME[0], accuracy="low") + 1.0
        )

    def test_an_unknown_keyword_is_reported(self):
        with pytest.raises(TypeError):
            self.subclass().scan(CTIME[0], nonsense=1)

    def test_defaults_apply_under_the_caller(self):
        class Defaulted(qpoint2.QPoint):
            @qpoint2.qp_settings(accuracy="low")
            def which(self):
                return self.get_param("accuracy")

        q = Defaulted()
        assert q.which() == "low"
        assert q.which(accuracy="high") == "high"
        assert q.get("accuracy") == qp(qpoint2).get("accuracy")


class TestLengthOneArraysKeepTheirAxis:
    """
    A deliberate divergence, and the one place the packages disagree on
    shape rather than on numbers.

    qpoint2 decides from how the arguments arrived: a true scalar (or a
    0-d array) degrades, an array keeps its axis even at length one. That
    cannot be read off the output, which is length one either way, so
    qpoint -- which decides from the output size -- cannot tell `f(1.0)`
    from `f(np.array([1.0]))` and collapses both.

    qpoint is not consistent about it either, which is the other reason to
    diverge: it keeps the axis for azel2radec, bore2pix and radec2azel and
    drops it for the six below. The values are identical throughout; only
    the shape differs.
    """

    NAMES = ["gmst", "lmst", "dipole", "bore2radec", "det_offset", "radec2gal"]
    ONE = np.array([CTIME[0]])

    def calls(self, mod):
        q = qp(mod)
        qb = np.atleast_2d(
            q.azel2bore(45.0, 45.0, None, None, LON[0], LAT[0], CTIME[0])
        )
        q_off = q.det_offset(0.0, 0.0, 0.0)
        return {
            "gmst": lambda: q.gmst(self.ONE),
            "lmst": lambda: q.lmst(self.ONE, LON[:1]),
            "dipole": lambda: q.dipole(self.ONE, RA[:1], DEC[:1]),
            "bore2radec": lambda: q.bore2radec(q_off, self.ONE, qb),
            "det_offset": lambda: q.det_offset(
                np.array([1.0]), np.array([2.0]), np.array([3.0])
            ),
            # copies: radec2gal rotates in place, and these are slices of
            # the module-level arrays every other test reads
            "radec2gal": lambda: q.radec2gal(
                RA[:1].copy(), DEC[:1].copy(), PA[:1].copy()
            ),
        }

    @pytest.mark.parametrize("name", NAMES)
    def test_qpoint2_keeps_the_axis(self, name):
        got = self.calls(qpoint2)[name]()
        for part in got if isinstance(got, tuple) else (got,):
            arr = np.asarray(part)
            assert arr.ndim >= 1 and arr.shape[0] == 1, (name, arr.shape)

    @pytest.mark.parametrize("name", NAMES)
    def test_the_divergence_is_real(self, name):
        """Guards the list above: qpoint really does collapse each of these."""
        ref = self.calls(qpoint)[name]()
        first = ref[0] if isinstance(ref, tuple) else ref
        assert np.asarray(first).shape in ((), (4,)), (name, np.shape(first))

    @pytest.mark.parametrize("name", NAMES)
    def test_only_the_shape_differs(self, name):
        """The numbers are still bit-identical, which is the point."""
        ref, got = self.calls(qpoint)[name](), self.calls(qpoint2)[name]()
        ref = ref if isinstance(ref, tuple) else (ref,)
        got = got if isinstance(got, tuple) else (got,)
        for r, g in zip(ref, got):
            assert_identical(np.ravel(r), np.ravel(g), name)

    @pytest.mark.parametrize("mod", IMPLS)
    @pytest.mark.parametrize("name", ["gmst", "lmst", "dipole"])
    def test_true_scalars_still_agree(self, mod, name):
        """Where the input is genuinely scalar, the two still match exactly."""
        args = {
            "gmst": (CTIME[0],),
            "lmst": (CTIME[0], LON[0]),
            "dipole": (CTIME[0], RA[0], DEC[0]),
        }[name]
        ref = getattr(qp(qpoint), name)(*args)
        got = getattr(qp(mod), name)(*args)
        assert np.ndim(ref) == np.ndim(got) == 0, name
        assert_identical(ref, got, name)


@pytest.mark.parametrize("mod", IMPLS)
class TestQuaternionConstruction:
    @pytest.mark.parametrize(
        "args",
        [
            pytest.param((1.0, 2.0, 3.0), id="scalar"),
            pytest.param(
                (np.linspace(-5, 5, 7), np.linspace(0, 3, 7), np.linspace(-90, 90, 7)),
                id="vector",
            ),
            pytest.param((1.0, np.linspace(0, 3, 7), 0.0), id="broadcast"),
        ],
    )
    def test_det_offset(self, mod, args):
        ref = np.asarray(qp(qpoint).det_offset(*args))
        got = np.asarray(qp(mod).det_offset(*args))
        assert_identical(ref, got, "det_offset")

    @pytest.mark.parametrize(
        "theta", [pytest.param(22.5, id="scalar"), pytest.param(HWP, id="vector")]
    )
    def test_hwp_quat(self, mod, theta):
        ref = np.asarray(qp(qpoint).hwp_quat(theta))
        got = np.asarray(qp(mod).hwp_quat(theta))
        assert_identical(ref, got, "hwp_quat")


@pytest.mark.parametrize("mod", IMPLS)
class TestOutputShapes:
    """True scalars degrade to scalars; an array keeps its axis."""

    def test_single_quat_is_flat(self, mod):
        assert np.asarray(qp(mod).det_offset(1.0, 2.0, 3.0)).shape == (4,)

    def test_boresight_keeps_leading_axis(self, mod):
        """azel2bore is the exception: it stays (1, 4) for a single sample."""
        qb = qp(mod).azel2bore(45.0, 45.0, None, None, LON[0], LAT[0], CTIME[0])
        assert np.asarray(qb).shape == (1, 4)

    def test_true_scalars_degrade(self, mod):
        q = qp(mod)
        qb = q.azel2bore(45.0, 45.0, None, None, LON[0], LAT[0], CTIME[0])
        out = q.bore2radec(q.det_offset(0.0, 0.0, 0.0), CTIME[0], np.asarray(qb)[0])
        assert all(np.ndim(v) == 0 for v in out)

    def test_a_length_one_array_keeps_its_axis(self, mod):
        """
        The distinction qpoint cannot make: it reads the output size, which
        is one either way. qpoint2 reads how the arguments arrived.
        """
        q = qp(mod)
        qb = np.atleast_2d(
            q.azel2bore(45.0, 45.0, None, None, LON[0], LAT[0], CTIME[0])
        )
        out = q.bore2radec(q.det_offset(0.0, 0.0, 0.0), CTIME[:1], qb)
        assert all(np.asarray(v).shape == (1,) for v in out)


class TestDecPolAxesIndependent:
    """
    qpoint2 treats the dec and pol output axes as independent, so it also
    offers (sindec, pa), which the C library has no entry point for.
    """

    def test_new_combination_agrees_with_the_others(self):
        q, qb = bore(qpoint2)
        q_off = q.det_offset(1.0, 2.0, 3.0)

        ra_a, sindec, pa_a = q.bore2radec(q_off, CTIME, qb, sindec=True, return_pa=True)
        q.reset_rates()
        ra_b, dec, pa_b = q.bore2radec(q_off, CTIME, qb, return_pa=True)

        assert identical(ra_a, ra_b)
        assert identical(pa_a, pa_b)
        assert np.allclose(sindec, np.sin(np.deg2rad(dec)), rtol=0, atol=1e-15)

    def test_c_library_rejects_it(self):
        q, qb = bore(qpoint)
        with pytest.raises(Exception):
            q.bore2radec(
                q.det_offset(1.0, 2.0, 3.0),
                CTIME,
                qb,
                sindec=True,
                return_pa=True,
            )


class TestZeroCopyContract:
    """
    qpoint2 never copies sample data. Arguments go straight from Python to
    the binding layer, which validates and rejects rather than converting.
    qpoint is not held to this: its ctypes layer copies freely.
    """

    @pytest.mark.parametrize(
        "bad, exc",
        [
            pytest.param(np.ones(3, dtype=np.float32), TypeError, id="wrong-dtype"),
            pytest.param(np.ones(3, dtype=np.int64), TypeError, id="int-dtype"),
            pytest.param(np.ones((3, 2))[:, 0], ValueError, id="non-contiguous"),
            pytest.param(np.ones((2, 3)), ValueError, id="wrong-rank"),
        ],
    )
    def test_rejects_rather_than_converts(self, bad, exc):
        """An existing array is never copied behind the caller's back."""
        ok = np.ones(3)
        with pytest.raises(exc):
            qp(qpoint2).det_offset(bad, ok, ok)

    def test_error_names_the_argument(self):
        with pytest.raises(TypeError, match="delta_el"):
            qp(qpoint2).det_offset(np.ones(3), np.ones(3, dtype=np.float32), 0.0)

    def test_length_mismatch_is_rejected(self):
        with pytest.raises(ValueError, match="length"):
            qp(qpoint2).det_offset(np.ones(3), np.ones(5), 0.0)


@pytest.mark.parametrize("mod", IMPLS)
class TestSequenceInputs:
    """
    Lists and tuples have no buffer to alias, so materializing one is the
    only option and is allowed -- unlike an array, which is never copied.
    """

    def test_list_matches_array(self, mod):
        q = qp(mod)
        ref = np.asarray(q.det_offset(np.array([0.0, 1.0, -1.0]), 0.0, 0.0))
        got = np.asarray(q.det_offset([0.0, 1.0, -1.0], 0.0, 0.0))
        assert identical(ref, got)

    def test_tuple_matches_array(self, mod):
        q = qp(mod)
        ref = np.asarray(q.det_offset(np.array([0.0, 1.0]), 0.0, 0.0))
        got = np.asarray(q.det_offset((0.0, 1.0), 0.0, 0.0))
        assert identical(ref, got)

    def test_int_list_is_converted(self, mod):
        """Unlike an int array, an int list has no buffer to preserve."""
        q = qp(mod)
        ref = np.asarray(q.det_offset(np.array([1.0, 2.0]), 0.0, 0.0))
        got = np.asarray(q.det_offset([1, 2], 0.0, 0.0))
        assert identical(ref, got)

    def test_nested_list_as_quaternion(self, mod):
        q = qp(mod)
        quats = [[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0]]
        ref = tuple(q.quat2radecpa(np.array(quats)))
        got = tuple(q.quat2radecpa(quats))
        assert_identical(ref, got, "quat2radecpa from list")

    def test_list_q_off(self, mod):
        out = []
        for q_off in (np.array([1.0, 0.0, 0.0, 0.0]), [1.0, 0.0, 0.0, 0.0]):
            q, qb = bore(mod)
            out.append(tuple(q.bore2radec(q_off, CTIME, qb)))
        assert_identical(out[0], out[1], "bore2radec with list q_off")


class TestScalarBroadcasting:
    """
    Scalars stand in for full-length columns and are resolved in C++ with a
    stride-0 view, so no broadcast array is ever materialized. The results
    must match passing the value repeated.
    """

    def test_scalar_matches_repeated(self):
        q = qp(qpoint2)
        a = np.asarray(q.det_offset(1.0, np.linspace(0, 3, 7), 0.0))
        b = np.asarray(q.det_offset(np.full(7, 1.0), np.linspace(0, 3, 7), np.zeros(7)))
        assert identical(a, b)

    def test_length_one_array_broadcasts(self):
        q = qp(qpoint2)
        a = np.asarray(q.det_offset(np.array([1.0]), np.linspace(0, 3, 7), 0.0))
        b = np.asarray(q.det_offset(1.0, np.linspace(0, 3, 7), 0.0))
        assert identical(a, b)

    def test_all_scalar_gives_one_sample(self):
        assert np.asarray(qp(qpoint2).det_offset(1.0, 2.0, 3.0)).shape == (4,)

    def test_scalar_ctime_matches_repeated(self):
        """The case that used to allocate a full-length array per scalar."""
        out = []
        for ct in (CTIME[0], np.full(N, CTIME[0])):
            q, qb = bore(qpoint2)
            q.reset_rates()
            out.append(tuple(q.bore2radec(q.det_offset(1.0, 2.0, 3.0), ct, qb)))
        assert_identical(out[0], out[1], "scalar ctime")

    def test_matches_qpoint_with_scalar_args(self):
        """qpoint broadcasts in Python; qpoint2 in C++. Same answer."""
        out = []
        for m in (qpoint, qpoint2):
            q = qp(m)
            out.append(
                np.asarray(q.azel2bore(AZ, EL[0], None, None, LON[0], LAT[0], CTIME))
            )
        assert identical(out[0], out[1])


# ---------------------------------------------------------------------------
# Mapmaking
#
# Bit-exact comparisons pin num_threads=1. The threaded reduction merges
# thread-local maps in whatever order threads reach the critical section, so
# it is not reproducible even against itself -- see TestThreadedReduction.
# ---------------------------------------------------------------------------

NSIDE_MAP = 16
NPIX_MAP = 12 * NSIDE_MAP * NSIDE_MAP
NS = 600
NDET = 4

MAP_CT = 1418662800.0 + np.arange(NS) / 10.0
MAP_AZ = np.linspace(0.0, 360.0, NS) % 360
MAP_EL = np.full(NS, 45.0) + 5 * np.sin(np.arange(NS) / 50.0)
MAP_LON = np.full(NS, 165.7)
MAP_LAT = np.full(NS, -77.6)

_mrng = np.random.default_rng(0)
TOD = _mrng.normal(size=(NDET, NS))
OFF = (
    _mrng.uniform(-3, 3, NDET),
    _mrng.uniform(-3, 3, NDET),
    _mrng.uniform(0, 180, NDET),
)
TOD2 = _mrng.normal(size=(2 * NDET, NS))
OFF2 = (
    _mrng.uniform(-3, 3, 2 * NDET),
    _mrng.uniform(-3, 3, 2 * NDET),
    _mrng.uniform(0, 180, 2 * NDET),
)
SOURCE_MAP = _mrng.normal(size=(3, NPIX_MAP))
DET_WEIGHTS = np.abs(np.random.default_rng(1).normal(size=(NDET, NS))) + 0.1
# the differencing kernel pairs the two halves, so it needs both
DET_WEIGHTS2 = np.abs(np.random.default_rng(2).normal(size=(2 * NDET, NS))) + 0.1


def qmap(mod, nthreads=1, **kwargs):
    """A QMap with pointing already initialized."""
    qm = mod.QMap(num_threads=nthreads, **kwargs)
    qb = qm.azel2bore(MAP_AZ, MAP_EL, None, None, MAP_LON, MAP_LAT, MAP_CT)
    qm.init_point(qb, ctime=MAP_CT)
    return qm, np.asarray(qm.det_offset(*OFF))


def hit_pixels():
    """Every pixel touched by any detector, for the partial-map tests."""
    qm = qpoint.QMap(nside=NSIDE_MAP)
    qb = qm.azel2bore(MAP_AZ, MAP_EL, None, None, MAP_LON, MAP_LAT, MAP_CT)
    hits = []
    for i in range(NDET):
        q = qpoint.QMap(nside=NSIDE_MAP)
        off = np.asarray(q.det_offset(OFF[0][i], OFF[1][i], OFF[2][i]))
        hits.append(np.asarray(q.bore2pix(off, MAP_CT, qb, nside=NSIDE_MAP)[0]))
    return np.unique(np.concatenate(hits)).astype(np.int64)


PARTIAL_PIX = hit_pixels()


@pytest.mark.parametrize("mod", IMPLS)
class TestTod2Map:
    @pytest.mark.parametrize("pol", [True, False])
    def test_from_tod(self, mod, pol):
        out = []
        for m in (qpoint, mod):
            qm, off = qmap(m, nside=NSIDE_MAP, pol=pol)
            out.append(tuple(qm.from_tod(off, tod=TOD.copy())))
        assert_identical(out[0], out[1], "from_tod")

    def test_vpol(self, mod):
        out = []
        for m in (qpoint, mod):
            qm, off = qmap(m, nside=NSIDE_MAP, pol=True, vpol=True)
            out.append(tuple(qm.from_tod(off, tod=TOD.copy())))
        assert_identical(out[0], out[1], "from_tod vpol")

    def test_vpol_weights(self, mod):
        """
        The V term picks up the per-sample weight separately from T, Q
        and U, so vpol has to be crossed with weights to reach it.
        """
        out = []
        for m in (qpoint, mod):
            qm, off = qmap(m, nside=NSIDE_MAP, pol=True, vpol=True)
            out.append(tuple(qm.from_tod(off, tod=TOD.copy(), weights=DET_WEIGHTS)))
        assert_identical(out[0], out[1], "from_tod vpol weighted")

    def test_flags(self, mod):
        flag = np.zeros((NDET, NS), dtype=np.uint8)
        flag[:, ::7] = 1
        out = []
        for m in (qpoint, mod):
            qm, off = qmap(m, nside=NSIDE_MAP)
            out.append(tuple(qm.from_tod(off, tod=TOD.copy(), flag=flag)))
        assert_identical(out[0], out[1], "from_tod flagged")

    def test_weights(self, mod):
        out = []
        for m in (qpoint, mod):
            qm, off = qmap(m, nside=NSIDE_MAP)
            out.append(tuple(qm.from_tod(off, tod=TOD.copy(), weights=DET_WEIGHTS)))
        assert_identical(out[0], out[1], "from_tod weighted")

    def test_hwp(self, mod):
        out = []
        for m in (qpoint, mod):
            qm, off = qmap(m, nside=NSIDE_MAP)
            qm.init_point(q_hwp=np.asarray(qm.hwp_quat(np.linspace(0, 180, NS))))
            out.append(tuple(qm.from_tod(off, tod=TOD.copy())))
        assert_identical(out[0], out[1], "from_tod hwp")

    def test_count_hits_only(self, mod):
        out = []
        for m in (qpoint, mod):
            qm, off = qmap(m, nside=NSIDE_MAP)
            out.append(np.asarray(qm.from_tod(off)))
        assert_identical(out[0], out[1], "from_tod hits only")

    def test_partial_map(self, mod):
        out = []
        for m in (qpoint, mod):
            qm, off = qmap(m)
            qm.init_dest(nside=NSIDE_MAP, pol=True, pixels=PARTIAL_PIX)
            res = tuple(qm.from_tod(off, tod=TOD.copy()))
            assert np.asarray(res[0]).shape == (3, len(PARTIAL_PIX))
            out.append(res)
        assert_identical(out[0], out[1], "from_tod partial")

    def test_diff_pairs(self, mod):
        out = []
        for m in (qpoint, mod):
            qm = m.QMap(nside=NSIDE_MAP, pol=True, num_threads=1)
            qb = qm.azel2bore(MAP_AZ, MAP_EL, None, None, MAP_LON, MAP_LAT, MAP_CT)
            qm.init_point(qb, ctime=MAP_CT)
            off = np.asarray(qm.det_offset(*OFF2))
            out.append(tuple(qm.from_tod(off, tod=TOD2.copy(), do_diff=True)))
        assert_identical(out[0], out[1], "from_tod do_diff")

    def test_diff_pairs_hwp(self, mod):
        """
        A HWP through the differencing kernel, which has its own copy of
        the per-detector rotation: it calls bore2det_hwp once per
        detector of the pair, and testing the HWP and the differencing
        separately reaches neither call.
        """
        out = []
        for m in (qpoint, mod):
            qm = m.QMap(nside=NSIDE_MAP, pol=True, num_threads=1)
            qb = qm.azel2bore(MAP_AZ, MAP_EL, None, None, MAP_LON, MAP_LAT, MAP_CT)
            qm.init_point(
                qb, ctime=MAP_CT, q_hwp=np.asarray(qm.hwp_quat(np.linspace(0, 180, NS)))
            )
            off = np.asarray(qm.det_offset(*OFF2))
            out.append(tuple(qm.from_tod(off, tod=TOD2.copy(), do_diff=True)))
        assert_identical(out[0], out[1], "from_tod do_diff hwp")

    def test_diff_pairs_vpol(self, mod):
        """
        Differencing into a T,Q,U,V map. The VPOL arms of both switches
        in tod2map1_diff -- the V row of vec, and the ten-row proj --
        are reached only by crossing vpol with the differencing, which
        the vpol and do_diff tests each miss on their own.
        """
        out = []
        for m in (qpoint, mod):
            qm = m.QMap(nside=NSIDE_MAP, pol=True, vpol=True, num_threads=1)
            qb = qm.azel2bore(MAP_AZ, MAP_EL, None, None, MAP_LON, MAP_LAT, MAP_CT)
            qm.init_point(qb, ctime=MAP_CT)
            off = np.asarray(qm.det_offset(*OFF2))
            out.append(tuple(qm.from_tod(off, tod=TOD2.copy(), do_diff=True)))
        assert_identical(out[0], out[1], "from_tod do_diff vpol")

    @pytest.mark.parametrize(
        "kwargs",
        [
            pytest.param({}, id="plain"),
            pytest.param({"weights": DET_WEIGHTS2}, id="weights"),
        ],
    )
    def test_diff_pairs_partial(self, mod, kwargs):
        """
        The differencing kernel on a partial map. It repixelizes both
        detectors of a pair separately, so there are two lookups and two
        missing-pixel branches, and neither the full-sky differencing
        test nor the non-differenced partial test reaches them.

        PARTIAL_PIX covers where OFF points, not OFF2, so the pair falls
        outside it -- which is the point: error_missing=False takes the
        skip branch for each detector, and the default raises instead.
        """
        out = []
        for m in (qpoint, mod):
            # no nside here: that initializes dest, and the partial map
            # has to be installed by init_dest instead
            qm = m.QMap(num_threads=1, error_missing=False)
            qb = qm.azel2bore(MAP_AZ, MAP_EL, None, None, MAP_LON, MAP_LAT, MAP_CT)
            qm.init_point(qb, ctime=MAP_CT)
            qm.init_dest(nside=NSIDE_MAP, pol=True, pixels=PARTIAL_PIX)
            off = np.asarray(qm.det_offset(*OFF2))
            out.append(tuple(qm.from_tod(off, tod=TOD2.copy(), do_diff=True, **kwargs)))
        assert_identical(out[0], out[1], "from_tod do_diff partial")

    def test_diff_pairs_partial_raises_by_default(self, mod):
        """The other half of that branch: both packages refuse."""
        for m in (qpoint, mod):
            qm = m.QMap(num_threads=1)
            qb = qm.azel2bore(MAP_AZ, MAP_EL, None, None, MAP_LON, MAP_LAT, MAP_CT)
            qm.init_point(qb, ctime=MAP_CT)
            qm.init_dest(nside=NSIDE_MAP, pol=True, pixels=PARTIAL_PIX)
            off = np.asarray(qm.det_offset(*OFF2))
            with pytest.raises(RuntimeError, match="out of bounds"):
                qm.from_tod(off, tod=TOD2.copy(), do_diff=True)


@pytest.mark.parametrize("mod", IMPLS)
class TestMap2Tod:
    @pytest.mark.parametrize(
        "options",
        [
            pytest.param({}, id="defaults"),
            pytest.param({"interp_pix": True}, id="interp"),
            pytest.param({"pix_order": "nest"}, id="nest"),
            pytest.param({"fast_pix": True}, id="fast-pix"),
        ],
    )
    def test_to_tod(self, mod, options):
        out = []
        for m in (qpoint, mod):
            qm, off = qmap(m, **options)
            qm.init_source(SOURCE_MAP, pol=True)
            out.append(np.asarray(qm.to_tod(off)))
        assert_identical(out[0], out[1], "to_tod")

    @pytest.mark.parametrize(
        "options",
        [
            pytest.param({}, id="defaults"),
            pytest.param({"interp_pix": True}, id="interp"),
        ],
    )
    def test_to_tod_hwp(self, mod, options):
        """
        Scanning a map with a HWP. from_tod had a HWP test and to_tod did
        not, so map2tod1's bore2det_hwp branch never ran.
        """
        out = []
        for m in (qpoint, mod):
            qm, off = qmap(m, **options)
            qm.init_point(q_hwp=np.asarray(qm.hwp_quat(np.linspace(0, 180, NS))))
            qm.init_source(SOURCE_MAP, pol=True)
            out.append(np.asarray(qm.to_tod(off)))
        assert_identical(out[0], out[1], "to_tod hwp")

    @pytest.mark.parametrize(
        "options",
        [
            pytest.param({}, id="defaults"),
            pytest.param(
                {"interp_pix": True, "error_missing": False, "interp_missing": True},
                id="interp-missing",
            ),
            pytest.param(
                {"interp_pix": True, "error_missing": False, "nan_missing": True},
                id="nan-missing",
            ),
        ],
    )
    def test_to_tod_partial(self, mod, options):
        """
        Scanning a partial source map. map2tod1 repixelizes through the
        hash, and with interp_pix it then has to decide what to do about
        neighbours falling outside the map -- reweight them away under
        interp_missing, or write NaN under nan_missing. The full-sky
        tests reach none of it.
        """
        out = []
        for m in (qpoint, mod):
            qm, off = qmap(m, **options)
            qm.init_source(
                SOURCE_MAP[:, PARTIAL_PIX],
                pol=True,
                pixels=PARTIAL_PIX,
                nside=NSIDE_MAP,
            )
            out.append(np.asarray(qm.to_tod(off)))
        assert_identical(out[0], out[1], "to_tod partial")

    def test_to_tod_partial_interp_raises_by_default(self, mod):
        """
        Without those flags it is an error in both packages: the
        interpolation wants neighbours the partial map does not have.
        """
        for m in (qpoint, mod):
            qm, off = qmap(m, interp_pix=True)
            qm.init_source(
                SOURCE_MAP[:, PARTIAL_PIX],
                pol=True,
                pixels=PARTIAL_PIX,
                nside=NSIDE_MAP,
            )
            with pytest.raises(RuntimeError, match="out of bounds"):
                qm.to_tod(off)

    def test_to_tod_flagged(self, mod):
        """Flagged samples are skipped on the scanning side too."""
        flag = np.zeros((NDET, NS), dtype=np.uint8)
        flag[:, ::7] = 1
        out = []
        for m in (qpoint, mod):
            qm, off = qmap(m)
            qm.init_source(SOURCE_MAP, pol=True)
            out.append(np.asarray(qm.to_tod(off, flag=flag)))
        assert_identical(out[0], out[1], "to_tod flagged")

    def test_roundtrip_recovers_a_constant(self, mod):
        """A constant T map scanned and rebinned must come back constant."""
        qm, off = qmap(mod, pol=False)
        qm.init_source(np.full(NPIX_MAP, 3.0), pol=False)
        tod = np.asarray(qm.to_tod(off))
        assert np.allclose(tod, 3.0)


@pytest.mark.parametrize("mod", IMPLS)
class TestMapInit:
    """
    init_dest and init_source decide map shapes, modes and nside. The
    shapes are written out rather than compared against another package,
    since the only other one left keeps its map state somewhere else.
    """

    @staticmethod
    def installed(qm):
        """
        Shapes of the maps init_dest actually installed. qpoint2's
        init_dest returns nothing, so the installed maps are the only
        thing to compare.
        """
        out = []
        if qm._dest.has_vec():
            out.append(np.asarray(qm._dest.get_vec()).shape)
        if qm._dest.has_proj():
            out.append(np.asarray(qm._dest.get_proj()).shape)
        return tuple(out)

    @staticmethod
    def make(mod, **kwargs):
        qm = mod.QMap()
        qm.init_dest(**kwargs)
        return qm

    @pytest.mark.parametrize(
        "kwargs, shapes",
        [
            pytest.param({"nside": 8}, ((3, 768), (6, 768)), id="pol-default"),
            pytest.param(
                {"nside": 8, "pol": False}, ((1, 768), (1, 768)), id="temperature"
            ),
            pytest.param(
                {"nside": 8, "pol": True, "vpol": True},
                ((4, 768), (10, 768)),
                id="vpol",
            ),
            pytest.param({"nside": 8, "vec": False}, ((6, 768),), id="proj-only"),
            pytest.param({"nside": 8, "proj": False}, ((3, 768),), id="vec-only"),
            pytest.param({}, ((3, 786432), (6, 786432)), id="nside-defaults-to-256"),
        ],
    )
    def test_default_allocation(self, mod, kwargs, shapes):
        assert self.installed(self.make(mod, **kwargs)) == shapes

    def test_nside_and_mode_inferred_from_vec(self, mod):
        for nrow, want_pol, want_vpol in [
            (1, False, False),
            (3, True, False),
            (4, True, True),
        ]:
            vec = np.zeros((nrow, 12 * 8 * 8))
            qm = mod.QMap()
            qm.init_dest(vec=vec)
            assert qm.dest_is_pol() is want_pol, nrow
            assert qm.dest_is_vpol() is want_vpol, nrow

    def test_both_false_raises(self, mod):
        with pytest.raises(ValueError, match="vec or proj"):
            mod.QMap().init_dest(nside=8, vec=False, proj=False)

    def test_second_init_raises(self, mod):
        qm = mod.QMap(nside=8)
        with pytest.raises(RuntimeError, match="already initialized"):
            qm.init_dest(nside=8)

    def test_reset_allows_reinit(self, mod):
        qm = mod.QMap(nside=8)
        qm.init_dest(nside=16, reset=True)
        assert qm.dest_is_init()
        assert self.installed(qm)[0] == (3, 12 * 16 * 16)

    def test_update_without_maps_zeros_the_dest(self, mod):
        """
        init_dest(update=True) with no vec/proj given resets the
        accumulators in place, which is how a scan chunk starts over.
        """
        qm = mod.QMap(nside=8)
        np.asarray(qm._dest.get_vec())[:] = 1.0
        np.asarray(qm._dest.get_proj())[:] = 1.0
        qm.init_dest(nside=8, update=True)
        assert self.installed(qm) == ((3, 12 * 8 * 8), (6, 12 * 8 * 8))
        assert np.allclose(np.asarray(qm._dest.get_vec()), 0)
        assert np.allclose(np.asarray(qm._dest.get_proj()), 0)

    def test_partial_requires_nside(self, mod):
        pix = np.arange(50, dtype=np.int64)
        with pytest.raises(ValueError, match="nside"):
            mod.QMap().init_dest(pixels=pix)

    def test_partial_shapes(self, mod):
        pix = np.arange(50, dtype=np.int64)
        assert self.installed(self.make(mod, nside=8, pixels=pix)) == ((3, 50), (6, 50))

    def test_partial_map_with_few_pixels_keeps_orientation(self, mod):
        """
        A (3, 2) map covering two pixels is taller than it is wide, so a
        "transpose if taller than wide" rule would flip it into three bogus
        pixels.
        """
        pix = np.array([0, 1], dtype=np.int64)
        qm = self.make(mod, nside=8, pol=True, pixels=pix)
        assert self.installed(qm) == ((3, 2), (6, 2))

    def test_transposed_input_raises(self, mod):
        """
        A full-sky map handed over as (npix, nrow) used to be flipped into
        shape. Nothing in the API produces one, and guessing is what read a
        small partial map as a pile of bogus pixels, so it is an error --
        of npix, since 3 columns are not a whole sky.
        """
        npix = 12 * 8 * 8
        with pytest.raises(ValueError):
            self.make(mod, vec=np.zeros((npix, 3)))

    def test_transposed_partial_input_raises(self, mod):
        """
        A partial map cannot be judged by its shape, so the pixel list is
        what catches it. The two packages word the complaint differently.
        """
        pix = np.arange(6, dtype=np.int64)
        with pytest.raises(ValueError):
            self.make(mod, nside=8, pixels=pix, vec=np.zeros((6, 3)))

    @pytest.mark.parametrize(
        "shape",
        [pytest.param((12 * 8 * 8,), id="1d"), pytest.param((3, 12 * 8 * 8), id="2d")],
    )
    def test_source_accepts_1d_and_2d(self, mod, shape):
        m = np.zeros(shape)
        qm = mod.QMap()
        qm.init_source(m, pol=len(shape) > 1)
        assert qm.source_is_init()

    def test_source_second_init_raises(self, mod):
        qm = mod.QMap()
        qm.init_source(np.zeros((3, 12 * 8 * 8)))
        with pytest.raises(RuntimeError, match="already initialized"):
            qm.init_source(np.zeros((3, 12 * 8 * 8)))

    def test_source_reset_and_update(self, mod):
        npix = 12 * 8 * 8
        qm = mod.QMap()
        qm.init_source(np.zeros((3, npix)))
        qm.init_source(np.ones((3, npix)), reset=True)
        assert qm.source_is_init()
        qm.init_source(np.full((3, npix), 2.0), update=True)
        assert qm.source_is_init()

    def test_source_partial_requires_nside(self, mod):
        with pytest.raises(ValueError, match="nside"):
            mod.QMap().init_source(
                np.zeros((3, 50)), pixels=np.arange(50, dtype=np.int64)
            )

    def test_copy_does_not_alias_the_input(self, mod):
        npix = 12 * 8 * 8
        vec = np.zeros((3, npix))
        qm = self.make(mod, vec=vec, copy=True)
        assert not np.shares_memory(np.asarray(qm._dest.get_vec()), vec)

    def test_without_copy_the_input_is_aliased(self, mod):
        """The default keeps the caller's buffer, which is the point."""
        npix = 12 * 8 * 8
        vec = np.require(np.zeros((3, npix)), float, ["A", "C"])
        qm = self.make(mod, vec=vec)
        assert np.shares_memory(np.asarray(qm._dest.get_vec()), vec)

    def test_dest_mode_queries_need_init(self, mod):
        qm = mod.QMap()
        for fn in ("dest_is_pol", "dest_is_vpol"):
            with pytest.raises(RuntimeError, match="not initialized"):
                getattr(qm, fn)()
        for fn in ("source_is_pol", "source_is_vpol"):
            with pytest.raises(RuntimeError, match="not initialized"):
                getattr(qm, fn)()


@pytest.mark.parametrize("mod", IMPLS)
class TestPolarizationReporting:
    """
    is_pol and is_vpol answer what the map contains, which is not the
    question the accumulation kernels ask.

    The kernels branch on `at_least(vec_mode, Pol)`, an ordering test that
    mirrors the C's `>=` on the same enum and has to keep doing so. But
    that enum is ordered T, Pol, VPol, D1, D1Pol, D2, D2Pol, which puts the
    *unpolarized* derivative modes above Pol -- so read as a predicate it
    calls a T-plus-derivatives map polarized. Only Pol, VPol, D1Pol and
    D2Pol carry Q and U.

    The other half is that a dest map can be proj-only, with no vec mode to
    read at all, so the proj has to answer when the vec cannot.

    Both are checked against `qpoint`, which gets these right from a
    different direction: it lists the polarized modes outright, and its
    dest query consults proj_mode too.
    """

    NPIX = 12 * 8 * 8

    def z(self, nrow):
        return np.zeros((nrow, self.NPIX))

    @pytest.mark.parametrize(
        "nrow,kwargs",
        [
            (1, {}),
            (3, {}),
            (3, {"pol": False}),  # D1: three rows, not polarized
            (4, {"vpol": True}),
            (6, {}),  # D2: six rows, not polarized
            (9, {}),  # D1Pol: polarized
            (18, {}),  # D2Pol: polarized
        ],
    )
    def test_source(self, mod, nrow, kwargs):
        got = []
        for m in (qpoint, mod):
            qm = m.QMap()
            qm.init_source(self.z(nrow), **kwargs)
            got.append((qm.source_is_pol(), qm.source_is_vpol()))
        assert got[0] == got[1], f"{nrow} rows {kwargs}"

    @pytest.mark.parametrize(
        "kwargs",
        [
            {},
            {"pol": False},
            {"vpol": True},
            {"vec": 1},
            {"vec": 3},
            {"vec": 4},
            {"vec": False, "proj": 1},
            {"vec": False, "proj": 6},
            {"vec": False, "proj": 10},
        ],
    )
    def test_dest(self, mod, kwargs):
        kwargs = {
            k: self.z(v) if k in ("vec", "proj") and v is not False else v
            for k, v in kwargs.items()
        }
        got = []
        for m in (qpoint, mod):
            qm = m.QMap()
            qm.init_dest(nside=8, **kwargs)
            got.append((qm.dest_is_pol(), qm.dest_is_vpol()))
        assert got[0] == got[1], str(kwargs.keys())

    def test_a_derivative_map_is_not_polarized(self, mod):
        """
        The case the ordering got wrong, stated directly so the reason
        survives even if the parametrized comparisons above are reworked.
        """
        qm = mod.QMap()
        qm.init_source(self.z(3), pol=False)
        assert qm.source_is_pol() is False
        qm = mod.QMap()
        qm.init_source(self.z(6))
        assert qm.source_is_pol() is False

    def test_a_proj_only_dest_still_reports_its_mode(self, mod):
        """The vec is switched off, so the proj is the only thing to read."""
        for nproj, want_pol, want_vpol in [
            (1, False, False),
            (6, True, False),
            (10, True, True),
        ]:
            qm = mod.QMap()
            qm.init_dest(nside=8, vec=False, proj=self.z(nproj))
            assert qm.dest_is_pol() is want_pol, nproj
            assert qm.dest_is_vpol() is want_vpol, nproj


class TestInitDestReturnValue:
    """
    qpoint2 drops init_dest's return value; the maps are read back with
    get_vec/get_proj, or via from_tod. qpoint still returns them, so this
    is a deliberate divergence rather than an oversight.
    """

    def test_qpoint2_returns_nothing(self):
        assert qpoint2.QMap().init_dest(nside=8) is None

    def test_qpoint_still_returns_maps(self):
        assert qpoint.QMap().init_dest(nside=8) is not None


@pytest.mark.parametrize("mod", IMPLS)
class TestUpdateComponents:
    """
    vec and proj each take three forms, at init and at update alike: a map
    installs it, None asks for one of zeros, and False says the map has no
    such component.

    An update may not change the row count of a component that is already
    installed. update_vec used to re-decide the vec mode from pol/vpol and
    the proj mode from the row count, while num_vec/num_proj kept indexing
    the old shape -- so a 1-row replacement for a TQU dest left
    qp_reshape_map building three row pointers into a one-row buffer, and a
    partly filled accumulator changed meaning halfway through a scan.
    """

    npix = 12 * 8 * 8

    @staticmethod
    def shapes(qm):
        vec = np.asarray(qm._dest.get_vec()).shape if qm._dest.has_vec() else None
        proj = np.asarray(qm._dest.get_proj()).shape if qm._dest.has_proj() else None
        return vec, proj

    @pytest.mark.parametrize("kw", ["vec", "proj"])
    def test_rejects_a_new_row_count(self, mod, kw):
        qm = mod.QMap(nside=8)
        with pytest.raises(ValueError, match="does not match"):
            qm.init_dest(nside=8, update=True, **{kw: np.ones((1, self.npix))})

    def test_accepts_a_matching_map(self, mod):
        qm = mod.QMap(nside=8)
        qm.init_dest(nside=8, update=True, vec=np.ones((3, self.npix)))
        assert np.allclose(np.asarray(qm._dest.get_vec()), 1.0)

    def test_zeros_by_replacing_the_buffer(self, mod):
        """
        A component left out comes back as a fresh map of zeros rather than
        being zeroed where it sits, so whatever the caller is still holding
        keeps its data.
        """
        qm = mod.QMap(nside=8)
        vec, proj = qm._dest.get_vec(), qm._dest.get_proj()
        vec[:] = 1.0
        proj[:] = 1.0
        qm.init_dest(nside=8, update=True)
        for old, new in [(vec, qm._dest.get_vec()), (proj, qm._dest.get_proj())]:
            assert not np.shares_memory(old, new)
            assert not np.asarray(new).any()
            assert old.all()

    @pytest.mark.parametrize("kw", ["vec", "proj"])
    def test_false_disables_a_component(self, mod, kw):
        qm = mod.QMap(nside=8)
        qm.init_dest(nside=8, update=True, **{kw: False})
        assert self.shapes(qm)[kw == "proj"] is None
        assert self.shapes(qm)[kw != "proj"] is not None

    @pytest.mark.parametrize("kw", ["vec", "proj"])
    def test_update_leaves_a_disabled_component_off(self, mod, kw):
        """
        init_dest's None means "not supplied", so zeroing the accumulators
        of a proj-only dest does not grow a vec back.
        """
        qm = mod.QMap()
        qm.init_dest(nside=8, **{kw: False})
        qm.init_dest(nside=8, update=True)
        assert self.shapes(qm)[kw == "proj"] is None

    @pytest.mark.parametrize("kw", ["vec", "proj"])
    def test_none_re_enables_at_the_mode_shape(self, mod, kw):
        """
        The binding's own None does ask for a fresh component. pol/vpol are
        remembered from setup, so one that was switched off comes back at
        the shape it would have had, zeroed.
        """
        want = {"vec": (3, self.npix), "proj": (6, self.npix)}[kw]
        qm = mod.QMap()
        qm.init_dest(nside=8, **{kw: False})
        getattr(qm._dest, "update_" + kw)(None)
        assert self.shapes(qm)[kw == "proj"] == want
        assert not np.asarray(getattr(qm._dest, "get_" + kw)()).any()

    def test_disabling_both_is_rejected(self, mod):
        qm = mod.QMap(nside=8)
        with pytest.raises(ValueError, match="vec or proj"):
            qm.init_dest(nside=8, update=True, vec=False, proj=False)

    @pytest.mark.parametrize("kw", ["vec", "proj"])
    def test_disabling_the_last_component_is_rejected(self, mod, kw):
        other = {"vec": "proj", "proj": "vec"}[kw]
        qm = mod.QMap()
        qm.init_dest(nside=8, **{other: False})
        with pytest.raises(ValueError, match="vec or proj"):
            qm.init_dest(nside=8, update=True, **{kw: False})

    def test_a_cycled_component_still_accumulates(self, mod):
        """
        In the C a component owns a table of row pointers into its
        buffer, sized for its row count, so switching it off has to free
        that table -- qp_reshape_map only builds a new one when the old is
        gone, and would otherwise leave the extra rows pointing off the
        end of the buffer.
        """
        ref, off = qmap(mod)
        ref.init_dest(nside=NSIDE_MAP, pol=True)
        want = [np.asarray(x) for x in ref.from_tod(off, tod=TOD.copy())]

        qm, off = qmap(mod)
        qm.init_dest(nside=NSIDE_MAP, pol=True)
        for kw in ("vec", "proj"):
            getattr(qm._dest, "update_" + kw)(False)
            getattr(qm._dest, "update_" + kw)(None)
        got = [np.asarray(x) for x in qm.from_tod(off, tod=TOD.copy())]
        assert want[0].any()
        assert all(np.array_equal(a, b) for a, b in zip(want, got))

    def test_a_default_follows_the_supplied_component(self, mod):
        """
        A supplied proj sizes the default vec, so the two always describe
        the same number of map components. Both used to size the default
        from pol, which for a 1-row proj gave a 3-row vec -- shapes the
        accumulation kernels cannot use together.
        """
        qm = mod.QMap()
        qm.init_dest(proj=np.zeros((1, self.npix)))
        assert self.shapes(qm) == ((1, self.npix), (1, self.npix))


@pytest.mark.parametrize("mod", IMPLS)
class TestDerivativeMapModes:
    """
    Row count fixes the vec mode, and the mode fixes which rows map2tod
    reads. Getting the two out of step reads off the end of the map, which
    is what qp_num_maps did for D2 and D1_POL.

    D2   = T with 1st and 2nd derivatives      -> 6 rows
    D1POL = (T, Q, U) with 1st derivatives     -> 9 rows
    """

    def _tod(self, mod, source):
        qm = mod.QMap(num_threads=1)
        qb = qm.azel2bore(MAP_AZ, MAP_EL, None, None, MAP_LON, MAP_LAT, MAP_CT)
        qm.init_point(qb, ctime=MAP_CT)
        qm.init_source(source, pol=True, vpol=len(source) == 4)
        return np.asarray(qm.to_tod(np.asarray(qm.det_offset(0.0, 0.0, 0.0))))

    @pytest.mark.parametrize("nrow", [1, 3, 4, 6, 9, 18])
    def test_every_row_is_read(self, mod, nrow):
        """
        Perturbing any row must change the tod. A row that never matters
        means the mode was inferred too small; the converse -- reading past
        the end -- shows up as a mismatch against the reference.
        """
        rng = np.random.default_rng(1)
        base = rng.normal(size=(nrow, NPIX_MAP))
        ref = self._tod(mod, base.copy())
        for row in range(nrow):
            bumped = base.copy()
            bumped[row] += 100.0
            assert not np.allclose(
                ref, self._tod(mod, bumped)
            ), "row {} of a {}-row map is never read".format(row, nrow)

    @pytest.mark.parametrize("nrow", [1, 3, 4, 6, 9, 18])
    def test_matches_reference(self, mod, nrow):
        rng = np.random.default_rng(1)
        m = rng.normal(size=(nrow, NPIX_MAP))
        assert_identical(self._tod(qpoint, m), self._tod(mod, m), "to_tod")


@pytest.mark.parametrize("mod", IMPLS)
class TestMapSolve:
    def test_solve_map(self, mod):
        out = []
        for m in (qpoint, mod):
            qm, off = qmap(m, nside=NSIDE_MAP)
            qm.from_tod(off, tod=TOD.copy())
            out.append(np.asarray(qm.solve_map()))
        assert_identical(out[0], out[1], "solve_map")

    def test_proj_cond(self, mod):
        out = []
        for m in (qpoint, mod):
            qm, off = qmap(m, nside=NSIDE_MAP)
            qm.from_tod(off, tod=TOD.copy())
            out.append(np.asarray(qm.proj_cond()))
        assert_identical(out[0], out[1], "proj_cond")


@pytest.mark.parametrize("mod", IMPLS)
class TestMapErrors:
    def test_missing_pixel_raises(self, mod):
        """A partial map that does not cover the scan is an error by default."""
        qm, off = qmap(mod)
        qm.init_dest(nside=NSIDE_MAP, pol=True, pixels=np.array([0, 1], dtype=np.int64))
        with pytest.raises(RuntimeError, match="out of bounds"):
            qm.from_tod(off, tod=TOD.copy())

    def test_missing_pixel_skipped_when_allowed(self, mod):
        qm, off = qmap(mod, error_missing=False)
        qm.init_dest(nside=NSIDE_MAP, pol=True, pixels=np.array([0, 1], dtype=np.int64))
        vec, proj = qm.from_tod(off, tod=TOD.copy())
        assert np.asarray(proj)[0].sum() == 0

    def test_missing_pixel_raises_when_threaded(self, mod):
        """An exception raised inside the OpenMP region must still surface."""
        qm, off = qmap(mod, nthreads=8)
        qm.init_dest(nside=NSIDE_MAP, pol=True, pixels=np.array([0, 1], dtype=np.int64))
        with pytest.raises(RuntimeError, match="out of bounds"):
            qm.from_tod(off, tod=TOD.copy())


HAS_OPENMP = qpoint2._libqpoint2.HAS_OPENMP


class TestThreadedReduction:
    """
    With OpenMP, tod2map merges thread-local maps in arrival order, so the
    result is not bit-reproducible -- in qpoint either. Hits are integer
    counts and stay exact regardless.

    OpenMP is off on macOS: linking a runtime collides with the one healpy
    bundles and segfaults whichever loads second. These tests still run
    there, they just compare a serial result against itself.
    """

    def _run(self, mod, nthreads):
        qm, off = qmap(mod, nthreads=nthreads, nside=NSIDE_MAP)
        return [np.asarray(x) for x in qm.from_tod(off, tod=TOD.copy())]

    def test_serial_is_bit_identical_to_qpoint(self):
        a, b = self._run(qpoint, 1), self._run(qpoint2, 1)
        assert all(identical(x, y) for x, y in zip(a, b))

    def test_threaded_matches_serial_to_rounding(self):
        serial, threaded = self._run(qpoint2, 1), self._run(qpoint2, 8)
        for x, y in zip(serial, threaded):
            assert np.allclose(x, y, rtol=0, atol=1e-9)

    def test_hits_are_exact_under_threading(self):
        serial, threaded = self._run(qpoint2, 1), self._run(qpoint2, 8)
        assert identical(serial[1][0], threaded[1][0])

    def test_threading_is_a_noop_without_openmp(self):
        if HAS_OPENMP:
            pytest.skip("OpenMP is enabled; the reduction is not reproducible")
        serial, threaded = self._run(qpoint2, 1), self._run(qpoint2, 8)
        assert all(identical(x, y) for x, y in zip(serial, threaded))

    def test_qpoint_is_equally_nonreproducible(self):
        """Any divergence is inherited, not introduced by the C++ port."""
        if not HAS_OPENMP:
            pytest.skip("OpenMP disabled; nothing to diverge")
        s1, t1 = self._run(qpoint, 1), self._run(qpoint, 8)
        s3, t3 = self._run(qpoint2, 1), self._run(qpoint2, 8)
        d1 = max(np.max(np.abs(x - y)) for x, y in zip(s1, t1))
        d3 = max(np.max(np.abs(x - y)) for x, y in zip(s3, t3))
        assert d1 > 0 and d3 > 0
        assert np.isclose(d1, d3, rtol=10)


GROUPS = ("rates", "options", "weather", "params")


def flat_params(mod):
    """
    Every parameter of a freshly built QPoint, as one flat dict. Both
    packages group them, so this flattens the groups to compare values
    without the grouping getting in the way; test_same_grouping is what
    checks the grouping itself.

    thread_num is dropped: qpoint2 has no such parameter, because in the C
    it only ever fed debug printing.
    """
    out = {}
    for group in GROUPS:
        out.update(qp(mod).get(group))
    out.pop("thread_num", None)
    return out


def qpoint_defaults():
    return flat_params(qpoint)


@pytest.mark.parametrize("mod", IMPLS)
class TestStandaloneRefraction:
    """
    tools.refraction, which takes the weather as arguments rather than
    from the stored state, so there is nothing to initialize and no rate
    cache to carry over.
    """

    EL = np.linspace(5.0, 85.0, N)
    TEMP = np.linspace(-40.0, 30.0, N)
    PRESS = np.linspace(500.0, 1050.0, N)
    HUM = np.linspace(0.0, 1.0, N)
    FREQ = np.linspace(30.0, 300.0, N)

    @pytest.mark.parametrize(
        "args",
        [
            pytest.param((45.0, 20.0, 1013.25, 0.5, 150.0), id="scalar"),
            pytest.param((EL, TEMP, PRESS, HUM, FREQ), id="vector"),
            pytest.param((EL, 20.0, 1013.25, 0.5, 150.0), id="scalar-weather"),
            # el > 90 folds back through the zenith, a branch of its own
            pytest.param(
                (np.linspace(85.0, 175.0, N), 20.0, 1013.25, 0.5, 150.0),
                id="past-zenith",
            ),
            # a vacuum, where the correction collapses to nearly nothing
            pytest.param((EL, 20.0, 0.0, 0.0, 150.0), id="no-atmosphere"),
        ],
    )
    def test_refraction(self, mod, args):
        ref = qpoint.tools.refraction(*args)
        got = mod.tools.refraction(*args)
        assert_identical(ref, got, "refraction")

    def test_default_frequency_agrees(self, mod):
        """The 150 GHz default is the Python layer's in both packages."""
        ref = qpoint.tools.refraction(self.EL, 20.0, 1013.25, 0.5)
        got = mod.tools.refraction(self.EL, 20.0, 1013.25, 0.5)
        assert_identical(ref, got, "refraction default freq")


@pytest.mark.parametrize("mod", IMPLS)
class TestParameterAPI:
    def test_same_parameter_set(self, mod):
        assert set(qpoint_defaults()) == set(flat_params(mod))

    def test_same_grouping(self, mod):
        """
        get() with no arguments returns the groups, not a flat dict, and
        puts the same parameter in the same group as qpoint does -- down
        to the order the groups come out in, since callers iterate it.
        """
        ref, got = qp(qpoint).get(), qp(mod).get()
        assert list(ref) == list(got) == list(GROUPS)
        for group in GROUPS:
            expected = set(ref[group]) - {"thread_num"}
            assert expected == set(got[group]), group

    def test_a_group_among_keys_stays_nested(self, mod):
        """
        Asking for a group alongside a plain key nests the group under its
        own name rather than merging it in, which is what qpoint does and
        what keeps the result unambiguous.
        """
        ref = qp(qpoint).get("weather", "accuracy")
        got = qp(mod).get("weather", "accuracy")
        assert set(got) == {"weather", "accuracy"}
        assert ref == got

    def test_same_defaults(self, mod):
        assert qpoint_defaults() == flat_params(mod)

    def test_same_default_types(self, mod):
        """
        Value equality is not enough: bool is a subclass of int, so True == 1
        and a flag/count mix-up compares equal.
        """
        ref, got = qpoint_defaults(), flat_params(mod)
        mismatched = {
            k: (type(ref[k]).__name__, type(got[k]).__name__)
            for k in ref
            if type(ref[k]) is not type(got[k])
        }
        assert not mismatched

    def test_num_threads_is_a_count(self, mod):
        """Not a flag -- this is the case the equality check above misses."""
        val = qp(mod).get("num_threads")
        assert type(val) is int
        assert val == qp(qpoint).get("options")["num_threads"]

    @pytest.mark.parametrize(
        "key, val",
        [
            ("accuracy", "low"),
            ("polconv", "iau"),
            ("pix_order", "nest"),
            ("fast_math", True),
            ("rate_npb", "never"),
            ("rate_aaber", 12.5),
            ("temperature", 5.0),
            ("dut1", 0.3),
        ],
    )
    def test_roundtrip(self, mod, key, val):
        """Setting a parameter and reading it back agrees with qpoint."""
        a, b = qp(qpoint), qp(mod)
        a.set(**{key: val})
        b.set(**{key: val})
        assert a.get(key) == b.get(key)

    def test_set_rejects_an_unknown_key(self, mod):
        """
        A divergence: qpoint drops a keyword it does not recognize, so a
        typo silently does nothing. qpoint2's parameter names are a fixed
        set and set() holds to it, which is also what lets the per-call
        keywords be split from a method's own arguments by name.
        """
        qp(qpoint).set(definitely_not_a_parameter=1)  # ignored, no raise
        with pytest.raises(KeyError):
            qp(mod).set(definitely_not_a_parameter=1)

    def test_a_method_rejects_an_unknown_keyword(self, mod):
        """The same holds through the per-call override path."""
        qp(qpoint).gmst(CTIME[0], definitely_not_a_parameter=1)
        with pytest.raises(TypeError):
            qp(mod).gmst(CTIME[0], definitely_not_a_parameter=1)

    def test_settings_raises_before_it_changes_anything(self, mod):
        q = qp(mod)
        before = q.get("accuracy")
        with pytest.raises(KeyError):
            with q.settings(accuracy="low", definitely_not_a_parameter=1):
                pass
        assert q.get("accuracy") == before

    def test_get_raises_on_unknown_key(self, mod):
        with pytest.raises(KeyError):
            qp(mod).get("definitely_not_a_parameter")

    def test_settings_restores(self, mod):
        q = qp(mod, accuracy="low")
        with q.settings(accuracy="high"):
            assert q.get("accuracy") == "high"
        assert q.get("accuracy") == "low"
