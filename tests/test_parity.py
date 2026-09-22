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
