"""Tests for qpoint.QPoint coordinate conversion class."""

import numpy as np
import pytest
import qpoint
from qpoint._libqpoint import libqp

# Reference observing parameters
CTIME = 1418662800.0  # 2014-12-15 12:00 UTC
LON = 165.7  # McMurdo Station longitude
LAT = -77.6  # McMurdo Station latitude
N = 20
AZ = 100.0 + np.linspace(0, 10, N)
EL = 32.0 * np.ones(N)
CTIMES = CTIME + np.arange(N, dtype=float)


qpoint2 = pytest.importorskip("qpoint2")

# The two packages implement the same API over different cores, so every
# behavioural test here runs against both. Where they are meant to differ
# the test says so; tests/test_parity.py is what pins them to the same
# numbers, and this is what pins each to the documented behaviour.
IMPLS = [qpoint, qpoint2]


@pytest.fixture(params=IMPLS, ids=lambda m: m.__name__)
def mod(request):
    return request.param


@pytest.fixture
def qp():
    return qpoint.QPoint(mean_aber=True, accuracy="low")


# ---------------------------------------------------------------------------
# Initialization
# ---------------------------------------------------------------------------


class TestInit:
    def test_default_init(self, mod):
        q = mod.QPoint()
        assert q is not None

    def test_kwargs_in_init(self, mod):
        q = mod.QPoint(accuracy="low", fast_math=True, mean_aber=True)
        assert q.get("accuracy") == "low"
        assert q.get("fast_math") is True
        assert q.get("mean_aber") is True

    def test_del(self, mod):
        q = mod.QPoint()
        del q  # should not raise


# ---------------------------------------------------------------------------
# set / get
# ---------------------------------------------------------------------------


class TestSetGet:
    def test_get_all_is_grouped(self, qp):
        """
        The parameters come back under 'rates', 'options', 'weather' and
        'params' rather than as one flat dict, so a caller can hand a
        whole group to set().
        """
        state = qp.get()
        assert list(state) == ["rates", "options", "weather", "params"]
        assert all(isinstance(group, dict) for group in state.values())
        assert "rate_npb" in state["rates"]
        assert "accuracy" in state["options"]
        # the parameters live in the groups, not beside them
        assert "accuracy" not in state

    def test_get_all_round_trips_through_set(self, qp):
        """Each group is accepted back by set(), which is the point of it."""
        state = qp.get()
        for group in state.values():
            qp.set(**group)
        assert qp.get() == state

    def test_get_group_rates(self, qp):
        rates = qp.get("rates")
        assert isinstance(rates, dict)
        assert "rate_lonlat" in rates

    def test_get_group_options(self, qp):
        opts = qp.get("options")
        assert isinstance(opts, dict)
        assert "accuracy" in opts

    def test_get_multiple_keys(self, qp):
        result = qp.get("accuracy", "polconv")
        assert isinstance(result, dict)
        assert "accuracy" in result and "polconv" in result

    def test_get_single_key(self, qp):
        val = qp.get("accuracy")
        assert val in ("high", "low")

    def test_set_get_accuracy(self, qp):
        qp.set(accuracy="low")
        assert qp.get("accuracy") == "low"
        qp.set(accuracy="high")
        assert qp.get("accuracy") == "high"

    def test_set_get_polconv(self, qp):
        qp.set(polconv="iau")
        assert qp.get("polconv") == "iau"
        qp.set(polconv="cosmo")
        assert qp.get("polconv") == "cosmo"

    def test_set_get_pix_order(self, qp):
        qp.set(pix_order="nest")
        assert qp.get("pix_order") == "nest"
        qp.set(pix_order="ring")
        assert qp.get("pix_order") == "ring"

    @pytest.mark.parametrize(
        "opt",
        [
            "mean_aber",
            "fast_aber",
            "fast_math",
            "interp_pix",
            "fast_pix",
            "error_missing",
            "nan_missing",
            "interp_missing",
        ],
    )
    def test_set_get_bool_option(self, qp, opt):
        qp.set(**{opt: True})
        assert qp.get(opt) is True
        qp.set(**{opt: False})
        assert qp.get(opt) is False

    @pytest.mark.parametrize(
        "rate",
        [
            "rate_lonlat",
            "rate_npb",
            "rate_erot",
            "rate_daber",
            "rate_aaber",
            "rate_wobble",
            "rate_dut1",
            "rate_ref",
        ],
    )
    def test_set_get_rate_strings(self, qp, rate):
        for val in ("always", "once", "never"):
            qp.set(**{rate: val})
            assert qp.get(rate) == val

    @pytest.mark.parametrize(
        "rate",
        ["rate_lonlat", "rate_npb", "rate_erot"],
    )
    def test_set_get_rate_float(self, qp, rate):
        qp.set(**{rate: 60.0})
        assert qp.get(rate) == pytest.approx(60.0)

    def test_set_get_weather(self, qp):
        qp.set(temperature=20.0, pressure=1013.25, humidity=0.5, frequency=150.0)
        assert qp.get("temperature") == pytest.approx(20.0)
        assert qp.get("pressure") == pytest.approx(1013.25)
        assert qp.get("humidity") == pytest.approx(0.5)
        assert qp.get("frequency") == pytest.approx(150.0)

    def test_set_get_dut1(self, qp):
        qp.set(dut1=0.1)
        assert qp.get("dut1") == pytest.approx(0.1)

    def test_set_get_ref_delta(self, qp):
        qp.set(ref_delta=0.05)
        assert qp.get("ref_delta") == pytest.approx(0.05)

    def test_set_unknown_key_ignored(self, qp):
        qp.set(nonexistent_key=123)  # should not raise

    def test_get_unknown_key_raises(self, qp):
        with pytest.raises(KeyError):
            qp.get("nonexistent_key")


# ---------------------------------------------------------------------------
# reset_rates / reset_inv_rates
# ---------------------------------------------------------------------------


class TestBulletinA:
    """
    IERS Bulletin A loading. The file path had two faults that concealed
    each other: the columns came out rotated, and numpy's unpack=True
    returns strided rows that ctypes would not take -- so the rotation was
    never reached.
    """

    COLUMNS = ["mjd", "dut1", "x", "y"]
    # distinct constants, so a rotation is visible rather than plausible
    VALUES = {"dut1": 0.11, "x": 0.22, "y": 0.33}
    MJD0 = 57000
    NDAY = 40

    def _write(self, tmp_path, order):
        col = {
            "mjd": np.arange(self.MJD0, self.MJD0 + self.NDAY, dtype=float),
            **{k: np.full(self.NDAY, v) for k, v in self.VALUES.items()},
        }
        path = tmp_path / "bulletin.txt"
        np.savetxt(path, np.column_stack([col[c] for c in order]), fmt="%.6f")
        return str(path)

    def test_round_trip(self, mod, tmp_path):
        q = mod.QPoint()
        path = self._write(tmp_path, self.COLUMNS)
        mjd, dut1, x, y = q.load_bulletin_a(path)
        assert mjd[0] == self.MJD0
        assert np.allclose(dut1, self.VALUES["dut1"])
        assert np.allclose(x, self.VALUES["x"])
        assert np.allclose(y, self.VALUES["y"])

    def test_stored_values_are_not_rotated(self, mod, tmp_path):
        q = mod.QPoint()
        q.load_bulletin_a(self._write(tmp_path, self.COLUMNS))
        got = q.get_bulletin_a(self.MJD0 + 10)
        assert np.allclose(
            got, [self.VALUES["dut1"], self.VALUES["x"], self.VALUES["y"]]
        )

    def test_a_reordered_file(self, mod, tmp_path):
        """What the columns argument is for."""
        order = ["x", "mjd", "y", "dut1"]
        q = mod.QPoint()
        q.load_bulletin_a(self._write(tmp_path, order), columns=order)
        got = q.get_bulletin_a(self.MJD0 + 10)
        assert np.allclose(
            got, [self.VALUES["dut1"], self.VALUES["x"], self.VALUES["y"]]
        )

    def test_missing_columns_raise(self, mod, tmp_path):
        q = mod.QPoint()
        path = self._write(tmp_path, self.COLUMNS)
        with pytest.raises(KeyError):
            q.load_bulletin_a(path, columns=["mjd", "dut1", "x"])


# Derived from the package rather than written out, so a rate added later
# is covered without anyone remembering to add it here.
FORWARD_RATES = sorted(
    k for k in qpoint.QPoint().get("rates") if not k.endswith("_inv")
)


class TestResetRates:
    """
    reset_rates has to put back every correction, not most of them.

    qp_reset_rates enumerates them one at a time, so it can be incomplete, and
    rate_defl was: reset_rates carried the previous chunk's sun position into
    the next one, where the documented use is to call it at the start of each
    chunk. The rates are read off the package, so this keeps holding as new
    ones appear.
    """

    N = 5
    AZ = np.linspace(10, 340, N)
    EL = np.linspace(30, 80, N)
    LONS = np.full(N, LON)
    LATS = np.full(N, LAT)
    T1 = np.full(N, CTIME)
    # half a year on: every slowly-varying correction has moved a long way
    T2 = T1 + 180 * 86400.0

    def radec(self, ctime, warm=None, **kwargs):
        q = qpoint.QPoint(mean_aber=True, **kwargs)
        args = (0.0, 0.0, 0.0, self.AZ, self.EL, None, None, self.LONS, self.LATS)
        if warm is not None:
            q.azel2radec(*args, warm)
            q.reset_rates()
        return np.asarray(q.azel2radec(*args, ctime))

    @pytest.mark.parametrize("rate", FORWARD_RATES)
    def test_reset_leaves_it_as_good_as_new(self, rate):
        """
        With one rate pinned to 'once' the correction is computed at the
        first sample and frozen, so a stale cache is visible: after
        reset_rates the answer has to match a freshly built QPoint.
        """
        fresh = self.radec(self.T2, **{rate: "once"})
        reused = self.radec(self.T2, warm=self.T1, **{rate: "once"})
        assert np.array_equal(fresh, reused), rate

    @pytest.mark.parametrize("rate", ["rate_npb", "rate_defl"])
    def test_and_the_comparison_can_fail(self, rate):
        """
        The teeth: without the reset these two really do carry the stale
        correction forward, so the assertion above is not vacuous.
        """
        q = qpoint.QPoint(mean_aber=True, **{rate: "once"})
        args = (0.0, 0.0, 0.0, self.AZ, self.EL, None, None, self.LONS, self.LATS)
        q.azel2radec(*args, self.T1)
        stale = np.asarray(q.azel2radec(*args, self.T2))
        assert not np.array_equal(self.radec(self.T2, **{rate: "once"}), stale)

    def test_reset_inv_rates_runs(self, qp):
        qp.reset_inv_rates()


# ---------------------------------------------------------------------------
# det_offset
# ---------------------------------------------------------------------------


class TestDetOffset:
    def test_zero_offset_shape(self, qp):
        q = qp.det_offset(0.0, 0.0, 0.0)
        assert q.shape == (4,)

    def test_zero_offset_is_unit(self, qp):
        q = qp.det_offset(0.0, 0.0, 0.0)
        assert np.isclose(np.linalg.norm(q), 1.0)

    def test_zero_offset_is_identity(self, qp):
        q = qp.det_offset(0.0, 0.0, 0.0)
        assert np.allclose(np.abs(q), [1, 0, 0, 0])

    def test_nonzero_offset_unit(self, qp):
        q = qp.det_offset(1.0, -1.0, 22.5)
        assert np.isclose(np.linalg.norm(q), 1.0)

    def test_vector_offset_shape(self, qp):
        q = qp.det_offset([0.0, 1.0, -1.0], [0.0, 0.0, 1.0], [0.0, 0.0, 0.0])
        assert q.shape == (3, 4)

    def test_vector_offset_unit(self, qp):
        q = qp.det_offset([0.0, 1.0, -1.0], [0.0, 0.0, 1.0], [0.0, 0.0, 0.0])
        norms = np.linalg.norm(q, axis=1)
        assert np.allclose(norms, 1.0)


# ---------------------------------------------------------------------------
# hwp_quat
# ---------------------------------------------------------------------------


class TestHwpQuat:
    def test_scalar_shape(self, qp):
        q = qp.hwp_quat(0.0)
        assert q.shape == (4,)

    def test_scalar_unit(self, qp):
        q = qp.hwp_quat(0.0)
        assert np.isclose(np.linalg.norm(q), 1.0)

    def test_zero_angle_identity(self, qp):
        q = qp.hwp_quat(0.0)
        assert np.allclose(np.abs(q), [1, 0, 0, 0])

    def test_vector_shape(self, qp):
        theta = np.linspace(0, 180, 10)
        q = qp.hwp_quat(theta)
        assert q.shape == (10, 4)

    def test_vector_unit(self, qp):
        theta = np.linspace(0, 180, 10)
        q = qp.hwp_quat(theta)
        norms = np.linalg.norm(q, axis=1)
        assert np.allclose(norms, 1.0)


# ---------------------------------------------------------------------------
# azel2bore
# ---------------------------------------------------------------------------


class TestAzel2Bore:
    def test_shape(self, qp):
        q = qp.azel2bore(AZ, EL, None, None, LON, LAT, CTIMES)
        assert q.shape == (N, 4)

    def test_unit_quaternions(self, qp):
        q = qp.azel2bore(AZ, EL, None, None, LON, LAT, CTIMES)
        norms = np.linalg.norm(q, axis=1)
        assert np.allclose(norms, 1.0, atol=1e-10)

    def test_with_pitch_roll(self, qp):
        pitch = np.zeros(N)
        roll = np.zeros(N)
        q = qp.azel2bore(AZ, EL, pitch, roll, LON, LAT, CTIMES)
        assert q.shape == (N, 4)
        norms = np.linalg.norm(q, axis=1)
        assert np.allclose(norms, 1.0, atol=1e-10)

    def test_pitch_roll_none_matches_zeros(self, qp):
        qp.reset_rates()
        q1 = qp.azel2bore(AZ, EL, None, None, LON, LAT, CTIMES)
        qp.reset_rates()
        q2 = qp.azel2bore(AZ, EL, np.zeros(N), np.zeros(N), LON, LAT, CTIMES)
        assert np.allclose(q1, q2, atol=1e-10)

    def test_returns_unit_quaternions_repeated_call(self, qp):
        # Ensure repeated calls are stable and always produce unit quaternions
        for _ in range(3):
            qp.reset_rates()
            q = qp.azel2bore(AZ, EL, None, None, LON, LAT, CTIMES)
            norms = np.linalg.norm(q, axis=1)
            assert np.allclose(norms, 1.0, atol=1e-10)


# ---------------------------------------------------------------------------
# azelpsi2bore
# ---------------------------------------------------------------------------


class TestAzelpsi2Bore:
    def test_shape(self, qp):
        psi = np.zeros(N)
        q = qp.azelpsi2bore(AZ, EL, psi, None, None, LON, LAT, CTIMES)
        assert q.shape == (N, 4)

    def test_unit_quaternions(self, qp):
        psi = np.zeros(N)
        q = qp.azelpsi2bore(AZ, EL, psi, None, None, LON, LAT, CTIMES)
        norms = np.linalg.norm(q, axis=1)
        assert np.allclose(norms, 1.0, atol=1e-10)

    def test_zero_psi_matches_azel2bore(self, qp):
        psi = np.zeros(N)
        qp.reset_rates()
        q1 = qp.azel2bore(AZ, EL, None, None, LON, LAT, CTIMES)
        qp.reset_rates()
        q2 = qp.azelpsi2bore(AZ, EL, psi, None, None, LON, LAT, CTIMES)
        assert np.allclose(q1, q2, atol=1e-10)


# ---------------------------------------------------------------------------
# bore2radec
# ---------------------------------------------------------------------------


class TestBore2Radec:
    def _make_bore(self, qp):
        return qp.azel2bore(AZ, EL, None, None, LON, LAT, CTIMES)

    def _make_qoff(self, qp):
        return qp.det_offset(0.0, 0.0, 0.0)

    def test_sindec_with_return_pa_raises(self, qp):
        """The C has no entry point taking both, so the combination is refused."""
        q_bore = self._make_bore(qp)
        q_off = self._make_qoff(qp)
        with pytest.raises(ValueError):
            qp.bore2radec(q_off, CTIMES, q_bore, sindec=True, return_pa=True)

    def test_true_scalars_give_scalars(self, qp):
        """Scalar in, scalar out -- nothing supplied an axis to keep."""
        q_bore = qp.azel2bore(AZ[0], EL[0], None, None, LON, LAT, CTIMES[0])
        q_off = qp.det_offset(0.0, 0.0, 0.0)
        ra, dec, sin2psi, cos2psi = qp.bore2radec(q_off, CTIMES[0], q_bore[0])
        assert np.isscalar(ra) or np.asarray(ra).ndim == 0

    def test_shape_default(self, qp):
        q_bore = self._make_bore(qp)
        q_off = self._make_qoff(qp)
        ra, dec, sin2psi, cos2psi = qp.bore2radec(q_off, CTIMES, q_bore)
        assert ra.shape == (N,)
        assert dec.shape == (N,)
        assert sin2psi.shape == (N,)
        assert cos2psi.shape == (N,)

    def test_ra_range(self, qp):
        q_bore = self._make_bore(qp)
        q_off = self._make_qoff(qp)
        ra, dec, _, _ = qp.bore2radec(q_off, CTIMES, q_bore)
        # RA is returned in degrees; may be in (-180, 180] or (0, 360)
        assert np.all(ra >= -180) and np.all(ra < 360)

    def test_dec_range(self, qp):
        q_bore = self._make_bore(qp)
        q_off = self._make_qoff(qp)
        _, dec, _, _ = qp.bore2radec(q_off, CTIMES, q_bore)
        assert np.all(dec >= -90) and np.all(dec <= 90)

    def test_sincos_bounded(self, qp):
        q_bore = self._make_bore(qp)
        q_off = self._make_qoff(qp)
        _, _, sin2psi, cos2psi = qp.bore2radec(q_off, CTIMES, q_bore)
        assert np.all(np.abs(sin2psi) <= 1.0 + 1e-10)
        assert np.all(np.abs(cos2psi) <= 1.0 + 1e-10)

    def test_return_pa_shape(self, qp):
        q_bore = self._make_bore(qp)
        q_off = self._make_qoff(qp)
        result = qp.bore2radec(q_off, CTIMES, q_bore, return_pa=True)
        assert len(result) == 3
        ra, dec, pa = result
        assert pa.shape == (N,)

    def test_sindec_range(self, qp):
        q_bore = self._make_bore(qp)
        q_off = self._make_qoff(qp)
        _, sindec, _, _ = qp.bore2radec(q_off, CTIMES, q_bore, sindec=True)
        assert np.all(np.abs(sindec) <= 1.0 + 1e-10)

    def test_pa_consistent_with_sincos(self, qp):
        q_bore = self._make_bore(qp)
        q_off = self._make_qoff(qp)
        qp.reset_rates()
        ra1, dec1, sin2psi, cos2psi = qp.bore2radec(q_off, CTIMES, q_bore)
        qp.reset_rates()
        ra2, dec2, pa = qp.bore2radec(q_off, CTIMES, q_bore, return_pa=True)
        assert np.allclose(ra1, ra2, atol=1e-10)
        assert np.allclose(dec1, dec2, atol=1e-10)
        pa_rad = np.deg2rad(2.0 * pa)
        assert np.allclose(np.sin(pa_rad), sin2psi, atol=1e-8)
        assert np.allclose(np.cos(pa_rad), cos2psi, atol=1e-8)

    def test_with_hwp(self, qp):
        q_bore = self._make_bore(qp)
        q_off = self._make_qoff(qp)
        q_hwp = qp.hwp_quat(22.5 * np.ones(N))
        ra, dec, sin2psi, cos2psi = qp.bore2radec(q_off, CTIMES, q_bore, q_hwp=q_hwp)
        assert ra.shape == (N,)

    def test_with_hwp_return_pa(self, qp):
        q_bore = self._make_bore(qp)
        q_off = self._make_qoff(qp)
        q_hwp = qp.hwp_quat(22.5 * np.ones(N))
        ra, dec, pa = qp.bore2radec(q_off, CTIMES, q_bore, q_hwp=q_hwp, return_pa=True)
        assert pa.shape == (N,)

    def test_ctime_none_with_mean_aber(self):
        qp = qpoint.QPoint(mean_aber=True, accuracy="low")
        q_bore = qp.azel2bore(AZ, EL, None, None, LON, LAT, CTIMES)
        q_off = qp.det_offset(0.0, 0.0, 0.0)
        ra, dec, sin2psi, cos2psi = qp.bore2radec(q_off, None, q_bore)
        assert ra.shape == (N,)

    def test_ctime_none_without_mean_aber_raises(self):
        qp = qpoint.QPoint(mean_aber=False, accuracy="low")
        q_bore = qp.azel2bore(AZ, EL, None, None, LON, LAT, CTIMES)
        q_off = qp.det_offset(0.0, 0.0, 0.0)
        with pytest.raises(ValueError):
            qp.bore2radec(q_off, None, q_bore)

    def test_single_sample_scalar_output(self, qp):
        q_bore = qp.azel2bore(AZ[:1], EL[:1], None, None, LON, LAT, CTIMES[:1])
        q_off = qp.det_offset(0.0, 0.0, 0.0)
        ra, dec, sin2psi, cos2psi = qp.bore2radec(q_off, CTIMES[:1], q_bore)
        assert np.isscalar(ra) or ra.ndim == 0


# ---------------------------------------------------------------------------
# azel2radec
# ---------------------------------------------------------------------------


class TestAzel2Radec:
    def test_shape(self, qp):
        ra, dec, sin2psi, cos2psi = qp.azel2radec(
            0, 0, 0, AZ, EL, None, None, LON, LAT, CTIMES
        )
        assert ra.shape == (N,)
        assert dec.shape == (N,)
        assert sin2psi.shape == (N,)
        assert cos2psi.shape == (N,)

    def test_consistent_with_bore2radec(self, qp):
        """Zero-offset azel2radec must agree with azel2bore + bore2radec."""
        qp.reset_rates()
        ra1, dec1, s1, c1 = qp.azel2radec(0, 0, 0, AZ, EL, None, None, LON, LAT, CTIMES)
        qp.reset_rates()
        q_bore = qp.azel2bore(AZ, EL, None, None, LON, LAT, CTIMES)
        q_off = qp.det_offset(0.0, 0.0, 0.0)
        ra2, dec2, s2, c2 = qp.bore2radec(q_off, CTIMES, q_bore)
        assert np.allclose(ra1, ra2, atol=1e-10)
        assert np.allclose(dec1, dec2, atol=1e-10)

    def test_return_pa(self, qp):
        result = qp.azel2radec(
            0, 0, 0, AZ, EL, None, None, LON, LAT, CTIMES, return_pa=True
        )
        assert len(result) == 3

    def test_sindec(self, qp):
        ra, sindec, sin2psi, cos2psi = qp.azel2radec(
            0, 0, 0, AZ, EL, None, None, LON, LAT, CTIMES, sindec=True
        )
        assert np.all(np.abs(sindec) <= 1.0 + 1e-10)

    def test_with_hwp(self, qp):
        hwp = 22.5 * np.ones(N)
        result = qp.azel2radec(0, 0, 0, AZ, EL, None, None, LON, LAT, CTIMES, hwp=hwp)
        assert len(result) == 4

    def test_with_hwp_return_pa(self, qp):
        hwp = 22.5 * np.ones(N)
        result = qp.azel2radec(
            0, 0, 0, AZ, EL, None, None, LON, LAT, CTIMES, hwp=hwp, return_pa=True
        )
        assert len(result) == 3


# ---------------------------------------------------------------------------
# azelpsi2radec
# ---------------------------------------------------------------------------


class TestAzelpsi2Radec:
    def test_shape(self, qp):
        psi = np.zeros(N)
        ra, dec, sin2psi, cos2psi = qp.azelpsi2radec(
            0, 0, 0, AZ, EL, psi, None, None, LON, LAT, CTIMES
        )
        assert ra.shape == (N,)

    def test_zero_psi_matches_azel2radec(self, qp):
        psi = np.zeros(N)
        qp.reset_rates()
        ra1, dec1, s1, c1 = qp.azel2radec(0, 0, 0, AZ, EL, None, None, LON, LAT, CTIMES)
        qp.reset_rates()
        ra2, dec2, s2, c2 = qp.azelpsi2radec(
            0, 0, 0, AZ, EL, psi, None, None, LON, LAT, CTIMES
        )
        assert np.allclose(ra1, ra2, atol=1e-10)
        assert np.allclose(dec1, dec2, atol=1e-10)


# ---------------------------------------------------------------------------
# radec2azel round-trip
# ---------------------------------------------------------------------------


class TestRadec2Azel:
    def test_roundtrip(self, qp):
        """azel2radec(return_pa=True) -> radec2azel should recover az/el."""
        qp.reset_rates()
        ra, dec, pa = qp.azel2radec(
            0, 0, 0, AZ, EL, None, None, LON, LAT, CTIMES, return_pa=True
        )
        qp.reset_inv_rates()
        az2, el2, hpa = qp.radec2azel(ra, dec, pa, LON, LAT, CTIMES)
        assert np.allclose(az2 % 360, AZ % 360, atol=1e-4)
        assert np.allclose(el2, EL, atol=1e-4)

    def test_output_shapes(self, qp):
        qp.reset_rates()
        ra, dec, pa = qp.azel2radec(
            0, 0, 0, AZ, EL, None, None, LON, LAT, CTIMES, return_pa=True
        )
        az, el, hpa = qp.radec2azel(ra, dec, pa, LON, LAT, CTIMES)
        assert az.shape == (N,)
        assert el.shape == (N,)
        assert hpa.shape == (N,)


# ---------------------------------------------------------------------------
# bore2azel round-trip
# ---------------------------------------------------------------------------


class TestBore2Azel:
    def test_roundtrip(self, qp):
        """azel2bore -> bore2azel should recover az/el."""
        qp.reset_rates()
        q_bore = qp.azel2bore(AZ, EL, None, None, LON, LAT, CTIMES)
        qp.reset_inv_rates()
        az2, el2, pa2 = qp.bore2azel(q_bore, LON, LAT, CTIMES)
        assert np.allclose(az2 % 360, AZ % 360, atol=1e-4)
        assert np.allclose(el2, EL, atol=1e-4)

    def test_output_shapes(self, qp):
        q_bore = qp.azel2bore(AZ, EL, None, None, LON, LAT, CTIMES)
        az, el, pa = qp.bore2azel(q_bore, LON, LAT, CTIMES)
        assert az.shape == (N,)
        assert el.shape == (N,)
        assert pa.shape == (N,)


# ---------------------------------------------------------------------------
# radecpa2quat / quat2radecpa
# ---------------------------------------------------------------------------


class TestRadecpaQuat:
    def test_radecpa2quat_shape(self, qp):
        ra = np.array([10.0, 45.0, 180.0])
        dec = np.array([20.0, -30.0, 0.0])
        pa = np.array([0.0, 45.0, 90.0])
        q = qp.radecpa2quat(ra, dec, pa)
        assert q.shape == (3, 4)

    def test_radecpa2quat_unit(self, qp):
        ra = np.linspace(0, 360, 8, endpoint=False)
        dec = np.linspace(-80, 80, 8)
        pa = np.zeros(8)
        q = qp.radecpa2quat(ra, dec, pa)
        norms = np.linalg.norm(q, axis=1)
        assert np.allclose(norms, 1.0, atol=1e-10)

    def test_quat2radecpa_shape(self, qp):
        ra = np.array([10.0, 45.0, 180.0])
        dec = np.array([20.0, -30.0, 0.0])
        pa = np.array([0.0, 45.0, 90.0])
        q = qp.radecpa2quat(ra, dec, pa)
        ra2, dec2, pa2 = qp.quat2radecpa(q)
        assert ra2.shape == (3,)
        assert dec2.shape == (3,)
        assert pa2.shape == (3,)

    def test_roundtrip(self, qp):
        ra = np.array([10.0, 45.0, 90.0, 270.0])
        dec = np.array([20.0, -30.0, 0.0, 60.0])
        pa = np.array([0.0, 45.0, 90.0, 180.0])
        q = qp.radecpa2quat(ra, dec, pa)
        ra2, dec2, pa2 = qp.quat2radecpa(q)
        assert np.allclose(ra2 % 360, ra % 360, atol=1e-10)
        assert np.allclose(dec2, dec, atol=1e-10)
        assert np.allclose(pa2 % 360, pa % 360, atol=1e-10)

    def test_single_quat_scalar_output(self, qp):
        q = np.array([1.0, 0.0, 0.0, 0.0])
        ra, dec, pa = qp.quat2radecpa(q)
        assert np.isscalar(ra) or ra.ndim == 0


# ---------------------------------------------------------------------------
# gmst / lmst
# ---------------------------------------------------------------------------


class TestGmst:
    def test_scalar_output(self, qp):
        gmst = qp.gmst(CTIME)
        assert np.isscalar(gmst) or gmst.ndim == 0

    def test_scalar_range(self, qp):
        gmst = qp.gmst(CTIME)
        assert 0.0 <= float(gmst) < 24.0

    def test_vector_shape(self, qp):
        ctimes = CTIME + np.arange(10) * 3600.0
        gmst = qp.gmst(ctimes)
        assert gmst.shape == (10,)

    def test_vector_range(self, qp):
        ctimes = CTIME + np.arange(100) * 3600.0
        gmst = qp.gmst(ctimes)
        assert np.all(gmst >= 0) and np.all(gmst < 24)


class TestLmst:
    def test_scalar_output(self, qp):
        lmst = qp.lmst(CTIME, LON)
        assert np.isscalar(lmst) or lmst.ndim == 0

    def test_scalar_range(self, qp):
        lmst = qp.lmst(CTIME, LON)
        assert 0.0 <= float(lmst) < 24.0

    def test_vector_shape(self, qp):
        ctimes = CTIME + np.arange(10) * 3600.0
        lons = LON * np.ones(10)
        lmst = qp.lmst(ctimes, lons)
        assert lmst.shape == (10,)

    def test_lmst_vs_gmst(self, qp):
        """LMST and GMST should differ by the longitude / 15."""
        gmst = qp.gmst(CTIME)
        lmst = qp.lmst(CTIME, LON)
        diff = (lmst - gmst - LON / 15.0) % 24
        assert np.isclose(diff, 0.0, atol=1e-6) or np.isclose(diff, 24.0, atol=1e-6)


# ---------------------------------------------------------------------------
# Galactic rotation
# ---------------------------------------------------------------------------


class TestGalactic:
    def test_radec2gal_gal2radec_roundtrip(self, qp):
        ra_orig = np.array([10.0, 45.0, 180.0, 270.0])
        dec_orig = np.array([20.0, -30.0, 0.0, 60.0])
        l, b, pa_g = qp.radec2gal(ra_orig.copy(), dec_orig.copy())
        ra2, dec2, _ = qp.gal2radec(l, b, pa_g)
        assert np.allclose(ra2 % 360, ra_orig % 360, atol=1e-10)
        assert np.allclose(dec2, dec_orig, atol=1e-10)

    def test_rotate_coord_cg_gc_roundtrip(self, qp):
        ra_orig = np.array([10.0, 90.0])
        dec_orig = np.array([20.0, 0.0])
        pa_orig = np.array([0.0, 45.0])
        l, b, pa_g = qp.rotate_coord(
            ra_orig.copy(), dec_orig.copy(), pa_orig.copy(), coord=["C", "G"]
        )
        ra2, dec2, pa2 = qp.rotate_coord(l, b, pa_g, coord=["G", "C"])
        assert np.allclose(ra2 % 360, ra_orig % 360, atol=1e-10)
        assert np.allclose(dec2, dec_orig, atol=1e-10)

    def test_rotate_coord_sincos_roundtrip(self, qp):
        ra = np.array([10.0, 90.0])
        dec = np.array([20.0, 0.0])
        s2 = np.array([0.5, -0.5])
        c2 = np.sqrt(1 - s2**2)
        ra_c, dec_c, s2_c, c2_c = ra.copy(), dec.copy(), s2.copy(), c2.copy()
        l, b, sl, cl = qp.rotate_coord(
            ra_c, dec_c, sin2psi=s2_c, cos2psi=c2_c, coord=["C", "G"]
        )
        ra2, dec2, s2r, c2r = qp.rotate_coord(
            l, b, sin2psi=sl, cos2psi=cl, coord=["G", "C"]
        )
        assert np.allclose(ra2 % 360, ra % 360, atol=1e-10)
        assert np.allclose(dec2, dec, atol=1e-10)

    def test_rotate_quat_cg_gc_roundtrip(self, qp):
        ra = np.array([10.0, 90.0])
        dec = np.array([20.0, 0.0])
        pa = np.zeros(2)
        q_orig = qp.radecpa2quat(ra, dec, pa)
        q = q_orig.copy()
        qp.rotate_quat(q, coord=["C", "G"])
        qp.rotate_quat(q, coord=["G", "C"])
        assert np.allclose(q, q_orig, atol=1e-10)

    def test_rotate_quat_unsupported_raises(self, qp):
        q = np.eye(1, 4)
        with pytest.raises(ValueError):
            qp.rotate_quat(q, coord=["C", "X"])


# ---------------------------------------------------------------------------
# radec2pix / quat2pix
# ---------------------------------------------------------------------------


class TestRadec2Pix:
    def test_shape(self, qp):
        nside = 64
        ra = np.array([0.0, 90.0, 180.0, 270.0])
        dec = np.array([0.0, 30.0, -30.0, 60.0])
        pix = qp.radec2pix(ra, dec, nside=nside)
        assert pix.shape == (4,)

    def test_valid_range(self, qp):
        nside = 64
        ra = np.array([0.0, 90.0, 180.0, 270.0])
        dec = np.array([0.0, 30.0, -30.0, 60.0])
        pix = qp.radec2pix(ra, dec, nside=nside)
        npix = 12 * nside * nside
        assert np.all(pix >= 0) and np.all(pix < npix)

    def test_scalar_output(self, qp):
        pix = qp.radec2pix(0.0, 0.0, nside=64)
        assert np.isscalar(pix) or pix.ndim == 0


class TestQuat2Pix:
    def test_shape(self, qp):
        nside = 64
        ra = np.array([0.0, 90.0, 180.0])
        dec = np.array([0.0, 30.0, -30.0])
        pa = np.zeros(3)
        q = qp.radecpa2quat(ra, dec, pa)
        pix, sin2psi, cos2psi = qp.quat2pix(q, nside=nside)
        assert pix.shape == (3,)
        assert sin2psi.shape == (3,)
        assert cos2psi.shape == (3,)

    def test_consistent_with_radec2pix(self, qp):
        nside = 64
        ra = np.array([0.0, 90.0, 180.0, 270.0])
        dec = np.array([10.0, -10.0, 30.0, -30.0])
        pa = np.zeros(4)
        q = qp.radecpa2quat(ra, dec, pa)
        pix_q, _, _ = qp.quat2pix(q, nside=nside)
        pix_rd = qp.radec2pix(ra, dec, nside=nside)
        assert np.all(pix_q == pix_rd)


class TestQuat2PixPa:
    """quat2pixpa previously dropped nside from its ctypes call."""

    def _quats(self, qp):
        self.ra = np.array([0.0, 90.0, 180.0, 270.0])
        self.dec = np.array([10.0, -10.0, 30.0, -30.0])
        return qp.radecpa2quat(self.ra, self.dec, np.zeros(4))

    def test_shape(self, qp):
        pix, pa = qp.quat2pixpa(self._quats(qp), nside=64)
        assert pix.shape == (4,)
        assert pa.shape == (4,)

    def test_consistent_with_radec2pix(self, qp):
        q = self._quats(qp)
        pix, _ = qp.quat2pixpa(q, nside=64)
        assert np.all(pix == qp.radec2pix(self.ra, self.dec, nside=64))

    def test_nside_is_honored(self, qp):
        q = self._quats(qp)
        assert not np.array_equal(
            qp.quat2pixpa(q, nside=64)[0], qp.quat2pixpa(q, nside=256)[0]
        )

    def test_pa_matches_quat2radecpa(self, qp):
        q = self._quats(qp)
        _, pa = qp.quat2pixpa(q, nside=64)
        _, _, pa_ref = qp.quat2radecpa(q)
        assert np.allclose(pa, pa_ref)


# ---------------------------------------------------------------------------
# bore2pix
# ---------------------------------------------------------------------------


class TestBore2Pix:
    def _make_bore(self, qp):
        return qp.azel2bore(AZ, EL, None, None, LON, LAT, CTIMES)

    def _make_qoff(self, qp):
        return qp.det_offset(0.0, 0.0, 0.0)

    def test_shape_pol_true(self, qp):
        q_bore = self._make_bore(qp)
        q_off = self._make_qoff(qp)
        pix, sin2psi, cos2psi = qp.bore2pix(q_off, CTIMES, q_bore, nside=64)
        assert pix.shape == (N,)
        assert sin2psi.shape == (N,)
        assert cos2psi.shape == (N,)

    def test_valid_range(self, qp):
        nside = 64
        q_bore = self._make_bore(qp)
        q_off = self._make_qoff(qp)
        pix, _, _ = qp.bore2pix(q_off, CTIMES, q_bore, nside=nside)
        npix = 12 * nside * nside
        assert np.all(pix >= 0) and np.all(pix < npix)

    def test_pol_false(self, qp):
        q_bore = self._make_bore(qp)
        q_off = self._make_qoff(qp)
        pix = qp.bore2pix(q_off, CTIMES, q_bore, nside=64, pol=False)
        assert pix.shape == (N,)

    def test_return_pa(self, qp):
        q_bore = self._make_bore(qp)
        q_off = self._make_qoff(qp)
        pix, pa = qp.bore2pix(q_off, CTIMES, q_bore, nside=64, return_pa=True)
        assert pix.shape == (N,)
        assert pa.shape == (N,)

    def test_consistent_with_radec2pix(self, qp):
        nside = 64
        q_bore = self._make_bore(qp)
        q_off = self._make_qoff(qp)
        qp.reset_rates()
        pix1 = qp.bore2pix(q_off, CTIMES, q_bore, nside=nside, pol=False)
        qp.reset_rates()
        ra, dec, _, _ = qp.bore2radec(q_off, CTIMES, q_bore)
        pix2 = qp.radec2pix(ra, dec, nside=nside)
        assert np.all(pix1 == pix2)

    def test_with_hwp(self, qp):
        q_bore = self._make_bore(qp)
        q_off = self._make_qoff(qp)
        q_hwp = qp.hwp_quat(22.5 * np.ones(N))
        pix, sin2psi, cos2psi = qp.bore2pix(
            q_off, CTIMES, q_bore, q_hwp=q_hwp, nside=64
        )
        assert pix.shape == (N,)


# ---------------------------------------------------------------------------
# dipole / bore2dipole
# ---------------------------------------------------------------------------


class TestDipole:
    def test_shape(self, qp):
        ra = np.zeros(N)
        dec = np.zeros(N)
        d = qp.dipole(CTIMES, ra, dec)
        assert d.shape == (N,)

    def test_physical_range(self, qp):
        ra = np.zeros(N)
        dec = np.zeros(N)
        d = qp.dipole(CTIMES, ra, dec)
        # CMB dipole amplitude ~3.36 mK
        assert np.all(np.abs(d) < 0.01)

    def test_scalar_output(self, qp):
        d = qp.dipole(CTIME, 0.0, 0.0)
        assert np.isscalar(d) or d.ndim == 0


class TestBore2Dipole:
    def test_shape(self, qp):
        q_bore = qp.azel2bore(AZ, EL, None, None, LON, LAT, CTIMES)
        q_off = qp.det_offset(0.0, 0.0, 0.0)
        d = qp.bore2dipole(q_off, CTIMES, q_bore)
        assert d.shape == (N,)

    def test_physical_range(self, qp):
        q_bore = qp.azel2bore(AZ, EL, None, None, LON, LAT, CTIMES)
        q_off = qp.det_offset(0.0, 0.0, 0.0)
        d = qp.bore2dipole(q_off, CTIMES, q_bore)
        assert np.all(np.abs(d) < 0.01)


# ---------------------------------------------------------------------------
# omega2azelpsi
# ---------------------------------------------------------------------------


class TestOmega2Azelpsi:
    def test_output_shapes(self, qp):
        n = 50
        omega_x = np.zeros(n)
        omega_y = np.zeros(n)
        omega_z = np.zeros(n)
        az, el, psi = qp.omega2azelpsi(45.0, 30.0, 0.0, omega_x, omega_y, omega_z, 0.01)
        assert az.shape == (n,)
        assert el.shape == (n,)
        assert psi.shape == (n,)

    def test_zero_omega_constant_position(self, qp):
        n = 50
        omega_x = np.zeros(n)
        omega_y = np.zeros(n)
        omega_z = np.zeros(n)
        az, el, psi = qp.omega2azelpsi(45.0, 30.0, 0.0, omega_x, omega_y, omega_z, 0.01)
        assert np.allclose(az, 45.0, atol=1e-8)
        assert np.allclose(el, 30.0, atol=1e-8)
        assert np.allclose(psi, 0.0, atol=1e-8)


# ---------------------------------------------------------------------------
# bore_offset
# ---------------------------------------------------------------------------


class TestBoreOffset:
    def test_output_shape(self, qp):
        q_bore = qp.azel2bore(AZ, EL, None, None, LON, LAT, CTIMES)
        ang1 = np.zeros(N)
        q_out = qp.bore_offset(q_bore, ang1=ang1)
        assert q_out.shape == (N, 4)

    def test_output_unit(self, qp):
        q_bore = qp.azel2bore(AZ, EL, None, None, LON, LAT, CTIMES)
        ang1 = np.zeros(N)
        q_out = qp.bore_offset(q_bore, ang1=ang1)
        norms = np.linalg.norm(q_out, axis=1)
        assert np.allclose(norms, 1.0, atol=1e-10)

    def test_no_angle_raises(self, qp):
        q_bore = qp.azel2bore(AZ, EL, None, None, LON, LAT, CTIMES)
        with pytest.raises(ValueError):
            qp.bore_offset(q_bore)


# ---------------------------------------------------------------------------
# refraction method
# ---------------------------------------------------------------------------


class TestRefractionMethod:
    def test_set_delta_scalar(self, qp):
        qp.refraction(0.05)
        assert qp.get("ref_delta") == pytest.approx(0.05)

    def test_delta_kwarg(self, qp):
        delta = qp.refraction(delta=0.01)
        assert delta == pytest.approx(0.01)
        assert qp.get("ref_delta") == pytest.approx(0.01)

    def test_get_stored_delta(self, qp):
        qp.set(ref_delta=0.02)
        delta = qp.refraction()
        assert delta == pytest.approx(0.02)


class TestRefractionAtTheZenith:
    """
    The refraction correction reads sin(el) off the quaternion as
    q0^2 - q1^2 - q2^2 + q3^2, and q is unit only to within rounding, so at
    el = 90 that can land a hair above 1. asin of it is NaN, which went into
    ref_delta and from there into the boresight: 11664 of 172824 samples
    over a day of azimuth sweeps, before the clamp.

    Only the forward transform is pinned. The inverse reads the elevation
    back after undoing the aberration, never lands on exactly 1, and did not
    reproduce it in bore2azel or radec2azel. poly_asin clamps on its own, so
    fast_math never failed; it is checked for staying finite, not as a
    regression.
    """

    AZ = np.linspace(0, 360, 721)
    CT = np.full(AZ.size, CTIME)

    @pytest.mark.parametrize("fast_math", [False, True])
    def test_forward(self, qp, fast_math):
        q = qp.azel2bore(
            self.AZ,
            90.0,
            None,
            None,
            LON,
            LAT,
            self.CT,
            rate_ref="always",
            fast_math=fast_math,
        )
        assert np.isfinite(q).all()


# ---------------------------------------------------------------------------
# get_interp_val
# ---------------------------------------------------------------------------


class TestRotateMap:
    """rotate_map previously raised on every call, so none of this was covered."""

    NSIDE = 8
    NPIX = 12 * 8 * 8

    def test_runs_and_keeps_shape(self, qp):
        m = np.zeros((3, self.NPIX))
        m[0] = 1.0
        out = qp.rotate_map(m, coord=("C", "G"))
        assert out.shape == (3, self.NPIX)

    @pytest.mark.parametrize("coord", [("C", "G"), ("G", "C")])
    def test_constant_temperature_is_preserved(self, qp, coord):
        """Resampling a constant map must give the same constant back."""
        m = np.zeros((3, self.NPIX))
        m[0] = 2.5
        out = qp.rotate_map(m, coord=coord)
        assert np.allclose(out[0], 2.5)
        assert np.allclose(out[1], 0.0)
        assert np.allclose(out[2], 0.0)

    def test_polarized_intensity_is_preserved(self, qp):
        """
        Rotation mixes Q into U, but the polarized intensity at a pixel is
        invariant.
        """
        m = np.zeros((3, self.NPIX))
        m[1] = 0.6
        m[2] = 0.8
        out = qp.rotate_map(m, coord=("C", "G"), interp_pix=False)
        assert np.allclose(np.hypot(out[1], out[2]), 1.0)

    def test_rotation_actually_mixes_q_and_u(self, qp):
        m = np.zeros((3, self.NPIX))
        m[1] = 1.0
        out = qp.rotate_map(m, coord=("C", "G"), interp_pix=False)
        assert not np.allclose(out[2], 0.0)

    @pytest.mark.parametrize("nrow", [1, 2, 4])
    def test_wrong_row_count_raises(self, qp, nrow):
        """Fewer than three rows used to read off the end and segfault."""
        with pytest.raises(ValueError, match="3 rows"):
            qp.rotate_map(np.ones((nrow, self.NPIX)), coord=("C", "G"))


class TestGetInterpVal:
    def test_constant_map(self, qp):
        nside = 8
        npix = 12 * nside * nside
        m = np.ones(npix)
        ra = np.array([0.0, 90.0, 180.0])
        dec = np.array([0.0, 30.0, -30.0])
        val = qp.get_interp_val(m, ra, dec)
        assert np.allclose(val, 1.0, atol=1e-10)

    def test_scalar_output(self, qp):
        nside = 8
        npix = 12 * nside * nside
        m = np.ones(npix)
        val = qp.get_interp_val(m, 0.0, 0.0)
        assert np.isscalar(val) or val.ndim == 0
        assert np.isclose(float(val), 1.0, atol=1e-10)

    def test_single_map_squeezes(self, qp):
        npix = 12 * 8 * 8
        ra = np.array([0.0, 90.0, 180.0])
        dec = np.array([0.0, 30.0, -30.0])
        assert qp.get_interp_val(np.ones(npix), ra, dec).shape == (3,)

    def test_multi_map_returns_every_map(self, qp):
        """The return previously collapsed to the last map."""
        npix = 12 * 8 * 8
        maps = np.array([np.full(npix, 1.0), np.full(npix, 2.0), np.full(npix, 3.0)])
        ra = np.array([0.0, 90.0, 180.0])
        dec = np.array([0.0, 30.0, -30.0])

        val = qp.get_interp_val(maps, ra, dec)
        assert val.shape == (3, 3)
        for i, level in enumerate([1.0, 2.0, 3.0]):
            assert np.allclose(val[i], level, atol=1e-10)

    def test_multi_map_rows_match_individual_calls(self, qp):
        rng = np.random.default_rng(0)
        npix = 12 * 8 * 8
        maps = rng.normal(size=(3, npix))
        ra = np.array([10.0, 45.0, 200.0])
        dec = np.array([5.0, -20.0, 60.0])

        val = qp.get_interp_val(maps, ra, dec)
        for i in range(3):
            assert np.array_equal(val[i], qp.get_interp_val(maps[i], ra, dec))


# ---------------------------------------------------------------------------
# Broadcasting behavior
# ---------------------------------------------------------------------------
# Each class verifies that scalar inputs broadcast correctly against array
# inputs, and that broadcast results are numerically identical to passing
# explicitly repeated arrays.


class TestLmstBroadcast:
    """lmst(ctime, lon) broadcasts ctime and lon against each other."""

    def test_scalar_lon_array_ctime_shape(self, qp):
        ctimes = CTIME + np.arange(5) * 3600.0
        result = qp.lmst(ctimes, LON)
        assert result.shape == (5,)

    def test_array_lon_scalar_ctime_shape(self, qp):
        lons = np.array([0.0, 90.0, 180.0, 270.0])
        result = qp.lmst(CTIME, lons)
        assert result.shape == (4,)

    def test_scalar_lon_matches_repeated_lon(self, qp):
        ctimes = CTIME + np.arange(5) * 3600.0
        r1 = qp.lmst(ctimes, LON)
        r2 = qp.lmst(ctimes, np.full(5, LON))
        assert np.allclose(r1, r2)

    def test_scalar_ctime_matches_repeated_ctime(self, qp):
        lons = np.array([0.0, 90.0, 180.0, 270.0])
        r1 = qp.lmst(CTIME, lons)
        r2 = qp.lmst(np.full(4, CTIME), lons)
        assert np.allclose(r1, r2)


class TestDipoleBroadcast:
    """dipole(ctime, ra, dec) broadcasts ctime, ra, and dec against each other."""

    def test_scalar_radec_array_ctime_shape(self, qp):
        d = qp.dipole(CTIMES, 0.0, 0.0)
        assert d.shape == (N,)

    def test_scalar_radec_matches_repeated_radec(self, qp):
        d1 = qp.dipole(CTIMES, 0.0, 0.0)
        d2 = qp.dipole(CTIMES, np.zeros(N), np.zeros(N))
        assert np.allclose(d1, d2)

    def test_scalar_ctime_array_radec_shape(self, qp):
        ra = np.linspace(0, 360, 8, endpoint=False)
        dec = np.zeros(8)
        d = qp.dipole(CTIME, ra, dec)
        assert d.shape == (8,)

    def test_scalar_ctime_matches_repeated_ctime(self, qp):
        ra = np.linspace(0, 360, 8, endpoint=False)
        dec = np.zeros(8)
        d1 = qp.dipole(CTIME, ra, dec)
        d2 = qp.dipole(np.full(8, CTIME), ra, dec)
        assert np.allclose(d1, d2)


class TestDetOffsetBroadcast:
    """det_offset(delta_az, delta_el, delta_psi) broadcasts all three arguments."""

    def test_scalar_az_el_array_psi_shape(self, qp):
        psi = np.array([0.0, 45.0, 90.0, 135.0])
        q = qp.det_offset(0.0, 0.0, psi)
        assert q.shape == (4, 4)

    def test_scalar_az_el_matches_repeated(self, qp):
        psi = np.array([0.0, 45.0, 90.0, 135.0])
        q1 = qp.det_offset(0.0, 0.0, psi)
        q2 = qp.det_offset(np.zeros(4), np.zeros(4), psi)
        assert np.allclose(q1, q2)

    def test_scalar_psi_array_az_el_shape(self, qp):
        az = np.array([-1.0, 0.0, 1.0])
        el = np.array([0.0, 0.0, 0.5])
        q = qp.det_offset(az, el, 0.0)
        assert q.shape == (3, 4)

    def test_scalar_psi_matches_repeated(self, qp):
        az = np.array([-1.0, 0.0, 1.0])
        el = np.array([0.0, 0.0, 0.5])
        q1 = qp.det_offset(az, el, 0.0)
        q2 = qp.det_offset(az, el, np.zeros(3))
        assert np.allclose(q1, q2)


class TestAzel2BoreBroadcast:
    """azel2bore broadcasts az, el, lon, lat, ctime against each other."""

    def test_scalar_el_array_az_ctime_shape(self, qp):
        q = qp.azel2bore(AZ, EL[0], None, None, LON, LAT, CTIMES)
        assert q.shape == (N, 4)

    def test_scalar_el_matches_repeated_el(self, qp):
        qp.reset_rates()
        q1 = qp.azel2bore(AZ, EL[0], None, None, LON, LAT, CTIMES)
        qp.reset_rates()
        q2 = qp.azel2bore(AZ, np.full(N, EL[0]), None, None, LON, LAT, CTIMES)
        assert np.allclose(q1, q2, atol=1e-12)

    def test_scalar_ctime_array_az_el_shape(self, qp):
        q = qp.azel2bore(AZ, EL, None, None, LON, LAT, CTIME)
        assert q.shape == (N, 4)

    def test_scalar_ctime_matches_repeated_ctime(self, qp):
        qp.reset_rates()
        q1 = qp.azel2bore(AZ, EL, None, None, LON, LAT, CTIME)
        qp.reset_rates()
        q2 = qp.azel2bore(AZ, EL, None, None, LON, LAT, np.full(N, CTIME))
        assert np.allclose(q1, q2, atol=1e-12)

    def test_scalar_lonlat_matches_repeated_lonlat(self, qp):
        qp.reset_rates()
        q1 = qp.azel2bore(AZ, EL, None, None, LON, LAT, CTIMES)
        qp.reset_rates()
        q2 = qp.azel2bore(AZ, EL, None, None, np.full(N, LON), np.full(N, LAT), CTIMES)
        assert np.allclose(q1, q2, atol=1e-12)


class TestAzel2RadecBroadcast:
    """azel2radec broadcasts az, el, lon, lat, ctime against each other."""

    def test_scalar_el_array_az_ctime_shape(self, qp):
        ra, dec, s, c = qp.azel2radec(0, 0, 0, AZ, EL[0], None, None, LON, LAT, CTIMES)
        assert ra.shape == (N,)

    def test_scalar_el_matches_repeated_el(self, qp):
        qp.reset_rates()
        ra1, dec1, s1, c1 = qp.azel2radec(
            0, 0, 0, AZ, EL[0], None, None, LON, LAT, CTIMES
        )
        qp.reset_rates()
        ra2, dec2, s2, c2 = qp.azel2radec(
            0, 0, 0, AZ, np.full(N, EL[0]), None, None, LON, LAT, CTIMES
        )
        assert np.allclose(ra1, ra2, atol=1e-12)
        assert np.allclose(dec1, dec2, atol=1e-12)

    def test_scalar_ctime_array_az_el_shape(self, qp):
        ra, dec, s, c = qp.azel2radec(0, 0, 0, AZ, EL, None, None, LON, LAT, CTIME)
        assert ra.shape == (N,)

    def test_scalar_ctime_matches_repeated_ctime(self, qp):
        qp.reset_rates()
        ra1, dec1, s1, c1 = qp.azel2radec(0, 0, 0, AZ, EL, None, None, LON, LAT, CTIME)
        qp.reset_rates()
        ra2, dec2, s2, c2 = qp.azel2radec(
            0, 0, 0, AZ, EL, None, None, LON, LAT, np.full(N, CTIME)
        )
        assert np.allclose(ra1, ra2, atol=1e-12)
        assert np.allclose(dec1, dec2, atol=1e-12)

    def test_scalar_lonlat_matches_repeated_lonlat(self, qp):
        qp.reset_rates()
        ra1, dec1, s1, c1 = qp.azel2radec(0, 0, 0, AZ, EL, None, None, LON, LAT, CTIMES)
        qp.reset_rates()
        ra2, dec2, s2, c2 = qp.azel2radec(
            0, 0, 0, AZ, EL, None, None, np.full(N, LON), np.full(N, LAT), CTIMES
        )
        assert np.allclose(ra1, ra2, atol=1e-12)
        assert np.allclose(dec1, dec2, atol=1e-12)

    def test_single_sample_output_shape(self, qp):
        """A length-one input keeps its axis; it is an array, not a scalar."""
        ra, dec, s, c = qp.azel2radec(
            0, 0, 0, AZ[:1], EL[:1], None, None, LON, LAT, CTIMES[:1]
        )
        assert ra.shape == (1,)


class TestRadec2AzelBroadcast:
    """radec2azel broadcasts ra, dec, pa, lon, lat, ctime against each other."""

    def _make_radecpa(self, qp):
        qp.reset_rates()
        ra, dec, pa = qp.azel2radec(
            0, 0, 0, AZ, EL, None, None, LON, LAT, CTIMES, return_pa=True
        )
        return ra, dec, pa

    def test_scalar_lonlat_array_sky_shape(self, qp):
        ra, dec, pa = self._make_radecpa(qp)
        qp.reset_inv_rates()
        az, el, hpa = qp.radec2azel(ra, dec, pa, LON, LAT, CTIMES)
        assert az.shape == (N,)

    def test_scalar_lonlat_matches_repeated_lonlat(self, qp):
        ra, dec, pa = self._make_radecpa(qp)
        qp.reset_inv_rates()
        az1, el1, _ = qp.radec2azel(ra, dec, pa, LON, LAT, CTIMES)
        qp.reset_inv_rates()
        az2, el2, _ = qp.radec2azel(
            ra, dec, pa, np.full(N, LON), np.full(N, LAT), CTIMES
        )
        assert np.allclose(az1, az2, atol=1e-12)
        assert np.allclose(el1, el2, atol=1e-12)


class TestBore2AzelBroadcast:
    """bore2azel broadcasts lon, lat, ctime against the q_bore array."""

    def test_scalar_lonlat_array_bore_ctime_shape(self, qp):
        q_bore = qp.azel2bore(AZ, EL, None, None, LON, LAT, CTIMES)
        qp.reset_inv_rates()
        az, el, pa = qp.bore2azel(q_bore, LON, LAT, CTIMES)
        assert az.shape == (N,)

    def test_scalar_lonlat_matches_repeated_lonlat(self, qp):
        q_bore = qp.azel2bore(AZ, EL, None, None, LON, LAT, CTIMES)
        qp.reset_inv_rates()
        az1, el1, _ = qp.bore2azel(q_bore, LON, LAT, CTIMES)
        qp.reset_inv_rates()
        az2, el2, _ = qp.bore2azel(q_bore, np.full(N, LON), np.full(N, LAT), CTIMES)
        assert np.allclose(az1, az2, atol=1e-12)
        assert np.allclose(el1, el2, atol=1e-12)


class TestBoreOffsetBroadcast:
    """bore_offset broadcasts ang1/ang2/ang3 against the boresight array length."""

    def test_scalar_ang1_broadcasts_shape(self, qp):
        q_bore = qp.azel2bore(AZ, EL, None, None, LON, LAT, CTIMES)
        q_out = qp.bore_offset(q_bore.copy(), ang1=1.0)
        assert q_out.shape == (N, 4)

    def test_scalar_ang1_matches_repeated_ang1(self, qp):
        q_bore = qp.azel2bore(AZ, EL, None, None, LON, LAT, CTIMES)
        q1 = qp.bore_offset(q_bore.copy(), ang1=1.0)
        q2 = qp.bore_offset(q_bore.copy(), ang1=np.full(N, 1.0))
        assert np.allclose(q1, q2, atol=1e-12)

    def test_scalar_ang2_broadcasts_shape(self, qp):
        q_bore = qp.azel2bore(AZ, EL, None, None, LON, LAT, CTIMES)
        q_out = qp.bore_offset(q_bore.copy(), ang2=0.5)
        assert q_out.shape == (N, 4)

    def test_scalar_ang3_matches_repeated_ang3(self, qp):
        q_bore = qp.azel2bore(AZ, EL, None, None, LON, LAT, CTIMES)
        q1 = qp.bore_offset(q_bore.copy(), ang3=2.0)
        q2 = qp.bore_offset(q_bore.copy(), ang3=np.full(N, 2.0))
        assert np.allclose(q1, q2, atol=1e-12)


class TestRadecpaQuatBroadcast:
    """radecpa2quat broadcasts ra, dec, pa against each other."""

    def test_scalar_pa_array_radec_shape(self, qp):
        ra = np.array([0.0, 90.0, 180.0, 270.0])
        dec = np.zeros(4)
        q = qp.radecpa2quat(ra, dec, 0.0)
        assert q.shape == (4, 4)

    def test_scalar_pa_matches_repeated_pa(self, qp):
        ra = np.array([0.0, 90.0, 180.0, 270.0])
        dec = np.zeros(4)
        q1 = qp.radecpa2quat(ra, dec, 0.0)
        q2 = qp.radecpa2quat(ra, dec, np.zeros(4))
        assert np.allclose(q1, q2)

    def test_scalar_dec_array_ra_pa_shape(self, qp):
        ra = np.linspace(0, 360, 6, endpoint=False)
        pa = np.zeros(6)
        q = qp.radecpa2quat(ra, 0.0, pa)
        assert q.shape == (6, 4)

    def test_scalar_dec_matches_repeated_dec(self, qp):
        ra = np.linspace(0, 360, 6, endpoint=False)
        pa = np.zeros(6)
        q1 = qp.radecpa2quat(ra, 0.0, pa)
        q2 = qp.radecpa2quat(ra, np.zeros(6), pa)
        assert np.allclose(q1, q2)

    def test_scalar_inputs_squeeze_to_quat(self, qp):
        q = qp.radecpa2quat(0.0, 0.0, 0.0)
        assert q.shape == (4,)


# ---------------------------------------------------------------------------
# Option and rate behaviour, as opposed to the values round-tripping
# ---------------------------------------------------------------------------


class TestMeanAberIsRestored:
    """
    The azel*2radec* family forces mean aberration on for the duration of
    the call. Whatever the caller set has to survive it.
    """

    def test_azel2radec_puts_it_back(self, mod):
        q = mod.QPoint(mean_aber=False, accuracy="low")
        q.azel2radec(1.0, 2.0, 3.0, 45.0, 45.0, None, None, LON, LAT, CTIME)
        assert q.get("mean_aber") is False

    def test_azelpsi2radec_puts_it_back(self, mod):
        q = mod.QPoint(mean_aber=False, accuracy="low")
        q.azelpsi2radec(1.0, 2.0, 3.0, 45.0, 45.0, 10.0, None, None, LON, LAT, CTIME)
        assert q.get("mean_aber") is False

    def test_it_is_still_on_where_the_caller_asked_for_it(self, mod):
        q = mod.QPoint(mean_aber=True, accuracy="low")
        q.azel2radec(1.0, 2.0, 3.0, 45.0, 45.0, None, None, LON, LAT, CTIME)
        assert q.get("mean_aber") is True


class TestRateCaching:
    """
    A rate says how often a correction is recomputed. Only the values
    round-tripping is covered elsewhere; this is what the rates do.
    """

    # two samples 200 days apart, so a frozen correction is visible
    TIMES = CTIME + np.array([0.0, 200.0 * 86400.0])

    def radec(self, mod, **kwargs):
        q = mod.QPoint(mean_aber=True, **kwargs)
        az, el = np.array([10.0, 10.0]), np.array([45.0, 45.0])
        q_bore = q.azel2bore(az, el, None, None, LON, LAT, self.TIMES)
        ra, dec, _, _ = q.bore2radec(q.det_offset(0.0, 0.0, 0.0), self.TIMES, q_bore)
        return np.asarray(ra).copy()

    def test_never_differs_from_always(self, mod):
        assert not np.allclose(
            self.radec(mod, rate_npb="never"), self.radec(mod, rate_npb="always")
        )

    def test_once_is_computed_at_the_first_sample(self, mod):
        once, always = self.radec(mod, rate_npb="once"), self.radec(
            mod, rate_npb="always"
        )
        assert np.isclose(once[0], always[0])

    def test_once_is_then_frozen(self, mod):
        once, always = self.radec(mod, rate_npb="once"), self.radec(
            mod, rate_npb="always"
        )
        assert not np.isclose(once[1], always[1])


class TestFastPix:
    """
    fast_pix skips the angle round trip and takes the pixel from the
    pointing vector. Away from the poles the two agree exactly, and the
    polarization angle agrees everywhere: both packages form cos^2(b) from
    the quaternion and hand it to one shared body, so there is nothing
    left to differ.
    """

    @pytest.mark.parametrize("order", ["ring", "nest"])
    @pytest.mark.parametrize(
        "lo, hi, label",
        [
            (-90.0, 90.0, "whole sky"),
            (89.0, 90.0, "north cap"),
            (-90.0, -89.0, "south cap"),
            (89.9999, 90.0, "inside the degenerate branch"),
        ],
    )
    def test_the_polarization_angle_is_exact(self, order, lo, hi, label):
        """
        Bit-identical, not close: cos^2(b) used to come from the pointing vector
        as (1 - vec[2]**2) / 4, which cancels as vec[2] tends to +/-1 and is then
        divided by. The caps are the point -- the old error was 2e-10 within a
        degree of a pole, so a whole-sky sample alone would barely see it.
        """
        rng = np.random.default_rng(4)
        n = 4000
        ra = rng.uniform(0.0, 360.0, n)
        dec = rng.uniform(lo, hi, n)
        pa = rng.uniform(-180.0, 180.0, n)
        q = qpoint.QPoint(mean_aber=True, pix_order=order)
        quat = np.asarray(q.radecpa2quat(ra, dec, pa), dtype=float)
        _, sin_f, cos_f = q.quat2pix(quat, nside=256, fast_pix=True)
        _, sin_s, cos_s = q.quat2pix(quat, nside=256, fast_pix=False)
        assert np.array_equal(np.asarray(sin_f), np.asarray(sin_s)), label
        assert np.array_equal(np.asarray(cos_f), np.asarray(cos_s)), label

    @pytest.mark.parametrize("fast_math", [False, True])
    def test_the_position_angle_is_exact(self, fast_math):
        """The same for quat2pixpa, which has its own copy of the path."""
        rng = np.random.default_rng(5)
        n = 4000
        ra = rng.uniform(0.0, 360.0, n)
        dec = np.concatenate(
            [rng.uniform(-90.0, 90.0, n // 2), rng.uniform(89.9, 90.0, n - n // 2)]
        )
        pa = rng.uniform(-180.0, 180.0, n)
        q = qpoint.QPoint(mean_aber=True, fast_math=fast_math)
        quat = np.asarray(q.radecpa2quat(ra, dec, pa), dtype=float)
        _, pa_fast = q.quat2pixpa(quat, nside=256, fast_pix=True)
        _, _, pa_slow = q.quat2radecpa(quat)
        assert np.array_equal(np.asarray(pa_fast), np.asarray(pa_slow))

    def test_it_is_on_by_default(self):
        """
        It costs nothing in accuracy now, and is worth 26-34% of tod2map.
        The pixel differs from the angle path only inside 2.1e-8 rad of a
        pole, where the fast path is the correct one.
        """
        assert qpoint.QPoint().get("fast_pix") is True
        assert qpoint.QMap().get("fast_pix") is True


FP_NSIDE = 128
FP_N = 50
FP_CTIME = CTIME + np.arange(FP_N, dtype=float)
FP_RA = np.linspace(0.0, 350.0, FP_N)
FP_DEC = np.linspace(-80.0, 80.0, FP_N)
FP_PA = np.linspace(-170.0, 170.0, FP_N)
FP_AZ = np.linspace(0.0, 360.0, FP_N)
FP_EL = np.linspace(30.0, 70.0, FP_N)
FP_LON = np.full(FP_N, LON)
FP_LAT = np.full(FP_N, LAT)


def fp_bore(**options):
    """A QPoint and a boresight quaternion to go with it."""
    q = qpoint.QPoint(**options)
    qb = q.azel2bore(FP_AZ, FP_EL, None, None, FP_LON, FP_LAT, FP_CTIME)
    return q, np.asarray(qb)


def same(expected, actual, label):
    """
    Exact equality, NaNs in matching slots included: a tolerance would
    defeat the point of holding a fast path to a slow one.
    """
    if not isinstance(expected, tuple):
        expected, actual = (expected,), (actual,)
    assert len(expected) == len(actual), "{}: arity differs".format(label)
    for i, (e, a) in enumerate(zip(expected, actual)):
        e, a = np.asarray(e), np.asarray(a)
        assert e.shape == a.shape, "{}[{}] shape differs".format(label, i)
        eq = np.array_equal(e, a, equal_nan=np.issubdtype(e.dtype, np.floating))
        assert eq, "{}[{}] differs".format(label, i)


class TestFastPixAgainstTheAnglePath:
    """
    fast_pix takes the pixel from the pointing vector instead of going
    round through ra/dec, and has to agree with that two-step path
    exactly. Crossed with pol and the pixel ordering, which the fast path
    picks itself.
    """

    @pytest.mark.parametrize("order", ["ring", "nest"])
    @pytest.mark.parametrize(
        "kwargs",
        [pytest.param({}, id="pol"), pytest.param({"pol": False}, id="no-pol")],
    )
    def test_quat2pix(self, order, kwargs):
        """fast_pix=False is the two-step path, taken inside quat2pix."""
        q = qpoint.QPoint(pix_order=order)
        quat = q.radecpa2quat(FP_RA, FP_DEC, FP_PA)
        slow = q.quat2pix(quat, nside=FP_NSIDE, fast_pix=False, **kwargs)
        fast = q.quat2pix(quat, nside=FP_NSIDE, fast_pix=True, **kwargs)
        same(tuple(slow), tuple(fast), "quat2pix fast vs two-step")

    @pytest.mark.parametrize("order", ["ring", "nest"])
    def test_quat2pixpa_against_the_public_two_step(self, order):
        """Here the two steps are separately reachable, so spell them out."""
        q = qpoint.QPoint(pix_order=order)
        quat = q.radecpa2quat(FP_RA, FP_DEC, FP_PA)
        ra, dec, pa = q.quat2radecpa(quat)
        want = (np.asarray(q.radec2pix(ra, dec, nside=FP_NSIDE)), np.asarray(pa))
        got = tuple(q.quat2pixpa(quat, nside=FP_NSIDE, fast_pix=True))
        same(want, got, "quat2pixpa fast vs quat2radecpa+radec2pix")

    @pytest.mark.parametrize("order", ["ring", "nest"])
    @pytest.mark.parametrize(
        "kwargs",
        [
            pytest.param({}, id="pol"),
            pytest.param({"pol": False}, id="no-pol"),
            pytest.param({"return_pa": True}, id="pa"),
        ],
    )
    def test_bore2pix(self, order, kwargs):
        q, qb = fp_bore(pix_order=order)
        off = q.det_offset(1.0, 2.0, 3.0)
        slow = q.bore2pix(off, FP_CTIME, qb, nside=FP_NSIDE, fast_pix=False, **kwargs)
        fast = q.bore2pix(off, FP_CTIME, qb, nside=FP_NSIDE, fast_pix=True, **kwargs)
        same(tuple(slow), tuple(fast), "bore2pix fast vs two-step")

    @pytest.mark.parametrize("order", ["ring", "nest"])
    def test_near_the_poles(self, order):
        """
        The polarization angle is exact at the poles, where the old form lost
        2e-10 within a degree of one.

        The pixel is a separate matter: within about a microdegree of a pole the
        two routes land on adjacent pixels, which no choice of cos^2(b) changes,
        so the pixel is only required to match outside that sliver.
        """
        dec = np.array([90.0, -90.0, 89.999999, -89.999999, 89.99, -89.99, 89.9, 89.0])
        ra = np.linspace(0.0, 350.0, len(dec))
        pa = np.linspace(-170.0, 170.0, len(dec))
        q = qpoint.QPoint(pix_order=order)
        quat = q.radecpa2quat(ra, dec, pa)
        pix_s, sin_s, cos_s = q.quat2pix(quat, nside=FP_NSIDE, fast_pix=False)
        pix_f, sin_f, cos_f = q.quat2pix(quat, nside=FP_NSIDE, fast_pix=True)
        same((sin_s, cos_s), (sin_f, cos_f), "pol angle at the poles")
        settled = np.abs(dec) <= 89.99
        same(
            np.asarray(pix_s)[settled],
            np.asarray(pix_f)[settled],
            "pixel away from the pole sliver",
        )

    @pytest.mark.parametrize("order", ["ring", "nest"])
    def test_the_fast_path_is_the_better_one_at_the_pole(self, order):
        """
        Inside theta < 2.1e-8 rad the two disagree, and the slow path is the wrong
        one: its cos(theta) rounds to 1, throwing the azimuth away and dumping
        every direction into pixel 0. So no agreement is forced there. The
        azimuths sit inside the quadrants rather than on their boundaries, where
        the assignment is a tie-break.
        """
        ra = np.array(
            [
                15.0,
                45.0,
                75.0,
                105.0,
                135.0,
                165.0,
                195.0,
                225.0,
                255.0,
                285.0,
                315.0,
                345.0,
            ]
        )
        pa = np.full(len(ra), 17.0)
        q = qpoint.QPoint(pix_order=order)

        def pix(theta, fast):
            dec = np.full(len(ra), 90.0 - np.degrees(theta))
            quat = q.radecpa2quat(ra, dec, pa)
            return np.asarray(q.quat2pix(quat, nside=FP_NSIDE, fast_pix=fast)[0])

        # far enough out that no rounding is in play, so this is the truth
        want = pix(1e-3, False)
        assert len(np.unique(want)) == 4, "expected one pixel per quadrant"
        same(want, pix(1e-3, True), "quadrants away from the pole")

        for theta in (2e-8, 1e-8, 1e-12):
            same(want, pix(theta, True), "fast_pix in the cap")
            assert np.all(pix(theta, False) == want[0]), "slow path collapses"


class TestBulletinARange:
    def test_a_lookup_outside_the_table_returns_zeros(self, mod):
        """
        Every caller in the C ignores the error return and uses the
        values, so an out-of-range date has to leave them at zero rather
        than raise.
        """
        dut1, x, y = mod.QPoint().get_bulletin_a(20000.0)
        assert (dut1, x, y) == (0.0, 0.0, 0.0)


class TestPrintMemory:
    def test_it_prints_the_state(self, mod, capfd):
        """
        print_memory writes from the C, so the file descriptor has to be
        captured rather than sys.stdout. It flushes itself, which is what
        makes the output readable here rather than after the test.
        """
        mod.QPoint(accuracy="low").print_memory()
        out = capfd.readouterr().out
        # the two head their dumps differently
        assert "QPOINT MEMORY" in out or "qpoint2 Pointing" in out
        assert "accuracy" in out


class TestInverseRatesAreIndependent:
    """
    The inverse transform has its own update rates, and they work, in both
    packages.

    The two aberrations took an inv flag and ignored it when checking the
    rate, so rate_daber_inv and rate_aaber_inv did nothing and a pass in one
    direction moved the other's clock -- which is what the separate states
    exist to prevent.
    """

    N = 200

    def bore(self, mod):
        ct = CTIME + np.arange(self.N) / 10.0
        lon = np.full(self.N, LON)
        lat = np.full(self.N, LAT)
        q = mod.QPoint()
        q_bore = np.asarray(
            q.azel2bore(
                np.linspace(0, 90, self.N),
                np.full(self.N, 45.0),
                None,
                None,
                lon,
                lat,
                ct,
            )
        )
        return q_bore, lon, lat, ct

    def inverse_el(self, mod, **kwargs):
        q_bore, lon, lat, ct = self.bore(mod)
        _, el, _ = mod.QPoint(**kwargs).bore2azel(q_bore, lon, lat, ct)
        return np.asarray(el)

    @pytest.mark.parametrize(
        "rate, inv_rate, floor_arcsec",
        [
            # annual aberration is ~20 arcsec, diurnal a fraction of one
            ("rate_aaber", "rate_aaber_inv", 1.0),
            ("rate_daber", "rate_daber_inv", 0.01),
        ],
    )
    def test_the_inverse_rate_controls_the_inverse(
        self, mod, rate, inv_rate, floor_arcsec
    ):
        """
        With the forward rate off, only the inverse one can act -- so if it
        is being read at all, turning it on has to move the answer.
        """
        on = self.inverse_el(mod, **{rate: "never", inv_rate: "always"})
        off = self.inverse_el(mod, **{rate: "never", inv_rate: "never"})
        assert np.abs(on - off).max() * 3600.0 > floor_arcsec

    @pytest.mark.parametrize(
        "rate, inv_rate",
        [("rate_aaber", "rate_aaber_inv"), ("rate_daber", "rate_daber_inv")],
    )
    def test_the_forward_rate_no_longer_decides_it(self, mod, rate, inv_rate):
        """
        The other half: with the inverse rate off, the forward one must not
        be able to switch the correction on behind its back.
        """
        a = self.inverse_el(mod, **{rate: "always", inv_rate: "never"})
        b = self.inverse_el(mod, **{rate: "never", inv_rate: "never"})
        assert np.array_equal(a, b)

    def test_neither_direction_moves_the_other_s_clock(self, mod):
        """
        qp_check_update writes ctime_last, so sharing a state let a run in
        one direction leave the other due for an update at a time it had
        already covered. Independent states mean an inverse pass cannot
        change what a forward pass then produces.
        """
        q_bore, lon, lat, ct = self.bore(mod)
        az = np.linspace(0, 90, self.N)
        el = np.full(self.N, 45.0)

        clean = mod.QPoint()
        want = np.asarray(clean.azel2bore(az, el, None, None, lon, lat, ct))

        used = mod.QPoint()
        used.bore2azel(q_bore, lon, lat, ct)  # run the inverse first
        got = np.asarray(used.azel2bore(az, el, None, None, lon, lat, ct))
        assert np.array_equal(want, got)

    def test_the_round_trip_still_closes(self, mod):
        q_bore, lon, lat, ct = self.bore(mod)
        az = np.linspace(0, 90, self.N)
        q = mod.QPoint()
        got_az, got_el, _ = q.bore2azel(q_bore, lon, lat, ct)
        assert np.abs(np.asarray(got_az) % 360 - az % 360).max() < 1e-5
        assert np.abs(np.asarray(got_el) - 45.0).max() < 1e-5


class TestRateStateSentinel:
    """
    An un-updated rate state is marked with NAN, not with a negative time.

    qp_check_update used to read ctime_last <= 0 as "never updated", which is
    also a legitimate ctime, so a run at non-positive ctime recomputed every
    correction on every sample. 'once' was worse: its "already done" guard was
    ctime_last > 0, which a negative ctime never satisfies, so it recomputed
    every sample -- the opposite of what was asked for.
    """

    N = 60

    def bore(self, ctime, **kwargs):
        n = len(ctime)
        q = qpoint.QPoint(**kwargs)
        out = np.asarray(
            q.azel2bore(
                np.linspace(0, 90, n),
                np.full(n, 45.0),
                None,
                None,
                np.full(n, LON),
                np.full(n, LAT),
                np.ascontiguousarray(ctime),
            )
        )
        # keep q alive until after the read: its __del__ frees the struct
        return out, q._memory.contents.state_npb.ctime_last

    @pytest.mark.parametrize(
        "label, ctime",
        [
            ("positive", CTIME + np.arange(N) * 5.0),
            ("negative", -50000.0 + np.arange(N) * 5.0),
            ("straddling zero", -150.0 + np.arange(N) * 5.0),
        ],
        ids=lambda v: v if isinstance(v, str) else "",
    )
    def test_once_freezes_whatever_the_epoch(self, label, ctime):
        """
        The decisive one: at non-positive ctime `once` used to behave like
        `always`, because its guard could never fire.
        """
        once, _ = self.bore(ctime, rate_npb="once")
        always, _ = self.bore(ctime, rate_npb="always")
        ten, _ = self.bore(ctime, rate_npb=10)
        assert not np.array_equal(once, always)
        assert not np.array_equal(once, ten)

    def test_a_fresh_state_is_nan(self):
        q = qpoint.QPoint()
        assert np.isnan(q._memory.contents.state_npb.ctime_last)

    def test_a_reset_returns_it_to_nan(self):
        q = qpoint.QPoint()
        ct = CTIME + np.arange(self.N) * 5.0
        q.azel2bore(
            np.linspace(0, 90, self.N),
            np.full(self.N, 45.0),
            None,
            None,
            np.full(self.N, LON),
            np.full(self.N, LAT),
            np.ascontiguousarray(ct),
        )
        assert not np.isnan(q._memory.contents.state_npb.ctime_last)
        q.reset_rates()
        assert np.isnan(q._memory.contents.state_npb.ctime_last)

    def test_the_state_records_the_real_time(self):
        """Including a negative one, rather than being stuck at a sentinel."""
        ctime = -50000.0 + np.arange(self.N) * 5.0
        _, last = self.bore(ctime, rate_npb="once")
        assert last == ctime[0]


# ---------------------------------------------------------------------------
# the UTC -> UT1 day cache
# ---------------------------------------------------------------------------

# 2015-07-01 00:00:00 UTC, the instant after that year's leap second.
LEAP_2015 = 1435708800.0

# Every correction recomputed per sample, so the day cache is the only thing
# a run of samples reuses. Without this the rate caches dominate: with
# rate_npb at its default of 10 s, samples a second apart share one
# nutation series and a fresh QPoint per sample does not, which looks like
# a caching bug and is not one.
ALL_RATES = dict(
    rate_npb="always",
    rate_aaber="always",
    rate_daber="always",
    rate_lonlat="always",
    rate_erot="always",
    rate_wobble="always",
)

UT1_SPANS = {
    "ordinary day": CTIME + np.arange(0, 86400, 3600.0),
    "across midnight": CTIME + np.arange(43000, 44200, 120.0),
    "leap second day": LEAP_2015 - 86400.0 + np.arange(0, 86401, 3600.0),
    # the margin itself: the cache stops ten seconds short of midnight
    "margin edges": LEAP_2015
    - 86400.0
    + np.array([0.0, 9.0, 10.0, 11.0, 86389.0, 86390.0, 86400.0]),
    "leap second crossing": LEAP_2015 + np.arange(-5, 5, 0.5),
    "one second apart": CTIME + np.arange(0, 20, 1.0),
    "100 Hz": CTIME + np.arange(0, 2, 0.01),
    "backwards": CTIME + np.arange(86400, 0, -3600.0),
}


class TestUt1Caching:
    """
    UT1 - UTC is cached for the interior of a calendar day, which is exact:
    the offset steps only at a leap second, and those fall at midnight.

    The awkward part is that a leap-second UTC day is 86401 seconds long, so
    ERFA and ctime / 86400 disagree about the date for one second. The cache
    stays ten seconds clear of midnight rather than reason about it, and a day
    with a leap second in it gets no cache.

    Being exact is the claim, so every test here compares an accumulating run
    against the same samples computed one at a time on a fresh QPoint.
    """

    @staticmethod
    def warm(ctime, **kwargs):
        q = qpoint.QPoint(**dict(ALL_RATES, **kwargs))
        n = len(ctime)
        return np.asarray(
            q.azel2bore(
                np.linspace(0, 90, n),
                np.full(n, 45.0),
                None,
                None,
                np.full(n, LON),
                np.full(n, LAT),
                np.ascontiguousarray(ctime),
            )
        )

    @staticmethod
    def cold(ctime, **kwargs):
        n = len(ctime)
        az = np.linspace(0, 90, n)
        out = [
            np.asarray(
                qpoint.QPoint(**dict(ALL_RATES, **kwargs)).azel2bore(
                    az[i], 45.0, None, None, LON, LAT, ctime[i]
                )
            )
            for i in range(n)
        ]
        return np.asarray(out).reshape(n, 4)

    @pytest.mark.parametrize(
        "span", sorted(UT1_SPANS), ids=lambda s: s.replace(" ", "-")
    )
    def test_the_cache_changes_nothing(self, span):
        ctime = UT1_SPANS[span]
        assert np.array_equal(self.warm(ctime), self.cold(ctime))

    @pytest.mark.parametrize("accuracy", ["high", "low"])
    def test_in_both_accuracy_modes(self, accuracy):
        ctime = UT1_SPANS["leap second day"]
        assert np.array_equal(
            self.warm(ctime, accuracy=accuracy), self.cold(ctime, accuracy=accuracy)
        )

    def test_a_dut1_that_changes_daily(self):
        """
        The cache keys on dut1 as well as the day, so a bulletin that steps
        from one day to the next has to invalidate it rather than carry an
        offset across the boundary.
        """
        nday = 6
        mjd0 = 57200
        dut1 = np.ascontiguousarray(np.linspace(-0.4, 0.4, nday))
        zeros = np.ascontiguousarray(np.zeros(nday))

        def go(ctime, cold):
            def build():
                q = qpoint.QPoint(rate_dut1="always", **ALL_RATES)
                libqp.qp_set_iers_bulletin_a(
                    q._memory, mjd0, mjd0 + nday - 1, dut1, zeros, zeros
                )
                return q

            n = len(ctime)
            az = np.linspace(0, 90, n)
            if not cold:
                return np.asarray(
                    build().azel2bore(
                        az,
                        np.full(n, 45.0),
                        None,
                        None,
                        np.full(n, LON),
                        np.full(n, LAT),
                        np.ascontiguousarray(ctime),
                    )
                )
            return np.asarray(
                [
                    np.asarray(
                        build().azel2bore(az[i], 45.0, None, None, LON, LAT, ctime[i])
                    )
                    for i in range(n)
                ]
            ).reshape(n, 4)

        ctime = LEAP_2015 + np.arange(0, 5 * 86400.0, 7200.0)
        assert np.array_equal(go(ctime, False), go(ctime, True))

    def test_the_cache_starts_empty(self):
        """
        A fresh QPoint must not answer from an uninitialized window, which
        is what lo > hi is for.
        """
        q = qpoint.QPoint()
        cache = q._memory.contents.ut1_cache
        assert cache.lo > cache.hi


class TestErotSlerp:
    """
    rate_erot as an interval interpolates the earth rotation quaternion
    across a window rather than holding it, which is exact for a
    constant-rate turn about a fixed axis.

    The floor is the two-part JD: one ulp of jd[1] at this epoch moves ERA
    by 0.0047 mas, and the two paths agree to less than that.
    """

    # one ulp of jd[1] near 2014, as mas of earth rotation
    ULP_MAS = 0.0047

    @staticmethod
    def bore(ctime, **kwargs):
        q = qpoint.QPoint(**dict(ALL_RATES, **kwargs))
        n = len(ctime)
        return np.asarray(
            q.azel2bore(
                np.linspace(0, 90, n),
                np.full(n, 45.0),
                None,
                None,
                np.full(n, LON),
                np.full(n, LAT),
                np.ascontiguousarray(ctime),
            )
        )

    @staticmethod
    def worst_mas(a, b):
        """Worst rotation between two runs of boresight quaternions, in mas."""
        a = np.array(a, copy=True)
        a[np.sum(a * b, axis=1) < 0] *= -1  # q and -q are the same rotation
        d = np.linalg.norm(a - b, axis=1)
        return float(
            np.max(np.rad2deg(2 * np.arcsin(np.clip(d / 2, 0, 1))) * 3.6e6 * 2)
        )

    def test_a_sample_on_a_grid_point_is_exact(self):
        """
        The window covers floor(ctime/rate)*rate, so an integer ctime with a
        one-second window lands on an endpoint.
        """
        ctime = CTIME + np.arange(20, dtype=float)
        assert np.array_equal(
            self.bore(ctime, rate_erot=1.0), self.bore(ctime, rate_erot="always")
        )

    @pytest.mark.parametrize("rate", [0.01, 1.0, 100.0, 1000.0])
    def test_within_one_ulp_of_the_per_sample_path(self, rate):
        """
        Off the grid the two paths differ by less than the JD can resolve, at
        every interval.
        """
        ctime = CTIME + np.arange(2000) * 0.01
        got = self.worst_mas(
            self.bore(ctime, rate_erot=rate), self.bore(ctime, rate_erot="always")
        )
        assert got < self.ULP_MAS

    def test_a_window_spanning_a_leap_second_is_refused(self):
        """
        ERA is affine in ctime only where ctime -> UT1 is a constant offset.
        The window is left empty across a leap second and those samples are
        computed exactly.
        """
        q = qpoint.QPoint(**dict(ALL_RATES, rate_erot=1.0))
        q.azel2bore(0.0, 45.0, None, None, LON, LAT, LEAP_2015 - 0.5)
        cache = q._memory.contents.erot_cache
        assert cache.lo > cache.hi

    def test_across_a_leap_second(self):
        ctime = LEAP_2015 + np.arange(-5, 5, 0.25)
        assert np.isfinite(self.bore(ctime, rate_erot=1.0)).all()
        assert (
            self.worst_mas(
                self.bore(ctime, rate_erot=1.0), self.bore(ctime, rate_erot="always")
            )
            < self.ULP_MAS
        )

    def test_backwards_time(self):
        """A miss rebuilds, so time need not arrive in order."""
        ctime = CTIME + np.arange(200, 0, -1.0)
        assert (
            self.worst_mas(
                self.bore(ctime, rate_erot=1.0), self.bore(ctime, rate_erot="always")
            )
            < self.ULP_MAS
        )

    def test_an_interval_short_enough_to_collapse_the_window(self):
        """
        A window this short leaves the endpoints equal to within a double, so
        alpha is zero and interpolation returns the endpoint.
        """
        ctime = CTIME + np.arange(50) * 0.01
        got = self.bore(ctime, rate_erot=1e-4)
        assert np.isfinite(got).all()
        assert self.worst_mas(got, self.bore(ctime, rate_erot="always")) < self.ULP_MAS

    def test_always_does_not_build_a_window(self):
        """'always' takes the exact path, so nothing reaches the cache."""
        q = qpoint.QPoint(**dict(ALL_RATES, rate_erot="always"))
        q.azel2bore(0.0, 45.0, None, None, LON, LAT, CTIME)
        cache = q._memory.contents.erot_cache
        assert cache.lo > cache.hi

    def test_the_cache_starts_empty(self):
        """The QPoint has to outlive the read: its __del__ frees the memory."""
        q = qpoint.QPoint()
        cache = q._memory.contents.erot_cache
        assert cache.lo > cache.hi

    def test_the_ctypes_mirror_agrees_with_the_struct(self):
        """
        Reads every field of the window back, nested slerp included, against
        the compiled layout.

        alpha is the angle between the endpoint quaternions, which for a
        rotation about a fixed axis is half the angle turned: the earth's
        7.2921e-5 rad/s over the window.
        """
        rate = 4.0
        q = qpoint.QPoint(**dict(ALL_RATES, rate_erot=rate))
        q.azel2bore(0.0, 45.0, None, None, LON, LAT, CTIME + 0.5)
        cache = q._memory.contents.erot_cache

        assert cache.hi - cache.lo == pytest.approx(rate)
        assert cache.lo <= CTIME + 0.5 < cache.hi
        assert cache.lo % rate == 0.0  # aligned to an absolute grid
        assert cache.rate == rate
        assert cache.dut1 == 0.0
        assert cache.slerp.alpha == pytest.approx(7.2921159e-5 * rate / 2, rel=1e-6)
        assert cache.slerp.sin_alpha == pytest.approx(cache.slerp.alpha, rel=1e-9)
        for end in (cache.slerp.q0, cache.slerp.q1):
            assert sum(x * x for x in end) == pytest.approx(1.0)


class TestLowAccuracyKeepsDut1:
    """
    gmst and lmst apply dut1 in both accuracy modes.

    The low path used to hand UTC to eraGmst00 as UT1, discarding a loaded
    bulletin -- 6 arcsec at a typical dut1 -- and misreading a leap-second
    day. Cached for the day the conversion costs about 2 ns, so low now gives
    up only the TT conversion, worth 0.1 mas.
    """

    # One leap-second day, sampled through it so any drift is visible.
    CTIME = np.ascontiguousarray(LEAP_2015 - 86400.0 + np.arange(0.0, 86400.0, 600.0))
    DUT1 = 0.4
    MJD0 = 57200
    NDAY = 40

    def qp(self, accuracy, dut1):
        q = qpoint.QPoint(accuracy=accuracy, rate_dut1="always")
        d = np.ascontiguousarray(np.full(self.NDAY, float(dut1)))
        z = np.ascontiguousarray(np.zeros(self.NDAY))
        libqp.qp_set_iers_bulletin_a(
            q._memory, self.MJD0, self.MJD0 + self.NDAY - 1, d, z, z
        )
        return q

    @staticmethod
    def arcsec(a, b):
        """Separation in arcsec of earth rotation; gmst is in hours."""
        return np.abs(np.asarray(a) - np.asarray(b)) * 15.0 * 3600.0

    def gmst(self, accuracy, dut1=None):
        dut1 = self.DUT1 if dut1 is None else dut1
        return np.asarray(self.qp(accuracy, dut1).gmst(self.CTIME))

    def test_low_tracks_high(self):
        """Only the TT term is given up: a tenth of a milliarcsecond."""
        assert self.arcsec(self.gmst("low"), self.gmst("high")).max() < 1e-3

    def test_the_term_is_really_applied(self):
        """
        Against a QPoint with no bulletin at all, so the test cannot pass
        by the low path quietly doing nothing.
        """
        plain = np.asarray(qpoint.QPoint(accuracy="low").gmst(self.CTIME))
        assert self.arcsec(self.gmst("low"), plain).min() > 0.9 * self.DUT1 * 15.0

    def test_the_leap_second_day_needs_no_bulletin(self):
        """
        With dut1 == 0 the low path still has to agree, which it did not
        before: the day's 86401 seconds were mapped onto 86400.
        """
        d = self.arcsec(self.gmst("low", dut1=0.0), self.gmst("high", dut1=0.0))
        assert d.max() < 1e-3

    def test_lmst_follows_gmst(self):
        """lmst is gmst plus a longitude, so it inherits the whole thing."""
        out = [
            np.asarray(self.qp(acc, self.DUT1).lmst(self.CTIME, LON))
            for acc in ("low", "high")
        ]
        assert self.arcsec(*out).max() < 1e-3

    def test_high_accuracy_applies_it_too(self):
        """
        It always did -- this is the reference the low path is held to, so
        it is worth pinning that it is not itself ignoring the bulletin.
        """
        plain = np.asarray(qpoint.QPoint(accuracy="high").gmst(self.CTIME))
        assert self.arcsec(self.gmst("high"), plain).min() > 0.9 * self.DUT1 * 15.0


# ---------------------------------------------------------------------------
# qpoint2 only: qp_settings has no counterpart in qpoint, and the
# zero-copy contract is a qpoint2 guarantee. Neither compares the two
# packages, so neither belongs in test_parity.py.
# ---------------------------------------------------------------------------


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
        assert q.scan(CTIME, accuracy="low") != q.scan(CTIME)
        assert q.get("accuracy") == before

    def test_the_methods_own_arguments_still_reach_it(self):
        q = self.subclass()
        assert q.scan(CTIME, offset=1.0) == q.scan(CTIME) + 1.0
        assert q.scan(CTIME, offset=1.0, accuracy="low") == (
            q.scan(CTIME, accuracy="low") + 1.0
        )

    def test_an_unknown_keyword_is_reported(self):
        with pytest.raises(TypeError):
            self.subclass().scan(CTIME, nonsense=1)

    def test_defaults_apply_under_the_caller(self):
        class Defaulted(qpoint2.QPoint):
            @qpoint2.qp_settings(accuracy="low")
            def which(self):
                return self.get_param("accuracy")

        q = Defaulted()
        assert q.which() == "low"
        assert q.which(accuracy="high") == "high"
        assert q.get("accuracy") == qpoint2.QPoint().get("accuracy")


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
            qpoint2.QPoint().det_offset(bad, ok, ok)

    def test_error_names_the_argument(self):
        with pytest.raises(TypeError, match="delta_el"):
            qpoint2.QPoint().det_offset(np.ones(3), np.ones(3, dtype=np.float32), 0.0)

    def test_length_mismatch_is_rejected(self):
        with pytest.raises(ValueError, match="length"):
            qpoint2.QPoint().det_offset(np.ones(3), np.ones(5), 0.0)
