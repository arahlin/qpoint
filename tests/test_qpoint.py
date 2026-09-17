"""Tests for qpoint.QPoint coordinate conversion class."""

import numpy as np
import pytest
import qpoint

# Reference observing parameters
CTIME = 1418662800.0  # 2014-12-15 12:00 UTC
LON = 165.7  # McMurdo Station longitude
LAT = -77.6  # McMurdo Station latitude
N = 20
AZ = 100.0 + np.linspace(0, 10, N)
EL = 32.0 * np.ones(N)
CTIMES = CTIME + np.arange(N, dtype=float)


@pytest.fixture
def qp():
    return qpoint.QPoint(mean_aber=True, accuracy="low")


# ---------------------------------------------------------------------------
# Initialization
# ---------------------------------------------------------------------------


class TestInit:
    def test_default_init(self):
        q = qpoint.QPoint()
        assert q is not None

    def test_kwargs_in_init(self):
        q = qpoint.QPoint(accuracy="low", fast_math=True, mean_aber=True)
        assert q.get("accuracy") == "low"
        assert q.get("fast_math") is True
        assert q.get("mean_aber") is True

    def test_del(self):
        q = qpoint.QPoint()
        del q  # should not raise


# ---------------------------------------------------------------------------
# set / get
# ---------------------------------------------------------------------------


class TestSetGet:
    def test_get_all_returns_nested_dict(self, qp):
        state = qp.get()
        assert isinstance(state, dict)
        for key in ("rates", "options", "weather", "params"):
            assert key in state

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


class TestResetRates:
    def test_reset_rates(self, qp):
        qp.reset_rates()

    def test_reset_inv_rates(self, qp):
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

    def test_sindec_with_return_pa_raises(self, qp):
        q_bore = self._make_bore(qp)
        q_off = self._make_qoff(qp)
        with pytest.raises(ValueError):
            qp.bore2radec(q_off, CTIMES, q_bore, sindec=True, return_pa=True)

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


# ---------------------------------------------------------------------------
# get_interp_val
# ---------------------------------------------------------------------------


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

    def test_scalar_inputs_squeeze_to_quat(self, qp):
        q = qp.det_offset(0.0, 0.0, 0.0)
        assert q.shape == (4,)


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
        # azel2radec always returns arrays (no scalar squeeze for n=1).
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
