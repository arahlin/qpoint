"""
Timing for the paths a real run spends its time in.

These measure, they do not check: there is no threshold to fail against,
because the number that matters is the one on your machine compared with
the one you got yesterday. They skip unless `pytest --benchmark` is given.

Each still asserts that the call produced something of the right shape,
so a benchmark cannot quietly end up timing a no-op -- which is the usual
way a benchmark suite rots.

The workloads are sized like a chunk of a real observation rather than
like a unit test: ~17 minutes of one detector at 100 Hz for the pointing,
and 8 detectors over 200 seconds for the mapmaking.
"""

import numpy as np
import pytest
import qpoint

pytestmark = pytest.mark.benchmark

# ~17 minutes at 100 Hz, which is the order of a scan chunk
NPOINT = 100_000
# mapmaking: 8 detectors over 200 s, into a patch-sized map
NDET, NMAP, NSIDE = 8, 20_000, 128

LON, LAT = 165.7, -77.6
CTIME = 1418662800.0


def scan(n):
    """A raster-ish scan: azimuth sweeping, elevation stepping slowly."""
    ctime = CTIME + np.arange(n) / 100.0
    az = 130.0 + 50.0 * np.sin(np.linspace(0, 40 * np.pi, n))
    el = 45.0 + 0.5 * np.floor(np.linspace(0, 20, n))
    return az, el, ctime


@pytest.fixture
def pointing():
    q = qpoint.QPoint(mean_aber=True, accuracy="low", fast_math=True)
    az, el, ctime = scan(NPOINT)
    return q, az, el, ctime


class TestPointing:
    def test_azel2bore(self, bench, pointing):
        q, az, el, ctime = pointing
        out = bench(
            "azel2bore",
            lambda: q.azel2bore(az, el, None, None, LON, LAT, ctime),
            samples=NPOINT,
        )
        assert np.asarray(out).shape == (NPOINT, 4)

    def test_bore2radec(self, bench, pointing):
        q, az, el, ctime = pointing
        q_bore = q.azel2bore(az, el, None, None, LON, LAT, ctime)
        q_off = q.det_offset(1.0, 2.0, 30.0)
        out = bench(
            "bore2radec",
            lambda: q.bore2radec(q_off, ctime, q_bore),
            samples=NPOINT,
        )
        assert np.asarray(out[0]).shape == (NPOINT,)

    def test_bore2pix(self, bench, pointing):
        q, az, el, ctime = pointing
        q_bore = q.azel2bore(az, el, None, None, LON, LAT, ctime)
        q_off = q.det_offset(1.0, 2.0, 30.0)
        out = bench(
            "bore2pix",
            lambda: q.bore2pix(q_off, ctime, q_bore, nside=256),
            samples=NPOINT,
        )
        assert np.asarray(out[0]).shape == (NPOINT,)

    def test_azel2radec(self, bench, pointing):
        q, az, el, ctime = pointing
        out = bench(
            "azel2radec",
            lambda: q.azel2radec(1.0, 2.0, 30.0, az, el, None, None, LON, LAT, ctime),
            samples=NPOINT,
        )
        assert np.asarray(out[0]).shape == (NPOINT,)

    def test_radec2azel(self, bench, pointing):
        q, az, el, ctime = pointing
        ra, dec, _, _ = q.azel2radec(
            1.0, 2.0, 30.0, az, el, None, None, LON, LAT, ctime
        )
        pa = np.zeros(NPOINT)
        out = bench(
            "radec2azel",
            lambda: q.radec2azel(ra, dec, pa, LON, LAT, ctime),
            samples=NPOINT,
        )
        assert np.asarray(out[0]).shape == (NPOINT,)

    @pytest.mark.parametrize("fast_math", [False, True])
    def test_fast_math(self, bench, fast_math):
        """
        What the polynomial trig is worth, measured where it dominates.

        bore2radec is mostly the quaternion-to-angle conversion, so this
        is the path fast_math changes: about 1.4x here. azel2bore barely
        moves, its cost being the correction chain rather than the trig,
        which is worth knowing before reaching for the option.
        """
        q = qpoint.QPoint(mean_aber=True, accuracy="low", fast_math=fast_math)
        az, el, ctime = scan(NPOINT)
        q_bore = q.azel2bore(az, el, None, None, LON, LAT, ctime)
        q_off = q.det_offset(1.0, 2.0, 30.0)
        out = bench(
            f"bore2radec fast_math={str(fast_math).lower()}",
            lambda: q.bore2radec(q_off, ctime, q_bore),
            samples=NPOINT,
        )
        assert np.asarray(out[0]).shape == (NPOINT,)


def mapper(**kwargs):
    qm = qpoint.QMap(nside=NSIDE, pol=True, mean_aber=True, fast_math=True, **kwargs)
    az, el, ctime = scan(NMAP)
    qm.init_point(qm.azel2bore(az, el, None, None, LON, LAT, ctime), ctime=ctime)
    pol = np.arange(NDET) * 180.0 / NDET
    q_off = np.atleast_2d(qm.det_offset(np.zeros(NDET), np.zeros(NDET), pol))
    return qm, q_off


class TestMapmaking:
    def test_from_tod(self, bench):
        qm, q_off = mapper()
        tod = np.random.default_rng(0).normal(size=(NDET, NMAP))
        out = bench(
            "from_tod (tod2map)",
            lambda: qm.from_tod(q_off, tod=tod.copy()),
            samples=NDET * NMAP,
        )
        assert np.any(np.asarray(out[1]))

    def test_to_tod(self, bench):
        qm, q_off = mapper()
        source = np.random.default_rng(0).normal(size=(3, 12 * NSIDE * NSIDE))
        qm.init_source(source, pol=True)
        out = bench(
            "to_tod (map2tod)",
            lambda: qm.to_tod(q_off),
            samples=NDET * NMAP,
        )
        assert np.asarray(out).shape == (NDET, NMAP)

    def test_solve_and_condition(self, bench):
        """
        The post-processing, which is where a full-sky destination map
        costs far more than the scan that filled it.
        """
        qm, q_off = mapper()
        vec, proj = qm.from_tod(
            q_off, tod=np.random.default_rng(0).normal(size=(NDET, NMAP))
        )
        vec, proj = np.asarray(vec), np.asarray(proj)
        npix = proj.shape[1]

        cond = bench(
            "proj_cond",
            lambda: qm.proj_cond(proj=proj.copy()),
            samples=npix,
            unit="pixel",
        )
        assert np.asarray(cond).shape == (npix,)

        solved = bench(
            "solve_map",
            lambda: qm.solve_map(vec=vec.copy(), proj=proj.copy()),
            samples=npix,
            unit="pixel",
        )
        assert np.asarray(solved).shape == (3, npix)


class TestInterpolation:
    def test_get_interp_val(self, bench):
        q = qpoint.QPoint(mean_aber=True)
        rng = np.random.default_rng(0)
        nside = 128
        m = rng.normal(size=12 * nside * nside)
        ra = rng.uniform(0, 360, NMAP)
        dec = rng.uniform(-80, -20, NMAP)
        out = bench(
            "get_interp_val",
            lambda: q.get_interp_val(m, ra, dec),
            samples=NMAP,
        )
        assert np.asarray(out).shape == (NMAP,)
