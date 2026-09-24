"""
qpoint against astropy.

Everything else in this suite checks qpoint against itself: round trips,
invariants, and values that were correct when they were written down.
These check it against an independent implementation of the same
transforms, so a shared mistake has to be made twice to survive.

astropy is not a test dependency, so these skip where it is absent. They
also never reach the network: the IERS tables ship with astropy, and the
auto-download that would otherwise fetch newer ones is switched off for
the duration.

Each package gets its Earth orientation data through
`QPoint(update_iers=True)`, which reads astropy's own table and installs
it as a Bulletin A. So every side is working from the same numbers
by construction, and the comparison is of the transforms rather than of
what each believes the Earth was doing.
"""

import numpy as np
import pytest
import qpoint

pytest.importorskip("astropy")
qpoint2 = pytest.importorskip("qpoint2")

from astropy import units as u  # noqa: E402
from astropy.coordinates import (  # noqa: E402
    AltAz,
    EarthLocation,
    get_sun,
    ICRS,
    SkyCoord,
)
from astropy.time import Time  # noqa: E402
from astropy.utils import iers  # noqa: E402

LON, LAT = 165.7, -77.6  # McMurdo Station
N = 7
CTIME = 1418662800.0 + np.linspace(0, 6 * 3600, N)
AZ = np.linspace(10, 340, N)
EL = np.linspace(30, 80, N)

# Weather for the refraction comparison. 150 GHz is 2 mm, which is what
# astropy wants as a wavelength.
PRESSURE, TEMPERATURE, HUMIDITY, FREQUENCY = 1000.0, 0.0, 0.0, 150.0
OBSWL = (299792458.0 / (FREQUENCY * 1e9)) * u.m

# Room above what the comparisons actually show, which is 0.8 mas once
# make_qpoint has switched everything on.
#
# It used to be 14 mas, and that was a single term: the Sun's
# gravitational light deflection, which astropy removes to recover the
# catalogue position. qpoint can model it now (rate_defl), so what is
# left is the genuine disagreement between two implementations.
#
# A tenth of an arcsecond would be a hundred times looser than that.
# Five mas clears the worst comparison by about six and is tight enough
# that a term going missing shows up at once -- leaving deflection at its
# default of never puts azel2radec back at 14 mas and radec2azel at 43,
# both of which TestLightDeflection asserts.
TOL_ARCSEC = 0.005


# Both packages are checked against astropy, not only against each other.
IMPLS = [qpoint, qpoint2]


@pytest.fixture(params=IMPLS, ids=lambda m: m.__name__)
def mod(request):
    return request.param


@pytest.fixture(scope="module", autouse=True)
def no_download():
    """
    Keep astropy off the network for the duration of these tests.

    Autouse because update_iers reads the table during construction, so
    every QPoint built here has to be built inside it, not just the
    astropy references.
    """
    with iers.conf.set_temp("auto_download", False):
        yield


@pytest.fixture(scope="module")
def location():
    return EarthLocation(lon=LON * u.deg, lat=LAT * u.deg, height=0 * u.m)


@pytest.fixture(scope="module")
def obstime(no_download):
    return Time(CTIME, format="unix", scale="utc")


def make_qpoint(mod, inverse=True, **kwargs):
    """
    A QPoint using the same Earth orientation data as astropy.

    `update_iers` does the loading: it calls update_bulletin_a, which
    reads astropy's IERS table and installs it as a Bulletin A.
    Without it the two disagree by the size of those terms, which is
    seconds of arc -- real, but a comparison of Earth orientation data
    rather than of the transforms built on it.

    The dut1 and wobble rates default to never, so the table is
    otherwise read and then ignored, and the inverse transforms are
    switched separately from the forward ones.

    Solar light deflection defaults to never for a different reason --
    it is a speed trade rather than missing data -- but it has to be on
    for the same purpose, since astropy applies it and it is 14 mas.
    """
    rates = dict(rate_dut1="always", rate_wobble="always", rate_defl=100)
    if inverse:
        rates.update(
            rate_dut1_inv="always", rate_wobble_inv="always", rate_defl_inv=100
        )
    # caller wins, so a test can put one of these back to its default
    rates.update(kwargs)
    return mod.QPoint(update_iers=True, **rates)


def separation(ra, dec, reference):
    got = SkyCoord(ra=np.asarray(ra) * u.deg, dec=np.asarray(dec) * u.deg, frame="icrs")
    return got.separation(reference).to_value(u.arcsec)


class TestAzel2Radec:
    def test_matches_astropy(self, mod, location, obstime):
        reference = AltAz(
            az=AZ * u.deg,
            alt=EL * u.deg,
            obstime=obstime,
            location=location,
            pressure=0 * u.hPa,
        ).transform_to(ICRS())

        q = make_qpoint(mod, pressure=0)
        ra, dec = q.azel2radec(
            0.0, 0.0, 0.0, AZ, EL, None, None, LON, LAT, CTIME, return_pa=True
        )[:2]
        assert separation(ra, dec, reference).max() < TOL_ARCSEC

    def test_matches_astropy_through_refraction(self, mod, location, obstime):
        """
        Both bend the incoming ray for the same atmosphere. Elevations
        here stay above 30 degrees, where the two refraction models have
        no room to disagree about the horizon.
        """
        reference = AltAz(
            az=AZ * u.deg,
            alt=EL * u.deg,
            obstime=obstime,
            location=location,
            pressure=PRESSURE * u.hPa,
            temperature=TEMPERATURE * u.deg_C,
            relative_humidity=HUMIDITY,
            obswl=OBSWL,
        ).transform_to(ICRS())

        q = make_qpoint(
            mod,
            pressure=PRESSURE,
            temperature=TEMPERATURE,
            humidity=HUMIDITY,
            frequency=FREQUENCY,
            rate_ref="always",
        )
        ra, dec = q.azel2radec(
            0.0, 0.0, 0.0, AZ, EL, None, None, LON, LAT, CTIME, return_pa=True
        )[:2]
        assert separation(ra, dec, reference).max() < TOL_ARCSEC


class TestRadec2Azel:
    SKY = dict(ra=np.linspace(20, 300, N), dec=np.linspace(-80, -20, N))

    def reference(self, location, obstime):
        sky = SkyCoord(
            ra=self.SKY["ra"] * u.deg, dec=self.SKY["dec"] * u.deg, frame="icrs"
        )
        return sky.transform_to(
            AltAz(obstime=obstime, location=location, pressure=0 * u.hPa)
        )

    def horizon(self, az, el, location, obstime):
        return SkyCoord(
            az=np.asarray(az) * u.deg,
            alt=np.asarray(el) * u.deg,
            frame=AltAz(obstime=obstime, location=location),
        )

    def test_matches_astropy(self, mod, location, obstime):
        q = make_qpoint(mod, pressure=0)
        az, el = q.radec2azel(
            self.SKY["ra"], self.SKY["dec"], np.zeros(N), LON, LAT, CTIME
        )[:2]
        got = self.horizon(az, el, location, obstime)
        assert (
            got.separation(self.reference(location, obstime)).to_value(u.arcsec).max()
            < TOL_ARCSEC
        )

    def test_the_inverse_rates_are_switched_separately(self, mod, location, obstime):
        """
        Setting rate_dut1 and rate_wobble alone leaves the inverse
        transform without polar motion, which is a quarter of an
        arcsecond and easy to mistake for noise.
        """
        q = make_qpoint(mod, inverse=False, pressure=0)
        az, el = q.radec2azel(
            self.SKY["ra"], self.SKY["dec"], np.zeros(N), LON, LAT, CTIME
        )[:2]
        got = self.horizon(az, el, location, obstime)
        worse = (
            got.separation(self.reference(location, obstime)).to_value(u.arcsec).max()
        )
        assert worse > TOL_ARCSEC
        assert worse < 1.0


class TestWithoutEarthOrientationData:
    """
    What the Earth orientation terms are worth, which is also what a
    caller who never loads a bulletin gives up: arcseconds, not degrees.
    """

    def test_the_disagreement_is_seconds_of_arc(self, mod, location, obstime):
        reference = AltAz(
            az=AZ * u.deg,
            alt=EL * u.deg,
            obstime=obstime,
            location=location,
            pressure=0 * u.hPa,
        ).transform_to(ICRS())

        q = mod.QPoint(pressure=0)
        ra, dec = q.azel2radec(
            0.0, 0.0, 0.0, AZ, EL, None, None, LON, LAT, CTIME, return_pa=True
        )[:2]
        worst = separation(ra, dec, reference).max()
        assert worst > TOL_ARCSEC
        assert worst < 60.0


class TestGalacticRotation:
    """
    The galactic rotation against an outside implementation, which is the
    only check it has: test_qpoint round trips radec2gal through
    gal2radec, and a round trip survives a wrong rotation -- turning by
    the wrong angle and back still returns where it started.

    The reference is erfa.icrs2g and not astropy's Galactic frame,
    because those two disagree by 25 mas and astropy is the one departing
    from the canonical definition. The IAU 1958 system was defined
    against FK4 B1950, so using it from a modern frame means accounting
    for the E-terms of aberration, the B1950 equinox, and the frame bias
    between ICRS and J2000 mean place. Hipparcos supplies a rotation
    straight from ICRS with all of that folded in, as three angles its
    catalogue calls exact for canonical purposes; qpoint is built from
    those angles and so is ERFA. astropy instead reads the same three
    numbers as FK5 and goes ICRS -> FK5 -> Galactic, picking the frame
    bias back up.

    So the agreement with ERFA is exact and asserted as such, and the
    astropy gap is pinned separately as the convention difference it is.
    Neither is a tolerance absorbing an unexplained residual, and the
    second is there so nobody "corrects" qpoint into matching astropy.
    """

    # Agreement with ERFA is exact; this is room for the last bits of a
    # different route to the same rotation, not for a modelling gap.
    TOL_UAS = 1.0
    # What astropy's FK5-mediated definition costs, which varies over the
    # sky because the frame bias is a rotation.
    BIAS_MAS = (0.5, 30.0)

    SKY = dict(ra=np.linspace(0.0, 350.0, N), dec=np.linspace(-80.0, 80.0, N))

    def erfa_reference(self):
        """The canonical Hipparcos ICRS -> Galactic rotation."""
        erfa = pytest.importorskip("erfa")
        lon, lat = erfa.icrs2g(np.deg2rad(self.SKY["ra"]), np.deg2rad(self.SKY["dec"]))
        return SkyCoord(
            l=np.rad2deg(lon) * u.deg, b=np.rad2deg(lat) * u.deg, frame="galactic"
        )

    def test_radec2gal_matches_erfa(self, mod):
        # radec2gal is in-place by default, so these are copies rather
        # than the arrays the class holds.
        lon, lat = mod.QPoint().radec2gal(
            self.SKY["ra"].copy(), self.SKY["dec"].copy()
        )[:2]
        got = SkyCoord(
            l=np.asarray(lon) * u.deg, b=np.asarray(lat) * u.deg, frame="galactic"
        )
        assert (
            got.separation(self.erfa_reference()).to_value(u.uas).max() < self.TOL_UAS
        )

    def test_gal2radec_matches_erfa(self, mod):
        """The inverse, taken back to where it started."""
        gal = self.erfa_reference()
        ra, dec = mod.QPoint().gal2radec(
            gal.l.to_value(u.deg).copy(), gal.b.to_value(u.deg).copy()
        )[:2]
        want = SkyCoord(
            ra=self.SKY["ra"] * u.deg, dec=self.SKY["dec"] * u.deg, frame="icrs"
        )
        got = SkyCoord(
            ra=np.asarray(ra) * u.deg, dec=np.asarray(dec) * u.deg, frame="icrs"
        )
        assert got.separation(want).to_value(u.uas).max() < self.TOL_UAS

    def test_astropy_differs_by_the_frame_bias(self, mod):
        """
        Pinned so the 25 mas is on the record as a convention difference
        rather than found later and mistaken for an error in qpoint.
        """
        lon, lat = mod.QPoint().radec2gal(
            self.SKY["ra"].copy(), self.SKY["dec"].copy()
        )[:2]
        got = SkyCoord(
            l=np.asarray(lon) * u.deg, b=np.asarray(lat) * u.deg, frame="galactic"
        )
        astropy_gal = SkyCoord(
            ra=self.SKY["ra"] * u.deg, dec=self.SKY["dec"] * u.deg, frame="icrs"
        ).galactic
        gap = got.separation(astropy_gal).to_value(u.mas)
        lo, hi = self.BIAS_MAS
        assert lo < gap.max() < hi


class TestSiderealTime:
    """
    gmst and lmst against astropy, the other transform with no outside
    reference of its own.

    These need the Earth orientation data for the same reason the
    coordinate transforms do: sidereal time is Earth rotation, so
    dropping dut1 costs the whole of that term, which is seconds of arc.
    """

    def astropy_gmst(self):
        return (
            Time(CTIME, format="unix", scale="utc")
            .sidereal_time("mean", "greenwich")
            .to_value(u.hourangle)
        )

    @staticmethod
    def arcsec(got, ref):
        """Separation in arcsec of Earth rotation, wrapped at 24 hours."""
        return np.abs((np.asarray(got) - ref + 12.0) % 24.0 - 12.0) * 15.0 * 3600.0

    def test_gmst_matches_astropy(self, mod):
        q = make_qpoint(mod, inverse=False)
        assert self.arcsec(q.gmst(CTIME), self.astropy_gmst()).max() < TOL_ARCSEC

    def test_lmst_matches_astropy(self, mod):
        """
        Longitude enters here and nowhere else, so this is what pins its
        sign: a flipped one leaves gmst untouched.
        """
        ref = (
            Time(CTIME, format="unix", scale="utc")
            .sidereal_time("mean", longitude=LON * u.deg)
            .to_value(u.hourangle)
        )
        q = make_qpoint(mod, inverse=False)
        assert self.arcsec(q.lmst(CTIME, LON), ref).max() < TOL_ARCSEC

    def test_without_the_bulletin_it_is_off_by_dut1(self, mod):
        """
        The same point TestWithoutEarthOrientationData makes for the
        coordinate transforms: what the Earth orientation data is worth
        here is seconds of arc, and it is exactly the dut1 term.
        """
        worst = self.arcsec(mod.QPoint().gmst(CTIME), self.astropy_gmst()).max()
        assert worst > TOL_ARCSEC
        assert worst < 0.9 * 15.0  # leap seconds bound |dut1| < 0.9 s


class TestLightDeflection:
    """
    The sun bends the incoming ray, and both packages model it.

    This is what the tolerance above is really resting on. Before the term
    was added the two implementations disagreed by 14 mas going out and 43
    mas coming back, and all of it was this -- so switching it off has to
    reproduce exactly that. Showing the residual is small proves nothing
    on its own; showing it is small *because* of a term with the right
    size and the right functional form is the point.
    """

    def radec(self, mod, rate):
        q = make_qpoint(mod, rate_defl=rate, rate_defl_inv=rate, pressure=0)
        return q.azel2radec(
            0.0, 0.0, 0.0, AZ, EL, None, None, LON, LAT, CTIME, return_pa=True
        )[:2]

    def reference(self, location, obstime):
        return AltAz(
            az=AZ * u.deg,
            alt=EL * u.deg,
            obstime=obstime,
            location=location,
            pressure=0 * u.hPa,
        ).transform_to(ICRS())

    def test_off_by_default(self, mod):
        """
        It is opt-in: ~20 mas against about 10% of azel2bore. Callers who
        want it say so, and both packages agree on that.
        """
        assert mod.QPoint().get("rate_defl") == "never"
        assert mod.QPoint().get("rate_defl_inv") == "never"

    def test_forward_needs_it_to_agree_with_astropy(self, mod, location, obstime):
        ref = self.reference(location, obstime)
        with_it = separation(*self.radec(mod, 100), ref).max()
        without = separation(*self.radec(mod, "never"), ref).max()
        assert with_it < TOL_ARCSEC
        assert without > 2 * TOL_ARCSEC, "turning it off must actually matter"
        assert 0.010 < without < 0.020, "and by the ~14 mas it is worth here"

    def test_inverse_needs_it_too(self, mod, location, obstime):
        """
        The inverse carried the larger error of the two, 43 mas, because
        there the term is applied rather than removed.
        """
        sky = SkyCoord(
            ra=np.linspace(20, 300, N) * u.deg,
            dec=np.linspace(-80, -20, N) * u.deg,
            frame="icrs",
        )
        ref = sky.transform_to(
            AltAz(obstime=obstime, location=location, pressure=0 * u.hPa)
        )
        got = []
        for rate in (100, "never"):
            q = make_qpoint(mod, rate_defl=rate, rate_defl_inv=rate, pressure=0)
            az, el = q.radec2azel(
                sky.ra.deg, sky.dec.deg, np.zeros(N), LON, LAT, CTIME
            )[:2]
            here = SkyCoord(
                az=np.asarray(az) * u.deg,
                alt=np.asarray(el) * u.deg,
                frame=AltAz(obstime=obstime, location=location),
            )
            got.append(here.separation(ref).to_value(u.arcsec).max())
        assert got[0] < TOL_ARCSEC
        assert 0.035 < got[1] < 0.050

    def test_it_matches_erfa_exactly(self, mod, obstime):
        """
        The size and direction of the term, against ERFA's own routine.

        This is the reference rather than the textbook
        4.07 mas / tan(elongation / 2), because that needs the elongation
        and `get_sun` hands back a GCRS position *with distance*: compare
        it to an ICRS coordinate and astropy dutifully applies parallax,
        which for a body one au away swings the direction by ~80 degrees.
        eraLdsun takes the sun-to-observer vector directly and sidesteps
        all of it.
        """
        erfa = pytest.importorskip("erfa")
        on, off = self.radec(mod, 100), self.radec(mod, "never")
        jd1, jd2 = obstime.tt.jd1, obstime.tt.jd2

        def unit(ra, dec):
            ra, dec = np.deg2rad(np.asarray(ra, float)), np.deg2rad(
                np.asarray(dec, float)
            )
            return np.array(
                [np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)]
            )

        got = unit(*on)
        base = unit(*off)
        for i in range(N):
            pvh, _ = erfa.epv00(jd1[i], jd2[i])
            ph = np.asarray(pvh["p"])
            em = np.linalg.norm(ph)
            # eraLdsun bends a catalogue direction into the observed one.
            # `got` is the corrected (catalogue) position and `base` the
            # uncorrected one, so bending `got` has to reproduce `base`.
            # That pins the direction and the sign, not merely the size.
            want = erfa.ldsun(got[:, i], ph / em, em)
            # 2 asin(|a-b|/2): stable where arccos of a dot product is not
            sep = (
                np.rad2deg(2 * np.arcsin(np.linalg.norm(base[:, i] - want) / 2)) * 3.6e6
            )
            assert sep < 0.01, "sample {}: {} mas from eraLdsun".format(i, sep)
