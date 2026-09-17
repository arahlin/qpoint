"""Tests for qpoint.tools module."""

import numpy as np
import pytest
from qpoint.tools import refraction


class TestRefraction:
    def test_scalar_input(self):
        delta = refraction(45.0, 20.0, 1013.25, 0.5, 150.0)
        assert np.isscalar(delta) or delta.ndim == 0

    def test_positive_at_low_elevation(self):
        """Refraction lifts objects upward, so delta > 0 at low elevation."""
        delta = refraction(10.0, 20.0, 1013.25, 0.5, 150.0)
        assert float(delta) > 0

    def test_larger_at_lower_elevation(self):
        """Refraction increases toward the horizon."""
        delta_high = refraction(60.0, 20.0, 1013.25, 0.5, 150.0)
        delta_low = refraction(10.0, 20.0, 1013.25, 0.5, 150.0)
        assert delta_low > delta_high

    def test_physical_range_at_45deg(self):
        """At 45 deg elevation, refraction should be ~1 arcmin ~ 0.017 deg."""
        delta = refraction(45.0, 20.0, 1013.25, 0.0, 150.0)
        assert 0 < float(delta) < 0.1  # less than 6 arcmin

    def test_zero_pressure_small_refraction(self):
        """With zero pressure (vacuum), refraction should be very small."""
        delta = refraction(45.0, 20.0, 0.0, 0.0, 150.0)
        assert float(delta) < float(refraction(45.0, 20.0, 1013.25, 0.5, 150.0))

    def test_vector_input(self):
        el = np.array([10.0, 20.0, 30.0, 45.0, 60.0, 80.0])
        temp = 20.0 * np.ones(6)
        press = 1013.25 * np.ones(6)
        hum = 0.5 * np.ones(6)
        freq = 150.0 * np.ones(6)
        delta = refraction(el, temp, press, hum, freq)
        assert delta.shape == (6,)
        assert np.all(delta > 0)

    def test_monotone_with_elevation(self):
        """Refraction is monotonically decreasing with elevation."""
        el = np.linspace(5, 85, 20)
        delta = refraction(el, 20.0, 1013.25, 0.5, 150.0)
        assert np.all(np.diff(delta) < 0)

    def test_frequency_argument_accepted(self):
        """refraction accepts a frequency argument without error."""
        delta = refraction(45.0, 20.0, 1013.25, 0.5, 150.0)
        assert np.isfinite(float(delta))
