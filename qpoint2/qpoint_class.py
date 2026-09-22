import numpy as np

from . import _libqpoint2 as lib

RATE_PARAMS = [
    "daber",
    "lonlat",
    "wobble",
    "dut1",
    "erot",
    "npb",
    "aaber",
    "ref",
    "daber_inv",
    "lonlat_inv",
    "wobble_inv",
    "dut1_inv",
    "erot_inv",
    "npb_inv",
    "aaber_inv",
    "ref_inv",
]

OPTION_PARAMS = [
    "accuracy",
    "mean_aber",
    "fast_aber",
    "fast_math",
    "polconv",
    "pix_order",
    "interp_pix",
    "fast_pix",
    "error_missing",
    "nan_missing",
    "interp_missing",
    "num_threads",
]

WEATHER_PARAMS = ["temperature", "pressure", "humidity", "frequency"]

DOUBLE_PARAMS = ["ref_delta", "dut1"]

__all__ = ["QPoint"]


# ---------------------------------------------------------------------------
# Array preparation helpers
#
# This is the single coercion point. The C++ bindings validate strictly and
# never copy, so anything they are handed must already be 1-D, C-contiguous,
# aligned float64 -- np.broadcast_arrays alone does not give that, since it
# returns 0-stride views.
# ---------------------------------------------------------------------------


def _bc(*args):
    """Broadcast args (None→0) to a common shape as C-contiguous float64."""
    arrays = [np.asarray(a if a is not None else 0.0, dtype=float) for a in args]
    return [
        np.require(np.atleast_1d(a), float, ["C", "A"])
        for a in np.broadcast_arrays(*arrays)
    ]


class QPoint:
    """Quaternion-based telescope pointing."""

    # ---- Quaternion construction ----

    def det_offset(self, delta_az, delta_el, delta_psi):
        """
        Quaternion for the requested detector centroid offset from boresight.

        Arguments
        ---------
        delta_az : array_like
            Azimuthal centroid offset of the detector in degrees
        delta_el : array_like
            Elevation centroid offset of the detector in degrees
        delta_psi : array_like
            Polarization offset of the detector from vertical in degrees

        Returns
        -------
        q : array_like
            Detector centroid offset quaternion for each detector
        """
        return lib.det_offset(*_bc(delta_az, delta_el, delta_psi))

    def hwp_quat(self, theta):
        """
        Quaternion for rotation by 2*theta (physical HWP angle).

        Arguments
        ---------
        theta : array_like
            HWP physical angle in degrees

        Returns
        -------
        q : array_like
            Quaternion for each hwp angle
        """
        return lib.hwp_quat(*_bc(theta))
