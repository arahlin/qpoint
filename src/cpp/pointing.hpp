#pragma once

#include "quat.hpp"
#include "units.hpp"

namespace qp {

// Stateless quaternion constructors. These take no Pointing state, so they
// stay free functions -- qp_det_offset and qp_hwp_quat in src/qpoint.c.

// Detector offset quaternion from boresight-relative angles, in degrees.
inline Quat det_offset(double delta_az, double delta_el, double delta_psi) {
  Quat q = Quat::r3(-deg2rad(delta_psi));
  r2_mul(deg2rad(delta_el), q);
  r1_mul(-deg2rad(delta_az), q);
  return q;
}

// Waveplate quaternion from the physical HWP angle, in degrees. The
// polarization angle rotates by twice the physical angle.
inline Quat hwp_quat(double ang) { return Quat::r3(-2. * deg2rad(ang)); }

}  // namespace qp
