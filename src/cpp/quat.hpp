#pragma once

#include <cmath>
#include <type_traits>

#include "vec3.hpp"

namespace qp {

// Layout-compatible with quat_t (double[4]), so an (n,4) numpy buffer can be
// viewed as a Quat array with no copy.
//
// Indexed access rather than named w/x/y/z members: the ported code in
// quat2radec and quat2pix is dense with literal indices, and keeping the same
// notation makes it diffable against the C.
//
// Every expression below is transcribed term for term from src/quaternion.c.
// Do not reassociate -- bit-exact agreement with the C build depends on it.
struct Quat {
  double v[4];

  constexpr double &operator[](int i) { return v[i]; }
  constexpr double operator[](int i) const { return v[i]; }

  static constexpr Quat identity() { return {{1., 0., 0., 0.}}; }

  static Quat load(const double *p) { return {{p[0], p[1], p[2], p[3]}}; }
  void store(double *p) const {
    p[0] = v[0];
    p[1] = v[1];
    p[2] = v[2];
    p[3] = v[3];
  }

  static Quat r1(double angle) {
    const double a2 = 0.5 * angle;
    return {{std::cos(a2), std::sin(a2), 0., 0.}};
  }
  static Quat r2(double angle) {
    const double a2 = 0.5 * angle;
    return {{std::cos(a2), 0., std::sin(a2), 0.}};
  }
  static Quat r3(double angle) {
    const double a2 = 0.5 * angle;
    return {{std::cos(a2), 0., 0., std::sin(a2)}};
  }

  // Quaternion_rot returns 1 on a zero-norm axis and leaves q untouched, and
  // its only caller (qp_aberration) ignores that return -- so the C reads an
  // uninitialized quaternion on this path. Return identity instead.
  static Quat rot(double angle, const Vec3 &axis) {
    const double a2 = 0.5 * angle;
    const double s = std::sin(a2);
    const double norm = std::sqrt(axis[0] * axis[0] + axis[1] * axis[1] +
                                  axis[2] * axis[2]);
    if (norm <= 0.) return identity();
    return {{std::cos(a2), s * axis[0] / norm, s * axis[1] / norm,
             s * axis[2] / norm}};
  }

  double norm2() const {
    return v[0] * v[0] + v[1] * v[1] + v[2] * v[2] + v[3] * v[3];
  }
  Quat inv() const {
    const double n2 = norm2();
    return {{v[0] / n2, -v[1] / n2, -v[2] / n2, -v[3] / n2}};
  }

  void unit() {
    const double invnorm = 1. / std::sqrt(norm2());
    v[0] *= invnorm;
    v[1] *= invnorm;
    v[2] *= invnorm;
    v[3] *= invnorm;
  }

  // Columns of the rotation matrix. Caller must normalize first.
  Vec3 col3() const {
    const double a2 = v[0] * v[0], b2 = v[1] * v[1], c2 = v[2] * v[2],
                 d2 = v[3] * v[3];
    return {{2. * (v[1] * v[3] + v[0] * v[2]),
             2. * (v[2] * v[3] - v[0] * v[1]), a2 - b2 - c2 + d2}};
  }
};

inline Quat operator*(const Quat &a, const Quat &b) {
  return {{a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3],
           a[0] * b[1] + a[1] * b[0] + a[2] * b[3] - a[3] * b[2],
           a[0] * b[2] - a[1] * b[3] + a[2] * b[0] + a[3] * b[1],
           a[0] * b[3] + a[1] * b[2] - a[2] * b[1] + a[3] * b[0]}};
}

inline void mul_left(const Quat &a, Quat &q) { q = a * q; }
inline void mul_right(Quat &q, const Quat &a) { q = q * a; }

// Axis-rotation multiplies, kept as the C's hand-specialized forms rather
// than q = Quat::rN(angle) * q. They drop the terms that are identically
// zero, so the general product would introduce extra 0*x additions. The
// cos/sin are the same two values Quat::rN would have put in its lanes.
inline void r1_mul(double angle, Quat &q) {
  const double a2 = 0.5 * angle;
  const double c = std::cos(a2), s = std::sin(a2);
  const Quat b = q;
  q[0] = c * b[0] - s * b[1];
  q[1] = c * b[1] + s * b[0];
  q[2] = c * b[2] - s * b[3];
  q[3] = c * b[3] + s * b[2];
}

inline void r2_mul(double angle, Quat &q) {
  const double a2 = 0.5 * angle;
  const double c = std::cos(a2), s = std::sin(a2);
  const Quat b = q;
  q[0] = c * b[0] - s * b[2];
  q[1] = c * b[1] + s * b[3];
  q[2] = c * b[2] + s * b[0];
  q[3] = c * b[3] - s * b[1];
}

inline void r3_mul(double angle, Quat &q) {
  const double a2 = 0.5 * angle;
  const double c = std::cos(a2), s = std::sin(a2);
  const Quat b = q;
  q[0] = c * b[0] - s * b[3];
  q[1] = c * b[1] - s * b[2];
  q[2] = c * b[2] + s * b[1];
  q[3] = c * b[3] + s * b[0];
}

static_assert(std::is_trivially_copyable_v<Quat>);
static_assert(std::is_standard_layout_v<Quat>);
static_assert(sizeof(Quat) == 4 * sizeof(double));

}  // namespace qp
