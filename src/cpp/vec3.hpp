#pragma once

#include <cmath>
#include <type_traits>

namespace qp {

// Layout-compatible with vec3_t (double[3]) so it can be loaded from and
// stored to numpy buffers directly. Arithmetic mirrors src/vec3.h term for
// term; do not reassociate.
struct Vec3 {
  double v[3];

  constexpr double &operator[](int i) { return v[i]; }
  constexpr double operator[](int i) const { return v[i]; }

  constexpr double *data() { return v; }

  double norm2() const { return v[0] * v[0] + v[1] * v[1] + v[2] * v[2]; }
  double norm() const { return std::sqrt(norm2()); }
};

inline double dot(const Vec3 &a, const Vec3 &b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

inline Vec3 cross(const Vec3 &a, const Vec3 &b) {
  return {{a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
           a[0] * b[1] - a[1] * b[0]}};
}

static_assert(std::is_trivially_copyable_v<Vec3>);
static_assert(std::is_standard_layout_v<Vec3>);
static_assert(sizeof(Vec3) == 3 * sizeof(double));

}  // namespace qp
