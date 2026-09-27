#pragma once

#include <cmath>
#include <cstddef>
#include <memory>
#include <vector>

namespace qp {

struct IersValues {
  double dut1;
  double x;
  double y;
};

// IERS Bulletin A polar motion and UT1 corrections.
//
// The table is immutable once set and held by shared_ptr, so copying a
// Pointing (which OpenMP does once per thread) is a refcount bump rather
// than the malloc + memcpy of the whole table that qp_copy_iers_bulletin_a
// performs.
class BulletinA {
 public:
  // Stored as float, as qp_bulletina_entry_t does. Widening to double here
  // would change the interpolated results.
  struct Entry {
    float x;
    float y;
    float dut1;
  };

  BulletinA() = default;

  void set(int mjd_min, int mjd_max, const double *dut1, const double *x,
           const double *y) {
    auto t = std::make_shared<Table>();
    t->mjd_min = mjd_min;
    t->mjd_max = mjd_max;
    const std::size_t n = static_cast<std::size_t>(mjd_max - mjd_min + 1);
    t->entries.resize(n);
    for (std::size_t k = 0; k < n; ++k)
      t->entries[k] = {static_cast<float>(x[k]), static_cast<float>(y[k]),
                       static_cast<float>(dut1[k])};
    table_ = std::move(t);
  }

  // Out of range returns zeros rather than throwing: every caller in the C
  // ignores the error return and proceeds with the zeroed values, so any
  // user who has not loaded a bulletin relies on this.
  IersValues interp(double mjd) const {
    if (!table_) return {0., 0., 0.};
    const Table &t = *table_;
    if (!(t.mjd_min <= mjd && mjd < t.mjd_max)) return {0., 0., 0.};

    double mjd_floor;
    const double r = std::modf(mjd, &mjd_floor);
    const std::size_t k =
        static_cast<std::size_t>(static_cast<int>(mjd_floor) - t.mjd_min);
    const Entry a = t.entries[k];
    const Entry b = t.entries[k + 1];

    // A jump of more than half a second across one day is a leap second,
    // not a real UT1 excursion; remove it before interpolating.
    double leap = b.dut1 - a.dut1;
    if (leap > 0.5)
      leap = 1.;
    else if (leap < -0.5)
      leap = -1.;
    else
      leap = 0;

    return {(1 - r) * a.dut1 + r * (b.dut1 - leap),
            (1 - r) * a.x + r * b.x,
            (1 - r) * a.y + r * b.y};
  }

 private:
  struct Table {
    std::vector<Entry> entries;
    int mjd_min = 0;
    int mjd_max = 0;
  };

  std::shared_ptr<const Table> table_;
};

}  // namespace qp
