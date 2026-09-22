#pragma once

#include <cstddef>
#include <cstdint>
#include <unordered_map>

#include "pixinfo.hpp"
#include "quat.hpp"
#include "span.hpp"

namespace qp {

// Maps sky pixel numbers onto their index in a partial map.
//
// Replaces the hand-rolled open-hash table in qp_pixhash.c. That version
// allocated one bucket per pixel and grew each bucket with realloc on every
// collision; std::unordered_map does the same job.
//
// Duplicate pixels keep their first index, matching qp_repixelize, which
// scans a bucket in insertion order and returns the first match.
class PixHash {
 public:
  PixHash() = default;

  explicit PixHash(Span<const long> pix) {
    map_.reserve(pix.size());
    for (std::size_t i = 0; i < pix.size(); ++i)
      map_.try_emplace(pix[i], static_cast<long>(i));
  }

  // Index into the partial map, or -1 if the pixel is not in it.
  long repixelize(long pix) const {
    const auto it = map_.find(pix);
    return (it == map_.end()) ? -1 : it->second;
  }


 private:
  std::unordered_map<long, long> map_;
};

enum class VecMode {
  None = 0,
  Temp,    // T
  Pol,     // T, Q, U
  VPol,    // T, Q, U, V
  D1,      // T + first derivatives
  D1Pol,   // T, Q, U + first derivatives
  D2,      // T + second derivatives
  D2Pol,   // T, Q, U + second derivatives
};

enum class ProjMode {
  None = 0,
  Temp,  // hits
  Pol,   // 3x3 upper triangle
  VPol,  // 4x4 upper triangle
};

std::size_t num_vec(VecMode mode);
std::size_t num_proj(ProjMode mode);

// The inverse: the mode a row count implies, or None if no mode has that
// many rows. pol/vpol disambiguate the counts two modes share.
VecMode vec_mode_for(std::size_t nrow, bool pol, bool vpol);
ProjMode proj_mode_for(std::size_t nrow);

// An N-component map has N vec rows and N*(N+1)/2 proj rows. These convert
// between the two counts, returning 0 when there is no such pairing.
std::size_t num_proj_for_vec(std::size_t nvec);
std::size_t num_vec_for_proj(std::size_t nproj);

// Ordering comparisons: the accumulation kernels branch on "at least
// polarized", mirroring the C's integer enum comparisons.
inline bool at_least(VecMode m, VecMode ref) {
  return static_cast<int>(m) >= static_cast<int>(ref);
}
inline bool at_least(ProjMode m, ProjMode ref) {
  return static_cast<int>(m) >= static_cast<int>(ref);
}

// Whether a vec mode actually carries Q and U -- which is not
// at_least(m, Pol), because the enum above is ordered as the C's is so
// the kernels can compare with >=, and that ordering puts the
// *unpolarized* derivative modes D1 and D2 above Pol. Only these four
// are polarized. Use at_least for a kernel branch that mirrors the C,
// and this wherever the question is what the map contains.
inline bool is_polarized(VecMode m) {
  return m == VecMode::Pol || m == VecMode::VPol || m == VecMode::D1Pol ||
         m == VecMode::D2Pol;
}

struct Mueller {
  double v[4];
  constexpr double operator[](int i) const { return v[i]; }
};

// A detector. Every array is a non-owning view of a numpy buffer; an empty
// span means the array was not supplied, which replaces the C's *_init
// flags. tod is writable because map2tod accumulates into it.
//
// No defaults: DetArrWrap::setup is the only producer and assigns every
// field, so a default here would only ever be wrong in secret -- the
// mueller one was, claiming {1,1,1,1} against the real {1,1,0,1}.
struct Det {
  Quat q_off;
  double weight;
  double gain;
  Mueller mueller;

  Span<double> tod;
  Span<const std::uint8_t> flag;
  Span<const double> weights;
};

// Boresight pointing shared by every detector.
struct PointData {
  Span<const Quat> q_bore;
  Span<const double> ctime;
  Span<const Quat> q_hwp;

  std::size_t n() const { return q_bore.size(); }
  double time(std::size_t i) const { return ctime ? ctime[i] : 0.; }
};

// A map, stored as a contiguous (nrow, npix) block.
//
// The C carries both a 1-D buffer and a 2-D table of row pointers, kept in
// sync by qp_reshape_map. A row is just an offset, so the table and the
// reshape step both go away -- along with the bug class where num_vec was
// computed from the element count rather than the row count, letting a
// V-component write run off the end of vec and into proj.
struct Map {
  double *vec = nullptr;
  double *proj = nullptr;
  std::size_t nside = 0;
  std::size_t npix = 0;
  VecMode vec_mode = VecMode::None;
  ProjMode proj_mode = ProjMode::None;

  // A partial map is exactly one with a pixel hash, so there is no
  // separate flag to fall out of step with it.
  const PixHash *pixhash = nullptr;
  const PixInfo *pixinfo = nullptr;

  double *vrow(std::size_t i) const { return vec + i * npix; }
  double *prow(std::size_t i) const { return proj + i * npix; }

  bool partial() const { return pixhash != nullptr; }

  bool has_vec() const { return vec != nullptr && vec_mode != VecMode::None; }
  bool has_proj() const {
    return proj != nullptr && proj_mode != ProjMode::None;
  }
};

class Pointing;

// Accumulate one detector's timestream into the map.
void tod2map1(Pointing &mem, const Det &det, const PointData &pnt, Map &map);

// Accumulate a detector pair as a difference, for a polarization-differenced
// map. Detector i is paired with detector i + ndet/2.
void tod2map1_diff(Pointing &mem, const Det &det, const Det &pair,
                   const PointData &pnt, Map &map);

// Scan the map into one detector's timestream, accumulating with +=.
void map2tod1(Pointing &mem, const Det &det, const PointData &pnt,
              const Map &map);

// Merge a thread-local map into the shared one.
void add_map(Map &map, const Map &local);

}  // namespace qp
