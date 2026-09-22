#pragma once

#include <cstddef>
#include <cstdint>
#include <unordered_map>
#include <vector>

#include "error.hpp"
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

// A thread-local accumulator's index: pixel -> slot in a compact buffer,
// slots handed out in the order the thread first sees each pixel.
//
// This exists for the case where the destination map is far bigger than
// anything one thread can touch -- a full-sky map scanned over a fraction
// of a percent of the sky, where a private copy is 906 MB per thread at
// nside 1024 and threading comes out slower than serial. Compacting makes
// the private copy proportional to the coverage instead.
//
// It is deliberately not used for a map that is already a compact
// footprint, which is the usual case: there the accumulator is
// proportional to the coverage to begin with, and this would only add a
// lookup per sample.
class PixAccum {
 public:
  explicit PixAccum(std::size_t cap) : cap_(cap) {}

  // Slot for a pixel, assigning the next free one if it is new.
  long slot(long pix) {
    const auto it = slot_.find(pix);
    if (it != slot_.end()) return it->second;
    // cap_ bounds the distinct pixels a thread can reach, so this is a
    // guard against the bound being computed wrongly, not a real case.
    if (slot_.size() >= cap_) throw QpMapError("PixAccum: out of slots");
    const long s = static_cast<long>(slot_.size());
    slot_.emplace(pix, s);
    return s;
  }

  // Calls f(pixel, slot) for everything the thread touched.
  template <class Fn>
  void for_each(Fn &&f) const {
    for (const auto &kv : slot_) f(kv.first, kv.second);
  }

  std::size_t size() const { return slot_.size(); }

 private:
  std::unordered_map<long, long> slot_;
  std::size_t cap_;
};

// Records which pixels a thread-local accumulator wrote, so that merging
// it can visit those rather than the whole map.
//
// A scan covers a fraction of a percent of the sky at high nside, and the
// alternative is to read every element of every row looking for the
// nonzero ones -- 113M doubles per thread at nside 1024 to find 53k
// pixels. A bitmap is npix/8 bytes, and iterating it costs 600x less than
// the map scan it replaces.
class PixTouch {
 public:
  explicit PixTouch(std::size_t npix) : bits_((npix + 63) / 64, 0) {}

  void set(long pix) {
    const std::size_t p = static_cast<std::size_t>(pix);
    bits_[p >> 6] |= std::uint64_t(1) << (p & 63);
  }

  // Calls f(pix) for each recorded pixel, ascending.
  template <class Fn>
  void for_each(Fn &&f) const {
    for (std::size_t w = 0; w < bits_.size(); ++w) {
      std::uint64_t b = bits_[w];
      while (b) {
        f(static_cast<long>(w * 64 + __builtin_ctzll(b)));
        b &= b - 1;  // clear the lowest set bit
      }
    }
  }

 private:
  std::vector<std::uint64_t> bits_;
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

  // Both are set only on a thread-local accumulator, and never together:
  // touched records the pixels a directly-indexed accumulator wrote, and
  // accum both assigns and records the slots of a compact one. Null
  // everywhere else, including on the map the caller passed in.
  PixTouch *touched = nullptr;
  PixAccum *accum = nullptr;

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
// Merge local into map. Given the pixels local recorded while
// accumulating, only those are visited; without them the whole map is
// scanned to find the nonzero entries.
void add_map(Map &map, const Map &local, const PixTouch *touched = nullptr);

// Merge a compact accumulator, whose rows are indexed by slot rather than
// by pixel. The index knows which pixel each slot belongs to, so this is
// proportional to what the thread touched.
void add_map(Map &map, const Map &local, const PixAccum &accum);

}  // namespace qp
