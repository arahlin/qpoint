#pragma once

#include <vector>

namespace qp {

// HEALPix ring geometry and bilinear interpolation, ported from
// healpix-cxx's get_interpol, which the C library copied because chealpix
// has no equivalent.
//
// Rings are computed on first touch and cached, so a scan that only covers
// part of the sky never builds the rest of the table. That cache is mutated
// as a side effect of reading, which makes the lazy path unsafe to share
// across threads: qp_map2tod hands one map, and so one PixInfo, to every
// OpenMP worker, and with interp_pix enabled several of them can initialize
// the same ring at once. Call populate() before entering a parallel region
// and every later access is read-only.
class PixInfo {
 public:
  explicit PixInfo(long nside);

  // Build the whole ring table up front. Required before sharing this
  // object across threads; pointless otherwise.
  void populate();

  // Four neighbouring pixels and their weights, for sky coordinates in
  // degrees. Out-of-range coordinates leave pix and weight untouched.
  void interpol(double ra, double dec, bool nest, long pix[4],
                double weight[4]) const;

  // Bilinearly interpolate a map at sky coordinates, in degrees.
  double interp_val(const double *map, double ra, double dec,
                    bool nest) const;

 private:
  struct Ring {
    long startpix = 0;
    long ringpix = 0;
    double theta = 0.;
    int shifted = 0;
    bool init = false;
  };

  // Fills on first touch, which is why rings_ is mutable.
  const Ring &ring(long iring) const;
  Ring make_ring(long iring) const;

  bool get_interpol_ring(double theta, double phi, long pix[4],
                         double weight[4]) const;
  bool get_interpol_nest(double theta, double phi, long pix[4],
                         double weight[4]) const;

  long nside_;
  long npface_;
  long npix_;
  long ncap_;
  double fact2_;
  double fact1_;
  mutable std::vector<Ring> rings_;
};

}  // namespace qp
