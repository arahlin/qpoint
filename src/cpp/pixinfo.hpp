#pragma once

#include <vector>

namespace qp {

// HEALPix ring geometry and bilinear interpolation, ported from
// healpix-cxx's get_interpol.
//
// Rings are filled on first touch, so a partial scan never builds the whole
// table -- but that mutates the cache as a side effect of reading. Call
// populate() before a parallel region: map2tod hands one PixInfo to every
// thread.
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
