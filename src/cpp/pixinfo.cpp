#include "pixinfo.hpp"

#include <cmath>

#include "cshims.hpp"
#include "units.hpp"

namespace qp {

PixInfo::PixInfo(long nside)
    : nside_(nside),
      npface_(nside * nside),
      npix_(12 * npface_),
      ncap_((npface_ - nside_) << 1),
      fact2_(4. / static_cast<double>(npix_)),
      fact1_(static_cast<double>(nside_ << 1) * fact2_),
      rings_(static_cast<std::size_t>(4 * nside)) {}

void PixInfo::populate() {
  for (long i = 0; i < 4 * nside_; ++i) ring(i);
}

const PixInfo::Ring &PixInfo::ring(long iring) const {
  Ring &r = rings_[static_cast<std::size_t>(iring)];
  if (!r.init) r = make_ring(iring);
  return r;
}

PixInfo::Ring PixInfo::make_ring(long iring) const {
  Ring r;
  const long northring = (iring > 2 * nside_) ? (4 * nside_ - iring) : iring;

  if (northring < nside_) {
    const double tmp = static_cast<double>(northring) * northring * fact2_;
    const double costheta = 1 - tmp;
    const double sintheta = std::sqrt(tmp * (2 - tmp));
    r.theta = std::atan2(sintheta, costheta);
    r.ringpix = 4 * northring;
    r.shifted = 1;
    r.startpix = 2 * northring * (northring - 1);
  } else {
    r.theta = std::acos((2 * nside_ - northring) * fact1_);
    r.ringpix = 4 * nside_;
    r.shifted = (((northring - nside_) & 1) == 0);
    r.startpix = ncap_ + (northring - nside_) * r.ringpix;
  }

  if (northring != iring) {
    r.theta = kPi - r.theta;
    r.startpix = npix_ - r.startpix - r.ringpix;
  }

  r.init = true;
  return r;
}

bool PixInfo::get_interpol_ring(double theta, double phi, long pix[4],
                                double weight[4]) const {
  if (theta < 0 || theta > kPi) return false;

  const int nside = static_cast<int>(nside_);
  const double z = std::cos(theta);
  const double az = std::fabs(z);

  long ir1;
  if (az < 2. / 3.) {
    ir1 = static_cast<long>(nside * (2 - 1.5 * z));
  } else {
    ir1 = static_cast<long>(nside * std::sqrt(3 * (1 - az)));
    if (z <= 0) ir1 = 4 * nside - ir1 - 1;
  }
  const long ir2 = ir1 + 1;

  // The two bracketing rings are handled identically, into their own pair
  // of pix/weight slots; returns the ring's theta.
  auto fill = [&](long iring, long *pix2, double *wt2) {
    const Ring &r = ring(iring);
    const double dphi = kTwoPi / r.ringpix;
    const double tmp = phi / dphi - 0.5 * r.shifted;
    long i1 = static_cast<long>((tmp < 0) ? tmp - 1 : tmp);
    const double w1 = (phi - (i1 + 0.5 * r.shifted) * dphi) / dphi;
    if (i1 < 0) i1 += r.ringpix;
    long i2 = i1 + 1;
    if (i2 >= r.ringpix) i2 -= r.ringpix;
    pix2[0] = r.startpix + i1;
    pix2[1] = r.startpix + i2;
    wt2[0] = 1 - w1;
    wt2[1] = w1;
    return r.theta;
  };

  double theta1 = 0, theta2 = 0;
  if (ir1 > 0) theta1 = fill(ir1, pix, weight);
  if (ir2 < 4 * nside) theta2 = fill(ir2, pix + 2, weight + 2);

  if (ir1 == 0) {
    // north polar cap: fold the four top pixels together
    const double wtheta = theta / theta2;
    weight[2] *= wtheta;
    weight[3] *= wtheta;
    const double fac = (1 - wtheta) * 0.25;
    weight[0] = fac;
    weight[1] = fac;
    weight[2] += fac;
    weight[3] += fac;
    pix[0] = (pix[2] + 2) & 3;
    pix[1] = (pix[3] + 2) & 3;
  } else if (ir2 == 4 * nside) {
    const double wtheta = (theta - theta1) / (kPi - theta1);
    weight[0] *= (1 - wtheta);
    weight[1] *= (1 - wtheta);
    const double fac = wtheta * 0.25;
    weight[0] += fac;
    weight[1] += fac;
    weight[2] = fac;
    weight[3] = fac;
    pix[2] = ((pix[0] + 2) & 3) + npix_ - 4;
    pix[3] = ((pix[1] + 2) & 3) + npix_ - 4;
  } else {
    const double wtheta = (theta - theta1) / (theta2 - theta1);
    weight[0] *= (1 - wtheta);
    weight[1] *= (1 - wtheta);
    weight[2] *= wtheta;
    weight[3] *= wtheta;
  }

  return true;
}

bool PixInfo::get_interpol_nest(double theta, double phi, long pix[4],
                                double weight[4]) const {
  if (!get_interpol_ring(theta, phi, pix, weight)) return false;
  for (int i = 0; i < 4; ++i) ring2nest(nside_, pix[i], pix + i);
  return true;
}

void PixInfo::interpol(double ra, double dec, bool nest, long pix[4],
                       double weight[4]) const {
  const double theta = kPiHalf - deg2rad(dec);
  const double phi = deg2rad(ra);
  if (nest)
    get_interpol_nest(theta, phi, pix, weight);
  else
    get_interpol_ring(theta, phi, pix, weight);
}

double PixInfo::interp_val(const double *map, double ra, double dec,
                           bool nest) const {
  long pix[4];
  double weight[4];
  interpol(ra, dec, nest, pix, weight);

  double val = 0;
  for (int i = 0; i < 4; ++i) val += map[pix[i]] * weight[i];
  return val;
}

}  // namespace qp
