#include "cshims.hpp"
#include "pointing.hpp"

namespace qp {

long Pointing::radec2pix(double ra, double dec, int nside) const {
  long pix;
  const double theta = kPiHalf - deg2rad(dec);
  const double phi = deg2rad(ra);
  if (opt_.pix_order == 1)
    ang2pix_nest(nside, theta, phi, &pix);
  else
    ang2pix_ring(nside, theta, phi, &pix);
  return pix;
}

void Pointing::quat2pix(Quat q, int nside, PolOut pmode, long &pix,
                        double &p1, double *p2) const {
  if (!opt_.fast_pix) {
    double ra, dec;
    quat2radec(q, DecOut::Dec, pmode, ra, dec, p1, p2);
    pix = radec2pix(ra, dec, nside);
    return;
  }

  // fast_pix skips the angle round trip: the pixel comes straight from the
  // pointing vector, and the polarization angle from the quaternion.
  Vec3 vec = q.col3();
  if (opt_.pix_order == 1)
    vec2pix_nest(nside, vec.data(), &pix);
  else
    vec2pix_ring(nside, vec.data(), &pix);

  pol_out(q, (1 - vec[2] * vec[2]) / 4., vec[2] > 0, pmode, p1, p2);
}

void Pointing::pixel_offset(int nside, long pix, double ra, double dec,
                            double &dtheta, double &dphi) const {
  if (opt_.pix_order == 1)
    pix2ang_nest(nside, pix, &dtheta, &dphi);
  else
    pix2ang_ring(nside, pix, &dtheta, &dphi);
  dtheta = kPiHalf - deg2rad(dec) - dtheta;
  if (dtheta < -kPiHalf) dtheta += kPi;
  if (dtheta > kPiHalf) dtheta -= kPi;
  dphi = deg2rad(ra) - dphi;
  if (dphi < -kPi) dphi += kTwoPi;
  if (dphi > kPi) dphi -= kTwoPi;
}

void Pointing::init_gal() {
  if (gal_init_) return;
  // galactic pole, cf. sofa/g2icrs
  q_gal_ = radecpa2quat(192.85948, 27.12825, 90 + 32.93192);
  q_gal_inv_ = q_gal_.inv();
  gal_init_ = true;
}

void Pointing::radec2gal_quat(Quat &q) {
  init_gal();
  mul_left(q_gal_inv_, q);
}

void Pointing::gal2radec_quat(Quat &q) {
  init_gal();
  mul_left(q_gal_, q);
}

void Pointing::rotate_coord(double &ra, double &dec, double &sin2psi,
                            double &cos2psi, bool to_gal) {
  Quat q = radec2quat(ra, dec, sin2psi, cos2psi);
  if (to_gal)
    radec2gal_quat(q);
  else
    gal2radec_quat(q);
  quat2radec(q, DecOut::Dec, PolOut::SinCos, ra, dec, sin2psi, &cos2psi);
}

void Pointing::pix2radec(int nside, long pix, double &ra, double &dec) const {
  double theta, phi;
  if (opt_.pix_order == 1)
    pix2ang_nest(nside, pix, &theta, &phi);
  else
    pix2ang_ring(nside, pix, &theta, &phi);
  dec = rad2deg(kPiHalf - theta);
  ra = rad2deg(phi);
}

}  // namespace qp
