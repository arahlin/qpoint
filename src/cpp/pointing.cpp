#include "pointing.hpp"

#include <cfloat>
#include <cmath>
#include <string>
#include <utility>

#include "cshims.hpp"
#include "error.hpp"

namespace qp {

// ---------------------------------------------------------------------------
// Attitude quaternion and its inverse: the only free functions the
// bindings reach for. Everything else in this file is internal.
// ---------------------------------------------------------------------------

Quat azelpsi_quat(double az, double el, double psi, double pitch,
                  double roll) {
  Quat q = Quat::r3(kPi - deg2rad(psi));
  r2_mul(kPiHalf - deg2rad(el), q);
  r3_mul(-deg2rad(az), q);
  if (pitch != 0) r2_mul(-deg2rad(pitch), q);
  if (roll != 0) r1_mul(-deg2rad(roll), q);
  return q;
}

void quat_azelpsi(const Quat &q, double &az, double &el, double &psi) {
  const double w = q[0], x = q[1], y = q[2], z = q[3];

  const double sin_el_sq = x * x + y * y;
  const double cos_el_sq = w * w + z * z;

  const double s = std::atan2(z, w);
  const double d = std::atan2(x, y);

  az = s - d;
  psi = s + d;

  if (cos_el_sq > 1e-12)
    el = 2.0 * std::atan(std::sqrt(sin_el_sq / cos_el_sq));
  else
    el = (sin_el_sq < 0.5) ? 0.0 : kPi;

  el = kPiHalf - el;
  psi = kPi - psi;

  az = -rad2deg(az);
  el = rad2deg(el);
  psi = rad2deg(psi);

  if (psi > 180.0)
    psi -= 360.0;
  else if (psi < -180.0)
    psi += 360.0;

  if (std::fabs(az) < 1e-12) az = 0.0;
  if (std::fabs(el) < 1e-12) el = 0.0;
  if (std::fabs(psi) < 1e-12) psi = 0.0;
}

namespace {

// ---------------------------------------------------------------------------
// Time conversions
// ---------------------------------------------------------------------------

void ctime2jd(double ctime, double jd[2]) {
  jd[0] = kCtimeJdEpoch;
  jd[1] = secs2days(ctime);
}

void ctime2jdtt(double ctime, double jd_tt[2]) {
  double jd_utc[2], jd_tai[2];
  ctime2jd(ctime, jd_utc);
  eraUtctai(jd_utc[0], jd_utc[1], &jd_tai[0], &jd_tai[1]);
  eraTaitt(jd_tai[0], jd_tai[1], &jd_tt[0], &jd_tt[1]);
}

// ---------------------------------------------------------------------------
// Correction quaternions
// ---------------------------------------------------------------------------

// A position on the sphere with a rotation about it, in degrees except for
// ang, which is radians. The three public spellings below differ only in
// where that third angle comes from.
Quat radecang2quat(double ra, double dec, double ang) {
  Quat q = Quat::r3(kPi - ang);
  r2_mul(kPiHalf - deg2rad(dec), q);
  r3_mul(deg2rad(ra), q);
  return q;
}

// kPi - 0. is exactly kPi, so this is the ang = 0 case and nothing moves.
Quat lonlat_quat(double lon, double lat) {
  return radecang2quat(lon, lat, 0.);
}

Quat npb_quat(const double jd_tt[2], int accuracy) {
  double X, Y, s;
  if (accuracy == 0)
    eraXys06a(jd_tt[0], jd_tt[1], &X, &Y, &s);
  else
    eraXys00b(jd_tt[0], jd_tt[1], &X, &Y, &s);
  const double Z = std::sqrt(1.0 - X * X - Y * Y);
  const double E = std::atan2(Y, X);
  const double d = std::acos(Z);

  Quat q = Quat::r3(-E - s);
  r2_mul(d, q);
  r3_mul(E, q);
  return q;
}

Quat erot_quat(const double jd_ut1[2]) {
  return Quat::r3(eraEra00(jd_ut1[0], jd_ut1[1]));
}

Quat wobble_quat(const double jd_tt[2], double xp, double yp) {
  const double sprime = eraSp00(jd_tt[0], jd_tt[1]);
  Quat q = Quat::r1(-arcsec2rad(yp));
  r2_mul(-arcsec2rad(xp), q);
  r3_mul(sprime, q);
  return q;
}

Vec3 earth_orbital_beta(const double jd_tdb[2]) {
  double pvb[2][3];
  eraEpv00(jd_tdb[0], jd_tdb[1], pvb, pvb);
  Vec3 beta;
  for (int i = 0; i < 3; i++) beta[i] = pvb[1][i] / kCAud;
  return beta;
}

// qp_aberration: v = (R(q)*z) x beta, angle = |v|, qa = quat(-angle, v)
Quat aberration(const Quat &q, const Vec3 &beta, bool inv, bool fast) {
  const Vec3 u = q.col3();
  const Vec3 n = inv ? cross(beta, u) : cross(u, beta);
  if (fast) {
    // small angle approximation
    const double sa_2 = 0.5 * n.norm();
    return {{1. - 0.5 * sa_2 * sa_2, -0.5 * n[0], -0.5 * n[1], -0.5 * n[2]}};
  }
  return Quat::rot(-std::asin(n.norm()), n);
}

std::pair<Rate, bool> parse_rate(std::string_view name) {
  bool inv = false;
  if (name.size() > 4 && name.substr(name.size() - 4) == "_inv") {
    inv = true;
    name = name.substr(0, name.size() - 4);
  }
  for (std::size_t i = 0; i < kNumRates; ++i)
    if (kRates[i].name == name) return {static_cast<Rate>(i), inv};
  throw QpError("unknown rate: " + std::string(name));
}

}  // namespace

// ---------------------------------------------------------------------------
// Pointing: construction and parameters
// ---------------------------------------------------------------------------

double refraction(double el, double temp, double press, double hum,
                  double freq) {
  double A, B;
  eraRefco(press, temp, hum, kCMs * 1e-3 / freq, &A, &B);
  if (el > 90) el = 180 - el;
  const double tz = std::tan(kPiHalf - deg2rad(el));
  return rad2deg(tz * (A + B * tz * tz));
}

UpdateState &Pointing::state(Rate r, bool inv) {
  return (inv ? inv_ : fwd_)[static_cast<std::size_t>(r)];
}

const UpdateState &Pointing::state(Rate r, bool inv) const {
  return (inv ? inv_ : fwd_)[static_cast<std::size_t>(r)];
}

void Pointing::set_rate(std::string_view name, double rate) {
  auto [r, inv] = parse_rate(name);
  state(r, inv).set_rate(rate);
}

double Pointing::get_rate(std::string_view name) const {
  auto [r, inv] = parse_rate(name);
  return state(r, inv).rate();
}

void Pointing::reset_rate(std::string_view name) {
  auto [r, inv] = parse_rate(name);
  state(r, inv).reset();
}

void Pointing::reset_rates() {
  for (auto &s : fwd_) s.reset();
}

void Pointing::reset_inv_rates() {
  for (auto &s : inv_) s.reset();
}

void Pointing::set_opt(std::string_view name, int val) {
  for (const auto &d : kOptionParams) {
    if (d.name != name) continue;
    if (opt_.*(d.mem) != val) {
      opt_.*(d.mem) = val;
      if (d.reset != Rate::COUNT) state(d.reset, false).reset();
    }
    return;
  }
  throw QpError("unknown option: " + std::string(name));
}

int Pointing::get_opt(std::string_view name) const {
  for (const auto &d : kOptionParams)
    if (d.name == name) return opt_.*(d.mem);
  throw QpError("unknown option: " + std::string(name));
}

void Pointing::set_weather_param(std::string_view name, double val) {
  for (const auto &d : kWeatherParams) {
    if (d.name != name) continue;
    if (weather_.*(d.mem) != val) {
      weather_.*(d.mem) = val;
      state(Rate::ref, false).reset();
    }
    return;
  }
  throw QpError("unknown weather parameter: " + std::string(name));
}

double Pointing::get_weather_param(std::string_view name) const {
  for (const auto &d : kWeatherParams)
    if (d.name == name) return weather_.*(d.mem);
  throw QpError("unknown weather parameter: " + std::string(name));
}

void Pointing::set_double(std::string_view name, double val) {
  if (name == "ref_delta") {
    ref_delta_ = val;
    return;
  }
  if (name == "dut1") {
    if (val != dut1_) {
      dut1_ = val;
      state(Rate::erot, false).reset();
      state(Rate::wobble, false).reset();
    }
    return;
  }
  throw QpError("unknown parameter: " + std::string(name));
}

double Pointing::get_double(std::string_view name) const {
  if (name == "ref_delta") return ref_delta_;
  if (name == "dut1") return dut1_;
  throw QpError("unknown parameter: " + std::string(name));
}

// ---------------------------------------------------------------------------
// Sidereal time
// ---------------------------------------------------------------------------

// How far from midnight to stop trusting the cache, in days. A UTC day
// containing a leap second is 86401 seconds long, so ERFA still calls the
// instant 86400 seconds in "the previous day, fraction 1.0", while a
// jd = ctime / 86400 reading has already rolled over. The two notions of
// the date therefore disagree for exactly one second at the end of such a
// day. Staying ten seconds clear of midnight keeps the cache away from
// it, at the price of 20 seconds of exact calls per day.
constexpr double kUt1Margin = 10. / 86400.;

void Pointing::jdutc2jdut1(const double jd_utc[2], double jd_ut1[2]) {
  if (ut1_valid_.covers(jd_utc[1], jd_utc[0], dut1_)) {
    jd_ut1[0] = jd_utc[0] + ut1_off0_;
    jd_ut1[1] = jd_utc[1] + ut1_off1_;
    return;
  }

  // Always answer exactly; the cache is only ever an accelerator for the
  // samples that follow.
  eraUtcut1(jd_utc[0], jd_utc[1], dut1_, &jd_ut1[0], &jd_ut1[1]);

  // Rebuild it for the interior of this day. UT1 - UTC is dut1 minus the
  // TAI - UTC step, and that step only moves at a leap second, which
  // falls at midnight -- so if the offset agrees at both ends of the
  // interior it is constant across the whole of it. When it does not, the
  // day has a leap second in it and gets no cache at all.
  const double day = std::floor(jd_utc[0] + jd_utc[1] + 0.5);
  const double lo = day - 0.5 - jd_utc[0] + kUt1Margin;
  const double hi = day + 0.5 - jd_utc[0] - kUt1Margin;
  double a[2], b[2];
  eraUtcut1(jd_utc[0], lo, dut1_, &a[0], &a[1]);
  eraUtcut1(jd_utc[0], hi, dut1_, &b[0], &b[1]);

  if (a[0] - jd_utc[0] == b[0] - jd_utc[0] && a[1] - lo == b[1] - hi) {
    ut1_off0_ = a[0] - jd_utc[0];
    ut1_off1_ = a[1] - lo;
    ut1_valid_.set(lo, hi, jd_utc[0], dut1_);
  } else {
    ut1_valid_.clear();  // the day has a leap second in it
  }
}

double Pointing::gmst(double ctime) {
  double jd_utc[2];
  ctime2jd(ctime, jd_utc);

  // UTC -> UT1 happens in both accuracy modes, unlike the C, whose 'low'
  // path hands UTC to eraGmst00 for UT1 and so discards dut1 entirely --
  // silently undoing a bulletin the caller went to the trouble of loading.
  // The cache is what makes keeping it affordable: the offset is fixed
  // within the day, so 'low' pays about 2 ns a sample for the term that
  // dominates this calculation, and the transforms apply it in both modes
  // already (update_erot reaches the same conversion with no accuracy gate).
  if (state(Rate::dut1, false).check(ctime))
    dut1_ = bulletin_.interp(jd2mjd(jd_utc[0]) + jd_utc[1]).dut1;
  double jd_ut1[2];
  jdutc2jdut1(jd_utc, jd_ut1);

  double g;
  if (opt_.accuracy == 0) {
    double jd_tt[2];
    ctime2jdtt(ctime, jd_tt);
    g = eraGmst00(jd_ut1[0], jd_ut1[1], jd_tt[0], jd_tt[1]);
  } else {
    // 'low' gives up only the TT conversion, which is the whole of the
    // remaining cost. GMST reads TT through the precession polynomial
    // alone, so the 69 s error that represents is worth 0.1 mas.
    g = eraGmst00(jd_ut1[0], jd_ut1[1], jd_ut1[0], jd_ut1[1]);
  }
  return std::fmod(rad2deg(g) / 15.0, 24.);
}

double Pointing::lmst(double ctime, double lon) {
  return std::fmod(gmst(ctime) + lon / 15.0, 24.);
}

// ---------------------------------------------------------------------------
// Corrections
// ---------------------------------------------------------------------------

double Pointing::update_ref(const Quat &q) {
  const double s = q[0] * q[0] - q[1] * q[1] - q[2] * q[2] + q[3] * q[3];
  const double el = rad2deg(opt_.fast_math ? poly_asin(s) : std::asin(s));
  // qp_update_ref stores this as a side effect, and ref_delta is a readable
  // parameter, so the write is part of the contract.
  ref_delta_ = refraction(el, weather_.temperature, weather_.pressure,
                          weather_.humidity, weather_.frequency);
  return ref_delta_;
}

void Pointing::apply_refraction(double ctime, Quat &q, bool inv) {
  UpdateState &st = state(Rate::ref, inv);
  Quat &q_ref = inv ? q_ref_inv_ : q_ref_;

  if (st.check(ctime)) {
    double delta = update_ref(q);
    if (inv) delta *= -1;
    q_ref = Quat::r2(-deg2rad(delta));
  }
  if (st.should_apply()) mul_right(q, q_ref);
}

// Both aberration corrections read the forward rate state even when applied
// in the inverse direction; that is what the C does, and the *_inv rate
// parameters for daber and aaber are consequently inert.
void Pointing::apply_diurnal_aberration(double ctime, double lat, Quat &q,
                                        bool inv) {
  UpdateState &st = state(Rate::daber, false);
  if (st.check(ctime)) {
    const double clat =
        opt_.fast_math ? poly_cos(deg2rad(lat)) : std::cos(deg2rad(lat));
    beta_rot_[0] = beta_rot_[2] = 0;
    beta_rot_[1] = -clat * kDAberRad;
  }
  if (st.should_apply()) {
    const Quat q_aber = aberration(q, beta_rot_, inv, opt_.fast_aber);
    mul_left(q_aber, q);
  }
}

void Pointing::apply_annual_aberration(double ctime, Quat &q, bool inv) {
  UpdateState &st = state(Rate::aaber, false);
  if (st.check(ctime)) {
    double jd_tt[2];
    ctime2jdtt(ctime, jd_tt);
    beta_earth_ = earth_orbital_beta(jd_tt);
  }
  if (st.should_apply()) {
    const Quat q_aber = aberration(q, beta_earth_, inv, opt_.fast_aber);
    mul_left(q_aber, q);
  }
}

// ---------------------------------------------------------------------------
// Horizon -> equatorial
// ---------------------------------------------------------------------------

void Pointing::azelpsi2quat(double az, double el, double psi, double pitch,
                            double roll, double lon, double lat, double ctime,
                            Quat &q) {
  double jd_utc[2], jd_tt[2] = {0, 0}, jd_ut1[2];
  ctime2jd(ctime, jd_utc);

  // elevations through zenith fold over the horizon
  if (el > 90) {
    el = 180 - el;
    az += 180;
    psi -= 180;
  }

  // psi is deferred so that refraction can be inserted before it
  mul_left(azelpsi_quat(az, el, 0, pitch, roll), q);

  // NB: per-detector refraction is not fully implemented; applied here as a
  // mean correction, and right-applied unlike everything below
  apply_refraction(ctime, q, false);

  if (psi != 0) mul_right(q, Quat::r3(-deg2rad(psi)));

  apply_diurnal_aberration(ctime, lat, q, false);

  // rotate to ITRS
  if (state(Rate::lonlat, false).check(ctime)) q_lonlat_ = lonlat_quat(lon, lat);
  if (state(Rate::lonlat, false).should_apply()) mul_left(q_lonlat_, q);

  // polar motion, or just dut1 from the IERS bulletin
  const double mjd_utc = jd2mjd(jd_utc[0]) + jd_utc[1];
  if (state(Rate::wobble, false).check(ctime)) {
    const IersValues v = bulletin_.interp(mjd_utc);
    dut1_ = v.dut1;
    ctime2jdtt(ctime, jd_tt);
    q_wobble_ = wobble_quat(jd_tt, v.x, v.y);
  } else if (state(Rate::dut1, false).check(ctime)) {
    dut1_ = bulletin_.interp(mjd_utc).dut1;
  }
  if (state(Rate::wobble, false).should_apply()) mul_left(q_wobble_, q);

  if (state(Rate::erot, false).check(ctime)) {
    jdutc2jdut1(jd_utc, jd_ut1);
    q_erot_ = erot_quat(jd_ut1);
  }
  if (state(Rate::erot, false).should_apply()) mul_left(q_erot_, q);

  if (state(Rate::npb, false).check(ctime)) {
    if (jd_tt[0] == 0) ctime2jdtt(ctime, jd_tt);
    q_npb_ = npb_quat(jd_tt, opt_.accuracy);
  }
  if (state(Rate::npb, false).should_apply()) mul_left(q_npb_, q);

  // ~20 arcsec max
  if (opt_.mean_aber) apply_annual_aberration(ctime, q, false);
}

// ---------------------------------------------------------------------------
// Equatorial -> horizon
// ---------------------------------------------------------------------------

void Pointing::quat2azel(const Quat &q_in, double lon, double lat,
                         double ctime, double &az, double &el, double &pa) {
  double jd_utc[2], jd_tt[2] = {0, 0}, jd_ut1[2];
  Quat q = q_in;
  ctime2jd(ctime, jd_utc);

  apply_annual_aberration(ctime, q, true);

  if (state(Rate::npb, true).check(ctime)) {
    ctime2jdtt(ctime, jd_tt);
    q_npb_inv_ = npb_quat(jd_tt, opt_.accuracy).inv();
  }
  if (state(Rate::npb, true).should_apply()) mul_left(q_npb_inv_, q);

  const double mjd_utc = jd2mjd(jd_utc[0]) + jd_utc[1];
  if (state(Rate::wobble, true).check(ctime)) {
    const IersValues v = bulletin_.interp(mjd_utc);
    dut1_ = v.dut1;
    if (jd_tt[0] == 0) ctime2jdtt(ctime, jd_tt);
    q_wobble_inv_ = wobble_quat(jd_tt, v.x, v.y).inv();
  } else if (state(Rate::dut1, false).check(ctime)) {
    // NB: the forward dut1 state, not the inverse one -- as in the C
    dut1_ = bulletin_.interp(mjd_utc).dut1;
  }

  if (state(Rate::erot, true).check(ctime)) {
    jdutc2jdut1(jd_utc, jd_ut1);
    q_erot_inv_ = erot_quat(jd_ut1).inv();
  }
  if (state(Rate::erot, true).should_apply()) mul_left(q_erot_inv_, q);

  if (state(Rate::wobble, true).should_apply()) mul_left(q_wobble_inv_, q);

  if (state(Rate::lonlat, true).check(ctime))
    q_lonlat_inv_ = lonlat_quat(lon, lat).inv();
  if (state(Rate::lonlat, true).should_apply()) mul_left(q_lonlat_inv_, q);

  apply_refraction(ctime, q, true);
  apply_diurnal_aberration(ctime, lat, q, true);

  quat2radec(q, DecOut::Dec, PolOut::PA, az, el, pa, nullptr);
  az *= -1;
}

// ---------------------------------------------------------------------------
// Quaternion <-> sky coordinates
// ---------------------------------------------------------------------------

void Pointing::quat2radec(Quat q, DecOut dmode, PolOut pmode,
                          double &ra, double &dec, double &p1,
                          double *p2) const {
  // ZYZ euler angles; factors of two have been redistributed
  const double q00p33 = q[0] * q[0] + q[3] * q[3];
  const double q11p22 = q[1] * q[1] + q[2] * q[2];
  const double cosb2 = q00p33 * q11p22;
  const double sinb = q00p33 - q11p22;
  const double sinb_2 = 0.5 * sinb;

  if (cosb2 < DBL_EPSILON) {
    ra = 0;
    dec = (dmode == DecOut::Dec) ? (sinb_2 > 0 ? 90. : -90.) : sinb;
  } else {
    const double sina_2 = q[2] * q[3] - q[0] * q[1];
    const double cosa_2 = q[0] * q[2] + q[1] * q[3];

    if (opt_.fast_math)
      ra = rad2deg(poly_atan2(sina_2, cosa_2));
    else
      ra = rad2deg(std::atan2(sina_2, cosa_2));

    if (dmode == DecOut::Dec) {
      const double cosb_2 = std::sqrt(cosb2);
      dec = rad2deg(opt_.fast_math ? poly_atan2(sinb_2, cosb_2)
                                   : std::atan2(sinb_2, cosb_2));
    } else {
      dec = sinb;
    }
  }

  pol_out(q, cosb2, sinb_2 > 0, pmode, p1, p2);
}

// The polarization half of the quaternion -> sky conversion, shared with
// quat2pix: the two differ only in how cosb2 and the pole test are
// derived, and every term here keeps the operand order it had in both, so
// the results are bit-identical to the two copies this replaced.
void Pointing::pol_out(Quat q, double cosb2, bool north, PolOut pmode,
                       double &p1, double *p2) const {
  double sing, cosg, norm = 0.;

  if (cosb2 < DBL_EPSILON) {
    if (north) {
      cosg = q[3] * q[3] - q[0] * q[0];
      sing = 2 * q[0] * q[3];
    } else {
      cosg = q[1] * q[1] - q[2] * q[2];
      sing = 2 * q[1] * q[2];
    }
    if (pmode == PolOut::SinCos) norm = 2. * cosg;
  } else {
    cosg = q[1] * q[3] - q[0] * q[2];
    sing = q[0] * q[1] + q[2] * q[3];
    if (pmode == PolOut::SinCos) norm = 2. * cosg / cosb2;
  }

  if (pmode == PolOut::SinCos) {
    p1 = norm * sing;
    *p2 = norm * cosg - 1.;
  } else {
    p1 = rad2deg(opt_.fast_math ? poly_atan2(sing, cosg)
                                : std::atan2(sing, cosg));
  }
}

Quat Pointing::radecpa2quat(double ra, double dec, double pa) const {
  return radecang2quat(ra, dec, deg2rad(pa));
}

Quat Pointing::radec2quat(double ra, double dec, double sin2psi,
                          double cos2psi) const {
  return radecang2quat(ra, dec, opt_.fast_math
                                    ? poly_atan2(sin2psi, cos2psi + 1)
                                    : std::atan2(sin2psi, cos2psi + 1));
}

// ---------------------------------------------------------------------------
// CMB dipole
// ---------------------------------------------------------------------------

namespace {
// Planck 2015 values, (l, b) = (264.00, 48.24)
constexpr double kDipoleRa = 167.923;
constexpr double kDipoleDec = -6.947;
}  // namespace

double Pointing::cdist2dipole(double cdist, double ctime) const {
  const double tcmb = 2.7255;             // Fixsen 2009
  const double beta = 3364.5e-6 / tcmb;   // Planck 2015
  const double vhelio = 0.00027;          // annual modulation
  const double dipole_epoch = 2451170;

  double out = tcmb * beta * (cdist + beta / 2. * (2 * cdist * cdist - 1));

  double jd[2];
  ctime2jd(ctime, jd);
  const double delta = (jd[1] + (jd[0] - dipole_epoch)) / 365.25;
  out += vhelio * (opt_.fast_math ? poly_cos(2 * kPi * delta)
                                  : std::cos(2 * kPi * delta));
  return out;
}

void Pointing::init_dipole() {
  if (dipole_init_) return;
  v_dipole_ = radecpa2quat(kDipoleRa, kDipoleDec, 0.).col3();
  dipole_init_ = true;
}

double Pointing::quat2dipole(double ctime, const Quat &q) {
  init_dipole();
  return cdist2dipole(dot(v_dipole_, q.col3()), ctime);
}

double Pointing::dipole(double ctime, double ra, double dec) const {
  const double dipole_phi = deg2rad(kDipoleRa);
  const double dipole_theta = kPi / 2 - deg2rad(kDipoleDec);
  const double sdtheta = std::sin(dipole_theta);
  const double cdtheta = std::cos(dipole_theta);

  const double theta = kPi / 2 - deg2rad(dec);
  const double phi = deg2rad(ra);

  double stheta, ctheta, cdphi;
  if (opt_.fast_math) {
    stheta = poly_sin(theta);
    ctheta = poly_cos(theta);
    cdphi = poly_cos(dipole_phi - phi);
  } else {
    stheta = std::sin(theta);
    ctheta = std::cos(theta);
    cdphi = std::cos(dipole_phi - phi);
  }

  return cdist2dipole(cdtheta * ctheta + sdtheta * stheta * cdphi, ctime);
}

Quat Pointing::bore2det(const Quat &q_off, double ctime, const Quat &q_bore) {
  Quat q_det = q_off;
  mul_left(q_bore, q_det);
  if (!opt_.mean_aber) apply_annual_aberration(ctime, q_det, false);
  return q_det;
}

Quat Pointing::bore2det_hwp(const Quat &q_off, double ctime,
                            const Quat &q_bore, const Quat &q_hwp) {
  Quat q_det = bore2det(q_off, ctime, q_bore);
  mul_right(q_det, q_hwp);
  return q_det;
}

}  // namespace qp
