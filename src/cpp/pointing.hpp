#pragma once

#include <array>
#include <string_view>

#include "bulletin_a.hpp"
#include "quat.hpp"
#include "state.hpp"
#include "units.hpp"
#include "vec3.hpp"

namespace qp {

// Stateless quaternion constructors -- qp_det_offset and qp_hwp_quat.

// Detector offset quaternion from boresight-relative angles, in degrees.
inline Quat det_offset(double delta_az, double delta_el, double delta_psi) {
  Quat q = Quat::r3(-deg2rad(delta_psi));
  r2_mul(deg2rad(delta_el), q);
  r1_mul(-deg2rad(delta_az), q);
  return q;
}

// Waveplate quaternion from the physical HWP angle, in degrees. The
// polarization angle rotates by twice the physical angle.
inline Quat hwp_quat(double ang) { return Quat::r3(-2. * deg2rad(ang)); }

// Attitude quaternion from horizon angles, in degrees.
Quat azelpsi_quat(double az, double el, double psi, double pitch, double roll);

// Its inverse. Used by the gyro integration, which needs the round trip.
void quat_azelpsi(const Quat &q, double &az, double &el, double &psi);

// Atmospheric refraction in degrees, from the elevation and the weather.
// Stateless, so it is a free function: Pointing::update_ref is the one that
// reads the stored weather and caches the result.
double refraction(double el, double temp, double press, double hum,
                  double freq);

// The time conversions and the individual correction quaternions are
// internal to pointing.cpp -- nothing outside it ever wanted one.

// The two output axes of the quat -> sky-coordinate conversion. They are
// independent, so all four combinations are available -- the C library only
// provides three, lacking (SinDec, PA).
enum class DecOut { Dec, SinDec };
enum class PolOut { SinCos, PA };

// Replaces qp_memory_t. Holds the correction rate caches, the cached
// correction quaternions, weather, options and the IERS table.
//
// Copyable, and deliberately cheap to copy: OpenMP takes one copy per
// thread, and the IERS table is shared rather than duplicated.
class Pointing {
 public:
  Pointing() = default;

  // ---- rate states ----
  void set_rate(std::string_view name, double rate);
  double get_rate(std::string_view name) const;
  void reset_rate(std::string_view name);
  void reset_rates();
  void reset_inv_rates();

  // ---- scalar parameters ----
  void set_opt(std::string_view name, int val);
  int get_opt(std::string_view name) const;
  void set_weather_param(std::string_view name, double val);
  double get_weather_param(std::string_view name) const;
  void set_double(std::string_view name, double val);
  double get_double(std::string_view name) const;

  const Options &opt() const { return opt_; }
  BulletinA &bulletin() { return bulletin_; }

  // ---- scalar workers ----
  double gmst(double ctime);
  double lmst(double ctime, double lon);

  // Accumulates into q rather than initializing it: the azel2radec family
  // pre-seeds q with the detector offset, and depends on this.
  void azelpsi2quat(double az, double el, double psi, double pitch,
                    double roll, double lon, double lat, double ctime,
                    Quat &q);

  void quat2azel(const Quat &q_in, double lon, double lat, double ctime,
                 double &az, double &el, double &pa);

  // Takes its Quat by value, not by const reference, and so do the other
  // workers a binding loop hands a freshly computed orientation to. A Quat
  // is four doubles, which AAPCS passes and returns in v0-v3; a reference
  // parameter instead forces a caller holding one in registers -- which
  // every loop is, having just built it -- to spill it to the stack purely
  // to have an address to pass, and the callee to load it back. That round
  // trip is on the loop's critical path and costs about 8% of bore2radec.
  // The copy it replaces is free. Don't restore the const Quat &.
  //
  // p2 is only written in PolOut::SinCos mode, so PA-mode callers pass
  // nullptr rather than inventing somewhere for it to go.
  void quat2radec(Quat q, DecOut dmode, PolOut pmode, double &ra,
                  double &dec, double &p1, double *p2) const;

  Quat radecpa2quat(double ra, double dec, double pa) const;
  Quat radec2quat(double ra, double dec, double sin2psi,
                  double cos2psi) const;

  Quat bore2det(const Quat &q_off, double ctime, const Quat &q_bore);
  Quat bore2det_hwp(const Quat &q_off, double ctime, const Quat &q_bore,
                    const Quat &q_hwp);

  double update_ref(const Quat &q);

  // ---- pixelization ----
  long radec2pix(double ra, double dec, int nside) const;

  // Combines qp_quat2pix and qp_quat2pixpa; p2 is untouched for PolOut::PA.
  // By value, for the reason quat2radec is.
  void quat2pix(Quat q, int nside, PolOut pmode, long &pix, double &p1,
                double *p2) const;

  void pixel_offset(int nside, long pix, double ra, double dec,
                    double &dtheta, double &dphi) const;

  // ---- galactic rotation ----
  void radec2gal_quat(Quat &q);
  void gal2radec_quat(Quat &q);

  // Rotate a sky position and its polarization basis in place.
  void rotate_coord(double &ra, double &dec, double &sin2psi, double &cos2psi,
                    bool to_gal);

  // Sky coordinates of a pixel centre, in degrees.
  void pix2radec(int nside, long pix, double &ra, double &dec) const;

  // ---- CMB dipole ----
  double dipole(double ctime, double ra, double dec) const;
  double quat2dipole(double ctime, const Quat &q);

 private:
  UpdateState &state(Rate r, bool inv);
  const UpdateState &state(Rate r, bool inv) const;

  // UTC -> UT1, with the conversion cached for the day. The offset it
  // applies is fixed within a calendar day while the Earth rotation angle
  // it feeds changes every sample, and eraUtcut1 is 52 ns against
  // eraEra00's 9.6, so recomputing it per sample dominates the cost of
  // the whole correction.
  void jdutc2jdut1(const double jd_utc[2], double jd_ut1[2]);

  void apply_refraction(double ctime, Quat &q, bool inv);
  void apply_diurnal_aberration(double ctime, double lat, Quat &q, bool inv);
  void apply_annual_aberration(double ctime, Quat &q, bool inv);

  // The polarization half of quat2radec, shared with quat2pix. By value
  // for the reason quat2radec is: it is handed a Quat the caller already
  // holds in registers.
  void pol_out(Quat q, double cosb2, bool north, PolOut pmode,
               double &p1, double *p2) const;

  void init_gal();
  void init_dipole();
  double cdist2dipole(double cdist, double ctime) const;

  std::array<UpdateState, kNumRates> fwd_ = kInitialRateStates;
  std::array<UpdateState, kNumRates> inv_ = kInitialRateStates;

  Options opt_;
  Weather weather_;
  BulletinA bulletin_;

  Quat q_lonlat_{}, q_lonlat_inv_{};
  Quat q_wobble_{}, q_wobble_inv_{};
  Quat q_npb_{}, q_npb_inv_{};
  Quat q_erot_{}, q_erot_inv_{};
  Quat q_ref_{}, q_ref_inv_{};

  Quat q_gal_{}, q_gal_inv_{};
  bool gal_init_ = false;

  Vec3 v_dipole_{};
  bool dipole_init_ = false;

  Vec3 beta_earth_{}, beta_rot_{};

  double ref_delta_ = 0.;
  double dut1_ = 0.;

  // The cached UT1 - UTC, valid while jd_utc[1] is inside the window at
  // this jd_utc[0] and this dut1_. Starts empty.
  ValidWindow ut1_valid_;
  double ut1_off0_ = 0., ut1_off1_ = 0.;
};

}  // namespace qp
