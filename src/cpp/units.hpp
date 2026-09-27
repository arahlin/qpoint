#pragma once

namespace qp {

// Defined unconditionally rather than behind #ifndef M_PI, as qpoint.h does.
// That header only works because -std=c99 sets __STRICT_ANSI__ and hides
// math.h's M_PI, which in turn lets its M_TWOPI definition be reached; under
// a GNU dialect or a libc that defines M_PI regardless, M_TWOPI silently
// vanishes while qp_pixel.c still uses it.
inline constexpr double kPi = 3.14159265358979323846;
inline constexpr double kTwoPi = 6.28318530717958647692;
inline constexpr double kPiHalf = 1.57079632679489661923;

inline constexpr double kD2R = kPi / 180.;
inline constexpr double kR2D = 180. / kPi;
inline constexpr double kAs2R = kPi / (180. * 3600.);

constexpr double deg2rad(double deg) { return deg * kD2R; }
constexpr double rad2deg(double rad) { return rad * kR2D; }
constexpr double arcsec2rad(double sec) { return sec * kAs2R; }

constexpr double secs2days(double s) { return s / 86400.; }
constexpr double jd2mjd(double jd) { return jd - 2400000.5; }

// JD for ctime = 0
inline constexpr double kCtimeJdEpoch = 2440587.5;
// diurnal aberration constant, radians (-0.3191 arcsec)
inline constexpr double kDAberRad = 1.54716541e-06;
// speed of light, AU/day
inline constexpr double kCAud = 173.14463269999999;
// speed of light, m/s
inline constexpr double kCMs = 299792458.0;

}  // namespace qp
