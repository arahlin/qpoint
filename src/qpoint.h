#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>

  /* *************************************************************************
     Types and parameters
     ********************************************************************** */

  /* Quaternion */
  typedef double quat_t[4];

  /* Mueller element vector */
  typedef double mueller_t[4];

  /* 3-vector */
  typedef double vec3_t[3];

  /* state structure for keeping track of transformation updates */
  typedef struct {
    double update_rate; // period in seconds
    double ctime_last;  // time of last update
  } qp_state_t;

  /* structure for storing refraction data */
  typedef struct {
    double temperature; // temperature, C
    double pressure;    // pressure, mbar
    double humidity;    // humidity, fraction
    double frequency;   // frequency, ghz
  } qp_weather_t;

  /* structures for storing Bulletin A data (for wobble correction) */
  typedef struct {
    float x;
    float y;
    float dut1;
  } qp_bulletina_entry_t;

  typedef struct {
    qp_bulletina_entry_t *entries;
    int mjd_min;
    int mjd_max;
  } qp_bulletina_t;

  /* parameter structure for storing corrections computed at variable rates */
  typedef struct qp_memory_t {
    int init;

    // update state
    qp_state_t state_daber;     // diurnal aberration
    qp_state_t state_lonlat;    // lat/lon
    qp_state_t state_wobble;    // polar motion
    qp_state_t state_dut1;      // ut1 correction
    qp_state_t state_erot;      // earth's rotation
    qp_state_t state_npb;       // nutation, precession, frame bias
    qp_state_t state_aaber;     // annual aberration
    qp_state_t state_ref;       // refraction

    // update inverse state
    qp_state_t state_daber_inv;     // diurnal aberration
    qp_state_t state_lonlat_inv;    // lat/lon
    qp_state_t state_wobble_inv;    // polar motion
    qp_state_t state_dut1_inv;      // ut1 correction
    qp_state_t state_erot_inv;      // earth's rotation
    qp_state_t state_npb_inv;       // nutation, precession, frame bias
    qp_state_t state_aaber_inv;     // annual aberration
    qp_state_t state_ref_inv;       // refraction

    // state data
    qp_weather_t weather;     // weather
    double ref_delta;         // refraction correction, deg
    quat_t q_ref;             // refraction quaternion
    quat_t q_ref_inv;         // inverse refraction quaternion
    double dut1;              // UT1 correction
    quat_t q_lonlat;          // lonlat quaternion
    quat_t q_lonlat_inv;      // inverse lonlat quaternion
    quat_t q_wobble;          // wobble quaternion
    quat_t q_wobble_inv;      // inverse wobble quaternion
    quat_t q_npb;             // nutation etc quaternion
    quat_t q_npb_inv;         // inverse nutation etc quaternion
    quat_t q_erot;            // earth's rotation quaternion
    quat_t q_erot_inv;        // inverse earth's rotation quaternion
    quat_t q_gal;             // galactic coordinates
    quat_t q_gal_inv;         // inverse of q_gal
    int gal_init;             // q_gal* initialized?
    vec3_t v_dipole;          // dipole direction
    int dipole_init;          // q_dipole initialized?
    vec3_t beta_earth;        // earth orbital velocity
    vec3_t beta_rot;          // earth rotational velocity
    qp_bulletina_t bulletinA; // bulletin A data

    // options
    int accuracy;          // 0=full accuracy, 1=low accuracy
    int mean_aber;         // 0=per-detector aberration, 1=mean
    int fast_math;         // 0=regular trig, 1=polynomial trig approximations
    int polconv;           // polarization convention (0=healpix,1=IAU)
    int pix_order;         // pixel ordering (1=nest, 0=ring)
    int interp_pix;        // interpolate between pixels in map2tod (1=yes, 0=no)
    int fast_pix;          // use vec2pix instead of ang2pix in binners
    int error_missing;     // raise an error when reading/writing missing pixels
    int nan_missing;       // set missing samples to NaN (used if !error_missing)
    int interp_missing;    // drop missing neighbors when interp_pix=1
    int num_threads;       // number of parallel threads
    int thread_num;        // current thread number

    // error handling
    int error_code;
    char *error_string;
  } qp_memory_t;

  /* parameter initialization */
  qp_memory_t * qp_init_memory(void);
  void qp_free_memory(qp_memory_t *mem);
  qp_memory_t * qp_copy_memory(qp_memory_t *memsrc);

  /* common update rates */
  extern const int QP_DO_ALWAYS;
  extern const int QP_DO_ONCE;
  extern const int QP_DO_NEVER;

  /* Set correction rates for each state, in seconds; control accuracy and speed
     Use above macros to allow states to be applied always, once or never. */
  void qp_set_rates(qp_memory_t *mem,
		    double daber_rate,
		    double lonlat_rate,
		    double wobble_rate,
		    double dut1_rate,
		    double erot_rate,
		    double npb_rate,
		    double aaber_rate,
		    double ref_rate);

  /* Set correction rates for each inverse state, in seconds; control accuracy and speed
     Use above macros to allow states to be applied always, once or never. */
  void qp_set_inv_rates(qp_memory_t *mem,
			double daber_rate,
			double lonlat_rate,
			double wobble_rate,
			double dut1_rate,
			double erot_rate,
			double npb_rate,
			double aaber_rate,
			double ref_rate);

  /* reset counters so that all corrections are recalculated at next sample */
  void qp_reset_rates(qp_memory_t *mem);
  void qp_reset_inv_rates(qp_memory_t *mem);

  /* Check whether a correction needs to be updated */
  int qp_check_update(qp_state_t *state, double ctime);

  /* check whether a correction needs to be applied */
  int qp_check_apply(qp_state_t *state);

  /* print stuff */
  void qp_print_vec3(const char *tag, vec3_t v);
  void qp_print_quat(const char *tag, quat_t q);
  void qp_print_state(const char *tag, qp_state_t *state);
  void qp_print_weather(qp_weather_t *w);
  void qp_print_vec3_mp(int thread, const char *tag, vec3_t v);
  void qp_print_quat_mp(int thread, const char *tag, quat_t q);
  void qp_print_state_mp(int thread, const char *tag, qp_state_t *state);
  void qp_print_weather_mp(int thread, qp_weather_t *w);
  void qp_print_memory(qp_memory_t *mem);

  /* per-correction functions */
#define RATEFUNC(state)                                         \
  void qp_set_rate_##state(qp_memory_t *mem, double rate);      \
  void qp_reset_rate_##state(qp_memory_t *mem);                 \
  double qp_get_rate_##state(qp_memory_t *mem);
  RATEFUNC(daber)
  RATEFUNC(lonlat)
  RATEFUNC(wobble)
  RATEFUNC(dut1)
  RATEFUNC(erot)
  RATEFUNC(npb)
  RATEFUNC(aaber)
  RATEFUNC(ref)
  RATEFUNC(daber_inv)
  RATEFUNC(lonlat_inv)
  RATEFUNC(wobble_inv)
  RATEFUNC(dut1_inv)
  RATEFUNC(erot_inv)
  RATEFUNC(npb_inv)
  RATEFUNC(aaber_inv)
  RATEFUNC(ref_inv)

  /* per-option functions */
  void qp_set_options(qp_memory_t *mem,
		      int accuracy,
		      int mean_aber,
		      int fast_math,
		      int polconv,
		      int pix_order,
                      int interp_pix,
                      int fast_pix,
                      int error_missing,
                      int nan_missing,
                      int interp_missing,
		      int num_threads);

#define OPTIONFUNC(opt)                                 \
  void qp_set_opt_##opt(qp_memory_t *mem, int val);     \
  int qp_get_opt_##opt(qp_memory_t *mem);
  OPTIONFUNC(accuracy);
  OPTIONFUNC(mean_aber);
  OPTIONFUNC(fast_math);
  OPTIONFUNC(polconv);
  OPTIONFUNC(pix_order);
  OPTIONFUNC(interp_pix);
  OPTIONFUNC(fast_pix);
  OPTIONFUNC(error_missing);
  OPTIONFUNC(nan_missing);
  OPTIONFUNC(interp_missing);
#ifndef ENABLE_LITE
  OPTIONFUNC(num_threads);
  OPTIONFUNC(thread_num);
#endif

  /* Set weather data */
  void qp_set_weather(qp_memory_t *mem, double temperature, double pressure,
		      double humidity, double frequency);

#define WEATHFUNC(param)                                        \
  void qp_set_weather_##param(qp_memory_t *mem, double val);    \
  double qp_get_weather_##param(qp_memory_t *mem);
  WEATHFUNC(temperature)
  WEATHFUNC(pressure)
  WEATHFUNC(humidity)
  WEATHFUNC(frequency)

#define DOUBLEFUNC(param)                               \
  void qp_set_##param(qp_memory_t *mem, double val);    \
  double qp_get_##param(qp_memory_t *mem);
  DOUBLEFUNC(ref_delta)
  DOUBLEFUNC(dut1)

  /* *************************************************************************
     Exception handling
     ********************************************************************** */

  typedef enum {
    QP_ERROR = 1,
    QP_ERROR_INIT,
    QP_ERROR_POINT,
    QP_ERROR_MAP
  } qp_error_codes;

  int qp_get_error_code(qp_memory_t *mem);
  char *qp_get_error_string(qp_memory_t *mem);
  void qp_set_error(qp_memory_t *mem, int error_code, const char *error_string);
  int qp_check_error(qp_memory_t *mem, int condition, int error_code,
                     const char *error_string);

  /* *************************************************************************
     Utility functions
     ********************************************************************** */

  /* diurnal aberration constant (radians) */
#define D_ABER_RAD 1.54716541e-06 // -0.3191 arcsec
  /* speed of light, AU/day */
#define C_AUD 173.14463269999999
  /* speed of light, m/s */
#define C_MS 299792458.0

  /* Return interpolated values from IERS Bulletin A */
  int qp_get_iers_bulletin_a( qp_memory_t *mem, double mjd,
                              double *dut1, double *x, double *y );
  /* Set IERS Bulletin A */
  int qp_set_iers_bulletin_a( qp_memory_t *mem, int mjd_min_, int mjd_max_,
                              double *dut1, double *x, double *y );
  /* Copy IERS Bulletin A */
  int qp_copy_iers_bulletin_a( qp_memory_t *memdest, qp_memory_t *memsrc );

  /* Time conversion */
#define CTIME_JD_EPOCH 2440587.5 /* JD for ctime = 0 */
  void ctime2jd(double ctime, double jd[2]);
  double jd2ctime(double jd[2]);
  void ctime2jdtt(double ctime, double jd_tt[2]);
  void jdutc2jdut1(double jd_utc[2], double dut1, double jd_ut1[2]);
  double ctime2gmst(double ctime, double dut1, int accuracy);
  static inline double secs2days( double s ) { return s/86400.; }
  static inline double days2secs( double d ) { return d*86400.; }
  static inline double jd2mjd( double jd ) { return jd - 2400000.5; }
  static inline double mjd2jd( double mjd ) { return mjd + 2400000.5; }

  /* Unit conversions */
#ifndef M_PI
#define M_PI		3.14159265358979323846	// pi
#define M_TWOPI         6.28318530717958647692  // 2*pi
#define M_PI_2		1.57079632679489661923	// pi/2
#endif
  static const double d2r = M_PI/180.;
  static const double r2d = 180./M_PI;
  static inline double deg2rad( double deg ) { return deg*d2r; }
  static inline double rad2deg( double rad ) { return rad*r2d; }
  static const double as2r = M_PI/(180.*3600.);
  static const double r2as = 3600.*180./M_PI;
  static inline double arcsec2rad( double sec ) { return sec*as2r; }
  static inline double rad2arcsec( double rad ) { return rad*r2as; }

  /* *************************************************************************
     Intermediate rotations and corrections
     ********************************************************************** */

  /* Calculate aberration correction quaternion
     v = (R(q)*z) x beta, angle = |v|, qa = quat(-angle,v) */
  void qp_aberration(quat_t q, vec3_t beta, quat_t qa, int inv);

  /* Calculate earth orbital velocity vector as fraction of speed of light */
  void qp_earth_orbital_beta(double jd_tdb[2], vec3_t beta);

  /* Apply annual aberration correction to given quaternion */
  void qp_apply_annual_aberration(qp_memory_t *mem, double ctime, quat_t q, int inv);

  /* Apply diurnal aberration correction to given quaternion */
  void qp_apply_diurnal_aberration(qp_memory_t *mem, double ctime, double lat,
                                   quat_t q, int inv);

  /* Calculate nutation/precession/bias correction quaternion
     use (faster) truncated series if accuracy > 0*/
  void qp_npb_quat(double jd_tt[2], quat_t q, int accuracy);

  /* Calculate ERA quaternion */
  void qp_erot_quat(double jd_ut1[2], quat_t q);

  /* Calcuate wobble correction quaternion */
  void qp_wobble_quat(double jd_tt[2], double xp, double yp, quat_t q);

  /* Calculate gondola orientation quaternion */
  void qp_azel_quat(double az, double el, double pitch, double roll, quat_t q);

  /* Calculate gondola orientation quaternion, accounting for FPU boresight rotation */
  void qp_azelpsi_quat(double az, double el, double psi, double pitch, double roll, quat_t q);

  /* Calculate longitude/latitude quaternion */
  void qp_lonlat_quat(double lon, double lat, quat_t q);

  /* Calculate Greenwich mean sidereal time in hours */
  double qp_gmst(qp_memory_t *mem, double ctime);

  /* Calculate Greenwich mean sidereal time in hours */
  void qp_gmstn(qp_memory_t *mem, double *ctime, double *gmst, int n);

  /* Calculate local mean sidereal time in hours */
  double qp_lmst(qp_memory_t *mem, double ctime, double lon);

  /* Calculate local mean sidereal time in hours */
  void qp_lmstn(qp_memory_t *mem, double *ctime, double *lon, double *lmst,
		int n);

  /* Calculate dipole amplitude(s) from quaternion */
  double qp_quat2dipole(qp_memory_t *mem, double ctime, quat_t q);
  void qp_bore2dipole(qp_memory_t *mem, quat_t q_off, double *ctime,
                        quat_t *q_bore, double *dipole, int n);

  /* Calculate dipole amplitude from ra/dec */
  double qp_dipole(qp_memory_t *mem, double ctime, double ra, double dec);

  /* Calculate dipole amplitudes from ra/dec */
  void qp_dipolen(qp_memory_t *mem, double *ctime, double *ra, double *dec,
                  double *dipole, int n);

  /* Calculate waveplate quaternion, given _physical_ HWP angle */
  void qp_hwp_quat(double ang, quat_t q);

  /* Calculate waveplate quaternions */
  void qp_hwp_quatn(double *ang, quat_t *q, int n);

  /* Calculate atmospheric refraction */
  double qp_refraction(double el, double temp, double press, double hum,
                       double freq);

  /* Update atmospheric refraction using stored parameters */
  double qp_update_ref(qp_memory_t *mem, quat_t q);

  /* Apply refraction correction */
  void qp_apply_refraction(qp_memory_t *mem, double ctime, quat_t q, int inv);

  /* *************************************************************************
     Output functions
     ********************************************************************** */

  /* Compute boresight quaternion for a single gondola orientation. */
  void qp_azel2quat(qp_memory_t *mem, double az, double el, double pitch,
		    double roll, double lon, double lat, double ctime,
		    quat_t q);

  /* Compute boresight quaternions for n gondola orientations. */
  void qp_azel2bore(qp_memory_t *mem, double *az, double *el, double *pitch,
		    double *roll, double *lon, double *lat, double *ctime,
		    quat_t *q, int n);

  /* Compute boresight quaternion for a single gondola orientation. Accounts for FPU boresight rotation. */
  void qp_azelpsi2quat(qp_memory_t *mem, double az, double el, double psi, double pitch,
		    double roll, double lon, double lat, double ctime,
		    quat_t q);

  /* Compute boresight quaternions for n gondola orientations. Accounts for FPU boresight rotation. */
  void qp_azelpsi2bore(qp_memory_t *mem, double *az, double *el, double *psi, double *pitch,
		    double *roll, double *lon, double *lat, double *ctime,
		    quat_t *q, int n);

  /* Compute horizon coordinates for a given quaternion in equatorial coordinates */
  void qp_quat2azel(qp_memory_t *mem, quat_t q, double lon, double lat,
		    double ctime, double *az, double *el, double *hpa);

  /* Compute horizon coordinates for n sets boresight equatorial coordinates */
  void qp_bore2azel(qp_memory_t *mem, quat_t *q, double *lon, double *lat,
		    double *ctime, double *az, double *el, double *pa, int n);

  /* Compute detector offset quaternion. */
  void qp_det_offset(double delta_az, double delta_el, double delta_psi,
		     quat_t q);


  /* Compute detector offset quaternion. */
  void qp_det_offsetn(double *delta_az, double *delta_el, double *delta_psi,
		      quat_t *q, int n);

  /* Adjust boresight pointing by a variable quaternion */
  void qp_bore_offset(qp_memory_t *mem, quat_t *q_bore, double *ang1, double *ang2,
                      double *ang3, int n, int post);

  /* Compute ra/dec/pa from quaternion */
  void qp_quat2radecpa(qp_memory_t *mem, quat_t q, double *ra, double *dec,
                       double *pa);
  void qp_quat2radecpan(qp_memory_t *mem, quat_t *q, double *ra, double *dec,
                        double *pa, int n);

  /* Compute ra/dec and sin(2*psi)/cos(2*psi) for a given quaternion */
  void qp_quat2radec(qp_memory_t *mem, quat_t q, double *ra, double *dec,
		     double *sin2psi, double *cos2psi);

  /* Compute quaternion given ra/dec and sin(2*psi)/cos(2*psi) */
  void qp_radec2quat(qp_memory_t *mem, double ra, double dec,
                     double sin2psi, double cos2psi, quat_t q);

  /* Compute quaternion given ra/dec/pa */
  void qp_radecpa2quat(qp_memory_t *mem, double ra, double dec,
                       double pa, quat_t q);

  /* Compute quaternions given ra/dec/pa */
  void qp_radecpa2quatn(qp_memory_t *mem, double *ra, double *dec,
                        double *pa, quat_t *q, int n);

  /* Compute ra/sin(dec) and sin(2*psi)/cos(2*psi) for a given quaternion */
  void qp_quat2rasindec(qp_memory_t *mem, quat_t q, double *ra, double *sindec,
			double *sin2psi, double *cos2psi);

  /* Calculate the detector quaternion from the boresight and offset. */
  void qp_bore2det(qp_memory_t *mem, quat_t q_off, double ctime, quat_t q_bore,
		   quat_t q_det);

  /* Calculate the detector quaternion from the boresight, offset and HWP angle. */
  void qp_bore2det_hwp(qp_memory_t *mem, quat_t q_off, double ctime, quat_t q_bore,
		       quat_t q_hwp, quat_t q_det);

  /* Calculate ra/dec and sin(2*psi)/cos(2*psi) for a given detector offset,
     from an array of boresight quaternions. */
  void qp_bore2radec(qp_memory_t *mem, quat_t q_off, double *ctime, quat_t *q_bore,
		     double *ra, double *dec, double *sin2psi, double *cos2psi,
		     int n);

  /* Calculate ra/dec and sin(2*psi)/cos(2*psi) for a given detector offset,
     from an array of boresight and waveplate quaternions. */
  void qp_bore2radec_hwp(qp_memory_t *mem, quat_t q_off, double *ctime,
			 quat_t *q_bore, quat_t *q_hwp, double *ra, double *dec,
			 double *sin2psi, double *cos2psi, int n);

  /* Calculate ra/sin(dec) and sin(2*psi)/cos(2*psi) for a given detector offset. */
  void qp_bore2rasindec(qp_memory_t *mem, quat_t q_off, double *ctime, quat_t *q_bore,
			double *ra, double *sindec, double *sin2psi, double *cos2psi,
			int n);

  /* Calculate ra/sin(dec) and sin(2*psi)/cos(2*psi) for a given detector offset. */
  void qp_bore2rasindec_hwp(qp_memory_t *mem, quat_t q_off, double *ctime,
			    quat_t *q_bore, quat_t *q_hwp, double *ra, double *sindec,
			    double *sin2psi, double *cos2psi, int n);

  /* Calculate ra/dec/pa for a given detector offset,
     from an array of boresight quaternions. */
  void qp_bore2radecpa(qp_memory_t *mem, quat_t q_off, double *ctime, quat_t *q_bore,
                       double *ra, double *dec, double *pa, int n);

  /* Calculate ra/dec/pa for a given detector offset,
     from an array of boresight and waveplate quaternions. */
  void qp_bore2radecpa_hwp(qp_memory_t *mem, quat_t q_off, double *ctime,
                           quat_t *q_bore, quat_t *q_hwp, double *ra, double *dec,
                           double *pa, int n);

  /* Calculate ra/dec and sin(2*psi)/cos(2*psi) for a given detector offset,
     from a set of boresight orientations. */
  void qp_azelpsi2radec(qp_memory_t *mem,
		     double delta_az, double delta_el, double delta_psi,
		     double *az, double *el, double *psi, double *pitch, double *roll,
		     double *lon, double *lat, double *ctime,
		     double *ra, double *dec, double *sin2psi, double *cos2psi,
		     int n);

  /* Calculate ra/dec and sin(2*psi)/cos(2*psi) for a given detector offset,
     from a set of boresight orientations. */
  void qp_azel2radec(qp_memory_t *mem,
		     double delta_az, double delta_el, double delta_psi,
		     double *az, double *el, double *pitch, double *roll,
		     double *lon, double *lat, double *ctime,
		     double *ra, double *dec, double *sin2psi, double *cos2psi,
		     int n);

  /* Calculate ra/dec/pa for a given detector offset,
     from a set of boresight orientations. */
  void qp_azelpsi2radecpa(qp_memory_t *mem,
		       double delta_az, double delta_el, double delta_psi,
		       double *az, double *el, double *psi, double *pitch, double *roll,
		       double *lon, double *lat, double *ctime,
		       double *ra, double *dec, double *pa, int n);

  /* Calculate ra/dec/pa for a given detector offset,
     from a set of boresight orientations. */
  void qp_azel2radecpa(qp_memory_t *mem,
		       double delta_az, double delta_el, double delta_psi,
		       double *az, double *el, double *pitch, double *roll,
		       double *lon, double *lat, double *ctime,
		       double *ra, double *dec, double *pa, int n);

  /* Calculate az/el/pa for a set of equatorial coordinates */
  void qp_radec2azel(qp_memory_t *mem,
		     double *ra, double *dec, double *pa, double *lon,
		     double *lat, double *ctime, double *az, double *el,
		     double *hpa, int n);

  /* Calculate ra/dec and sin(2*psi)/cos(2*psi) for a given detector offset,
     from a set of boresight orientations and waveplate angles. */
  void qp_azelpsi2radec_hwp(qp_memory_t *mem,
			 double delta_az, double delta_el, double delta_psi,
			 double *az, double *el, double *psi, double *pitch, double *roll,
			 double *lon, double *lat, double *ctime, double *hwp,
			 double *ra, double *dec, double *sin2psi, double *cos2psi,
			 int n);

  /* Calculate ra/dec and sin(2*psi)/cos(2*psi) for a given detector offset,
     from a set of boresight orientations and waveplate angles. */
  void qp_azel2radec_hwp(qp_memory_t *mem,
			 double delta_az, double delta_el, double delta_psi,
			 double *az, double *el, double *pitch, double *roll,
			 double *lon, double *lat, double *ctime, double *hwp,
			 double *ra, double *dec, double *sin2psi, double *cos2psi,
			 int n);

  /* Calculate ra/dec and sin(2*psi)/cos(2*psi) for a given detector offset,
     from a set of boresight orientations and waveplate angles. */
  void qp_azelpsi2radecpa_hwp(qp_memory_t *mem,
			   double delta_az, double delta_el, double delta_psi,
			   double *az, double *el, double *psi, double *pitch, double *roll,
			   double *lon, double *lat, double *ctime, double *hwp,
			   double *ra, double *dec, double *pa, int n);

  /* Calculate ra/dec and sin(2*psi)/cos(2*psi) for a given detector offset,
     from a set of boresight orientations and waveplate angles. */
  void qp_azel2radecpa_hwp(qp_memory_t *mem,
			   double delta_az, double delta_el, double delta_psi,
			   double *az, double *el, double *pitch, double *roll,
			   double *lon, double *lat, double *ctime, double *hwp,
			   double *ra, double *dec, double *pa, int n);

  /* Calculate ra/sin(dec) and sin(2*psi)/cos(2*psi) for a given detector offset,
     from a set of boresight orientations.  */
  void qp_azelpsi2rasindec(qp_memory_t *mem,
			double delta_az, double delta_el, double delta_psi,
			double *az, double *el, double *psi, double *pitch, double *roll,
			double *lon, double *lat, double *ctime,
			double *ra, double *sindec, double *sin2psi, double *cos2psi,
			int n);

  /* Calculate ra/sin(dec) and sin(2*psi)/cos(2*psi) for a given detector offset,
     from a set of boresight orientations.  */
  void qp_azel2rasindec(qp_memory_t *mem,
			double delta_az, double delta_el, double delta_psi,
			double *az, double *el, double *pitch, double *roll,
			double *lon, double *lat, double *ctime,
			double *ra, double *sindec, double *sin2psi, double *cos2psi,
			int n);

  /* Calculate ra/sin(dec) and sin(2*psi)/cos(2*psi) for a given detector offset,
     from a set of boresight orientations.  */
  void qp_azel2rasindec_hwp(qp_memory_t *mem,
			    double delta_az, double delta_el, double delta_psi,
			    double *az, double *el, double *pitch, double *roll,
			    double *lon, double *lat, double *ctime, double *hwp,
			    double *ra, double *sindec, double *sin2psi,
			    double *cos2psi, int n);

  /* Calculate ra/sin(dec) and sin(2*psi)/cos(2*psi) for a given detector offset,
     from a set of boresight orientations.  */
  void qp_azelpsi2rasindec_hwp(qp_memory_t *mem,
			    double delta_az, double delta_el, double delta_psi,
			    double *az, double *el, double *psi, double *pitch, double *roll,
			    double *lon, double *lat, double *ctime, double *hwp,
			    double *ra, double *sindec, double *sin2psi,
			    double *cos2psi, int n);

#ifndef ENABLE_LITE

  /* *************************************************************************
     Pixelization
     ********************************************************************** */

  /* Ordering */
#define QP_ORDER_NEST 1
#define QP_ORDER_RING 0

  /* Compute healpix pixel number for given nside and ra/dec */
  long qp_radec2pix(qp_memory_t *mem, double ra, double dec, int nside);

  /* Compute healpix pixel number for given nside and ra/dec */
  void qp_radec2pixn(qp_memory_t *mem, double *ra, double *dec, int nside,
                     long *pix, int n);

  /* Compute pixel number and pol angle for given nside and quaternion */
  void qp_quat2pix(qp_memory_t *mem, quat_t q, int nside, long *pix,
                   double *sin2psi, double *cos2psi);
  void qp_quat2pixpa(qp_memory_t *mem, quat_t q, int nside, long *pix,
                     double *pa);

  /* Compute pixel numbers and pol angles for given nside and quaternions */
  void qp_quat2pixn(qp_memory_t *mem, quat_t *q, int nside, long *pix,
                    double *sin2psi, double *cos2psi, int n);
  void qp_quat2pixpan(qp_memory_t *mem, quat_t *q, int nside, long *pix,
                      double *pa, int n);

  /* Rotate from celestial to galactic coordinates */
  void qp_radec2gal_quat(qp_memory_t *mem, quat_t q);
  void qp_radec2gal_quatn(qp_memory_t *mem, quat_t *q, int n);
  void qp_radec2gal(qp_memory_t *mem, double *ra, double *dec,
                    double *sin2psi, double *cos2psi);
  void qp_radec2galn(qp_memory_t *mem, double *ra, double *dec,
                     double *sin2psi, double *cos2psi, int n);
  void qp_radecpa2gal(qp_memory_t *mem, double *ra, double *dec, double *pa);
  void qp_radecpa2galn(qp_memory_t *mem, double *ra, double *dec,
                       double *pa, int n);

  /* Rotate from galactic to celestial coordinates */
  void qp_gal2radec_quat(qp_memory_t *mem, quat_t q);
  void qp_gal2radec_quatn(qp_memory_t *mem, quat_t *q, int n);
  void qp_gal2radec(qp_memory_t *mem, double *ra, double *dec,
                    double *sin2psi, double *cos2psi);
  void qp_gal2radecn(qp_memory_t *mem, double *ra, double *dec,
                     double *sin2psi, double *cos2psi, int n);
  void qp_gal2radecpa(qp_memory_t *mem, double *ra, double *dec, double *pa);
  void qp_gal2radecpan(qp_memory_t *mem, double *ra, double *dec,
                       double *pa, int n);

  /* Rotate a TQU map from one coordinate system to another */
  void qp_rotate_map(qp_memory_t *mem, int nside,
                     double **map_in, const char coord_in,
                     double **map_out, const char coord_out);

  /* Compute pix/pol timestreams for given boresight timestream and detector
     offset */
  void qp_bore2pix(qp_memory_t *mem, quat_t q_off, double *ctime, quat_t *q_bore,
                   int nside, long *pix, double *sin2psi, double *cos2psi, int n);
  void qp_bore2pixpa(qp_memory_t *mem, quat_t q_off, double *ctime, quat_t *q_bore,
                     int nside, long *pix, double *pa, int n);

  /* Compute pix/pol timestreams for given boresight timestream,
     waveplate timestream and detector offset */
  void qp_bore2pix_hwp(qp_memory_t *mem, quat_t q_off, double *ctime,
                       quat_t *q_bore, quat_t *q_hwp, int nside, long *pix,
                       double *sin2psi, double *cos2psi, int n);
  void qp_bore2pixpa_hwp(qp_memory_t *mem, quat_t q_off, double *ctime,
                         quat_t *q_bore, quat_t *q_hwp, int nside, long *pix,
                         double *pa, int n);

  /* *************************************************************************
     Mapmaking and Projection
     ********************************************************************** */

  typedef enum {
    QP_STRUCT_INIT = 0x01,     // structure initialized
    QP_STRUCT_MALLOC = 0x02,   // malloc'd structure
    QP_ARR_INIT_PTR = 0x04,    // pointer to array
    QP_ARR_MALLOC_1D = 0x08,   // top dimension malloc'd
    QP_ARR_MALLOC_2D = 0x10,   // second dimension malloc'd
  } qp_init_mode;

  typedef struct {
    int init;        // initialized?
    quat_t q_off;    // offset quaternion
    double weight;   // det weight
    double gain;     // det gain
    mueller_t mueller;  // mueller matrix parameters

    size_t n;        // samples

    int tod_init;    // tod initialized?
    double *tod;     // tod array

    int flag_init;   // flag initialized?
    uint8_t *flag;   // flag array

    int weights_init;   // weight tod initialized?
    double *weights;    // weight tod array
  } qp_det_t;

  typedef struct {
    int init;        // initialized?
    size_t n;        // number of dets
    int arr_init;    // array init?
    size_t diff;     //do differencing?
    qp_det_t *arr;   // array

  } qp_detarr_t;

  typedef struct {
    int init;        // initialized?
    size_t n;        // number of elements

    int q_bore_init; // q_bore init?
    quat_t *q_bore;  // boresight quaternion

    int ctime_init;  // ctime array set?
    double *ctime;   // ctime array

    int q_hwp_init;  // q_hwp array set?
    quat_t *q_hwp;   // hwp quaternion
  } qp_point_t;

  /* map type enum */
  typedef enum {
    QP_VEC_NONE = 0,  // no map
    QP_VEC_TEMP,      // unpolarized (T-only)
    QP_VEC_POL,       // polarized
    QP_VEC_VPOL,      // polarized + V-pol
    QP_VEC_D1,        // unpolarized + 1st derivs
    QP_VEC_D1_POL,    // polarized + 1st derivs
    QP_VEC_D2,        // unpolarized + 2nd derivs
    QP_VEC_D2_POL     // polarized + 2nd derivs
  } qp_vec_mode;

  /* projection type enum */
  typedef enum {
    QP_PROJ_NONE = 0, // no projection
    QP_PROJ_TEMP,     // unpolarized (hits-only)
    QP_PROJ_POL,      // polarized
    QP_PROJ_VPOL      // polarized + V-pol
  } qp_proj_mode;

  typedef struct {
    int init;
    long idx;
    long startpix;
    long ringpix;
    double theta;
    int shifted;
  } qp_ring_t;

  typedef struct {
    int init;
    int nside;
    long npix;
    long npface;
    long ncap;
    double fact1;
    double fact2;
    int rings_init;
    qp_ring_t *rings;
  } qp_pixinfo_t;

  typedef struct {
    long key;
    long index;
  } qp_pix_pair_t;

  typedef struct {
    size_t count;
    qp_pix_pair_t *pairs;
  } qp_pix_bucket_t;

  typedef struct {
    int init;
    size_t count;
    qp_pix_bucket_t *buckets;
  } qp_pixhash_t;

  typedef struct {
    int init;                // initialized?
    int partial;             // partial map?
    size_t nside;            // map nside
    size_t npix;             // map npix

    int pixinfo_init;        // pix info initialized?
    qp_pixinfo_t *pixinfo;   // pixel info structure

    int pixhash_init;        // pix hash initialized?
    qp_pixhash_t *pixhash;   // repixelization hash table

    size_t num_vec;          // number of map columns
    qp_vec_mode vec_mode;    // map mode
    int vec1d_init;          // vec1d initialized?
    double *vec1d;           // 1d map array
    int vec_init;            // vec initialized?
    double **vec;            // 2d map array

    size_t num_proj;         // number of projection columns
    qp_proj_mode proj_mode;  // projection mode
    int proj1d_init;         // proj1d initialized?
    double *proj1d;          // 1d proj array
    int proj_init;           // proj initialized?
    double **proj;           // projection array
  } qp_map_t;

  /* initialize detectors */
  qp_det_t * qp_init_det(quat_t q_off, double weight, double gain, mueller_t mueller);
  qp_det_t * qp_default_det(void);
  void qp_init_det_tod(qp_det_t *det, size_t n);
  void qp_init_det_tod_from_array(qp_det_t *det, double *tod, size_t n, int copy);
  void qp_init_det_flag(qp_det_t *det, size_t n);
  void qp_init_det_flag_from_array(qp_det_t *det, uint8_t *flag, size_t n, int copy);
  void qp_init_det_weights(qp_det_t *det, size_t n);
  void qp_init_det_weights_from_array(qp_det_t *det, double *weights, size_t n,
                                      int copy);
  void qp_free_det(qp_det_t *det);
  qp_detarr_t * qp_init_detarr(quat_t *q_off, double *weight, double *gain,
                               mueller_t *mueller, size_t n);
  void qp_init_detarr_tod(qp_detarr_t *dets, size_t n);
  void qp_init_detarr_tod_from_array(qp_detarr_t *dets, double **tod,
                                     size_t n, int copy);
  void qp_init_detarr_flag(qp_detarr_t *dets, size_t n);
  void qp_init_detarr_flag_from_array(qp_detarr_t *dets, uint8_t **flag,
                                      size_t n, int copy);
  void qp_init_detarr_weights(qp_detarr_t *dets, size_t n);
  void qp_init_detarr_weights_from_array(qp_detarr_t *dets, double **weights,
                                         size_t n, int copy);
  void qp_free_detarr(qp_detarr_t *dets);

  /* initialize pointing */
  qp_point_t * qp_init_point(size_t n, int time, int pol);
  qp_point_t *qp_init_point_from_arrays(quat_t *q_bore, double *ctime, quat_t *q_hwp,
                                        size_t n, int copy);
  void qp_free_point(qp_point_t *pnt);

  /* map repixelization */
  qp_pixhash_t * qp_init_pixhash(long *pix, size_t npix);
  qp_pixhash_t * qp_copy_pixhash(qp_pixhash_t *pixhash);
  void qp_free_pixhash(qp_pixhash_t *pixhash);
  long qp_repixelize(qp_pixhash_t *pixhash, long pix);
  int qp_init_map_pixhash(qp_map_t *map, long *pix, size_t npix);

  /* initialize maps */
  qp_map_t * qp_init_map(size_t nside, size_t npix, qp_vec_mode vec_mode,
                         qp_proj_mode proj_mode);
  qp_map_t * qp_init_map_from_arrays(double **vec, double **proj, size_t nside,
                                     size_t npix, qp_vec_mode vec_mode,
                                     qp_proj_mode proj_mode, int copy);
  qp_map_t * qp_init_map_from_map(qp_map_t *map, int blank, int copy);
  void qp_free_map(qp_map_t *map);

  /* tod -> map */
  int qp_add_map(qp_memory_t *mem, qp_map_t *map, qp_map_t *maploc);
  int qp_tod2map1(qp_memory_t *mem, qp_det_t *det, qp_point_t *pnt, qp_map_t *map);
  int qp_tod2map(qp_memory_t *mem, qp_detarr_t *dets, qp_point_t *pnt,
                 qp_map_t *map);

  /* Compute theta/phi offset from pixel center.  For interpolating maps. */
  void qp_pixel_offset(qp_memory_t *mem, int nside, long pix, double ra,
                       double dec, double *dtheta, double *dphi);

  /* C implementation of healpix-cxx/healpy get_interp_val method */
  qp_pixinfo_t * qp_init_pixinfo(size_t nside, int populate);
  void qp_free_pixinfo(qp_pixinfo_t *pixinfo);
  int qp_init_map_pixinfo(qp_map_t *map);
  void qp_get_interpol(qp_memory_t *mem, qp_pixinfo_t *pixinfo, double ra,
                       double dec, long pix[4], double weight[4]);
  double qp_get_interp_val(qp_memory_t *mem, qp_pixinfo_t *pixinfo, double *map,
                           double ra, double dec);
  void qp_get_interp_valn(qp_memory_t *mem, int nside, double *map,
                          double *ra, double *dec, double *val, int n);

  /* map -> tod */
  int qp_map2tod1(qp_memory_t *mem, qp_det_t *det, qp_point_t *pnt, qp_map_t *map);
  int qp_map2tod(qp_memory_t *mem, qp_detarr_t *dets, qp_point_t *pnt,
                 qp_map_t *map);

#endif // ENABLE_LITE

#ifdef __cplusplus
}
#endif
