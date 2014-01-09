#pragma once

#ifdef __cplusplus
extern "C" {
#endif

  /* Quaternion */
  typedef double quat_t[4];
  
  /* 3-vector */
  typedef double vec3_t[3];
  
  /* ************************************************************************* 
     Intermediate rotations and corrections
     ********************************************************************** */
  
  /* Calculate aberration correction quaternion
     v = (R(q)*z) x beta, angle = |v|, qa = quat(-angle,v) */
  void qp_aberration(quat_t q, vec3_t beta, quat_t qa);
  
  /* Calculate earth orbital velocity vector as fraction of speed of light */
  void qp_earth_orbital_beta(double jd_tdb[2], vec3_t beta);
  
  /* Apply annual aberration correction to given quaternion */
  void qp_apply_annual_aberration(double ctime, quat_t q);
  
  /* Calculate nutation/precession/bias correction quaternion
     use (faster) truncated series if accuracy > 0*/
  void qp_npb_quat(double jd_tt[2], quat_t q, int accuracy);
  
  /* Calculate ERA quaternion */
  void qp_erot_quat(double jd_ut1[2], quat_t q);
  
  /* Calcuate wobble correction quaternion */
  void qp_wobble_quat(double mjd_utc, double *dut1, quat_t q);
  
  /* Calculate gondola orientation quaternion */
  void qp_azel_quat(double az, double el, double pitch, double roll, quat_t q);
  
  /* Calculate longitude/latitude quaternion */
  void qp_lonlat_quat(double lon, double lat, quat_t q);

  /* ************************************************************************* 
     Internal parameter settings
     ********************************************************************** */
  
  /* step structure for keeping track of transformation updates */
  typedef struct qp_step_t {
    double update_rate; // period in seconds
    double ctime_last; // time of last update
  } qp_step_t;
  
  /* parameter structure for storing corrections computed at variable rates */
  typedef struct qp_params_t {
    int initialized;
    qp_step_t s_daber;  // diurnal aberration
    qp_step_t s_lonlat; // lat/lon
    qp_step_t s_wobble; // polar motion
    qp_step_t s_dut1;   // ut1 correction
    qp_step_t s_erot;   // earth's rotation
    qp_step_t s_npb;    // nutation, precession, frame bias    
    qp_step_t s_aaber;  // annual aberration
    int accuracy;       // 0=full accuracy, 1=low accuracy
    int mean_aber;      // 0=per-detector aberration, 1=mean
    int fast_math;      // 0=regular trig, 1=polynomial trig approximations
    int polconv;        // polarization convention (0=healpix,1=IAU)
  } qp_params_t;
  extern qp_params_t qp_params;
  
  /* parameter initialization */
  void qp_init_params(void);
  
  /* common update rates */
#define QP_DO_ALWAYS 0
#define QP_DO_ONCE   -1
#define QP_DO_NEVER  -999
  
  /* Set correction rates for each step, in seconds; control accuracy and speed
     Use above macros to allow steps to be applied always, once or never. */
  void qp_set_params(double daber_rate,
		     double lonlat_rate,
		     double wobble_rate,
		     double dut1_rate,
		     double erot_rate,
		     double npb_rate,
		     double aaber_rate,
		     int accuracy,
		     int mean_aber,
		     int fast_math,
		     int polconv);
  
  /* reset counters so that all steps are recalculated at next sample */
  void qp_reset_all(void);
  
  /* Check whether a step needs to be updated */
  int qp_check_update(qp_step_t *step, double ctime);
  
  /* check whether a step needs to be applied */
  int qp_check_apply(qp_step_t *step);
  
  /* print quaternion to screen */
  void qp_print_debug(const char *tag, quat_t q);
  
  /* per-parameter functions */
#define PARAMFUNC(param)	\
  void qp_set_##param(int val);	\
  int qp_get_##param(void);
  PARAMFUNC(accuracy);
  PARAMFUNC(mean_aber);
  PARAMFUNC(fast_math);
  PARAMFUNC(polconv);
  
  /* per-step functions */
#define STEPFUNC(step)		    \
  void qp_set_##step(double rate);  \
  void qp_reset_##step(void);	    \
  double qp_get_##step(void);
  STEPFUNC(daber)
  STEPFUNC(lonlat)
  STEPFUNC(wobble)
  STEPFUNC(dut1)
  STEPFUNC(erot)
  STEPFUNC(npb)
  STEPFUNC(aaber)
  
  /* ************************************************************************* 
     Utility functions
     ********************************************************************** */
  
  /* diurnal aberration constant (radians) */
#define D_ABER_RAD 1.430408829156e-06 // -0.295043 arcsec
  /* speed of light, AU/day */
#define C_AUD 173.14463269999999 

  /* Return interpolated values from IERS Bulletin A */
  int get_iers_bulletin_a( double mjd, double *dut1, double *x, double *y );
  
  /* Time conversion */
#define CTIME_JD_EPOCH 2440587.5 /* JD for ctime = 0 */
  void ctime2jd(double ctime, double jd[2]);
  double jd2ctime(double jd[2]);
  void ctime2jdtt(double ctime, double jd_tt[2]);
  void jdutc2jdut1(double jd_utc[2], double dut1, double jd_ut1[2]);
  static inline double secs2days( double s ) { return s/86400.; }
  static inline double days2secs( double d ) { return d*86400.; }
  static inline double jd2mjd( double jd ) { return jd - 2400000.5; }
  static inline double mjd2jd( double mjd ) { return mjd + 2400000.5; }
  
  /* Unit conversions */
#ifndef M_PI
#define M_PI		3.14159265358979323846	// pi
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
     Output functions
     ********************************************************************** */
  
  /* Compute boresight quaternion for a single gondola orientation.
     Use low-accuracy NPB correction if accuracy > 0.
     Apply a mean annual aberration correction to the boresight quaternion
     if mean_aber > 0 (otherwise apply to each detector later) */
  void qp_azel2quat(double az, double el, double pitch, double roll,
		    double lon, double lat, double ctime, quat_t q);
  
  /* Compute boresight quaternions for n gondola orientations. */
  void qp_azel2bore(double *az, double *el, double *pitch, double *roll,
		    double *lon, double *lat, double *ctime, quat_t *q, int n);
  
  /* Compute detector offset quaternion. */
  void qp_det_offset(double delta_az, double delta_el, double delta_psi,
		     quat_t q);
  
  /* Compute ra/dec and sin(2*psi)/cos(2*psi) for a given quaternion */
  void qp_quat2radec(quat_t q, double *ra, double *dec, double *sin2psi,
		     double *cos2psi);
  
  /* Compute ra/sin(dec) and sin(2*psi)/cos(2*psi) for a given quaternion */
  void qp_quat2rasindec(quat_t q, double *ra, double *sindec, double *sin2psi,
			double *cos2psi);
  
  /* Calculate the detector quaternion from the boresight and offset. */
  void qp_bore2det(quat_t q_off, double ctime, quat_t q_bore, quat_t q_det);
  
  /* Calculate ra/dec and sin(2*psi)/cos(2*psi) for a given detector offset,
     from an array of boresight quaternions. */
  void qp_bore2radec(double delta_az, double delta_el, double delta_psi,
		     double *ctime, quat_t *q_bore,
		     double *ra, double *dec, double *sin2psi, double *cos2psi,
		     int n);
  
  /* Calculate ra/sin(dec) and sin(2*psi)/cos(2*psi) for a given detector offset. */
  void qp_bore2rasindec(double delta_az, double delta_el, double delta_psi,
			double *ctime, quat_t *q_bore,
			double *ra, double *sindec, double *sin2psi, double *cos2psi,
			int n);
  
  /* Calculate ra/dec and sin(2*psi)/cos(2*psi) for a given detector offset,
     from a set of boresight orientations. */
  void qp_azel2radec(double delta_az, double delta_el, double delta_psi,
		     double *az, double *el, double *pitch, double *roll,
		     double *lon, double *lat, double *ctime, 
		     double *ra, double *dec, double *sin2psi, double *cos2psi,
		     int n);
  
  /* Calculate ra/sin(dec) and sin(2*psi)/cos(2*psi) for a given detector offset,
     from a set of boresight orientations.  */
  void qp_azel2rasindec(double delta_az, double delta_el, double delta_psi,
			double *az, double *el, double *pitch, double *roll,
			double *lon, double *lat, double *ctime, 
			double *ra, double *sindec, double *sin2psi, double *cos2psi,
			int n);
  
#ifdef __cplusplus
}
#endif
