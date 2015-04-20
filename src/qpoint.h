#pragma once

#ifdef __cplusplus
extern "C" {
#endif

  /* Quaternion */
  typedef double quat_t[4];
  
  /* 3-vector */
  typedef double vec3_t[3];
  
  /* ************************************************************************* 
     Internal parameter settings
     ********************************************************************** */
  
  /* state structure for keeping track of transformation updates */
  typedef struct {
    double update_rate; // period in seconds
    double ctime_last;  // time of last update
  } qp_state_t;
  
  /* structure for storing refraction data */
  typedef struct {
    double height;      // height, m
    double temperature; // temperature, C
    double pressure;    // pressure, mbar
    double humidity;    // humidity, fraction
    double frequency;   // frequency, ghz
    double lapse_rate;  // tropospheric lapse rate, K/m
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
    int initialized;
    
    // update state
    qp_state_t state_daber;     // diurnal aberration
    qp_state_t state_lonlat;    // lat/lon
    qp_state_t state_wobble;    // polar motion
    qp_state_t state_dut1;      // ut1 correction
    qp_state_t state_erot;      // earth's rotation
    qp_state_t state_npb;       // nutation, precession, frame bias    
    qp_state_t state_aaber;     // annual aberration
    qp_state_t state_ref;       // refraction
    
    // state data
    qp_weather_t weather;     // weather
    double ref_tol;           // refraction tolerance, rad
    double ref_delta;         // refraction correction, deg
    quat_t q_ref;             // refraction quaternion
    double dut1;              // UT1 correction
    quat_t q_lonlat;          // lonlat quaternion
    quat_t q_wobble;          // wobble quaternion
    quat_t q_npb;             // nutation etc quaternion
    quat_t q_erot;            // earth's rotation quaternion
    quat_t q_gal;             // galactic coordinates
    quat_t q_gal_inv;         // inverse of q_gal
    int gal_init;             // q_gal* initialized?
    vec3_t beta_earth;        // earth orbital velocity
    vec3_t beta_rot;          // earth rotational velocity
    qp_bulletina_t bulletinA; // bulletin A data
    
    // options
    int accuracy;          // 0=full accuracy, 1=low accuracy
    int mean_aber;         // 0=per-detector aberration, 1=mean
    int fast_math;         // 0=regular trig, 1=polynomial trig approximations
    int polconv;           // polarization convention (0=healpix,1=IAU)
    int pair_dets;         // pair A/B detectors for bore2map (1=True, 0=False)
    int pix_order;         // pixel ordering (1=nest, 0=ring)
    int fast_pix;          // use vec2pix instead of ang2pix in binners
    int num_threads;       // number of parallel threads
    int thread_num;        // current thread number
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
  
  /* reset counters so that all corrections are recalculated at next sample */
  void qp_reset_rates(qp_memory_t *mem);
  
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
  
  /* per-option functions */
  void qp_set_options(qp_memory_t *mem,
		      int accuracy,
		      int mean_aber,
		      int fast_math,
		      int polconv,
		      int pair_dets,
		      int pix_order,
                      int fast_pix,
		      int num_threads);

#define OPTIONFUNC(opt)                                 \
  void qp_set_opt_##opt(qp_memory_t *mem, int val);     \
  int qp_get_opt_##opt(qp_memory_t *mem);
  OPTIONFUNC(accuracy);
  OPTIONFUNC(mean_aber);
  OPTIONFUNC(fast_math);
  OPTIONFUNC(polconv);
  OPTIONFUNC(pair_dets);
  OPTIONFUNC(pix_order);
  OPTIONFUNC(fast_pix);
  OPTIONFUNC(num_threads);
  OPTIONFUNC(thread_num);
  
  /* Set weather data */
  void qp_set_weather(qp_memory_t *mem,
		      double height, double temperature, double pressure,
		      double humidity, double frequency, double lapse_rate);
  
#define WEATHFUNC(param)                                        \
  void qp_set_weather_##param(qp_memory_t *mem, double val);    \
  double qp_get_weather_##param(qp_memory_t *mem);
  WEATHFUNC(height)
  WEATHFUNC(temperature)
  WEATHFUNC(pressure)
  WEATHFUNC(humidity)
  WEATHFUNC(frequency)
  WEATHFUNC(lapse_rate)

#define DOUBLEFUNC(param)                               \
  void qp_set_##param(qp_memory_t *mem, double val);    \
  double qp_get_##param(qp_memory_t *mem);
  DOUBLEFUNC(ref_tol)
  DOUBLEFUNC(ref_delta)
  DOUBLEFUNC(dut1)

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
  int get_iers_bulletin_a( qp_memory_t *mem, double mjd,
			   double *dut1, double *x, double *y );
  /* Set IERS Bulletin A */
  int set_iers_bulletin_a( qp_memory_t *mem, int mjd_min_, int mjd_max_,
			   double *dut1, double *x, double *y );
  /* Copy IERS Bulletin A */
  int copy_iers_bulletin_a( qp_memory_t *memdest, qp_memory_t *memsrc );

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
  void qp_aberration(quat_t q, vec3_t beta, quat_t qa);
  
  /* Calculate earth orbital velocity vector as fraction of speed of light */
  void qp_earth_orbital_beta(double jd_tdb[2], vec3_t beta);
  
  /* Apply annual aberration correction to given quaternion */
  void qp_apply_annual_aberration(qp_memory_t *mem, double ctime, quat_t q);
  
  /* Apply diurnal aberration correction to given quaternion */
  void qp_apply_diurnal_aberration(qp_memory_t *mem, double ctime, double lat,
                                   quat_t q);

  /* Calculate nutation/precession/bias correction quaternion
     use (faster) truncated series if accuracy > 0*/
  void qp_npb_quat(double jd_tt[2], quat_t q, int accuracy);
  
  /* Calculate ERA quaternion */
  void qp_erot_quat(double jd_ut1[2], quat_t q);
  
  /* Calcuate wobble correction quaternion */
  void qp_wobble_quat(double jd_tt[2], double xp, double yp, quat_t q);
  
  /* Calculate gondola orientation quaternion */
  void qp_azel_quat(double az, double el, double pitch, double roll, quat_t q);
  
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
  
  /* Calculate dipole amplitude */
  double qp_dipole(qp_memory_t *mem, double ctime, double ra, double dec);

  /* Calculate dipole amplitudes */
  void qp_dipolen(qp_memory_t *mem, double *ctime, double *ra, double *dec,
                  double *dipole, int n);

  /* Calculate waveplate quaternion, given _physical_ HWP angle */
  void qp_hwp_quat(double ang, quat_t q);
  
  /* Calculate waveplate quaternions */
  void qp_hwp_quatn(double *ang, quat_t *q, int n);
  
  /* Calculate atmospheric refraction */
  double qp_refraction(double el, double lat, double height, double temp,
		       double press, double hum, double freq, double lapse,
		       double tol);
  
  /* Update atmospheric refraction using stored parameters */
  double qp_update_ref(qp_memory_t *mem, quat_t q, double lat);

  /* Apply refraction correction */
  void qp_apply_refraction(qp_memory_t *mem, double ctime, double lat, quat_t q);

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
  
  /* Compute detector offset quaternion. */
  void qp_det_offset(double delta_az, double delta_el, double delta_psi,
		     quat_t q);

  
  /* Compute detector offset quaternion. */
  void qp_det_offsetn(double *delta_az, double *delta_el, double *delta_psi,
		      quat_t *q, int n);
  
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
  
  /* Calculate ra/dec and sin(2*psi)/cos(2*psi) for a given detector offset,
     from a set of boresight orientations. */
  void qp_azel2radec(qp_memory_t *mem,
		     double delta_az, double delta_el, double delta_psi,
		     double *az, double *el, double *pitch, double *roll,
		     double *lon, double *lat, double *ctime, 
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
  
  /* ************************************************************************* 
     Pixelization
     ********************************************************************** */
  
  /* Pixel, contains (hits, p01, p02, p11, p12, p22) */
  typedef double vec6_t[6];
  
  /* map with first derivatives */
  typedef double vec9_t[9];

  /* map with first and second derivatives */
  typedef double vec18_t[18];

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

  /* Compute pixel numbers and pol angles for given nside and quaternions */
  void qp_quat2pixn(qp_memory_t *mem, quat_t *q, int nside, long *pix,
                    double *sin2psi, double *cos2psi, int n);

  /* Rotate from celestial to galactic coordinates */
  void qp_radec2gal(qp_memory_t *mem, double *ra, double *dec,
                    double *sin2psi, double *cos2psi);
  void qp_radec2galn(qp_memory_t *mem, double *ra, double *dec,
                     double *sin2psi, double *cos2psi, int n);

  /* Rotate from galactic to celestial coordinates */
  void qp_gal2radec(qp_memory_t *mem, double *ra, double *dec,
                    double *sin2psi, double *cos2psi);
  void qp_gal2radecn(qp_memory_t *mem, double *ra, double *dec,
                     double *sin2psi, double *cos2psi, int n);

  /* Rotate a TQU map from one coordinate system to another */
  void qp_rotate_map(qp_memory_t *mem, int nside,
                     vec3_t *map_in, const char coord_in,
                     vec3_t *map_out, const char coord_out);

  /* Compute pix/pol timestreams for given boresight timestream and detector
     offset */
  void qp_bore2pix(qp_memory_t *mem, quat_t q_off, double *ctime, quat_t *q_bore,
                   int nside, long *pix, double *sin2psi, double *cos2psi, int n);

  /* Compute pix/pol timestreams for given boresight timestream,
     waveplate timestream and detector offset */
  void qp_bore2pix_hwp(qp_memory_t *mem, quat_t q_off, double *ctime,
                       quat_t *q_bore, quat_t *q_hwp, int nside, long *pix,
                       double *sin2psi, double *cos2psi, int n);
  
  /* Compute pointing matrix map for given boresight timestream and detector
     offset. pmap is a npix-x-6 array containing (hits, p01, p02, p11, p12, p22) */
  void qp_tod2map_pnt_single(qp_memory_t *mem, quat_t q_off,
                             double *ctime, quat_t *q_bore, int n,
                             vec6_t *pmap, int nside);
  
  /* Compute pointing matrix map for given boresight timestream and detector
     offset. pmap is a npix-x-6 array containing (hits, p01, p02, p11, p12, p22) */
  void qp_tod2map_pnt_single_hwp(qp_memory_t *mem, quat_t q_off,
                                 double *ctime, quat_t *q_bore, quat_t *q_hwp, int n,
                                 vec6_t *pmap, int nside);
  
  /* Compute pointing matrix map for given boresight timestream and detector
     offset for both A and B polarizations.
     pmap is a npix-x-6 array containing (hits, p01, p02, p11, p12, p22) */
  void qp_tod2map_pnt_pair(qp_memory_t *mem, quat_t q_off,
                           double *ctime, quat_t *q_bore, int n,
                           vec6_t *pmap, int nside);
  
  /* Compute pointing matrix map for given boresight timestream and detector
     offset for both A and B polarizations.
     pmap is a npix-x-6 array containing (hits, p01, p02, p11, p12, p22) */
  void qp_tod2map_pnt_pair_hwp(qp_memory_t *mem, quat_t q_off,
                               double *ctime, quat_t *q_bore, quat_t *q_hwp, int n,
                               vec6_t *pmap, int nside);
  
  /* Compute pointing matrix map for given boresight timestream and many detector
     offsets. pmap is a npix-x-6 array containing (hits, p01, p02, p11, p12, p22) */
  void qp_tod2map_pnt(qp_memory_t *mem, quat_t *q_off, int ndet,
                      double *ctime, quat_t *q_bore, int n,
                      vec6_t *pmap, int nside);
  
  /* Compute pointing matrix map for given boresight timestream and many detector
     offsets. pmap is a npix-x-6 array containing (hits, p01, p02, p11, p12, p22) */
  void qp_tod2map_pnt_hwp(qp_memory_t *mem, quat_t *q_off, int ndet,
                          double *ctime, quat_t *q_bore, quat_t *q_hwp, int n,
                          vec6_t *pmap, int nside);

  /* Compute hits map for given boresight timestream and many detector
     offsets. pmap is a npix array containing (hits) */
  void qp_tod2map_pnt_nopol(qp_memory_t *mem, quat_t *q_off, int ndet,
                            double *ctime, quat_t *q_bore, int n,
                            double *pmap, int nside);
  
  /* Compute signal map for given boresight timestream, signal timestream,
     and detector offset.
     smap is a npix-x-3 array containing (d, d*cos(2 psi), d*sin(2 psi)). */
  void qp_tod2map_sig_single(qp_memory_t *mem, quat_t q_off,
                             double *ctime, quat_t *q_bore, double *tod, int n,
                             vec3_t *smap, int nside);

  /* Compute signal map for given boresight timestream, hwp timestream,
     signal timestream and detector offset.
     smap is a npix-x-3 array containing (d, d*cos(2 psi), d*sin(2 psi)). */
  void qp_tod2map_sig_single_hwp(qp_memory_t *mem, quat_t q_off,
                                 double *ctime, quat_t *q_bore, quat_t *q_hwp,
                                 double *tod, int n, vec3_t *smap, int nside);

  /* Compute signal map for given boresight timestream, many signal timestreams,
     and many detector offsets.
     smap is a npix-x-3 array containing (d, d*cos(2 psi), d*sin(2 psi)). */
  void qp_tod2map_sig(qp_memory_t *mem, quat_t *q_off, int ndet,
                      double *ctime, quat_t *q_bore, double **tod, int n,
                      vec3_t *smap, int nside);

  /* Compute signal map for given boresight timestream, hwp timestream,
     many signal timestreams, and many detector offsets.
     smap is a npix-x-3 array containing (d, d*cos(2 psi), d*sin(2 psi)). */
  void qp_tod2map_sig_hwp(qp_memory_t *mem, quat_t *q_off, int ndet,
                          double *ctime, quat_t *q_bore, quat_t *q_hwp,
                          double **tod, int n, vec3_t *smap, int nside);

  /* Compute signal map for given boresight timestream, many signal timestreams,
     and many detector offsets.
     smap is a npix array. */
  void qp_tod2map_sig_nopol(qp_memory_t *mem, quat_t *q_off, int ndet,
                            double *ctime, quat_t *q_bore, double **tod, int n,
                            double *smap, int nside);

  /* Compute signal and pointing matrix maps for given boresight timestream,
     signal timestream, and detector offset.
     smap is a npix-x-3 array containing (d, d*cos(2 psi), d*sin(2 psi)).
     pmap is a npix-x-6 array containing (hits, p01, p02, p11, p12, p22) */
  void qp_tod2map_sigpnt_single(qp_memory_t *mem, quat_t q_off,
                                double *ctime, quat_t *q_bore, double *tod, int n,
                                vec3_t *smap, vec6_t *pmap, int nside);

  /* Compute signal and pointing matrix maps for given boresight timestream,
     hwp timestream, signal timestream, and detector offset.
     smap is a npix-x-3 array containing (d, d*cos(2 psi), d*sin(2 psi)).
     pmap is a npix-x-6 array containing (hits, p01, p02, p11, p12, p22) */
  void qp_tod2map_sigpnt_single_hwp(qp_memory_t *mem, quat_t q_off,
                                    double *ctime, quat_t *q_bore, quat_t *q_hwp,
                                    double *tod, int n, vec3_t *smap,
                                    vec6_t *pmap, int nside);

  /* Compute signal and pointing matrix maps for given boresight timestream,
     many signal timestreams, and many detector offsets.
     smap is a npix-x-3 array containing (d, d*cos(2 psi), d*sin(2 psi)).
     pmap is a npix-x-6 array containing (hits, p01, p02, p11, p12, p22) */
  void qp_tod2map_sigpnt(qp_memory_t *mem, quat_t *q_off, int ndet,
                         double *ctime, quat_t *q_bore, double **tod, int n,
                         vec3_t *smap, vec6_t *pmap, int nside);

  /* Compute signal and pointing matrix maps for given boresight timestream,
     hwp timestream, many signal timestreams, and many detector offsets.
     smap is a npix-x-3 array containing (d, d*cos(2 psi), d*sin(2 psi)).
     pmap is a npix-x-6 array containing (hits, p01, p02, p11, p12, p22) */
  void qp_tod2map_sigpnt_hwp(qp_memory_t *mem, quat_t *q_off, int ndet,
                             double *ctime, quat_t *q_bore, quat_t *q_hwp,
                             double **tod, int n, vec3_t *smap, vec6_t *pmap,
                             int nside);

  /* Compute signal and hits maps for given boresight timestream,
     many signal timestreams, and many detector offsets.
     smap is a npix array. pmap is a npix array containing (hits) */
  void qp_tod2map_sigpnt_nopol(qp_memory_t *mem, quat_t *q_off, int ndet,
                               double *ctime, quat_t *q_bore, double **tod, int n,
                               double *smap, double *pmap, int nside);

  /* Compute signal timestream given an input map, boresight pointing
     and detector offset. smap is a npix-x-3 array containing (T,Q,U) maps. */
  void qp_map2tod_single(qp_memory_t *mem, quat_t q_off,
                         double *ctime, quat_t *q_bore, vec3_t *smap,
                         int nside, double *tod, int n);

  /* Compute signal timestream given an input map, boresight pointing,
     HWP and detector offset. smap is a npix-x-3 array containing (T,Q,U) maps. */
  void qp_map2tod_single_hwp(qp_memory_t *mem, quat_t q_off,
                             double *ctime, quat_t *q_bore, quat_t *q_hwp,
                             vec3_t *smap, int nside, double *tod, int n);

  /* Compute signal timestream given an input map, boresight pointing
     and detector offset. smap is a npix-x-9 array containing (T,Q,U) maps
     and their first derivatives */
  void qp_map2tod_der1_single(qp_memory_t *mem, quat_t q_off,
                              double *ctime, quat_t *q_bore, vec9_t *smap,
                              int nside, double *tod, int n);

  /* Compute signal timestream given an input map, boresight pointing
     and detector offset. smap is a npix-x-9 array containing (T,Q,U) maps
     and their first derivatives. */
  void qp_map2tod_der1_single_hwp(qp_memory_t *mem, quat_t q_off,
                                  double *ctime, quat_t *q_bore, quat_t *q_hwp,
                                  vec9_t *smap, int nside, double *tod, int n);

  /* Compute signal timestream given an input map, boresight pointing
     and detector offset. smap is a npix-x-18 array containing (T,Q,U) maps
     and their first and second derivatives. */
  void qp_map2tod_der2_single(qp_memory_t *mem, quat_t q_off,
                              double *ctime, quat_t *q_bore, vec18_t *smap,
                              int nside, double *tod, int n);

  /* Compute signal timestream given an input map, boresight pointing
     and detector offset. smap is a npix-x-18 array containing (T,Q,U) maps
     and their first and second derivatives. */
  void qp_map2tod_der2_single_hwp(qp_memory_t *mem, quat_t q_off,
                                  double *ctime, quat_t *q_bore, quat_t *q_hwp,
                                  vec18_t *smap, int nside, double *tod, int n);

  /* Compute signal timestream given an input map, boresight pointing
     and detector offsets. smap is a npix-x-3 array containing (T,Q,U) maps. */
  void qp_map2tod(qp_memory_t *mem, quat_t *q_off, int ndet,
                  double *ctime, quat_t *q_bore, vec3_t *smap, int nside,
                  double **tod, int n);

  /* Compute signal timestream given an input map, boresight pointing,
     hwp timestream, and detector offsets.
     smap is a npix-x-3 array containing (T,Q,U) maps. */
  void qp_map2tod_hwp(qp_memory_t *mem, quat_t *q_off, int ndet,
                      double *ctime, quat_t *q_bore, quat_t *q_hwp,
                      vec3_t *smap, int nside, double **tod, int n);

  /* Compute signal timestream given an input map, boresight pointing
     and detector offset. smap is a npix-x-9 array containing (T,Q,U) maps
     and their first derivatives. */
  void qp_map2tod_der1(qp_memory_t *mem, quat_t *q_off, int ndet,
                       double *ctime, quat_t *q_bore, vec9_t *smap,
                       int nside, double **tod, int n);
  
  /* Compute signal timestream given an input map, boresight pointing
     hwp and detector offset. smap is a npix-x-9 array containing (T,Q,U) maps
     and their first derivatives. */
  void qp_map2tod_der1_hwp(qp_memory_t *mem, quat_t *q_off, int ndet,
                           double *ctime, quat_t *q_bore, quat_t *q_hwp,
                           vec9_t *smap, int nside, double **tod, int n);

  /* Compute signal timestream given an input map, boresight pointing
     and detector offset. smap is a npix-x-3 array containing (T,Q,U) maps
     and their first and second derivatives. */
  void qp_map2tod_der2(qp_memory_t *mem, quat_t *q_off, int ndet,
                       double *ctime, quat_t *q_bore, vec18_t *smap,
                       int nside, double **tod, int n);

  /* Compute signal timestream given an input map, boresight pointing,
     hwp and detector offset. smap is a npix-x-3 array containing (T,Q,U) maps
     and their first and second derivatives. */
  void qp_map2tod_der2_hwp(qp_memory_t *mem, quat_t *q_off, int ndet,
                           double *ctime, quat_t *q_bore, quat_t *q_hwp,
                           vec18_t *smap, int nside, double **tod, int n);

  /* Compute signal timestream given an input map, boresight pointing
     and detector offsets. smap is a npix array. */
  void qp_map2tod_nopol(qp_memory_t *mem, quat_t *q_off, int ndet,
                        double *ctime, quat_t *q_bore, double *smap, int nside,
                        double **tod, int n);

  /* Compute signal timestream given an input map, boresight pointing
     and detector offset. smap is a npix-x-3 array containing a map and its
     first derivatives. */
  void qp_map2tod_der1_nopol(qp_memory_t *mem, quat_t *q_off, int ndet,
                             double *ctime, quat_t *q_bore, vec3_t *smap,
                             int nside, double **tod, int n);

  /* Compute signal timestream given an input map, boresight pointing
     and detector offset. smap is a npix-x-6 array containing a map and
     its first and second derivatives. */
  void qp_map2tod_der2_nopol(qp_memory_t *mem, quat_t *q_off, int ndet,
                             double *ctime, quat_t *q_bore, vec6_t *smap,
                             int nside, double **tod, int n);

#ifdef __cplusplus
}
#endif
