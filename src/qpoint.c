#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include "sofa.h"
#include "slarefro.h"
#include "qpoint.h"
#include "fast_math.h"
#include "vec3.h"
#include "quaternion.h"

void ctime2jd(double ctime, double jd[2]) {
  jd[0] = CTIME_JD_EPOCH;
  jd[1] = secs2days(ctime);
}

double jd2ctime(double jd[2]) {
  double d = (jd[0] - CTIME_JD_EPOCH) + jd[1];
  return days2secs(d);
}

void ctime2jdtt(double ctime, double jd_tt[2]) {
  double jd_utc[2], jd_tai[2];
  int stat;
  
  ctime2jd(ctime, jd_utc);
  stat = iauUtctai(jd_utc[0], jd_utc[1], &jd_tai[0], &jd_tai[1]);
  assert(stat == 0);
  stat = iauTaitt(jd_tai[0], jd_tai[1], &jd_tt[0], &jd_tt[1]);
  assert(stat == 0);  
}

void jdutc2jdut1(double jd_utc[2], double dut1, double jd_ut1[2]) {
  int stat;
  stat = iauUtcut1(jd_utc[0], jd_utc[1], dut1, &jd_ut1[0], &jd_ut1[1]);
  assert(stat == 0);
}

void qp_aberration(quat_t q, vec3_t beta, quat_t qa) {
  vec3_t u,n;
  Quaternion_to_matrix_col3(q,u);
  vec3_cross_product(n, u, beta);
  // small angle approximation!
  double sa_2 = 0.5*vec3_norm(n);
  qa[0] = 1. - 0.5*sa_2*sa_2;
  qa[1] = -0.5*n[0];
  qa[2] = -0.5*n[1];
  qa[3] = -0.5*n[2];
};

void qp_earth_orbital_beta(double jd_tdb[2], vec3_t beta) {
  double pvb[2][3];
  iauEpv00(jd_tdb[0], jd_tdb[1], pvb, pvb);
  for (int i = 0; i < 3; i++) beta[i] = pvb[1][i]/C_AUD;
}

void qp_lonlat_quat(double lon, double lat, quat_t q) {
  Quaternion_r3(q, M_PI);
  Quaternion_r2_mul(M_PI_2 - deg2rad(lat), q);
  Quaternion_r3_mul(deg2rad(lon), q);
}

void qp_azel_quat(double az, double el, double pitch, double roll, quat_t q) {
  Quaternion_r3(q, M_PI);
  Quaternion_r2_mul(M_PI_2 - deg2rad(el), q);
  Quaternion_r3_mul(-deg2rad(az), q);
  Quaternion_r2_mul(-deg2rad(pitch), q);
  Quaternion_r1_mul(-deg2rad(roll), q);
}

void qp_npb_quat(double jd_tt[2], quat_t q, int accuracy) {
  double X,Y,s,Z,E,d;
  if (accuracy == 0)
    iauXys06a(jd_tt[0],jd_tt[1], &X, &Y, &s);
  else
    iauXys00b(jd_tt[0],jd_tt[1], &X, &Y, &s);
  Z = sqrt(1.0 - X*X - Y*Y);
  E = atan2(Y, X);
  d = acos(Z);
  
  Quaternion_r3(q, -E - s);
  Quaternion_r2_mul(d, q);
  Quaternion_r3_mul(E, q);
}

void qp_erot_quat(double jd_ut1[2], quat_t q) {
  double theta = iauEra00(jd_ut1[0], jd_ut1[1]);
  Quaternion_r3(q, theta);
}

void qp_wobble_quat(double xp, double yp, quat_t q) {
  double sprime = 0;
  
  Quaternion_r1(q, -arcsec2rad(yp));
  Quaternion_r2_mul(-arcsec2rad(xp), q);
  Quaternion_r3_mul(sprime, q);
}

/* Calculate atmospheric refraction */
double qp_refraction(double el, double lat, double height, double temp,
		     double press, double hum, double freq, double lapse,
		     double tol) {
  double ref;
  slaf_refro(M_PI_2 - deg2rad(el),
	     height, temp + 273.15, // temperature, K
	     press, hum, C_MS * 1e-3 / freq, // wavelength, um
	     deg2rad(lat), lapse, tol, // precision, radians
	     &ref);
  return rad2deg(ref);
}

double qp_update_ref(qp_memory_t *mem, double el, double lat) {
  qp_weather_t *W = &mem->weather;
  double ref = qp_refraction(el, W->height, W->temperature,
			     W->pressure, W->humidity, W->frequency,
			     lat, W->lapse_rate, mem->ref_tol);
  mem->ref_delta = ref;
  return ref;
}

void qp_apply_annual_aberration(qp_memory_t *mem, double ctime, quat_t q) {
  quat_t q_aber;
  double jd_tt[2];
  
  if (qp_check_update(&mem->state_aaber, ctime)) {
    ctime2jdtt(ctime, jd_tt);
    qp_earth_orbital_beta(jd_tt, mem->beta_earth);
  }
  if (qp_check_apply(&mem->state_aaber)) {
    qp_aberration(q, mem->beta_earth, q_aber);
    Quaternion_mul_left(q_aber, q);
#ifdef DEBUG
    qp_print_debug("aaber", q_aber);
    qp_print_debug("state aaber", q);
#endif
  }
}

void qp_azel2quat(qp_memory_t *mem, double az, double el, double pitch,
		  double roll, double lon, double lat, double ctime,
		  quat_t q) {
  
  double jd_utc[2], jd_tt[2], jd_ut1[2], mjd_utc;
  double x,y;
  double delta;
  quat_t q_step;
  static const vec3_t beta_diurnal = {0, -D_ABER_RAD, 0};
  
  // deal with times
  ctime2jd(ctime, jd_utc);
  
#ifdef DEBUG
  printf("ctime %f, jd_utc %f %f\n", ctime, jd_utc[0], jd_utc[1]);
#endif
  
#ifdef DEBUG
  qp_print_debug("state init", q);
#endif
  
  // apply refraction correction
  // NB: if rate is never, correction can still be calculated externally
  // by the user with qp_update_ref()
  if (qp_check_apply(&mem->state_ref))
    delta = qp_update_ref(mem, el, lat);
  else
    delta = mem->ref_delta;
#ifdef DEBUG
  printf("Refraction: %.6g\n", delta);
#endif
  
  // apply boresight rotation
  qp_azel_quat(az, el-delta, pitch, roll, q);
#ifdef DEBUG
  qp_print_debug("azel", q);
#endif
  
  // apply diurnal aberration
  // TODO propagate this to aberration step?
  if (qp_check_apply(&mem->state_daber)) {
    qp_aberration(q, (double *)beta_diurnal, q_step);
    Quaternion_mul_left(q_step, q);
#ifdef DEBUG
    qp_print_debug("daber", q_step);
    qp_print_debug("state daber", q);
#endif
  }
  
  // rotate to ITRS (by lon/lat)
  if (qp_check_update(&mem->state_lonlat, ctime)) {
    qp_lonlat_quat(lon, lat, mem->q_lonlat);
#ifdef DEBUG
    qp_print_debug("lonlat", mem->q_lonlat);
#endif
  }
  if (qp_check_apply(&mem->state_lonlat)) {
    Quaternion_mul_left(mem->q_lonlat, q);
#ifdef DEBUG
    qp_print_debug("state lonlat", q);
#endif
  }
  
  // apply wobble correction (polar motion)
  // or get dut1 from IERS bulletin
  mjd_utc = jd2mjd(jd_utc[0]) + jd_utc[1];
  if (qp_check_update(&mem->state_wobble, ctime)) {
    get_iers_bulletin_a(mem, mjd_utc, &mem->dut1, &x, &y);
    qp_wobble_quat(x, y, mem->q_wobble);
#ifdef DEBUG
    qp_print_debug("wobble", mem->q_wobble);
#endif
  } else if (qp_check_update(&mem->state_dut1, ctime))
    get_iers_bulletin_a(mem, mjd_utc, &mem->dut1, &x, &y);
  if (qp_check_apply(&mem->state_wobble)) {
    Quaternion_mul_left(mem->q_wobble, q);
#ifdef DEBUG
    qp_print_debug("state wobble", q);
#endif
  }
  
  // apply earth rotation
  if (qp_check_update(&mem->state_erot, ctime)) {
    // get ut1
    jdutc2jdut1(jd_utc, mem->dut1, jd_ut1);
    qp_erot_quat(jd_ut1, mem->q_erot);
#ifdef DEBUG
    qp_print_debug("erot", mem->q_erot);
#endif
  }
  if (qp_check_apply(&mem->state_erot)) {
    Quaternion_mul_left(mem->q_erot, q);
#ifdef DEBUG
    qp_print_debug("state erot", q);
#endif
  }
  
  // apply nutation/precession/frame bias correction
  if (qp_check_update(&mem->state_npb, ctime)) {
    ctime2jdtt(ctime, jd_tt);
    qp_npb_quat(jd_tt, mem->q_npb, mem->accuracy);
#ifdef DEBUG
    qp_print_debug("npb", mem->q_npb);
#endif
  }
  if (qp_check_apply(&mem->state_npb)) {
    Quaternion_mul_left(mem->q_npb, q);
#ifdef DEBUG
    qp_print_debug("state npb", q);
#endif
  }
  
  // apply annual aberration
  // ~20 arcsec max
  if (mem->mean_aber)
    qp_apply_annual_aberration(mem, ctime, q);
  
#ifdef DEBUG
  qp_print_debug("state final", q);
#endif
}

void qp_azel2bore(qp_memory_t *mem, double *az, double *el, double *pitch,
		  double *roll, double *lon, double *lat, double *ctime, 
		  quat_t *q, int n) {
  for (int i=0; i<n; i++)
    qp_azel2quat(mem, az[i], el[i], pitch[i], roll[i], lon[i], lat[i],
		 ctime[i], q[i]);
}

void qp_hwp_quat(double ang, quat_t q) {
  Quaternion_r3(q, -deg2rad(ang));
#ifdef DEBUG
  qp_print_debug("hwp", q);
#endif
}

void qp_hwp_quatn(double *ang, quat_t *q, int n) {
  for (int ii=0; ii<n; ii++)
    qp_hwp_quat(ang[ii], q[ii]);
}

void qp_det_offset(double delta_az, double delta_el, double delta_psi, quat_t q) {
  Quaternion_r3(q, -deg2rad(delta_psi));
  Quaternion_r2_mul(deg2rad(delta_el), q);
  Quaternion_r1_mul(-deg2rad(delta_az), q);
#ifdef DEBUG
  qp_print_debug("offset", q);
#endif
}

void qp_det_offsetn(double *delta_az, double *delta_el, double *delta_psi, quat_t *q,
		    int n) {
  for (int ii=0; ii<n; ii++)
    qp_det_offset(delta_az[ii], delta_el[ii], delta_psi[ii], q[ii]);
}

void qp_bore2det(qp_memory_t *mem, quat_t q_off, double ctime, quat_t q_bore,
		 quat_t q_det) {
  Quaternion_copy(q_det,q_off);
  Quaternion_mul_left(q_bore, q_det);
  
  if (!mem->mean_aber)
    qp_apply_annual_aberration(mem, ctime, q_det);
}

void qp_bore2det_hwp(qp_memory_t *mem, quat_t q_off, double ctime, quat_t q_bore,
		     quat_t q_hwp, quat_t q_det) {

  qp_bore2det(mem, q_off, ctime, q_bore, q_det);
  Quaternion_mul_right(q_det, q_hwp);
}

void qp_quat2rasindec(qp_memory_t *mem, quat_t q, double *ra, double *sindec,
		      double *sin2psi, double *cos2psi) {
  
  // ZYZ euler angles
  
  double q01 = q[0]*q[1];
  double q02 = q[0]*q[2];
  double q13 = q[1]*q[3];
  double q23 = q[2]*q[3];
  double q00p33 = q[0]*q[0] + q[3]*q[3];
  double q11p22 = q[1]*q[1] + q[2]*q[2];
  // NB: factors of two have been redistributed...
  double sina_2 = q23 - q01;
  double cosa_2 = q02 + q13;
  double icosb2 = 1./(q00p33*q11p22);
  double sinb = q00p33 - q11p22;
  double sing_cb = q01 + q23;
  double cosg_cb = q02 - q13;
  
  // cosmo (healpix) or IAU polarization convention?
  if (!mem->polconv) sing_cb = -sing_cb;
  
  if (mem->fast_math)
    *ra = rad2deg(poly_atan2(sina_2, cosa_2));
  else
    *ra = rad2deg(atan2(sina_2, cosa_2));
  *sindec = sinb;
  
  *sin2psi = 2.*sing_cb*cosg_cb*icosb2;
  *cos2psi = 2.*cosg_cb*cosg_cb*icosb2 - 1.;
}

void qp_quat2radec(qp_memory_t *mem, quat_t q, double *ra, double *dec, 
		   double *sin2psi, double *cos2psi) {
  
  // ZYZ euler angles
  
  double q01 = q[0]*q[1];
  double q02 = q[0]*q[2];
  double q13 = q[1]*q[3];
  double q23 = q[2]*q[3];
  double q00p33 = q[0]*q[0] + q[3]*q[3];
  double q11p22 = q[1]*q[1] + q[2]*q[2];
  // NB: factors of 2 have been redistributed...
  double sina_2 = q23 - q01;
  double cosa_2 = q02 + q13;
  double cosb2 = q00p33*q11p22;
  double cosb_2 = sqrt(cosb2); // hypot(sina,cosa);
  double icosb2 = 1./cosb2;
  double sinb_2 = 0.5*(q00p33 - q11p22);
  double sing_cb = q01 + q23;
  double cosg_cb = q02 - q13;
  
  // cosmo (healpix) or IAU polarization convention?
  if (!mem->polconv) sing_cb = -sing_cb;
  
  if (mem->fast_math) {
    *ra = rad2deg(poly_atan2(sina_2, cosa_2));
    *dec = rad2deg(poly_atan2(sinb_2, cosb_2));
  } else {
    *ra = rad2deg(atan2(sina_2, cosa_2));
    *dec = rad2deg(atan2(sinb_2, cosb_2));
  }
  *sin2psi = 2.*sing_cb*cosg_cb*icosb2;
  *cos2psi = 2.*cosg_cb*cosg_cb*icosb2 - 1.;
}

void qp_bore2radec(qp_memory_t *mem, quat_t q_off, double *ctime, quat_t *q_bore,
		   double *ra, double *dec, double *sin2psi, 
		   double *cos2psi, int n) {
  quat_t q;
  
  for (int i=0; i<n; i++) {
    qp_bore2det(mem, q_off, ctime[i], q_bore[i], q);
    qp_quat2radec(mem, q, ra+i, dec+i, sin2psi+i, cos2psi+i);
  }
}

void qp_bore2radec_hwp(qp_memory_t *mem, quat_t q_off, double *ctime, quat_t *q_bore,
		       quat_t *q_hwp, double *ra, double *dec, double *sin2psi, 
		       double *cos2psi, int n) {
  quat_t q;
  
  for (int i=0; i<n; i++) {
    qp_bore2det_hwp(mem, q_off, ctime[i], q_bore[i], q_hwp[i], q);
    qp_quat2radec(mem, q, ra+i, dec+i, sin2psi+i, cos2psi+i);
  }
}

void qp_bore2rasindec(qp_memory_t *mem, quat_t q_off, double *ctime, quat_t *q_bore,
		      double *ra, double *sindec, double *sin2psi, 
		      double *cos2psi, int n) {
  quat_t q;

  for (int i=0; i<n; i++) {
    qp_bore2det(mem, q_off, ctime[i], q_bore[i], q);
    qp_quat2rasindec(mem, q, ra+i, sindec+i, sin2psi+i, cos2psi+i);
  }
}

void qp_bore2rasindec_hwp(qp_memory_t *mem, quat_t q_off, double *ctime, quat_t *q_bore,
			  quat_t *q_hwp, double *ra, double *sindec, double *sin2psi, 
			  double *cos2psi, int n) {
  quat_t q;
  
  for (int i=0; i<n; i++) {
    qp_bore2det_hwp(mem, q_off, ctime[i], q_bore[i], q_hwp[i], q);
    qp_quat2rasindec(mem, q, ra+i, sindec+i, sin2psi+i, cos2psi+i);
  }
}

// all input and output angles are in degrees!
void qp_azel2radec(qp_memory_t *mem,
		   double delta_az, double delta_el, double delta_psi,
		   double *az, double *el, double *pitch, double *roll,
		   double *lon, double *lat, double *ctime,
		   double *ra, double *dec, double *sin2psi, 
		   double *cos2psi, int n) {
  quat_t q_det, q_off, q_bore;
  
  qp_det_offset(delta_az, delta_el, delta_psi, q_off);
  
  for (int i=0; i<n; i++) {
    qp_azel2quat(mem, az[i], el[i], pitch[i], roll[i], lon[i], lat[i], ctime[i],
		 q_bore);
    qp_bore2det(mem, q_off, ctime[i], q_bore, q_det);
    qp_quat2radec(mem, q_det, &ra[i], &dec[i], &sin2psi[i], &cos2psi[i]);
  }
}

// all input and output angles are in degrees!
void qp_azel2radec_hwp(qp_memory_t *mem,
		       double delta_az, double delta_el, double delta_psi,
		       double *az, double *el, double *pitch, double *roll,
		       double *lon, double *lat, double *ctime, double *hwp,
		       double *ra, double *dec, double *sin2psi, 
		       double *cos2psi, int n) {
  quat_t q_det, q_off, q_bore, q_hwp;
  
  qp_det_offset(delta_az, delta_el, delta_psi, q_off);
  
  for (int i=0; i<n; i++) {
    qp_azel2quat(mem, az[i], el[i], pitch[i], roll[i], lon[i], lat[i], ctime[i],
		 q_bore);
    qp_hwp_quat(hwp[i], q_hwp);
    qp_bore2det_hwp(mem, q_off, ctime[i], q_bore, q_hwp, q_det);
    qp_quat2radec(mem, q_det, &ra[i], &dec[i], &sin2psi[i], &cos2psi[i]);
  }
}

// all input and output angles are in degrees!
void qp_azel2rasindec(qp_memory_t *mem,
		      double delta_az, double delta_el, double delta_psi,
		      double *az, double *el, double *pitch, double *roll,
		      double *lon, double *lat, double *ctime,
		      double *ra, double *sindec, double *sin2psi, 
		      double *cos2psi, int n) {
  quat_t q_det, q_off, q_bore;
  
  qp_det_offset(delta_az, delta_el, delta_psi, q_off);
  
  for (int i=0; i<n; i++) {
    qp_azel2quat(mem, az[i], el[i], pitch[i], roll[i], lon[i], lat[i], ctime[i],
		 q_bore);
    qp_bore2det(mem, q_off, ctime[i], q_bore, q_det);
    qp_quat2rasindec(mem, q_det, &ra[i], &sindec[i], &sin2psi[i], &cos2psi[i]);
  }
}

// all input and output angles are in degrees!
void qp_azel2rasindec_hwp(qp_memory_t *mem,
			  double delta_az, double delta_el, double delta_psi,
			  double *az, double *el, double *pitch, double *roll,
			  double *lon, double *lat, double *ctime, double *hwp,
			  double *ra, double *sindec, double *sin2psi, 
			  double *cos2psi, int n) {
  quat_t q_det, q_off, q_bore, q_hwp;
  
  qp_det_offset(delta_az, delta_el, delta_psi, q_off);
  
  for (int i=0; i<n; i++) {
    qp_azel2quat(mem, az[i], el[i], pitch[i], roll[i], lon[i], lat[i], ctime[i],
		 q_bore);
    qp_hwp_quat(hwp[i], q_hwp);
    qp_bore2det_hwp(mem, q_off, ctime[i], q_bore, q_hwp, q_det);
    qp_quat2rasindec(mem, q_det, &ra[i], &sindec[i], &sin2psi[i], &cos2psi[i]);
  }
}
