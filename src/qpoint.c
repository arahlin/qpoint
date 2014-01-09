#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include "sofa.h"
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

void qp_wobble_quat(double mjd_utc, double *dut1, quat_t q) {
  double xp, yp, sprime = 0;
  
  get_iers_bulletin_a(mjd_utc, dut1, &xp, &yp);
  
  Quaternion_r1(q, -arcsec2rad(yp));
  Quaternion_r2_mul(-arcsec2rad(xp), q);
  Quaternion_r3_mul(sprime, q);
}

void qp_apply_annual_aberration(double ctime, quat_t q) {
  quat_t q_aber;
  double jd_tt[2];
  static vec3_t beta_earth;
  
  if (qp_check_update(&qp_params.s_aaber, ctime)) {
    ctime2jdtt(ctime, jd_tt);
    qp_earth_orbital_beta(jd_tt, beta_earth);
  }
  if (qp_check_apply(&qp_params.s_aaber)) {
    qp_aberration(q, beta_earth, q_aber);
    Quaternion_mul_left(q_aber, q);
#ifdef DEBUG
    qp_print_debug("aaber", q_aber);
    qp_print_debug("state aaber", q);
#endif
  }
}

void qp_azel2quat(double az, double el, double pitch, double roll,
		  double lon, double lat, double ctime, quat_t q) {
  
  double jd_utc[2], jd_tt[2], jd_ut1[2], mjd_utc;
  double x,y;
  static double dut1;
  quat_t q_step;
  static quat_t q_lonlat, q_wobble, q_npb, q_erot;
  static const vec3_t beta_diurnal = {0, -D_ABER_RAD, 0};
  
  // deal with times
  ctime2jd(ctime, jd_utc);
  
#ifdef DEBUG
  printf("ctime %f, jd_utc %f %f\n", ctime, jd_utc[0], jd_utc[1]);
#endif
  
#ifdef DEBUG
  qp_print_debug("state init", q);
#endif
  
  // apply boresight rotation
  qp_azel_quat(az, el, pitch, roll, q);
#ifdef DEBUG
  qp_print_debug("azel", q);
#endif
  
  // apply diurnal aberration
  // TODO propagate this to aberration step?
  if (qp_check_apply(&qp_params.s_daber)) {
    qp_aberration(q, (double *)beta_diurnal, q_step);
    Quaternion_mul_left(q_step, q);
#ifdef DEBUG
    qp_print_debug("daber", q_step);
    qp_print_debug("state daber", q);
#endif
  }
  
  // rotate to ITRS (by lon/lat)
  if (qp_check_update(&qp_params.s_lonlat, ctime)) {
    qp_lonlat_quat(lon, lat, q_lonlat);
#ifdef DEBUG
    qp_print_debug("lonlat", q_lonlat);
#endif
  }
  if (qp_check_apply(&qp_params.s_lonlat)) {
    Quaternion_mul_left(q_lonlat, q);
#ifdef DEBUG
    qp_print_debug("state lonlat", q);
#endif
  }
  
  // apply wobble correction (polar motion)
  // or get dut1 from IERS bulletin
  mjd_utc = jd2mjd(jd_utc[0]) + jd_utc[1];
  if (qp_check_update(&qp_params.s_wobble, ctime)) {
    qp_wobble_quat(mjd_utc, &dut1, q_wobble);
#ifdef DEBUG
    qp_print_debug("wobble", q_wobble);
#endif
  } else if (qp_check_update(&qp_params.s_dut1, ctime))
    get_iers_bulletin_a(mjd_utc, &dut1, &x, &y);
  if (qp_check_apply(&qp_params.s_wobble)) {
    Quaternion_mul_left(q_wobble, q);
#ifdef DEBUG
    qp_print_debug("state wobble", q);
#endif
  }
  
  // apply earth rotation
  if (qp_check_update(&qp_params.s_erot, ctime)) {
    // get ut1
    jdutc2jdut1(jd_utc, dut1, jd_ut1);
    qp_erot_quat(jd_ut1, q_erot);
#ifdef DEBUG
    qp_print_debug("erot", q_erot);
#endif
  }
  if (qp_check_apply(&qp_params.s_erot)) {
    Quaternion_mul_left(q_erot, q);
#ifdef DEBUG
    qp_print_debug("state erot", q);
#endif
  }
  
  // apply nutation/precession/frame bias correction
  if (qp_check_update(&qp_params.s_npb, ctime)) {
    ctime2jdtt(ctime, jd_tt);
    qp_npb_quat(jd_tt, q_npb, qp_params.accuracy);
#ifdef DEBUG
    qp_print_debug("npb", q_npb);
#endif
  }
  if (qp_check_apply(&qp_params.s_npb)) {
    Quaternion_mul_left(q_npb, q);
#ifdef DEBUG
    qp_print_debug("state npb", q);
#endif
  }
  
  // apply annual aberration
  // ~20 arcsec max
  if (qp_params.mean_aber)
    qp_apply_annual_aberration(ctime, q);
  
#ifdef DEBUG
  qp_print_debug("state final", q);
#endif
}

void qp_azel2bore(double *az, double *el, double *pitch, double *roll,
		  double *lon, double *lat, double *ctime, quat_t *q, int n) {
  qp_init_params();
  for (int i=0; i<n; i++)
    qp_azel2quat(az[i], el[i], pitch[i], roll[i], lon[i], lat[i], ctime[i], q[i]);
}

void qp_det_offset(double delta_az, double delta_el, double delta_psi, quat_t q) {
  qp_init_params();
  Quaternion_r3(q, -deg2rad(delta_psi));
  Quaternion_r2_mul(deg2rad(delta_el), q);
  Quaternion_r1_mul(-deg2rad(delta_az), q);
#ifdef DEBUG
  qp_print_debug("offset", q);
#endif
}

void qp_bore2det(quat_t q_off, double ctime, quat_t q_bore, quat_t q_det) {
  Quaternion_copy(q_det,q_off);
  Quaternion_mul_left(q_bore, q_det);
  
  if (!qp_params.mean_aber)
    qp_apply_annual_aberration(ctime, q_det);
}

void qp_quat2rasindec(quat_t q, double *ra, double *sindec, double *sin2psi,
		      double *cos2psi) {
  
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
  if (!qp_params.polconv) sing_cb = -sing_cb;
  
  if (qp_params.fast_math)
    *ra = rad2deg(poly_atan2(sina_2, cosa_2));
  else
    *ra = rad2deg(atan2(sina_2, cosa_2));
  *sindec = sinb;
  
  *sin2psi = 2.*sing_cb*cosg_cb*icosb2;
  *cos2psi = 2.*cosg_cb*cosg_cb*icosb2 - 1.;
}

void qp_quat2radec(quat_t q, double *ra, double *dec, 
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
  if (!qp_params.polconv) sing_cb = -sing_cb;
  
  if (qp_params.fast_math) {
    *ra = rad2deg(poly_atan2(sina_2, cosa_2));
    *dec = rad2deg(poly_atan2(sinb_2, cosb_2));
  } else {
    *ra = rad2deg(atan2(sina_2, cosa_2));
    *dec = rad2deg(atan2(sinb_2, cosb_2));
  }
  *sin2psi = 2.*sing_cb*cosg_cb*icosb2;
  *cos2psi = 2.*cosg_cb*cosg_cb*icosb2 - 1.;
}

void qp_interp_radec(double *ra, double *dec, double *sin2psi, double *cos2psi,
		   double *t, int n) {
  int i;
  double ra1, ra2, dec1, dec2, dr, dd, ds, dc, icp12, isp12;
  double sp1 = sin2psi[0];
  double sp2 = sin2psi[n-1];
  double cp1 = cos2psi[0];
  double cp2 = cos2psi[n-1];
  int interps = (-0.5 < sp1 && sp1 < 0.5);
  if (interps)
    icp12 = 1.0/(cp1*cp1);
  else
    isp12 = 1.0/(sp1*sp1);
  
  dr = ra[n-1]-ra[0];
  ra1 = ra[0];
  ra2 = ra[n-1];
  if (dr > 180.) ra2 -= 360.;
  else if (dr < -180.) ra2 += 360.;
  
  dd = dec[n-1] - dec[0];
  dec1 = dec[0];
  dec2 = dec[n-1];
  if (dd > 180.) dec2 -= 360.;
  else if (dd < -180.) dec2 += 360.;

  for (i=1; i<n-1; i++) {
    ra[i] = (1-t[i])*ra1 + t[i]*ra2;
    if (ra[i] > 180.) ra[i] -= 360.;
    if (ra[i] < -180.) ra[i] += 360.;
    dec[i] = (1-t[i])*dec1 + t[i]*dec2;
    if (dec[i] > 180.) dec[i] -= 360.;
    if (dec[i] < -180.) dec[i] += 360.;
    if (interps) {
      sin2psi[i] = (1-t[i])*sp1 + t[i]*sp2;
      ds = sin2psi[i] - sp1;
      cos2psi[i] = cp1*sqrt(1. - ds*(2*sp1 + ds)*icp12);
    } else {
      cos2psi[i] = (1-t[i])*cp1 + t[i]*cp2;
      dc = cos2psi[i] - cp1;
      sin2psi[i] = sp1*sqrt(1. - dc*(2*cp1 + dc)*isp12);
    }
  }
}

void qp_bore2radec(double delta_az, double delta_el, double delta_psi,
		   double *ctime, quat_t *q_bore,
		   double *ra, double *dec, double *sin2psi, 
		   double *cos2psi, int n) {
  int i;
  quat_t q_off, q;
  double *t;
  int interp = 1; // disable for now
  
  if (interp<1) 
    interp=1;
  else if (interp>1) {
    t = malloc(sizeof(double)*(interp+1));
    for (i=0; i<=interp; i++)
      t[i] = i/(double)(interp);
  }
  
  qp_init_params();
  qp_det_offset(delta_az, delta_el, delta_psi, q_off);
  
  for (i=0; i<n; i+=interp) {
    qp_bore2det(q_off, ctime[i], q_bore[i], q);
    qp_quat2radec(q, ra+i, dec+i, sin2psi+i, cos2psi+i);
    if (interp>1 && i>0)
      qp_interp_radec(&ra[i-interp], &dec[i-interp], &sin2psi[i-interp],
		      &cos2psi[i-interp], t,interp+1);
    // TODO handle chunk ends!
  }
  
  if (interp>1) free(t);
}

void qp_bore2rasindec(double delta_az, double delta_el, double delta_psi,
		      double *ctime, quat_t *q_bore,
		      double *ra, double *sindec, double *sin2psi, 
		      double *cos2psi, int n) {
  int i;
  quat_t q_off, q;
  double *t;
  int interp = 1; //disable for now
  
  if (interp<1) 
    interp=1;
  else if (interp>1) {
    t = malloc(sizeof(double)*(interp+1));
    for (i=0; i<=interp; i++)
      t[i] = i/(double)(interp);
  }
  
  qp_init_params();
  qp_det_offset(delta_az, delta_el, delta_psi, q_off);
  
  for (i=0; i<n; i+=interp) {
    qp_bore2det(q_off, ctime[i], q_bore[i], q);
    qp_quat2rasindec(q, ra+i, sindec+i, sin2psi+i, cos2psi+i);
    if (interp>1 && i>0)
      qp_interp_radec(&ra[i-interp], &sindec[i-interp], &sin2psi[i-interp],
		      &cos2psi[i-interp], t,interp+1);
    // TODO handle chunk ends!
  }
  
  if (interp>1) free(t);
}

// all input and output angles are in degrees!
void qp_azel2radec(double delta_az, double delta_el, double delta_psi,
		   double *az, double *el, double *pitch, double *roll,
		   double *lon, double *lat, double *ctime,
		   double *ra, double *dec, double *sin2psi, 
		   double *cos2psi, int n) {
  quat_t q_det, q_off, q_bore;
  
  qp_init_params();
  qp_det_offset(delta_az, delta_el, delta_psi, q_off);
  
  for (int i=0; i<n; i++) {
    qp_azel2quat(az[i], el[i], pitch[i], roll[i], lon[i], lat[i], ctime[i],
		 q_bore);
    qp_bore2det(q_off, ctime[i], q_bore, q_det);
    qp_quat2radec(q_det, &ra[i], &dec[i], &sin2psi[i], &cos2psi[i]);
  }
}

// all input and output angles are in degrees!
void qp_azel2rasindec(double delta_az, double delta_el, double delta_psi,
		      double *az, double *el, double *pitch, double *roll,
		      double *lon, double *lat, double *ctime,
		      double *ra, double *sindec, double *sin2psi, 
		      double *cos2psi, int n) {
  quat_t q_det, q_off, q_bore;
  
  qp_init_params();
  qp_det_offset(delta_az, delta_el, delta_psi, q_off);
  
  for (int i=0; i<n; i++) {
    qp_azel2quat(az[i], el[i], pitch[i], roll[i], lon[i], lat[i], ctime[i],
		 q_bore);
    qp_bore2det(q_off, ctime[i], q_bore, q_det);
    qp_quat2rasindec(q_det, &ra[i], &sindec[i], &sin2psi[i], &cos2psi[i]);
  }
}
