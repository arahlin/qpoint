#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include "sofa.h"
#ifdef ENABLE_SLAREFRO
#include "slarefro.h"
#endif
#include "qpoint.h"
#include "fast_math.h"
#include "vec3.h"
#include "quaternion.h"

#define _unused(x) ((void)x)

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

  ctime2jd(ctime, jd_utc);
  iauUtctai(jd_utc[0], jd_utc[1], &jd_tai[0], &jd_tai[1]);
  iauTaitt(jd_tai[0], jd_tai[1], &jd_tt[0], &jd_tt[1]);
}

void jdutc2jdut1(double jd_utc[2], double dut1, double jd_ut1[2]) {
  iauUtcut1(jd_utc[0], jd_utc[1], dut1, &jd_ut1[0], &jd_ut1[1]);
}

double ctime2gmst(double ctime, double dut1, int accuracy) {
  double jd_utc[2], jd_ut1[2], jd_tt[2];

  ctime2jd(ctime, jd_utc);

  if (!accuracy) {
    jdutc2jdut1(jd_utc, dut1, jd_ut1);
    ctime2jdtt(ctime, jd_tt);
    return iauGmst00(jd_ut1[0], jd_ut1[1], jd_tt[0], jd_tt[1]);
  } else {
    return iauGmst00(jd_utc[0], jd_utc[1], jd_utc[0], jd_utc[1]);
  }
}

double qp_gmst(qp_memory_t *mem, double ctime) {
  double jd_utc[2];
  ctime2jd(ctime, jd_utc);
  double mjd_utc = jd2mjd(jd_utc[0]) + jd_utc[1];
  double x,y,gmst;

  if (mem->accuracy == 0) {
    if (qp_check_update(&mem->state_dut1, ctime)) {
      qp_get_iers_bulletin_a(mem, mjd_utc, &mem->dut1, &x, &y);
    }
    gmst = ctime2gmst(ctime, mem->dut1, mem->accuracy);
  } else {
    gmst = ctime2gmst(ctime, 0, mem->accuracy);
  }
  return fmod(rad2deg(gmst) / 15.0, 24.);
}

void qp_gmstn(qp_memory_t *mem, double *ctime, double *gmst, int n) {
  for (int ii=0; ii<n; ii++) {
    gmst[ii] = qp_gmst(mem, ctime[ii]);
  }
}

double qp_lmst(qp_memory_t *mem, double ctime, double lon) {
  double gmst = qp_gmst(mem, ctime);
  return fmod(gmst + lon / 15.0, 24.);
}

void qp_lmstn(qp_memory_t *mem, double *ctime, double *lon, double *lmst, int n) {
  for (int ii=0; ii<n; ii++) {
    lmst[ii] = qp_lmst(mem, ctime[ii], lon[ii]);
  }
}

#define DIPOLE_RA  167.987505
#define DIPOLE_DEC -7.22

void qp_dipole_init(qp_memory_t *mem) {
  if (mem->dipole_init)
    return;

  /* dipole angles */
  quat_t q;
  qp_radecpa2quat(mem, DIPOLE_RA, DIPOLE_DEC, 0., q);
  Quaternion_to_matrix_col3(q, mem->v_dipole);
  mem->dipole_init = 1;
}

double cdist2dipole(qp_memory_t *mem, double cdist, double ctime) {
  const double tcmb = 2.725;
  const double beta = 369e3 / C_MS;
  const double vhelio = 0.00027;
  const double dipole_epoch = 2451170;

  double out = tcmb * beta * (cdist + beta/2. * (2 * cdist * cdist - 1));

  double jd[2];
  ctime2jd(ctime, jd);
  double delta = (jd[1] + (jd[0] - dipole_epoch)) / 365.25;

  if (mem->fast_math) {
    out += vhelio * poly_cos(2 * M_PI * delta);
  } else {
    out += vhelio * cos(2 * M_PI * delta);
  }

  return out;
}

double qp_quat2dipole(qp_memory_t *mem, double ctime, quat_t q) {
  qp_dipole_init(mem);

  vec3_t v;
  Quaternion_to_matrix_col3(q, v);
  double cdist = vec3_dot_product(mem->v_dipole, v);
  return cdist2dipole(mem, cdist, ctime);
}

void qp_bore2dipole(qp_memory_t *mem, quat_t q_off, double *ctime,
                    quat_t *q_bore, double *dipole, int n) {
  quat_t q_det;

  qp_dipole_init(mem);

  for (int ii = 0; ii < n; ii++) {
    qp_bore2det(mem, q_off, ctime[ii], q_bore[ii], q_det);
    dipole[ii] = qp_quat2dipole(mem, ctime[ii], q_det);
  }
}

double qp_dipole(qp_memory_t *mem, double ctime, double ra, double dec) {

  const double dipole_phi = deg2rad(DIPOLE_RA);
  const double dipole_theta = M_PI/2 - deg2rad(DIPOLE_DEC);
  const double sdtheta = sin(dipole_theta);
  const double cdtheta = cos(dipole_theta);

  double theta = M_PI/2 - deg2rad(dec);
  double phi = deg2rad(ra);

  double stheta, ctheta, cdphi;
  if (mem->fast_math) {
    stheta = poly_sin(theta);
    ctheta = poly_cos(theta);
    cdphi = poly_cos(dipole_phi - phi);
  } else {
    stheta = sin(theta);
    ctheta = cos(theta);
    cdphi = cos(dipole_phi - phi);
  }

  double cdist = cdtheta * ctheta + sdtheta * stheta * cdphi;
  return cdist2dipole(mem, cdist, ctime);
}

void qp_dipolen(qp_memory_t *mem, double *ctime, double *ra, double *dec,
                double *dipole, int n) {
  for (int ii=0; ii<n; ii++) {
    dipole[ii] = qp_dipole(mem, ctime[ii], ra[ii], dec[ii]);
  }
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
}

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
  if (pitch != 0)
    Quaternion_r2_mul(-deg2rad(pitch), q);
  if (roll != 0)
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

void qp_wobble_quat(double jd_tt[2], double xp, double yp, quat_t q) {
  double sprime = iauSp00(jd_tt[0], jd_tt[1]);

  Quaternion_r1(q, -arcsec2rad(yp));
  Quaternion_r2_mul(-arcsec2rad(xp), q);
  Quaternion_r3_mul(sprime, q);
}

/* Calculate atmospheric refraction */
double qp_refraction(double el, double lat, double height, double temp,
		     double press, double hum, double freq, double lapse,
		     double tol) {
  double ref = 0;
#ifdef ENABLE_SLAREFRO
  slaf_refro(M_PI_2 - deg2rad(el),
	     height, temp + 273.15, // temperature, K
	     press, hum, C_MS * 1e-3 / freq, // wavelength, um
	     deg2rad(lat), lapse, tol, // precision, radians
	     &ref);
#else
  double A, B;
  iauRefco(press, temp, hum, C_MS * 1e-3 / freq, &A, &B);
  double tz = tan(M_PI_2 - deg2rad(el));
  ref = tz * (A + B * tz * tz);
#endif
  return rad2deg(ref);
}

double qp_update_ref(qp_memory_t *mem, quat_t q, double lat) {
  qp_weather_t *W = &mem->weather;
  double el;

  if (mem->fast_math)
    el = rad2deg(poly_asin(q[0]*q[0] - q[1]*q[1] - q[2]*q[2] + q[3]*q[3]));
  else
    el = rad2deg(asin(q[0]*q[0] - q[1]*q[1] - q[2]*q[2] + q[3]*q[3]));

  double ref = qp_refraction(el, lat, W->height, W->temperature,
                             W->pressure, W->humidity, W->frequency,
                             W->lapse_rate, mem->ref_tol);
  mem->ref_delta = ref;
  return ref;
}

void qp_apply_refraction(qp_memory_t *mem, double ctime, double lat, quat_t q) {
  if (qp_check_update(&mem->state_ref, ctime)) {
    double delta = qp_update_ref(mem, q, lat);
    quat_t q_delta;
    Quaternion_r2(q_delta, -deg2rad(delta));
    Quaternion_copy(mem->q_ref, q_delta);
#ifdef DEBUG
    qp_print_quat("ref", q_delta);
#endif
  }
  if (qp_check_apply(&mem->state_ref)) {
    Quaternion_mul_right(q, mem->q_ref);
#ifdef DEBUG
    qp_print_quat("state ref", q);
#endif
  }
}

void qp_apply_diurnal_aberration(qp_memory_t *mem, double ctime, double lat,
                                 quat_t q) {
  double clat;
  quat_t q_aber;

  if (qp_check_update(&mem->state_daber, ctime)) {
    if (mem->fast_math)
      clat = poly_cos(deg2rad(lat));
    else
      clat = cos(deg2rad(lat));
    mem->beta_rot[0] = mem->beta_rot[2] = 0;
    mem->beta_rot[1] = -clat * D_ABER_RAD;
  }
  if (qp_check_apply(&mem->state_daber)) {
    qp_aberration(q, (double *)mem->beta_rot, q_aber);
    Quaternion_mul_left(q_aber, q);
#ifdef DEBUG
    qp_print_quat("daber", q_aber);
    qp_print_quat("state daber", q);
#endif
  }
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
    qp_print_quat("aaber", q_aber);
    qp_print_quat("state aaber", q);
#endif
  }
}

void qp_azel2quat(qp_memory_t *mem, double az, double el, double pitch,
		  double roll, double lon, double lat, double ctime,
		  quat_t q) {

  double jd_utc[2], jd_tt[2] = {0,0}, jd_ut1[2], mjd_utc;
  double x,y;
  quat_t q_step;

  // deal with times
  ctime2jd(ctime, jd_utc);

#ifdef DEBUG
  qp_print_memory(mem);
#endif

#ifdef DEBUG
  printf("ctime %f, jd_utc %f %f\n", ctime, jd_utc[0], jd_utc[1]);
#endif

#ifdef DEBUG
  qp_print_quat("state init", q);
#endif

  // apply boresight rotation
  qp_azel_quat(az, el, pitch, roll, q_step);
  Quaternion_mul_left(q_step, q);
#ifdef DEBUG
  qp_print_quat("azel", q_step);
  qp_print_quat("state azel", q);
#endif

  // apply refraction correction
  // NB: per-detector refraction is not fully implemented!
  // can only be done properly using the azel2radec functions
  // otherwise, this is treated as a mean correction
  qp_apply_refraction(mem, ctime, lat, q);

  // apply diurnal aberration
  // NB: same issue as refraction
  // TODO propagate this to aberration step?
  qp_apply_diurnal_aberration(mem, ctime, lat, q);

  // rotate to ITRS (by lon/lat)
  if (qp_check_update(&mem->state_lonlat, ctime)) {
    qp_lonlat_quat(lon, lat, mem->q_lonlat);
#ifdef DEBUG
    qp_print_quat("lonlat", mem->q_lonlat);
#endif
  }
  if (qp_check_apply(&mem->state_lonlat)) {
    Quaternion_mul_left(mem->q_lonlat, q);
#ifdef DEBUG
    qp_print_quat("state lonlat", q);
#endif
  }

  // apply wobble correction (polar motion)
  // or get dut1 from IERS bulletin
  mjd_utc = jd2mjd(jd_utc[0]) + jd_utc[1];
  if (qp_check_update(&mem->state_wobble, ctime)) {
    qp_get_iers_bulletin_a(mem, mjd_utc, &mem->dut1, &x, &y);
    ctime2jdtt(ctime, jd_tt);
    qp_wobble_quat(jd_tt, x, y, mem->q_wobble);
#ifdef DEBUG
    qp_print_quat("wobble", mem->q_wobble);
#endif
  } else if (qp_check_update(&mem->state_dut1, ctime))
    qp_get_iers_bulletin_a(mem, mjd_utc, &mem->dut1, &x, &y);
  if (qp_check_apply(&mem->state_wobble)) {
    Quaternion_mul_left(mem->q_wobble, q);
#ifdef DEBUG
    qp_print_quat("state wobble", q);
#endif
  }

  // apply earth rotation
  if (qp_check_update(&mem->state_erot, ctime)) {
    // get ut1
    jdutc2jdut1(jd_utc, mem->dut1, jd_ut1);
    qp_erot_quat(jd_ut1, mem->q_erot);
#ifdef DEBUG
    qp_print_quat("erot", mem->q_erot);
#endif
  }
  if (qp_check_apply(&mem->state_erot)) {
    Quaternion_mul_left(mem->q_erot, q);
#ifdef DEBUG
    qp_print_quat("state erot", q);
#endif
  }

  // apply nutation/precession/frame bias correction
  if (qp_check_update(&mem->state_npb, ctime)) {
    if (jd_tt[0] == 0) ctime2jdtt(ctime, jd_tt);
    qp_npb_quat(jd_tt, mem->q_npb, mem->accuracy);
#ifdef DEBUG
    qp_print_quat("npb", mem->q_npb);
#endif
  }
  if (qp_check_apply(&mem->state_npb)) {
    Quaternion_mul_left(mem->q_npb, q);
#ifdef DEBUG
    qp_print_quat("state npb", q);
#endif
  }

  // apply annual aberration
  // ~20 arcsec max
  if (mem->mean_aber)
    qp_apply_annual_aberration(mem, ctime, q);

#ifdef DEBUG
  qp_print_quat("state final", q);
#endif
}

void qp_azel2bore(qp_memory_t *mem, double *az, double *el, double *pitch,
		  double *roll, double *lon, double *lat, double *ctime,
		  quat_t *q, int n) {
  for (int i=0; i<n; i++)
    qp_azel2quat(mem, az[i], el[i], (pitch == NULL) ? 0 : pitch[i],
                 (roll == NULL) ? 0 : roll[i], lon[i], lat[i], ctime[i], q[i]);
}

void qp_hwp_quat(double ang, quat_t q) {
  Quaternion_r3(q, -2.*deg2rad(ang));  // rotate psi by 2*theta!
#ifdef DEBUG
  qp_print_quat("hwp", q);
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
  qp_print_quat("offset", q);
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

  double q00p33 = q[0]*q[0] + q[3]*q[3];
  double q11p22 = q[1]*q[1] + q[2]*q[2];
  // NB: factors of two have been redistributed...
  double cosb2 = q00p33*q11p22;
  double sinb = q00p33 - q11p22;
  double norm, sing, cosg;

  if (cosb2 < DBL_EPSILON) {
    *ra = 0;

    if (sinb > 0) {
      cosg = q[0] * q[0] - q[3] * q[3];
      sing = 2 * q[0] * q[3];
    } else {
      cosg = q[2] * q[2] - q[1] * q[1];
      sing = 2 * q[1] * q[2];
    }
    norm = 2. * cosg;
  } else {
    double q01 = q[0]*q[1];
    double q02 = q[0]*q[2];
    double q13 = q[1]*q[3];
    double q23 = q[2]*q[3];
    double sina_2 = q23 - q01;
    double cosa_2 = q02 + q13;

    if (mem->fast_math)
      *ra = rad2deg(poly_atan2(sina_2, cosa_2));
    else
      *ra = rad2deg(atan2(sina_2, cosa_2));

    sing = q01 + q23;
    cosg = q02 - q13;
    norm = 2. * cosg / cosb2;
  }

  *sindec = sinb;

  // cosmo (healpix) or IAU polarization convention?
  if (!mem->polconv) sing = -sing;

  *sin2psi = norm * sing;
  *cos2psi = norm * cosg - 1.;
}

void qp_quat2radecpa(qp_memory_t *mem, quat_t q, double *ra, double *dec,
                     double *pa) {
  // ZYZ euler angles

  double q00p33 = q[0]*q[0] + q[3]*q[3];
  double q11p22 = q[1]*q[1] + q[2]*q[2];
  // NB: factors of 2 have been redistributed...
  double cosb2 = q00p33*q11p22;
  double sinb_2 = 0.5*(q00p33 - q11p22);
  double sing, cosg;

  if (cosb2 < DBL_EPSILON) {
    *ra = 0;

    if (sinb_2 > 0) {
      *dec = 90;
      cosg = q[0] * q[0] - q[3] * q[3];
      sing = 2 * q[0] * q[3];
    } else {
      *dec = -90;
      cosg = q[2] * q[2] - q[1] * q[1];
      sing = 2 * q[1] * q[2];
    }
  } else {
    double q01 = q[0]*q[1];
    double q02 = q[0]*q[2];
    double q13 = q[1]*q[3];
    double q23 = q[2]*q[3];
    double sina_2 = q23 - q01;
    double cosa_2 = q02 + q13;
    double cosb_2 = sqrt(cosb2); // hypot(sina,cosa);

    if (mem->fast_math) {
      *ra = rad2deg(poly_atan2(sina_2, cosa_2));
      *dec = rad2deg(poly_atan2(sinb_2, cosb_2));
    } else {
      *ra = rad2deg(atan2(sina_2, cosa_2));
      *dec = rad2deg(atan2(sinb_2, cosb_2));
    }

    sing = q01 + q23;
    cosg = q02 - q13;
  }

  // cosmo (healpix) or IAU polarization convention?
  if (!mem->polconv) cosg = -cosg;

  if (mem->fast_math) {
    *pa = rad2deg(poly_atan2(sing, cosg));
  } else {
    *pa = rad2deg(atan2(sing, cosg));
  }
}

void qp_quat2radecpan(qp_memory_t *mem, quat_t *q, double *ra, double *dec,
                     double *pa, int n) {
  for (int ii = 0; ii < n; ii++) {
    qp_quat2radecpa(mem, q[ii], ra+ii, dec+ii, pa+ii);
  }
}

void qp_quat2radec(qp_memory_t *mem, quat_t q, double *ra, double *dec,
		   double *sin2psi, double *cos2psi) {

  // ZYZ euler angles

  double q00p33 = q[0]*q[0] + q[3]*q[3];
  double q11p22 = q[1]*q[1] + q[2]*q[2];
  // NB: factors of 2 have been redistributed...
  double cosb2 = q00p33*q11p22;
  double sinb_2 = 0.5*(q00p33 - q11p22);
  double norm, sing, cosg;

  if (cosb2 < DBL_EPSILON) {
    *ra = 0;

    if (sinb_2 > 0) {
      *dec = 90;
      cosg = q[0] * q[0] - q[3] * q[3];
      sing = 2 * q[0] * q[3];
    } else {
      *dec = -90;
      cosg = q[2] * q[2] - q[1] * q[1];
      sing = 2 * q[1] * q[2];
    }
    norm = 2. * cosg;
  } else {
    double q01 = q[0]*q[1];
    double q02 = q[0]*q[2];
    double q13 = q[1]*q[3];
    double q23 = q[2]*q[3];
    double sina_2 = q23 - q01;
    double cosa_2 = q02 + q13;
    double cosb_2 = sqrt(cosb2); // hypot(sina,cosa);

    if (mem->fast_math) {
      *ra = rad2deg(poly_atan2(sina_2, cosa_2));
      *dec = rad2deg(poly_atan2(sinb_2, cosb_2));
    } else {
      *ra = rad2deg(atan2(sina_2, cosa_2));
      *dec = rad2deg(atan2(sinb_2, cosb_2));
    }

    sing = q01 + q23;
    cosg = q02 - q13;
    norm = 2. * cosg / cosb2;
  }

  // cosmo (healpix) or IAU polarization convention?
  if (!mem->polconv) sing = -sing;

  *sin2psi = norm * sing;
  *cos2psi = norm * cosg - 1.;
}

void qp_radec2quat(qp_memory_t *mem, double ra, double dec, double sin2psi,
                   double cos2psi, quat_t q) {
  double ang;
  if (!mem->polconv)
    sin2psi = -sin2psi;
  if (mem->fast_math)
    ang = poly_atan2(sin2psi, cos2psi + 1);
  else
    ang = atan2(sin2psi, cos2psi + 1);
  Quaternion_r3(q, ang);
  Quaternion_r2_mul(M_PI_2 - deg2rad(dec), q);
  Quaternion_r3_mul(deg2rad(ra), q);
}

void qp_radecpa2quat(qp_memory_t *mem, double ra, double dec, double pa, quat_t q) {
  if (!mem->polconv)
    pa = 180. - pa;
  Quaternion_r3(q, deg2rad(pa));
  Quaternion_r2_mul(M_PI_2 - deg2rad(dec), q);
  Quaternion_r3_mul(deg2rad(ra), q);
}

void qp_radecpa2quatn(qp_memory_t *mem, double *ra, double *dec, double *pa,
                      quat_t *q, int n) {
  for (int ii=0; ii < n; ii++) {
    qp_radecpa2quat(mem, ra[ii], dec[ii], pa[ii], q[ii]);
  }
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

// NB: for all azel2radec functions below:
// since the complete offset -> ra/dec operation is done in one go here,
// we should not ignore annual aberration.  in this case,
// enabling mean_aber actually applies the correction at the detector,
// since q_off is propagated all the way through.

// all input and output angles are in degrees!
void qp_azel2radec(qp_memory_t *mem,
		   double delta_az, double delta_el, double delta_psi,
		   double *az, double *el, double *pitch, double *roll,
		   double *lon, double *lat, double *ctime,
		   double *ra, double *dec, double *sin2psi,
		   double *cos2psi, int n) {
  quat_t q_det, q_off;
  int mean_aber = qp_get_opt_mean_aber(mem);
  qp_set_opt_mean_aber(mem, 1);

  qp_det_offset(delta_az, delta_el, delta_psi, q_off);

  for (int i=0; i<n; i++) {
    Quaternion_copy(q_det, q_off);
    qp_azel2quat(mem, az[i], el[i], (pitch == NULL) ? 0 : pitch[i],
                 (roll == NULL) ? 0 : roll[i], lon[i], lat[i], ctime[i],
                 q_det);
    qp_quat2radec(mem, q_det, &ra[i], &dec[i], &sin2psi[i], &cos2psi[i]);
  }

  qp_set_opt_mean_aber(mem, mean_aber);
}

// all input and output angles are in degrees!
void qp_azel2radec_hwp(qp_memory_t *mem,
		       double delta_az, double delta_el, double delta_psi,
		       double *az, double *el, double *pitch, double *roll,
		       double *lon, double *lat, double *ctime, double *hwp,
		       double *ra, double *dec, double *sin2psi,
		       double *cos2psi, int n) {
  quat_t q_det, q_off, q_hwp;
  int mean_aber = qp_get_opt_mean_aber(mem);
  qp_set_opt_mean_aber(mem, 1);

  qp_det_offset(delta_az, delta_el, delta_psi, q_off);

  for (int i=0; i<n; i++) {
    Quaternion_copy(q_det, q_off);
    qp_hwp_quat(hwp[i], q_hwp);
    Quaternion_mul_right(q_det, q_hwp);
    qp_azel2quat(mem, az[i], el[i], (pitch == NULL) ? 0 : pitch[i],
                 (roll == NULL) ? 0 : roll[i], lon[i], lat[i], ctime[i],
                 q_det);
    qp_quat2radec(mem, q_det, &ra[i], &dec[i], &sin2psi[i], &cos2psi[i]);
  }

  qp_set_opt_mean_aber(mem, mean_aber);
}

// all input and output angles are in degrees!
void qp_azel2rasindec(qp_memory_t *mem,
		      double delta_az, double delta_el, double delta_psi,
		      double *az, double *el, double *pitch, double *roll,
		      double *lon, double *lat, double *ctime,
		      double *ra, double *sindec, double *sin2psi,
		      double *cos2psi, int n) {
  quat_t q_det, q_off;
  int mean_aber = qp_get_opt_mean_aber(mem);
  qp_set_opt_mean_aber(mem, 1);

  qp_det_offset(delta_az, delta_el, delta_psi, q_off);

  for (int i=0; i<n; i++) {
    Quaternion_copy(q_det, q_off);
    qp_azel2quat(mem, az[i], el[i], (pitch == NULL) ? 0 : pitch[i],
                 (roll == NULL) ? 0 : roll[i], lon[i], lat[i], ctime[i],
                 q_det);
    qp_quat2rasindec(mem, q_det, &ra[i], &sindec[i], &sin2psi[i], &cos2psi[i]);
  }

  qp_set_opt_mean_aber(mem, mean_aber);
}

// all input and output angles are in degrees!
void qp_azel2rasindec_hwp(qp_memory_t *mem,
			  double delta_az, double delta_el, double delta_psi,
			  double *az, double *el, double *pitch, double *roll,
			  double *lon, double *lat, double *ctime, double *hwp,
			  double *ra, double *sindec, double *sin2psi,
			  double *cos2psi, int n) {
  quat_t q_det, q_off, q_hwp;
  int mean_aber = qp_get_opt_mean_aber(mem);
  qp_set_opt_mean_aber(mem, 1);

  qp_det_offset(delta_az, delta_el, delta_psi, q_off);

  for (int i=0; i<n; i++) {
    Quaternion_copy(q_det, q_off);
    qp_hwp_quat(hwp[i], q_hwp);
    Quaternion_mul_right(q_det, q_hwp);
    qp_azel2quat(mem, az[i], el[i], (pitch == NULL) ? 0 : pitch[i],
                 (roll == NULL) ? 0 : roll[i], lon[i], lat[i], ctime[i],
                 q_det);
    qp_quat2rasindec(mem, q_det, &ra[i], &sindec[i], &sin2psi[i], &cos2psi[i]);
  }

  qp_set_opt_mean_aber(mem, mean_aber);
}
