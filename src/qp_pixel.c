#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include "qpoint.h"
#include "fast_math.h"
#include "vec3.h"
#include "quaternion.h"
#include <chealpix.h>

/* Compute healpix pixel number for given nside and ra/dec */
long qp_radec2pix(qp_memory_t *mem, double ra, double dec, int nside) {
  long pix;
  if (mem->pix_order == QP_ORDER_NEST)
    ang2pix_nest(nside, M_PI_2 - deg2rad(dec), deg2rad(ra), &pix);
  else
    ang2pix_ring(nside, M_PI_2 - deg2rad(dec), deg2rad(ra), &pix);
  return pix;
}

void qp_radec2pixn(qp_memory_t *mem, double *ra, double *dec,
                   int nside, long *pix, int n) {
  for (int ii = 0; ii < n; ii++) {
    pix[ii] = qp_radec2pix(mem, ra[ii], dec[ii], nside);
  }
}

void qp_init_gal(qp_memory_t *mem) {
  if (mem->gal_init)
    return;

  /* galactic pole cf. sofa/g2icrs */
  double gp_ra = 192.85948;
  double gp_dec = 27.12825;
  double gp_pa = 32.93192+90;

  qp_radecpa2quat(mem, gp_ra, gp_dec, gp_pa, mem->q_gal);
  Quaternion_copy(mem->q_gal_inv, mem->q_gal);
  Quaternion_inv(mem->q_gal_inv);

  mem->gal_init = 1;
}

void qp_radec2gal(qp_memory_t *mem, double *ra, double *dec,
                  double *sin2psi, double *cos2psi) {
  quat_t q;
  qp_init_gal(mem);
  qp_radec2quat(mem, *ra, *dec, *sin2psi, *cos2psi, q);
  Quaternion_mul_left(mem->q_gal_inv, q);
  qp_quat2radec(mem, q, ra, dec, sin2psi, cos2psi);
}

void qp_radec2galn(qp_memory_t *mem, double *ra, double *dec,
                   double *sin2psi, double *cos2psi, int n) {
  for (int ii=0; ii<n; ii++) {
    qp_radec2gal(mem, ra+ii, dec+ii, sin2psi+ii, cos2psi+ii);
  }
}

void qp_gal2radec(qp_memory_t *mem, double *ra, double *dec,
                  double *sin2psi, double *cos2psi) {
  quat_t q;
  qp_init_gal(mem);
  qp_radec2quat(mem, *ra, *dec, *sin2psi, *cos2psi, q);
  Quaternion_mul_left(mem->q_gal, q);
  qp_quat2radec(mem, q, ra, dec, sin2psi, cos2psi);
}

void qp_gal2radecn(qp_memory_t *mem, double *ra, double *dec,
                   double *sin2psi, double *cos2psi, int n) {
  for (int ii=0; ii<n; ii++) {
    qp_gal2radec(mem, ra+ii, dec+ii, sin2psi+ii, cos2psi+ii);
  }
}

void qp_rotate_map(qp_memory_t *mem, int nside,
                   double **map_in, const char coord_in,
                   double **map_out, const char coord_out) {
  long npix = 12 * nside * nside;
  long pix;
  double ra, dec, sin2psi, cos2psi, norm;
  double t,q,u;

  /* check inputs */
  if (!(coord_in == 'C' || coord_in == 'G')) {
    return;
  }
  if (!(coord_out == 'C' || coord_out == 'G')) {
    return;
  }
  if (coord_in == coord_out) {
    return;
  }

  qp_init_gal(mem);

#pragma omp parallel for private(pix, ra, dec, sin2psi, cos2psi, norm, t, q, u)
  for (long ii=0; ii<npix; ii++) {
    /* ra/dec of output pixel */
    if (mem->pix_order == QP_ORDER_NEST)
      pix2ang_nest(nside, ii, &dec, &ra);
    else
      pix2ang_ring(nside, ii, &dec, &ra);
    dec = rad2deg(M_PI_2 - dec);
    ra = rad2deg(ra);
    sin2psi = 0;
    cos2psi = 1;

    /* find corresponding input pixel */
    if (coord_in == 'C' && coord_out == 'G') {
      qp_gal2radec(mem, &ra, &dec, &sin2psi, &cos2psi);
    } else if (coord_in == 'G' && coord_out == 'C') {
      qp_radec2gal(mem, &ra, &dec, &sin2psi, &cos2psi);
    }
    pix = qp_radec2pix(mem, ra, dec, nside);

    /* rotate input pixel to output pixel */
    t = map_in[0][pix];
    q = map_in[1][pix];
    u = map_in[2][pix];
    map_out[0][ii] = t;
    norm = sqrt(q*q + u*u);
    if (norm == 0) continue;
    cos2psi = q / norm; /* input pol */
    sin2psi = u / norm;
    if (coord_in == 'C' && coord_out == 'G') {
      qp_radec2gal(mem, &ra, &dec, &sin2psi, &cos2psi);
    } else if (coord_in == 'G' && coord_out == 'C') {
      qp_gal2radec(mem, &ra, &dec, &sin2psi, &cos2psi);
    }
    map_out[1][ii] = norm * cos2psi;
    map_out[2][ii] = norm * sin2psi;
  }
}

/* Compute pixel number and pol angle given nside and quaternion */
void qp_quat2pix(qp_memory_t *mem, quat_t q, int nside, long *pix,
                 double *sin2psi, double *cos2psi) {
  if (mem->fast_pix) {
    vec3_t vec;
    Quaternion_to_matrix_col3(q, vec);
    if (mem->pix_order == QP_ORDER_NEST)
      vec2pix_nest(nside, vec, pix);
    else
      vec2pix_ring(nside, vec, pix);

    double cosb2 = (1 - vec[2] * vec[2]) / 4.;
    double norm, cosg, sing;
    if (cosb2 < DBL_EPSILON) {
      if (vec[2] > 0) {
        cosg = q[0] * q[0] - q[3] * q[3];
        sing = 2 * q[0] * q[3];
      } else {
        cosg = q[2] * q[2] - q[1] * q[1];
        sing = 2 * q[1] * q[2];
      }
      norm = 2 * cosg;
    } else {
      cosg = q[0] * q[2] - q[1] * q[3];
      sing = q[0] * q[1] + q[2] * q[3];
      norm = 2. * cosg / cosb2;
    }
    if (!mem->polconv) sing = -sing;
    *sin2psi = norm * sing;
    *cos2psi = norm * cosg - 1;
  } else {
    double ra, dec;
    qp_quat2radec(mem, q, &ra, &dec, sin2psi, cos2psi);
    *pix = qp_radec2pix(mem, ra, dec, nside);
  }
}

void qp_quat2pixn(qp_memory_t *mem, quat_t *q, int nside, long *pix,
                  double *sin2psi, double *cos2psi, int n) {
  for (int ii = 0; ii < n; ii++) {
    qp_quat2pix(mem, q[ii], nside, pix+ii, sin2psi+ii, cos2psi+ii);
  }
}

void qp_bore2pix(qp_memory_t *mem, quat_t q_off, double *ctime, quat_t *q_bore,
                 int nside, long *pix, double *sin2psi, double *cos2psi, int n) {
  quat_t q;

  for (int ii = 0; ii < n; ii++) {
    qp_bore2det(mem, q_off, ctime[ii], q_bore[ii], q);
    qp_quat2pix(mem, q, nside, pix+ii, sin2psi+ii, cos2psi+ii);
  }
}

void qp_bore2pix_hwp(qp_memory_t *mem, quat_t q_off, double *ctime,
                     quat_t *q_bore, quat_t *q_hwp, int nside, long *pix,
                     double *sin2psi, double *cos2psi, int n) {
  quat_t q;

  for (int ii = 0; ii < n; ii++) {
    qp_bore2det_hwp(mem, q_off, ctime[ii], q_bore[ii], q_hwp[ii], q);
    qp_quat2pix(mem, q, nside, pix+ii, sin2psi+ii, cos2psi+ii);
  }
}
