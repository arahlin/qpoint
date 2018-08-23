#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include "qpoint.h"
#include "fast_math.h"
#include "vec3.h"
#include "quaternion.h"
#include "chealpix.h"

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
  double gp_pa = 90 + 32.93192;

  qp_radecpa2quat(mem, gp_ra, gp_dec, gp_pa, mem->q_gal);
  Quaternion_copy(mem->q_gal_inv, mem->q_gal);
  Quaternion_inv(mem->q_gal_inv);

  mem->gal_init = 1;
}

void qp_radec2gal_quat(qp_memory_t *mem, quat_t q) {
  qp_init_gal(mem);
  Quaternion_mul_left(mem->q_gal_inv, q);
}

void qp_gal2radec_quat(qp_memory_t *mem, quat_t q) {
  qp_init_gal(mem);
  Quaternion_mul_left(mem->q_gal, q);
}

void qp_radec2gal_quatn(qp_memory_t *mem, quat_t *q, int n) {
  for (int ii=0; ii<n; ii++) {
    qp_radec2gal_quat(mem, q[ii]);
  }
}

void qp_radec2gal(qp_memory_t *mem, double *ra, double *dec,
                  double *sin2psi, double *cos2psi) {
  quat_t q;
  qp_radec2quat(mem, *ra, *dec, *sin2psi, *cos2psi, q);
  qp_radec2gal_quat(mem, q);
  qp_quat2radec(mem, q, ra, dec, sin2psi, cos2psi);
}

void qp_radecpa2gal(qp_memory_t *mem, double *ra, double *dec,
                    double *pa) {
  quat_t q;
  qp_radecpa2quat(mem, *ra, *dec, *pa, q);
  qp_radec2gal_quat(mem, q);
  qp_quat2radecpa(mem, q, ra, dec, pa);
}

void qp_radec2galn(qp_memory_t *mem, double *ra, double *dec,
                   double *sin2psi, double *cos2psi, int n) {
  for (int ii=0; ii<n; ii++) {
    qp_radec2gal(mem, ra+ii, dec+ii, sin2psi+ii, cos2psi+ii);
  }
}

void qp_radecpa2galn(qp_memory_t *mem, double *ra, double *dec,
                     double *pa, int n) {
  for (int ii=0; ii<n; ii++) {
    qp_radecpa2gal(mem, ra+ii, dec+ii, pa+ii);
  }
}

void qp_gal2radec_quatn(qp_memory_t *mem, quat_t *q, int n) {
  for (int ii=0; ii<n; ii++) {
    qp_gal2radec_quat(mem, q[ii]);
  }
}

void qp_gal2radec(qp_memory_t *mem, double *ra, double *dec,
                    double *sin2psi, double *cos2psi) {
  quat_t q;
  qp_radec2quat(mem, *ra, *dec, *sin2psi, *cos2psi, q);
  qp_gal2radec_quat(mem, q);
  qp_quat2radec(mem, q, ra, dec, sin2psi, cos2psi);
}

void qp_gal2radecpa(qp_memory_t *mem, double *ra, double *dec,
                    double *pa) {
  quat_t q;
  qp_radecpa2quat(mem, *ra, *dec, *pa, q);
  qp_gal2radec_quat(mem, q);
  qp_quat2radecpa(mem, q, ra, dec, pa);
}

void qp_gal2radecn(qp_memory_t *mem, double *ra, double *dec,
                   double *sin2psi, double *cos2psi, int n) {
  for (int ii=0; ii<n; ii++) {
    qp_gal2radec(mem, ra+ii, dec+ii, sin2psi+ii, cos2psi+ii);
  }
}

void qp_gal2radecpan(qp_memory_t *mem, double *ra, double *dec,
                     double *pa, int n) {
  for (int ii=0; ii<n; ii++) {
    qp_gal2radecpa(mem, ra+ii, dec+ii, pa+ii);
  }
}

void qp_rotate_map(qp_memory_t *mem, int nside,
                   double **map_in, const char coord_in,
                   double **map_out, const char coord_out) {
  long npix = 12 * nside * nside;
  long pix, ipix[4];
  double ra, dec, sin2psi, cos2psi;
  double t, q, u, weight[4];

  /* check inputs */
  if (!(coord_in == 'C' || coord_in == 'G'))
    return;
  if (!(coord_out == 'C' || coord_out == 'G'))
    return;
  if (coord_in == coord_out)
    return;

  /* initialize */
  qp_init_gal(mem);
  qp_pixinfo_t *pixinfo = qp_init_pixinfo(nside, 1);

#pragma omp parallel for private(pix, ipix, weight, ra, dec, sin2psi, cos2psi, t, q, u)
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

    if (mem->interp_pix) {
      qp_get_interpol(mem, pixinfo, ra, dec, ipix, weight);
      /* rotate input pixel to output pixel */
      t = q = u = 0;
      for (int jj = 0; jj < 4; jj++) {
        t += map_in[0][ipix[jj]] * weight[jj];
        q += map_in[1][ipix[jj]] * weight[jj];
        u += map_in[2][ipix[jj]] * weight[jj];
      }
    } else {
      pix = qp_radec2pix(mem, ra, dec, nside);
      t = map_in[0][pix];
      q = map_in[1][pix];
      u = map_in[2][pix];
    }

    if (t == 0 && q == 0 && u == 0) continue;
    cos2psi = 1.;
    sin2psi = 0.;
    if (coord_in == 'C' && coord_out == 'G') {
      qp_radec2gal(mem, &ra, &dec, &sin2psi, &cos2psi);
    } else if (coord_in == 'G' && coord_out == 'C') {
      qp_gal2radec(mem, &ra, &dec, &sin2psi, &cos2psi);
    }
    map_out[0][ii] = t;
    map_out[1][ii] = q * cos2psi + u * sin2psi;
    map_out[2][ii] = u * cos2psi - q * sin2psi;
  }

  qp_free_pixinfo(pixinfo);
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
        cosg = q[3] * q[3] - q[0] * q[0];
        sing = 2 * q[0] * q[3];
      } else {
        cosg = q[1] * q[1] - q[2] * q[2];
        sing = 2 * q[1] * q[2];
      }
      norm = 2 * cosg;
    } else {
      cosg = q[1] * q[3] - q[0] * q[2];
      sing = q[0] * q[1] + q[2] * q[3];
      norm = 2. * cosg / cosb2;
    }
    *sin2psi = norm * sing;
    *cos2psi = norm * cosg - 1;
  } else {
    double ra, dec;
    qp_quat2radec(mem, q, &ra, &dec, sin2psi, cos2psi);
    *pix = qp_radec2pix(mem, ra, dec, nside);
  }
}

void qp_quat2pixpa(qp_memory_t *mem, quat_t q, int nside, long *pix, double *pa) {
  if (mem->fast_pix) {
    vec3_t vec;
    Quaternion_to_matrix_col3(q, vec);
    if (mem->pix_order == QP_ORDER_NEST)
      vec2pix_nest(nside, vec, pix);
    else
      vec2pix_ring(nside, vec, pix);

    double cosb2 = (1 - vec[2] * vec[2]) / 4.;
    double cosg, sing;
    if (cosb2 < DBL_EPSILON) {
      if (vec[2] > 0) {
        cosg = q[3] * q[3] - q[0] * q[0];
        sing = 2 * q[0] * q[3];
      } else {
        cosg = q[1] * q[1] - q[2] * q[2];
        sing = 2 * q[1] * q[2];
      }
    } else {
      cosg = q[1] * q[3] - q[0] * q[2];
      sing = q[0] * q[1] + q[2] * q[3];
    }

    if (mem->fast_math) {
      *pa = rad2deg(poly_atan2(sing, cosg));
    } else {
      *pa = rad2deg(atan2(sing, cosg));
    }
  } else {
    double ra, dec;
    qp_quat2radecpa(mem, q, &ra, &dec, pa);
    *pix = qp_radec2pix(mem, ra, dec, nside);
  }
}

void qp_quat2pixn(qp_memory_t *mem, quat_t *q, int nside, long *pix,
                  double *sin2psi, double *cos2psi, int n) {
  for (int ii = 0; ii < n; ii++) {
    qp_quat2pix(mem, q[ii], nside, pix+ii, sin2psi+ii, cos2psi+ii);
  }
}

void qp_quat2pixpan(qp_memory_t *mem, quat_t *q, int nside, long *pix,
                  double *pa, int n) {
  for (int ii = 0; ii < n; ii++) {
    qp_quat2pixpa(mem, q[ii], nside, pix+ii, pa+ii);
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

void qp_bore2pixpa(qp_memory_t *mem, quat_t q_off, double *ctime, quat_t *q_bore,
                   int nside, long *pix, double *pa, int n) {
  quat_t q;

  for (int ii = 0; ii < n; ii++) {
    qp_bore2det(mem, q_off, ctime[ii], q_bore[ii], q);
    qp_quat2pixpa(mem, q, nside, pix+ii, pa+ii);
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

void qp_bore2pixpa_hwp(qp_memory_t *mem, quat_t q_off, double *ctime,
                       quat_t *q_bore, quat_t *q_hwp, int nside, long *pix,
                       double *pa, int n) {
  quat_t q;

  for (int ii = 0; ii < n; ii++) {
    qp_bore2det_hwp(mem, q_off, ctime[ii], q_bore[ii], q_hwp[ii], q);
    qp_quat2pixpa(mem, q, nside, pix+ii, pa+ii);
  }
}

void qp_pixel_offset(qp_memory_t *mem, int nside, long pix,
                     double ra, double dec, double *dtheta,
                     double *dphi) {
  if (mem->pix_order == QP_ORDER_NEST)
    pix2ang_nest(nside, pix, dtheta, dphi);
  else
    pix2ang_ring(nside, pix, dtheta, dphi);
  *dtheta = M_PI_2 - deg2rad(dec) - *dtheta;
  if (*dtheta < -M_PI_2) *dtheta += M_PI;
  if (*dtheta > M_PI_2) *dtheta -= M_PI;
  *dphi = deg2rad(ra) - *dphi;
  if (*dphi < -M_PI) *dphi += M_TWOPI;
  if (*dphi > M_PI) *dphi -= M_TWOPI;
}

/* copied get_interpol from healpix-cxx, because there is no C equivalent. */

void get_ring_info2(qp_pixinfo_t *pixinfo, long iring, long *startpix,
                    long *ringpix, double *theta, int *shifted) {

  qp_ring_t *ring = pixinfo->rings + iring;

  if (!ring->init) {
    ring->idx = iring;

    long northring = \
      (iring > (2 * pixinfo->nside)) ? (4 * pixinfo->nside - iring) : iring;

    if (northring < pixinfo->nside) {
      double tmp = northring * northring * pixinfo->fact2;
      double costheta = 1 - tmp;
      double sintheta = sqrt(tmp * (2 - tmp));
      ring->theta = atan2(sintheta, costheta);
      ring->ringpix = 4 * northring;
      ring->shifted = 1;
      ring->startpix = 2 * northring * (northring - 1);
    } else {
      ring->theta = acos((2 * pixinfo->nside - northring) * pixinfo->fact1);
      ring->ringpix = 4 * pixinfo->nside;
      ring->shifted = (((northring - pixinfo->nside) & 1) == 0);
      ring->startpix = pixinfo->ncap + (northring - pixinfo->nside) * ring->ringpix;
    }

    if (northring != iring) {
      ring->theta = M_PI - ring->theta;
      ring->startpix = pixinfo->npix - ring->startpix - ring->ringpix;
    }

    ring->init = 1;
  }

  if (theta != NULL)
    *theta = ring->theta;
  if (ringpix != NULL)
    *ringpix = ring->ringpix;
  if (shifted != NULL)
    *shifted = ring->shifted;
  if (startpix != NULL)
    *startpix = ring->startpix;
}

qp_pixinfo_t * qp_init_pixinfo(size_t nside, int populate) {
  qp_pixinfo_t *pixinfo = malloc(sizeof(*pixinfo));
  pixinfo->nside = nside;
  pixinfo->npface = pixinfo->nside * pixinfo->nside;
  pixinfo->npix = 12 * pixinfo->npface;
  pixinfo->ncap = (pixinfo->npface - pixinfo->nside) << 1;
  pixinfo->fact2 = 4. / pixinfo->npix;
  pixinfo->fact1 = (pixinfo->nside << 1) * pixinfo->fact2;
  pixinfo->rings = calloc(4 * pixinfo->nside, sizeof(qp_ring_t));
  pixinfo->rings_init = QP_ARR_MALLOC_1D;
  pixinfo->init = QP_STRUCT_INIT | QP_STRUCT_MALLOC;

  if (populate) {
    for (int iring = 0; iring < 4 * pixinfo->nside; iring++)
      get_ring_info2(pixinfo, iring, NULL, NULL, NULL, NULL);
  }
  return pixinfo;
}

void qp_free_pixinfo(qp_pixinfo_t *pixinfo) {
  if (pixinfo->rings_init & QP_ARR_MALLOC_1D)
    free(pixinfo->rings);
  if (pixinfo->init & QP_STRUCT_MALLOC)
    free(pixinfo);
  else
    memset(pixinfo, 0, sizeof(*pixinfo));
}

int get_interpol_ring(qp_pixinfo_t *pixinfo, double theta, double phi,
                      long pix[4], double weight[4]) {
  if (theta < 0 || theta > M_PI)
    return -1;

  int nside = pixinfo->nside;
  long npix = pixinfo->npix;
  double z = cos(theta);
  double az = fabs(z);

  long ir1;
  if (az < 2. / 3.)
    ir1 = nside * (2 - 1.5*z);
  else {
    ir1 = nside * sqrt(3 * (1 - az));
    if (z <= 0)
      ir1 = 4 * nside - ir1 - 1;
  }
  long ir2 = ir1 + 1;

  double theta1=0, theta2=0, w1, tmp, dphi;
  long sp, nr;
  int shift;
  long i1, i2;

  if (ir1 > 0) {
    get_ring_info2(pixinfo, ir1, &sp, &nr, &theta1, &shift);
    dphi = M_TWOPI / nr;
    tmp = phi / dphi - 0.5 * shift;
    i1 = (tmp < 0) ? tmp - 1 : tmp;
    w1 = (phi - (i1 + 0.5 * shift) * dphi) / dphi;
    if (i1 < 0) i1 += nr;
    i2 = i1 + 1;
    if (i2 >= nr) i2 -= nr;
    pix[0] = sp + i1;
    pix[1] = sp + i2;
    weight[0] = 1 - w1;
    weight[1] = w1;
  }

  if (ir2 < (4 * nside)) {
    get_ring_info2(pixinfo, ir2, &sp, &nr, &theta2, &shift);
    dphi = M_TWOPI / nr;
    tmp = phi / dphi - 0.5 * shift;
    i1 = (tmp < 0) ? tmp - 1 : tmp;
    w1 = (phi - (i1 + 0.5 * shift) * dphi) / dphi;
    if (i1 < 0) i1 += nr;
    i2 = i1 + 1;
    if (i2 >= nr) i2 -= nr;
    pix[2] = sp + i1;
    pix[3] = sp + i2;
    weight[2] = 1 - w1;
    weight[3] = w1;
  }

  if (ir1 == 0) {
    double wtheta = theta / theta2;
    weight[2] *= wtheta;
    weight[3] *= wtheta;
    double fac = (1 - wtheta) * 0.25;
    weight[0] = fac;
    weight[1] = fac;
    weight[2] += fac;
    weight[3] += fac;
    pix[0] = (pix[2] + 2) & 3;
    pix[1] = (pix[3] + 2) & 3;
  } else if (ir2 == 4 * nside) {
    double wtheta = (theta - theta1)/(M_PI - theta1);
    weight[0] *= (1 - wtheta);
    weight[1] *= (1 - wtheta);
    double fac = wtheta * 0.25;
    weight[0] += fac;
    weight[1] += fac;
    weight[2] = fac;
    weight[3] = fac;
    pix[2] = ((pix[0] + 2) & 3) + npix - 4;
    pix[3] = ((pix[1] + 2) & 3) + npix - 4;
  } else {
    double wtheta = (theta - theta1) / (theta2 - theta1);
    weight[0] *= (1 - wtheta);
    weight[1] *= (1 - wtheta);
    weight[2] *= wtheta;
    weight[3] *= wtheta;
  }

  return 0;
}

int get_interpol_nest(qp_pixinfo_t *pixinfo, double theta, double phi,
                      long pix[4], double weight[4]) {
  if (get_interpol_ring(pixinfo, theta, phi, pix, weight))
    return -1;
  for (int ii = 0; ii < 4; ii++) {
    ring2nest(pixinfo->nside, pix[ii], pix + ii);
  }
  return 0;
}

double get_interp_val_ring(qp_pixinfo_t *pixinfo, double *map,
                           double theta, double phi) {
  long pix[4];
  double weight[4];
  get_interpol_ring(pixinfo, theta, phi, pix, weight);

  double val = 0;
  for (int ii = 0; ii < 4; ii++)
    val += map[pix[ii]] * weight[ii];
  return val;
}

double get_interp_val_nest(qp_pixinfo_t *pixinfo, double *map,
                           double theta, double phi) {
  long pix[4];
  double weight[4];
  get_interpol_nest(pixinfo, theta, phi, pix, weight);

  double val = 0;
  for (int ii = 0; ii < 4; ii++)
    val += map[pix[ii]] * weight[ii];
  return val;
}

void qp_get_interpol(qp_memory_t *mem, qp_pixinfo_t *pixinfo, double ra,
                     double dec, long pix[4], double weight[4]) {
  if (mem->pix_order == QP_ORDER_RING) {
    get_interpol_ring(pixinfo, M_PI_2 - deg2rad(dec), deg2rad(ra),
                      pix, weight);
  } else {
    get_interpol_nest(pixinfo, M_PI_2 - deg2rad(dec), deg2rad(ra),
                      pix, weight);
  }
}

double qp_get_interp_val(qp_memory_t *mem, qp_pixinfo_t *pixinfo, double *map,
                         double ra, double dec) {
  if (mem->pix_order == QP_ORDER_RING) {
    return get_interp_val_ring(pixinfo, map, M_PI_2 - deg2rad(dec),
                               deg2rad(ra));
  } else {
    return get_interp_val_nest(pixinfo, map, M_PI_2 - deg2rad(dec),
                               deg2rad(ra));
  }
}

void qp_get_interp_valn(qp_memory_t *mem, int nside, double *map, double *ra,
                        double *dec, double *val, int n) {
  qp_pixinfo_t *pixinfo = qp_init_pixinfo(nside, 0);

  for (int ii = 0; ii < n; ii++) {
    val[ii] = qp_get_interp_val(mem, pixinfo, map, ra[ii], dec[ii]);
  }

  qp_free_pixinfo(pixinfo);
}
