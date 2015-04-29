#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include "qpoint.h"
#include "fast_math.h"
#include "vec3.h"
#include "quaternion.h"
#include <chealpix.h>
#ifdef _OPENMP
#include <omp.h>
#endif

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
                   vec3_t *map_in, const char coord_in,
                   vec3_t *map_out, const char coord_out) {
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
    t = map_in[pix][0];
    q = map_in[pix][1];
    u = map_in[pix][2];
    map_out[ii][0] = t;
    norm = sqrt(q*q + u*u);
    if (norm == 0) continue;
    cos2psi = q / norm; /* input pol */
    sin2psi = u / norm;
    if (coord_in == 'C' && coord_out == 'G') {
      qp_radec2gal(mem, &ra, &dec, &sin2psi, &cos2psi);
    } else if (coord_in == 'G' && coord_out == 'C') {
      qp_gal2radec(mem, &ra, &dec, &sin2psi, &cos2psi);
    }
    map_out[ii][1] = norm * cos2psi;
    map_out[ii][2] = norm * sin2psi;
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

    double cosb2 = (1 - vec[2]*vec[2]) / 4.;
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

/* Compute pointing matrix map for given boresight timestream and detector
   offset. pmap is a 6-x-npix array containing (hits, p01, p02, p11, p12, p22) */
void qp_tod2map_pnt_single(qp_memory_t *mem, quat_t q_off,
                           double *ctime, quat_t *q_bore, int n,
                           vec6_t *pmap, int nside) {
  double sin2psi, cos2psi;
  long ipix;
  quat_t q;

  for (int ii=0; ii<n; ii++) {
    qp_bore2det(mem, q_off, ctime[ii], q_bore[ii], q);
    qp_quat2pix(mem, q, nside, &ipix, &sin2psi, &cos2psi);
    pmap[ipix][0] += 1;
    pmap[ipix][1] += cos2psi;
    pmap[ipix][2] += sin2psi;
    pmap[ipix][3] += cos2psi*cos2psi;
    pmap[ipix][4] += cos2psi*sin2psi;
    pmap[ipix][5] += sin2psi*sin2psi;
  }
}

void qp_tod2map_pnt_nopol_single(qp_memory_t *mem, quat_t q_off,
                                 double *ctime, quat_t *q_bore,
                                 int n, double *pmap, int nside) {
  double dum;
  long ipix;
  quat_t q;

  for (int ii=0; ii<n; ii++) {
    qp_bore2det(mem, q_off, ctime[ii], q_bore[ii], q);
    qp_quat2pix(mem, q, nside, &ipix, &dum, &dum);
    pmap[ipix] += 1;
  }
}

/* Compute signal map for given boresight timestream, signal timestream,
   and detector offset.
   smap is a npix-x-3 array containing (d, d*cos(2 psi), d*sin(2 psi)). */
void qp_tod2map_sig_single(qp_memory_t *mem, quat_t q_off,
                           double *ctime, quat_t *q_bore, double *tod, int n,
                           vec3_t *smap, int nside) {
  double sin2psi, cos2psi;
  long ipix;
  quat_t q;

  for (int ii=0; ii<n; ii++) {
    qp_bore2det(mem, q_off, ctime[ii], q_bore[ii], q);
    qp_quat2pix(mem, q, nside, &ipix, &sin2psi, &cos2psi);
    if (tod[ii] != tod[ii]) continue;
    smap[ipix][0] += tod[ii];
    smap[ipix][1] += tod[ii] * cos2psi;
    smap[ipix][2] += tod[ii] * sin2psi;
  }
}

void qp_tod2map_sig_nopol_single(qp_memory_t *mem, quat_t q_off,
                                 double *ctime, quat_t *q_bore, double *tod,
                                 int n, double *smap, int nside) {
  double dum;
  long ipix;
  quat_t q;

  for (int ii=0; ii<n; ii++) {
    qp_bore2det(mem, q_off, ctime[ii], q_bore[ii], q);
    qp_quat2pix(mem, q, nside, &ipix, &dum, &dum);
    if (tod[ii] != tod[ii]) continue;
    smap[ipix] += tod[ii];
  }
}

/* Compute signal and pointing matrix maps for given boresight timestream,
   signal timestream, and detector offset.
   smap is a npix-x-3 array containing (d, d*cos(2 psi), d*sin(2 psi)). */
void qp_tod2map_sigpnt_single(qp_memory_t *mem, quat_t q_off,
                              double *ctime, quat_t *q_bore, double *tod, int n,
                              vec3_t *smap, vec6_t *pmap, int nside) {
  double sin2psi, cos2psi;
  long ipix;
  quat_t q;

  for (int ii=0; ii<n; ii++) {
    qp_bore2det(mem, q_off, ctime[ii], q_bore[ii], q);
    qp_quat2pix(mem, q, nside, &ipix, &sin2psi, &cos2psi);
    if (tod[ii] != tod[ii]) continue;
    smap[ipix][0] += tod[ii];
    smap[ipix][1] += tod[ii] * cos2psi;
    smap[ipix][2] += tod[ii] * sin2psi;
    pmap[ipix][0] += 1;
    pmap[ipix][1] += cos2psi;
    pmap[ipix][2] += sin2psi;
    pmap[ipix][3] += cos2psi*cos2psi;
    pmap[ipix][4] += cos2psi*sin2psi;
    pmap[ipix][5] += sin2psi*sin2psi;
  }
}

void qp_tod2map_sigpnt_nopol_single(qp_memory_t *mem, quat_t q_off,
                                    double *ctime, quat_t *q_bore, double *tod,
                                    int n, double *smap, double *pmap, int nside) {
  double dum;
  long ipix;
  quat_t q;

  for (int ii=0; ii<n; ii++) {
    qp_bore2det(mem, q_off, ctime[ii], q_bore[ii], q);
    qp_quat2pix(mem, q, nside, &ipix, &dum, &dum);
    if (tod[ii] != tod[ii]) continue;
    pmap[ipix] += 1;
    smap[ipix] += tod[ii];
  }
}

/* Compute pointing matrix map for given boresight timestream and detector
   offset. pmap is a 6-x-npix array containing (hits, p01, p02, p11, p12, p22) */
void qp_tod2map_pnt_single_hwp(qp_memory_t *mem, quat_t q_off,
                               double *ctime, quat_t *q_bore, quat_t *q_hwp, int n,
                               vec6_t *pmap, int nside) {
  double sin2psi, cos2psi;
  long ipix;
  quat_t q;

  for (int ii=0; ii<n; ii++) {
    qp_bore2det_hwp(mem, q_off, ctime[ii], q_bore[ii], q_hwp[ii], q);
    qp_quat2pix(mem, q, nside, &ipix, &sin2psi, &cos2psi);
    pmap[ipix][0] += 1;
    pmap[ipix][1] += cos2psi;
    pmap[ipix][2] += sin2psi;
    pmap[ipix][3] += cos2psi*cos2psi;
    pmap[ipix][4] += cos2psi*sin2psi;
    pmap[ipix][5] += sin2psi*sin2psi;
  }
}

/* Compute signal map for given boresight timestream, hwp timestream,
   signal timestream, and detector offset.
   smap is a npix-x-3 array containing (d, d*cos(2 psi), d*sin(2 psi)). */
void qp_tod2map_sig_single_hwp(qp_memory_t *mem, quat_t q_off,
                               double *ctime, quat_t *q_bore, quat_t *q_hwp,
                               double *tod, int n, vec3_t *smap, int nside) {
  double sin2psi, cos2psi;
  long ipix;
  quat_t q;

  for (int ii=0; ii<n; ii++) {
    qp_bore2det_hwp(mem, q_off, ctime[ii], q_bore[ii], q_hwp[ii], q);
    qp_quat2pix(mem, q, nside, &ipix, &sin2psi, &cos2psi);
    if (tod[ii] != tod[ii]) continue;
    smap[ipix][0] += tod[ii];
    smap[ipix][1] += tod[ii] * cos2psi;
    smap[ipix][2] += tod[ii] * sin2psi;
  }
}

/* Compute signal and pointing matrix maps for given boresight timestream,
   hwp timestream, signal timestream, and detector offset.
   smap is a npix-x-3 array containing (d, d*cos(2 psi), d*sin(2 psi)). */
void qp_tod2map_sigpnt_single_hwp(qp_memory_t *mem, quat_t q_off,
                                  double *ctime, quat_t *q_bore, quat_t *q_hwp,
                                  double *tod, int n, vec3_t *smap,
                                  vec6_t *pmap, int nside) {
  double sin2psi, cos2psi;
  long ipix;
  quat_t q;

  for (int ii=0; ii<n; ii++) {
    qp_bore2det_hwp(mem, q_off, ctime[ii], q_bore[ii], q_hwp[ii], q);
    qp_quat2pix(mem, q, nside, &ipix, &sin2psi, &cos2psi);
    if (tod[ii] != tod[ii]) continue;
    smap[ipix][0] += tod[ii];
    smap[ipix][1] += tod[ii] * cos2psi;
    smap[ipix][2] += tod[ii] * sin2psi;
    pmap[ipix][0] += 1;
    pmap[ipix][1] += cos2psi;
    pmap[ipix][2] += sin2psi;
    pmap[ipix][3] += cos2psi*cos2psi;
    pmap[ipix][4] += cos2psi*sin2psi;
    pmap[ipix][5] += sin2psi*sin2psi;
  }
}

/* Compute pointing matrix map for given boresight timestream and detector
   pair. pmap is a 6-x-npix array containing (hits, p01, p02, p11, p12, p22) */
void qp_tod2map_pnt_pair(qp_memory_t *mem, quat_t q_off,
                         double *ctime, quat_t *q_bore, int n,
                         vec6_t *pmap, int nside) {
  double sin2psi, cos2psi;
  long ipix;
  quat_t q;

  for (int ii=0; ii<n; ii++) {
    qp_bore2det(mem, q_off, ctime[ii], q_bore[ii], q);
    qp_quat2pix(mem, q, nside, &ipix, &sin2psi, &cos2psi);
    pmap[ipix][0] += 2;
    pmap[ipix][3] += 2*cos2psi*cos2psi;
    pmap[ipix][4] += 2*cos2psi*sin2psi;
    pmap[ipix][5] += 2*sin2psi*sin2psi;
  }
}

/* Compute pointing matrix map for given boresight timestream and detector
   pair. pmap is a 6-x-npix array containing (hits, p01, p02, p11, p12, p22) */
void qp_tod2map_pnt_pair_hwp(qp_memory_t *mem, quat_t q_off,
                             double *ctime, quat_t *q_bore, quat_t *q_hwp, int n,
                             vec6_t *pmap, int nside) {
  double sin2psi, cos2psi;
  long ipix;
  quat_t q;

  for (int ii=0; ii<n; ii++) {
    qp_bore2det_hwp(mem, q_off, ctime[ii], q_bore[ii], q_hwp[ii], q);
    qp_quat2pix(mem, q, nside, &ipix, &sin2psi, &cos2psi);
    pmap[ipix][0] += 2;
    pmap[ipix][3] += 2*cos2psi*cos2psi;
    pmap[ipix][4] += 2*cos2psi*sin2psi;
    pmap[ipix][5] += 2*sin2psi*sin2psi;
  }
}

/* boilerplate */

void qp_reduce_smap(vec3_t *smap, vec3_t *smaploc, int npix) {
#pragma omp critical
  for (int ipix=0; ipix<npix; ipix++)
    for (int ii=0; ii<3; ii++)
      if (smaploc[ipix][ii] != 0)
        smap[ipix][ii] += smaploc[ipix][ii];
}

void qp_reduce_pmap(vec6_t *pmap, vec6_t *pmaploc, int npix) {
#pragma omp critical
  for (int ipix=0; ipix<npix; ipix++)
    for (int ii=0; ii<6; ii++)
      if (pmaploc[ipix][ii] != 0)
        pmap[ipix][ii] += pmaploc[ipix][ii];
}

void qp_reduce_map_nopol(double *map, double *maploc, int npix) {
#pragma omp critical
  for (int ipix=0; ipix<npix; ipix++)
    if (maploc[ipix] != 0)
      map[ipix] += maploc[ipix];
}

void qp_debug_mp(qp_memory_t *mem) {
#ifdef DEBUG
  qp_print_memory(mem);
#endif
}

#define INIT_MAP(type, x)                   \
  type * x##maploc;                         \
  if (nthreads > 1) {                       \
    x##maploc = calloc(npix, sizeof(type)); \
  } else {                                  \
    x##maploc = x##map;                     \
  }

#define INIT_MAP1(type, x)                             \
  qp_memory_t *memloc = qp_copy_memory(mem);           \
  const int nthreads = qp_get_opt_num_threads(memloc); \
  qp_debug_mp(memloc);                                 \
  INIT_MAP(type, x)

#define INIT_MAP2(type1, x1, type2, x2) \
  INIT_MAP1(type1, x1);                 \
  INIT_MAP(type2, x2)

#define INIT_PARALLELP INIT_MAP1(vec6_t, p)
#define INIT_PARALLELS INIT_MAP1(vec3_t, s)
#define INIT_PARALLELSP INIT_MAP2(vec3_t, s, vec6_t, p)
#define INIT_PARALLELP_NOPOL INIT_MAP1(double, p)
#define INIT_PARALLELS_NOPOL INIT_MAP1(double, s)
#define INIT_PARALLELSP_NOPOL INIT_MAP2(double, s, double, p)

#define END_MAP(x)                               \
  if (nthreads > 1) {                            \
    qp_reduce_##x##map(x##map, x##maploc, npix); \
    free(x##maploc);                             \
  }

#define END_MAP_NOPOL(x)                          \
  if (nthreads > 1) {                             \
    qp_reduce_map_nopol(x##map, x##maploc, npix); \
    free(x##maploc);                              \
  }

#define END_PARALLEL1(x) \
  END_MAP(x);            \
  qp_free_memory(memloc)

#define END_PARALLEL2(x1, x2) \
  END_MAP(x1);                \
  END_PARALLEL1(x2)

#define END_PARALLEL1N(x) \
  END_MAP_NOPOL(x); \
  qp_free_memory(memloc)

#define END_PARALLEL2N(x1, x2) \
  END_MAP_NOPOL(x1);           \
  END_PARALLEL1N(x2)

#define END_PARALLELP END_PARALLEL1(p)
#define END_PARALLELS END_PARALLEL1(s)
#define END_PARALLELSP END_PARALLEL2(s, p)
#define END_PARALLELP_NOPOL END_PARALLEL1N(p)
#define END_PARALLELS_NOPOL END_PARALLEL1N(s)
#define END_PARALLELSP_NOPOL END_PARALLEL2N(s, p)

#define INIT_MASTER                                                    \
  long npix = nside2npix(nside);                                       \
  int num_threads = ndet < mem->num_threads ? ndet : mem->num_threads; \
  omp_set_num_threads(num_threads)

#define INIT_MASTERP INIT_MASTER
#define INIT_MASTERS INIT_MASTER
#define INIT_MASTERSP INIT_MASTER
#define INIT_MASTERP_NOPOL INIT_MASTER
#define INIT_MASTERS_NOPOL INIT_MASTER
#define INIT_MASTERSP_NOPOL INIT_MASTER

#define END_MASTERP
#define END_MASTERS
#define END_MASTERSP
#define END_MASTERP_NOPOL
#define END_MASTERS_NOPOL
#define END_MASTERSP_NOPOL

/* Compute pointing matrix map for given boresight timestream and many detector
   offsets. pmap a 6-x-npix array containing (hits, p01, p02, p11, p12, p22).
   openMP-parallelized. */
void qp_tod2map_pnt(qp_memory_t *mem, quat_t *q_off, int ndet,
                    double *ctime, quat_t *q_bore, int n,
                    vec6_t *pmap, int nside) {

  INIT_MASTERP;

#pragma omp parallel
  {
    INIT_PARALLELP;

#pragma omp for
    for (int idet=0; idet<ndet; idet++) {
#ifdef DEBUG
      printf("thread %d, det %d\n", qp_get_opt_thread_num(memloc), idet);
      printf("offset %f %f %f %f\n", q_off[idet][0], q_off[idet][1],
             q_off[idet][2], q_off[idet][3]);
#endif
      if (memloc->pair_dets)
        qp_tod2map_pnt_pair(memloc, q_off[idet], ctime, q_bore, n, pmaploc, nside);
      else
        qp_tod2map_pnt_single(memloc, q_off[idet], ctime, q_bore, n, pmaploc, nside);
    }

    END_PARALLELP;
  }

  END_MASTERP;
}

/* Compute signal map for given boresight timestream, many signal timestreams,
   and many detector offsets.
   smap is a npix-x-3 array containing (d, d*cos(2 psi), d*sin(2 psi)).
   openMP-parallelized. */
void qp_tod2map_sig(qp_memory_t *mem, quat_t *q_off, int ndet,
                    double *ctime, quat_t *q_bore, double **tod, int n,
                    vec3_t *smap, int nside) {

  INIT_MASTERS;

#pragma omp parallel
  {
    INIT_PARALLELS;

#pragma omp for
    for (int idet=0; idet<ndet; idet++) {
      qp_tod2map_sig_single(memloc, q_off[idet], ctime, q_bore, tod[idet], n,
                            smaploc, nside);
    }

    END_PARALLELS;
  }

  END_MASTERS;
}

/* Compute signal and pointing matrix maps for given boresight timestream,
   many signal timestreams, and many detector offsets.
   smap is a npix-x-3 array containing (d, d*cos(2 psi), d*sin(2 psi)).
   pmap is a npix-x-6 array containing (hits, p01, p02, p11, p12, p22).
   openMP-parallelized. */
void qp_tod2map_sigpnt(qp_memory_t *mem, quat_t *q_off, int ndet,
                       double *ctime, quat_t *q_bore, double **tod, int n,
                       vec3_t *smap, vec6_t *pmap, int nside) {

  INIT_MASTERSP;

#pragma omp parallel
  {
    INIT_PARALLELSP;

#pragma omp for
    for (int idet=0; idet<ndet; idet++) {
      qp_tod2map_sigpnt_single(memloc, q_off[idet], ctime, q_bore, tod[idet], n,
                               smaploc, pmaploc, nside);
    }

    END_PARALLELSP;
  }

  END_MASTERSP;
}

void qp_tod2map_pnt_nopol(qp_memory_t *mem, quat_t *q_off, int ndet,
                          double *ctime, quat_t *q_bore, int n,
                          double *pmap, int nside) {

  INIT_MASTERP_NOPOL;

#pragma omp parallel
  {
    INIT_PARALLELP_NOPOL;

#pragma omp for
    for (int idet=0; idet<ndet; idet++) {
      qp_tod2map_pnt_nopol_single(memloc, q_off[idet], ctime, q_bore, n,
                                  pmaploc, nside);
    }

    END_PARALLELP_NOPOL;
  }

  END_MASTERP_NOPOL;
}

void qp_tod2map_sig_nopol(qp_memory_t *mem, quat_t *q_off, int ndet,
                          double *ctime, quat_t *q_bore, double **tod, int n,
                          double *smap, int nside) {

  INIT_MASTERS_NOPOL;

#pragma omp parallel
  {
    INIT_PARALLELS_NOPOL;

#pragma omp for
    for (int idet=0; idet<ndet; idet++) {
      qp_tod2map_sig_nopol_single(memloc, q_off[idet], ctime, q_bore, tod[idet], n,
                                  smaploc, nside);
    }

    END_PARALLELS_NOPOL;
  }

  END_MASTERS_NOPOL;
}

void qp_tod2map_sigpnt_nopol(qp_memory_t *mem, quat_t *q_off, int ndet,
                             double *ctime, quat_t *q_bore, double **tod, int n,
                             double *smap, double *pmap, int nside) {

  INIT_MASTERSP_NOPOL;

#pragma omp parallel
  {
    INIT_PARALLELSP_NOPOL;

#pragma omp for schedule(static, 1)
    for (int idet=0; idet<ndet; idet++) {
      qp_tod2map_sigpnt_nopol_single(memloc, q_off[idet], ctime, q_bore, tod[idet], n,
                                     smaploc, pmaploc, nside);
    }

    END_PARALLELSP_NOPOL;
  }

  END_MASTERSP_NOPOL;
}

/* Compute pointing matrix map for given boresight timestream and many detector
   offsets. pmap a 6-x-npix array containing (hits, p01, p02, p11, p12, p22).
   openMP-parallelized. */
void qp_tod2map_pnt_hwp(qp_memory_t *mem, quat_t *q_off, int ndet,
                        double *ctime, quat_t *q_bore, quat_t *q_hwp, int n,
                        vec6_t *pmap, int nside) {

  INIT_MASTERP;

#pragma omp parallel
  {
    INIT_PARALLELP;

#pragma omp for
    for (int idet=0; idet<ndet; idet++) {
#ifdef DEBUG
      printf("thread %d, det %d\n", omp_get_thread_num(), idet);
      printf("offset %f %f %f %f\n", q_off[idet][0], q_off[idet][1],
	     q_off[idet][2], q_off[idet][3]);
#endif
      if (memloc->pair_dets)
	qp_tod2map_pnt_pair_hwp(memloc, q_off[idet], ctime, q_bore, q_hwp, n,
                                pmaploc, nside);
      else
	qp_tod2map_pnt_single_hwp(memloc, q_off[idet], ctime, q_bore, q_hwp, n,
                                  pmaploc, nside);
    }

    END_PARALLELP;
  }

  END_MASTERP;
}

/* Compute signal map for given boresight timestream, hwp timestream,
   many signal timestreams, and many detector offsets.
   smap is a npix-x-3 array containing (d, d*cos(2 psi), d*sin(2 psi)).
   openMP-parallelized. */
void qp_tod2map_sig_hwp(qp_memory_t *mem, quat_t *q_off, int ndet,
                        double *ctime, quat_t *q_bore, quat_t *q_hwp,
                        double **tod, int n, vec3_t *smap, int nside) {

  INIT_MASTERS;

#pragma omp parallel
  {
    INIT_PARALLELS;

#pragma omp for
    for (int idet=0; idet<ndet; idet++) {
      qp_tod2map_sig_single_hwp(memloc, q_off[idet], ctime, q_bore, q_hwp,
                                tod[idet], n, smaploc, nside);
    }

    END_PARALLELS;
  }

  END_MASTERS;
}

/* Compute signal and pointing matrix maps for given boresight timestream,
   hwp timestrea, many signal timestreams, and many detector offsets.
   smap is a npix-x-3 array containing (d, d*cos(2 psi), d*sin(2 psi)).
   pmap is a npix-x-6 array containing (hits, p01, p02, p11, p12, p22).
   openMP-parallelized. */
void qp_tod2map_sigpnt_hwp(qp_memory_t *mem, quat_t *q_off, int ndet,
                           double *ctime, quat_t *q_bore, quat_t *q_hwp,
                           double **tod, int n, vec3_t *smap, vec6_t *pmap,
                           int nside) {

  INIT_MASTERSP;

#pragma omp parallel
  {
    INIT_PARALLELSP;

#pragma omp for
    for (int idet=0; idet<ndet; idet++) {
      qp_tod2map_sigpnt_single_hwp(memloc, q_off[idet], ctime, q_bore, q_hwp,
                                   tod[idet], n, smaploc, pmaploc, nside);
    }

    END_PARALLELSP;
  }

  END_MASTERSP;
}

double qp_m2d(double *map, double cpp, double spp) {
  return map[0] + map[1] * cpp + map[2] * spp;
}

/* Compute signal timestream given an input map, boresight pointing
   and detector offset. smap is a npix-x-3 array containing (T,Q,U) maps. */
void qp_map2tod_single(qp_memory_t *mem, quat_t q_off,
                       double *ctime, quat_t *q_bore, vec3_t *smap,
                       int nside, double *tod, int n) {

  double spp, cpp;
  long ipix;
  quat_t q;

  for (int ii=0; ii<n; ii++) {
    qp_bore2det(mem, q_off, ctime[ii], q_bore[ii], q);
    qp_quat2pix(mem, q, nside, &ipix, &spp, &cpp);
    tod[ii] = qp_m2d(smap[ipix], cpp, spp);
  }
}

/* Compute signal timestream given an input map, boresight pointing
   and detector offset. smap is a npix-x-3 array containing (T,Q,U) maps. */
void qp_map2tod_single_hwp(qp_memory_t *mem, quat_t q_off,
                           double *ctime, quat_t *q_bore, quat_t *q_hwp,
                           vec3_t *smap, int nside, double *tod, int n) {

  double spp, cpp;
  long ipix;
  quat_t q;

  for (int ii=0; ii<n; ii++) {
    qp_bore2det_hwp(mem, q_off, ctime[ii], q_bore[ii], q_hwp[ii], q);
    qp_quat2pix(mem, q, nside, &ipix, &spp, &cpp);
    tod[ii] = qp_m2d(smap[ipix], cpp, spp);
  }
}

/* Compute signal timestream given an input map, boresight pointing
   and detector offset. smap is a npix-x-3 array containing (T,Q,U) maps. */
void qp_map2tod_nopol_single(qp_memory_t *mem, quat_t q_off,
                             double *ctime, quat_t *q_bore, double *smap,
                             int nside, double *tod, int n) {

  double dum;
  long ipix;
  quat_t q;

  for (int ii=0; ii<n; ii++) {
    qp_bore2det(mem, q_off, ctime[ii], q_bore[ii], q);
    qp_quat2pix(mem, q, nside, &ipix, &dum, &dum);
    tod[ii] = smap[ipix];
  }
}

void qp_pixel_offset(qp_memory_t *mem, int nside, long pix, double ra, double dec,
                     double *dtheta, double *dphi) {
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

/* Compute signal timestream given an input map, boresight pointing
   and detector offset. smap is a npix-x-3 array containing (T,Q,U) maps. */
void qp_map2tod_der1_single(qp_memory_t *mem, quat_t q_off,
                            double *ctime, quat_t *q_bore, vec9_t *smap,
                            int nside, double *tod, int n) {

  double ra, dec, spp, cpp, dtheta, dphi;
  long ipix;
  quat_t q;
  double *m;

  for (int ii=0; ii<n; ii++) {
    qp_bore2det(mem, q_off, ctime[ii], q_bore[ii], q);
    qp_quat2radec(mem, q, &ra, &dec, &spp, &cpp);
    ipix = qp_radec2pix(mem, ra, dec, nside);
    qp_pixel_offset(mem, nside, ipix, ra, dec, &dtheta, &dphi);
    m = smap[ipix];
    tod[ii] = qp_m2d(m, cpp, spp)
      + dtheta * qp_m2d(m+3, cpp, spp)
      + dphi   * qp_m2d(m+6, cpp, spp);
  }
}

/* Compute signal timestream given an input map, boresight pointing
   and detector offset. smap is a npix-x-3 array containing (T,Q,U) maps. */
void qp_map2tod_der1_single_hwp(qp_memory_t *mem, quat_t q_off,
                                double *ctime, quat_t *q_bore, quat_t *q_hwp,
                                vec9_t *smap, int nside, double *tod,
                                int n) {

  double ra, dec, spp, cpp, dtheta, dphi;
  long ipix;
  quat_t q;
  double *m;

  for (int ii=0; ii<n; ii++) {
    qp_bore2det_hwp(mem, q_off, ctime[ii], q_bore[ii], q_hwp[ii], q);
    qp_quat2radec(mem, q, &ra, &dec, &spp, &cpp);
    ipix = qp_radec2pix(mem, ra, dec, nside);
    qp_pixel_offset(mem, nside, ipix, ra, dec, &dtheta, &dphi);
    m = smap[ipix];
    tod[ii] = qp_m2d(m, cpp, spp)
      + dtheta * qp_m2d(m+3, cpp, spp)
      + dphi   * qp_m2d(m+6, cpp, spp);
  }
}

/* Compute signal timestream given an input map, boresight pointing
   and detector offset. smap is a npix-x-3 array containing (T,Q,U) maps. */
void qp_map2tod_der1_nopol_single(qp_memory_t *mem, quat_t q_off,
                                  double *ctime, quat_t *q_bore, vec3_t *smap,
                                  int nside, double *tod, int n) {

  double ra, dec, dtheta, dphi, dum;
  long ipix;
  quat_t q;
  double *m;

  for (int ii=0; ii<n; ii++) {
    qp_bore2det(mem, q_off, ctime[ii], q_bore[ii], q);
    qp_quat2radec(mem, q, &ra, &dec, &dum, &dum);
    ipix = qp_radec2pix(mem, ra, dec, nside);
    qp_pixel_offset(mem, nside, ipix, ra, dec, &dtheta, &dphi);
    m = smap[ipix];
    tod[ii] = m[0]
      + dtheta * m[1]
      + dphi   * m[2];
  }
}

/* Compute signal timestream given an input map, boresight pointing
   and detector offset. smap is a npix-x-3 array containing (T,Q,U) maps. */
void qp_map2tod_der2_single(qp_memory_t *mem, quat_t q_off,
                            double *ctime, quat_t *q_bore, vec18_t *smap,
                            int nside, double *tod, int n) {

  double ra, dec, spp, cpp, dtheta, dphi;
  long ipix;
  quat_t q;
  double *m;

  for (int ii=0; ii<n; ii++) {
    qp_bore2det(mem, q_off, ctime[ii], q_bore[ii], q);
    qp_quat2radec(mem, q, &ra, &dec, &spp, &cpp);
    ipix = qp_radec2pix(mem, ra, dec, nside);
    qp_pixel_offset(mem, nside, ipix, ra, dec, &dtheta, &dphi);
    m = smap[ipix];
    tod[ii] = qp_m2d(m, cpp, spp)
      + dtheta *          qp_m2d(m+3,  cpp, spp)
      + dphi   *          qp_m2d(m+6,  cpp, spp)
      + dtheta * dtheta * qp_m2d(m+9,  cpp, spp) / 2.
      + dtheta * dphi   * qp_m2d(m+12, cpp, spp)
      + dphi   * dphi   * qp_m2d(m+15, cpp, spp) / 2.;
  }
}

/* Compute signal timestream given an input map, boresight pointing
   and detector offset. smap is a npix-x-3 array containing (T,Q,U) maps. */
void qp_map2tod_der2_single_hwp(qp_memory_t *mem, quat_t q_off,
                                double *ctime, quat_t *q_bore, quat_t *q_hwp,
                                vec18_t *smap, int nside, double *tod, int n) {

  double ra, dec, spp, cpp, dtheta, dphi;
  long ipix;
  quat_t q;
  double *m;

  for (int ii=0; ii<n; ii++) {
    qp_bore2det_hwp(mem, q_off, ctime[ii], q_bore[ii], q_hwp[ii], q);
    qp_quat2radec(mem, q, &ra, &dec, &spp, &cpp);
    ipix = qp_radec2pix(mem, ra, dec, nside);
    qp_pixel_offset(mem, nside, ipix, ra, dec, &dtheta, &dphi);
    m = smap[ipix];
    tod[ii] = qp_m2d(m, cpp, spp)
      + dtheta *          qp_m2d(m+3,  cpp, spp)
      + dphi   *          qp_m2d(m+6,  cpp, spp)
      + dtheta * dtheta * qp_m2d(m+9,  cpp, spp) / 2.
      + dtheta * dphi   * qp_m2d(m+12, cpp, spp)
      + dphi   * dphi   * qp_m2d(m+15, cpp, spp) / 2.;
  }
}

/* Compute signal timestream given an input map, boresight pointing
   and detector offset. smap is a npix array containing a T map. */
void qp_map2tod_der2_nopol_single(qp_memory_t *mem, quat_t q_off,
                                  double *ctime, quat_t *q_bore, vec6_t *smap,
                                  int nside, double *tod, int n) {

  double ra, dec, dtheta, dphi, dum;
  long ipix;
  quat_t q;
  double *m;

  for (int ii=0; ii<n; ii++) {
    qp_bore2det(mem, q_off, ctime[ii], q_bore[ii], q);
    qp_quat2radec(mem, q, &ra, &dec, &dum, &dum);
    ipix = qp_radec2pix(mem, ra, dec, nside);
    qp_pixel_offset(mem, nside, ipix, ra, dec, &dtheta, &dphi);
    m = smap[ipix];
    tod[ii] = m[0]
      + dtheta *          m[1]
      + dphi   *          m[2]
      + dtheta * dtheta * m[3] / 2.
      + dtheta * dphi   * m[4]
      + dphi   * dphi   * m[5] / 2.;
  }
}

/* Compute signal timestream given an input map, boresight pointing
   and detector offsets. smap is a npix-x-3 array containing (T,Q,U) maps. */
void qp_map2tod(qp_memory_t *mem, quat_t *q_off, int ndet,
                double *ctime, quat_t *q_bore, vec3_t *smap, int nside,
                double **tod, int n) {
#pragma omp parallel
  {
    // local copy of memory
    qp_memory_t *memloc = qp_copy_memory(mem);

#pragma omp for nowait
    for (int idet=0; idet<ndet; idet++) {
      qp_map2tod_single(memloc, q_off[idet], ctime, q_bore, smap, nside,
                        tod[idet], n);
    }
    qp_free_memory(memloc);
  }
}

/* Compute signal timestream given an input map, boresight pointing,
   hwp timestream, and detector offsets.
   smap is a npix-x-3 array containing (T,Q,U) maps. */
void qp_map2tod_hwp(qp_memory_t *mem, quat_t *q_off, int ndet,
                    double *ctime, quat_t *q_bore, quat_t *q_hwp,
                    vec3_t *smap, int nside, double **tod, int n) {
#pragma omp parallel
  {
    // local copy of memory
    qp_memory_t *memloc = qp_copy_memory(mem);

#pragma omp for nowait
    for (int idet=0; idet<ndet; idet++) {
      qp_map2tod_single_hwp(memloc, q_off[idet], ctime, q_bore, q_hwp,
                            smap, nside, tod[idet], n);
    }
    qp_free_memory(memloc);
  }
}

/* Compute signal timestream given an input map, boresight pointing
   and detector offsets. smap is a npix-x-3 array containing (T,Q,U) maps. */
void qp_map2tod_nopol(qp_memory_t *mem, quat_t *q_off, int ndet,
                      double *ctime, quat_t *q_bore, double *smap, int nside,
                      double **tod, int n) {
#pragma omp parallel
  {
    // local copy of memory
    qp_memory_t *memloc = qp_copy_memory(mem);

#pragma omp for nowait
    for (int idet=0; idet<ndet; idet++) {
      qp_map2tod_nopol_single(memloc, q_off[idet], ctime, q_bore, smap, nside,
                              tod[idet], n);
    }
    qp_free_memory(memloc);
  }
}

/* Compute signal timestream given an input map, boresight pointing
   and detector offsets. smap is a npix-x-3 array containing (T,Q,U) maps. */
void qp_map2tod_der1(qp_memory_t *mem, quat_t *q_off, int ndet,
                     double *ctime, quat_t *q_bore, vec9_t *smap,
                     int nside, double **tod, int n) {
#pragma omp parallel
  {
    // local copy of memory
    qp_memory_t *memloc = qp_copy_memory(mem);

#pragma omp for nowait
    for (int idet=0; idet<ndet; idet++) {
      qp_map2tod_der1_single(memloc, q_off[idet], ctime, q_bore, smap,
                             nside, tod[idet], n);
    }
    qp_free_memory(memloc);
  }
}

/* Compute signal timestream given an input map, boresight pointing,
   hwp timestream, and detector offsets.
   smap is a npix-x-3 array containing (T,Q,U) maps. */
void qp_map2tod_der1_hwp(qp_memory_t *mem, quat_t *q_off, int ndet,
                         double *ctime, quat_t *q_bore, quat_t *q_hwp,
                         vec9_t *smap, int nside, double **tod, int n) {
#pragma omp parallel
  {
    // local copy of memory
    qp_memory_t *memloc = qp_copy_memory(mem);

#pragma omp for nowait
    for (int idet=0; idet<ndet; idet++) {
      qp_map2tod_der1_single_hwp(memloc, q_off[idet], ctime, q_bore, q_hwp,
                                 smap, nside, tod[idet], n);
    }
    qp_free_memory(memloc);
  }
}

/* Compute signal timestream given an input map, boresight pointing
   and detector offsets. smap is a npix-x-3 array containing (T,Q,U) maps. */
void qp_map2tod_der1_nopol(qp_memory_t *mem, quat_t *q_off, int ndet,
                           double *ctime, quat_t *q_bore, vec3_t *smap,
                           int nside, double **tod, int n) {
#pragma omp parallel
  {
    // local copy of memory
    qp_memory_t *memloc = qp_copy_memory(mem);

#pragma omp for nowait
    for (int idet=0; idet<ndet; idet++) {
      qp_map2tod_der1_nopol_single(memloc, q_off[idet], ctime, q_bore, smap,
                                   nside, tod[idet], n);
    }
    qp_free_memory(memloc);
  }
}

/* Compute signal timestream given an input map, boresight pointing
   and detector offsets. smap is a npix-x-3 array containing (T,Q,U) maps. */
void qp_map2tod_der2(qp_memory_t *mem, quat_t *q_off, int ndet,
                     double *ctime, quat_t *q_bore, vec18_t *smap,
                     int nside, double **tod, int n) {
#pragma omp parallel
  {
    // local copy of memory
    qp_memory_t *memloc = qp_copy_memory(mem);

#pragma omp for nowait
    for (int idet=0; idet<ndet; idet++) {
      qp_map2tod_der2_single(memloc, q_off[idet], ctime, q_bore, smap,
                             nside, tod[idet], n);
    }
    qp_free_memory(memloc);
  }
}

/* Compute signal timestream given an input map, boresight pointing,
   hwp timestream, and detector offsets.
   smap is a npix-x-3 array containing (T,Q,U) maps. */
void qp_map2tod_der2_hwp(qp_memory_t *mem, quat_t *q_off, int ndet,
                         double *ctime, quat_t *q_bore, quat_t *q_hwp,
                         vec18_t *smap, int nside, double **tod, int n) {
#pragma omp parallel
  {
    // local copy of memory
    qp_memory_t *memloc = qp_copy_memory(mem);

#pragma omp for nowait
    for (int idet=0; idet<ndet; idet++) {
      qp_map2tod_der2_single_hwp(memloc, q_off[idet], ctime, q_bore, q_hwp,
                                 smap, nside, tod[idet], n);
    }
    qp_free_memory(memloc);
  }
}

/* Compute signal timestream given an input map, boresight pointing
   and detector offsets. smap is a npix-x-3 array containing (T,Q,U) maps. */
void qp_map2tod_der2_nopol(qp_memory_t *mem, quat_t *q_off, int ndet,
                           double *ctime, quat_t *q_bore, vec6_t *smap,
                           int nside, double **tod, int n) {
#pragma omp parallel
  {
    // local copy of memory
    qp_memory_t *memloc = qp_copy_memory(mem);

#pragma omp for nowait
    for (int idet=0; idet<ndet; idet++) {
      qp_map2tod_der2_nopol_single(memloc, q_off[idet], ctime, q_bore, smap,
                                   nside, tod[idet], n);
    }
    qp_free_memory(memloc);
  }
}

void qp_set_opt_num_threads(qp_memory_t *mem, int num_threads) {
  if (num_threads == 0) {
#ifdef _OPENMP
    num_threads = omp_get_num_procs();
#else
    num_threads = 1;
#endif
  }
  mem->num_threads = num_threads;
#ifdef _OPENMP
  omp_set_num_threads(num_threads);
#endif
}

int qp_get_opt_num_threads(qp_memory_t *mem) {
#ifdef _OPENMP
  if (omp_in_parallel())
    mem->num_threads = omp_get_num_threads();
#endif
  return mem->num_threads;
}

void qp_set_opt_thread_num(qp_memory_t *mem, int thread) {
#ifdef _OPENMP
  mem->thread_num = omp_get_thread_num();
#else
  mem->thread_num = thread;
#endif
}

int qp_get_opt_thread_num(qp_memory_t *mem) {
#ifdef _OPENMP
  qp_set_opt_thread_num(mem, 0);
#endif
  return mem->thread_num;
}
