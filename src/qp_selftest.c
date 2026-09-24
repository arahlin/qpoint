/* Exercise the C constructors and copy helpers from C.
 *
 * These are public API that neither Python package calls. The ctypes
 * layer builds qp_det_t, qp_detarr_t, qp_point_t and qp_map_t itself,
 * field by field through the struct mirrors, so nothing in the test
 * suite reaches the constructors and their coverage read zero -- not
 * because they are wrong, but because the suite cannot get to them.
 *
 * An external C consumer would use exactly these, so rather than
 * exclude them from the coverage report or delete them, this exercises
 * them the way such a consumer would, and the suite reaches all of it
 * through one binding.
 *
 * Each check is numbered in sequence. Returns 0 if they all pass, or
 * the number of the first that fails, with a description in msg.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "qpoint.h"

#define CHECK(cond, what)                                               \
  do {                                                                  \
    step++;                                                             \
    if (!(cond)) {                                                      \
      snprintf(msg, msglen, "check %d failed: %s", step, what);         \
      goto done;                                                        \
    }                                                                   \
  } while (0)

int qp_selftest(char *msg, size_t msglen) {
  const size_t n = 8, ndet = 3, nside = 8, npix = 12 * nside * nside;
  int step = 0, failed = 0;

  qp_det_t *det = NULL, *ddet = NULL;
  qp_detarr_t *dets = NULL;
  qp_point_t *pnt = NULL, *pnt2 = NULL;
  qp_map_t *map = NULL, *map2 = NULL, *map3 = NULL, *pmap = NULL;
  qp_pixhash_t *hash = NULL, *hash2 = NULL;
  qp_memory_t *mem = NULL, *memcpy_ = NULL;

  quat_t q_off = {1., 0., 0., 0.};
  mueller_t mueller = {1., 1., 0., 1.};
  double tod[8], weights[8], ctime[8];
  uint8_t flag[8];
  quat_t q_bore[8], q_hwp[8];
  double *vec[3], *proj[6];
  double vecbuf[3][12 * 8 * 8], projbuf[6][12 * 8 * 8];
  long pix[8];

  if (msglen) msg[0] = '\0';

  for (size_t i = 0; i < n; i++) {
    tod[i] = 0.5 * (double)i;
    weights[i] = 1. + (double)i;
    flag[i] = (i % 3 == 0);
    ctime[i] = 1418662800. + (double)i;
    for (int j = 0; j < 4; j++) {
      q_bore[i][j] = (j == 0) ? 1. : 0.;
      q_hwp[i][j] = (j == 0) ? 1. : 0.;
    }
    pix[i] = (long)(i * 5);
  }
  for (int r = 0; r < 3; r++) {
    vec[r] = vecbuf[r];
    for (size_t i = 0; i < npix; i++) vecbuf[r][i] = (double)(r + 1);
  }
  for (int r = 0; r < 6; r++) {
    proj[r] = projbuf[r];
    for (size_t i = 0; i < npix; i++) projbuf[r][i] = 1.;
  }

  /* ---- single detector ---- */
  det = qp_init_det(q_off, 2.5, 1.5, mueller);
  CHECK(det != NULL && det->init, "qp_init_det returned an uninitialized det");
  CHECK(det->weight == 2.5 && det->gain == 1.5, "det weight/gain not stored");
  CHECK(det->mueller[0] == 1. && det->mueller[2] == 0., "det mueller not stored");
  CHECK(!det->tod_init && !det->flag_init && !det->weights_init,
        "det arrays should start uninitialized");

  qp_init_det_tod(det, n);
  CHECK(det->tod_init && det->n == n && det->tod != NULL,
        "qp_init_det_tod did not allocate");
  CHECK(det->tod[0] == 0. && det->tod[n - 1] == 0.,
        "qp_init_det_tod should zero the array");

  qp_init_det_flag(det, n);
  CHECK(det->flag_init && det->flag != NULL, "qp_init_det_flag did not allocate");
  qp_init_det_weights(det, n);
  CHECK(det->weights_init && det->weights != NULL,
        "qp_init_det_weights did not allocate");
  qp_free_det(det);
  det = NULL;

  /* ---- and from arrays, copying or aliasing ---- */
  det = qp_init_det(q_off, 1., 1., mueller);
  qp_init_det_tod_from_array(det, tod, n, 0);
  CHECK(det->tod == tod, "copy=0 should alias the caller's tod");
  qp_init_det_flag_from_array(det, flag, n, 0);
  CHECK(det->flag == flag, "copy=0 should alias the caller's flag");
  qp_init_det_weights_from_array(det, weights, n, 0);
  CHECK(det->weights == weights, "copy=0 should alias the caller's weights");
  qp_free_det(det);

  det = qp_init_det(q_off, 1., 1., mueller);
  qp_init_det_tod_from_array(det, tod, n, 1);
  CHECK(det->tod != tod && det->tod[3] == tod[3], "copy=1 should copy the tod");
  qp_init_det_flag_from_array(det, flag, n, 1);
  CHECK(det->flag != flag && det->flag[3] == flag[3], "copy=1 should copy the flag");
  qp_init_det_weights_from_array(det, weights, n, 1);
  CHECK(det->weights != weights && det->weights[3] == weights[3],
        "copy=1 should copy the weights");
  qp_free_det(det);
  det = NULL;

  ddet = qp_default_det();
  CHECK(ddet != NULL && ddet->init, "qp_default_det returned an uninitialized det");
  CHECK(ddet->weight == 1. && ddet->gain == 1., "default det weight/gain");
  qp_free_det(ddet);
  ddet = NULL;

  /* ---- detector array ---- */
  {
    quat_t offs[3];
    double w[3] = {1., 2., 3.}, g[3] = {1., 1., 1.};
    mueller_t mu[3];
    double *todp[3], *wp[3];
    uint8_t *flagp[3];
    for (size_t d = 0; d < ndet; d++) {
      for (int j = 0; j < 4; j++) offs[d][j] = (j == 0) ? 1. : 0.;
      for (int j = 0; j < 4; j++) mu[d][j] = mueller[j];
      todp[d] = tod;
      wp[d] = weights;
      flagp[d] = flag;
    }
    dets = qp_init_detarr(offs, w, g, mu, ndet);
    CHECK(dets != NULL && dets->init && dets->n == ndet,
          "qp_init_detarr returned the wrong size");
    CHECK(dets->arr != NULL && dets->arr[1].weight == 2.,
          "qp_init_detarr did not fill the array");

    qp_init_detarr_tod(dets, n);
    CHECK(dets->arr[0].tod_init && dets->arr[2].tod != NULL,
          "qp_init_detarr_tod did not allocate every det");
    qp_init_detarr_flag(dets, n);
    CHECK(dets->arr[0].flag_init, "qp_init_detarr_flag did not allocate");
    qp_init_detarr_weights(dets, n);
    CHECK(dets->arr[0].weights_init, "qp_init_detarr_weights did not allocate");
    qp_free_detarr(dets);

    dets = qp_init_detarr(offs, w, g, mu, ndet);
    qp_init_detarr_tod_from_array(dets, todp, n, 0);
    CHECK(dets->arr[1].tod == tod, "detarr tod copy=0 should alias");
    qp_init_detarr_flag_from_array(dets, flagp, n, 0);
    CHECK(dets->arr[1].flag == flag, "detarr flag copy=0 should alias");
    qp_init_detarr_weights_from_array(dets, wp, n, 0);
    CHECK(dets->arr[1].weights == weights, "detarr weights copy=0 should alias");
    qp_free_detarr(dets);

    dets = qp_init_detarr(offs, w, g, mu, ndet);
    qp_init_detarr_tod_from_array(dets, todp, n, 1);
    CHECK(dets->arr[1].tod != tod && dets->arr[1].tod[2] == tod[2],
          "detarr tod copy=1 should copy");
    qp_free_detarr(dets);
    dets = NULL;
  }

  /* ---- pointing ---- */
  pnt = qp_init_point(n, 1, 1);
  CHECK(pnt != NULL && pnt->init && pnt->n == n, "qp_init_point wrong size");
  CHECK(pnt->q_bore_init && pnt->ctime_init && pnt->q_hwp_init,
        "qp_init_point(time=1, pol=1) should allocate all three");
  qp_free_point(pnt);

  pnt = qp_init_point(n, 0, 0);
  CHECK(pnt->q_bore_init && !pnt->ctime_init && !pnt->q_hwp_init,
        "qp_init_point(time=0, pol=0) should allocate only q_bore");
  qp_free_point(pnt);
  pnt = NULL;

  pnt2 = qp_init_point_from_arrays(q_bore, ctime, q_hwp, n, 0);
  CHECK(pnt2 != NULL && pnt2->q_bore == q_bore && pnt2->ctime == ctime,
        "qp_init_point_from_arrays copy=0 should alias");
  qp_free_point(pnt2);

  pnt2 = qp_init_point_from_arrays(q_bore, ctime, q_hwp, n, 1);
  CHECK(pnt2->ctime != ctime && pnt2->ctime[4] == ctime[4],
        "qp_init_point_from_arrays copy=1 should copy");
  qp_free_point(pnt2);
  pnt2 = NULL;

  /* ---- maps ---- */
  map = qp_init_map(nside, npix, QP_VEC_POL, QP_PROJ_POL);
  CHECK(map != NULL && map->init && map->npix == npix,
        "qp_init_map returned the wrong npix");
  CHECK(map->num_vec == 3 && map->num_proj == 6,
        "qp_init_map row counts for POL/POL");

  map2 = qp_init_map_from_arrays(vec, proj, nside, npix, QP_VEC_POL,
                                 QP_PROJ_POL, 0);
  CHECK(map2 != NULL && map2->vec[1] == vec[1],
        "qp_init_map_from_arrays copy=0 should alias");
  qp_free_map(map2);

  map2 = qp_init_map_from_arrays(vec, proj, nside, npix, QP_VEC_POL,
                                 QP_PROJ_POL, 1);
  CHECK(map2->vec[1] != vec[1] && map2->vec[1][7] == vec[1][7],
        "qp_init_map_from_arrays copy=1 should copy");

  map3 = qp_init_map_from_map(map2, 1, 0);
  CHECK(map3 != NULL && map3->npix == npix && map3->vec[0][7] == 0.,
        "qp_init_map_from_map(blank=1) should zero the copy");
  qp_free_map(map3);

  map3 = qp_init_map_from_map(map2, 0, 1);
  CHECK(map3->vec[1][7] == map2->vec[1][7],
        "qp_init_map_from_map(copy=1) should carry the values");
  qp_free_map(map3);
  map3 = NULL;

  /* ---- pixel hash ----
     qp_init_map_pixhash insists the pixel list is exactly as long as the
     map, so a partial map is the only place it fits: this one is n
     pixels wide rather than a full sky. */
  pmap = qp_init_map(nside, n, QP_VEC_TEMP, QP_PROJ_TEMP);
  CHECK(pmap != NULL && pmap->npix == n, "partial qp_init_map wrong npix");
  CHECK(qp_init_map_pixhash(pmap, pix, n) == 0, "qp_init_map_pixhash failed");
  CHECK(pmap->pixhash_init, "map should carry a pixhash after init");
  CHECK(qp_repixelize(pmap->pixhash, pix[2]) == 2, "qp_repixelize wrong index");
  CHECK(qp_init_map_pixhash(pmap, pix, n - 1) != 0,
        "a pixel list shorter than the map should be refused");

  hash = qp_init_pixhash(pix, n);
  CHECK(hash != NULL, "qp_init_pixhash returned NULL");
  hash2 = qp_copy_pixhash(hash);
  CHECK(hash2 != NULL && hash2 != hash, "qp_copy_pixhash should allocate a copy");
  CHECK(qp_repixelize(hash2, pix[3]) == 3, "the copied hash lost an entry");
  CHECK(qp_repixelize(hash2, 999999) < 0, "an absent pixel should be negative");
  qp_free_pixhash(hash2);
  hash2 = NULL;
  qp_free_pixhash(hash);
  hash = NULL;

  /* ---- memory copy, which carries the bulletin table ---- */
  mem = qp_init_memory();
  CHECK(mem != NULL && mem->init, "qp_init_memory failed");
  {
    /* the table stores floats, so these are chosen to survive the
       narrowing exactly and keep the comparison below an equality */
    double dut1[3] = {0.125, 0.25, 0.375}, x[3] = {0.0625, 0.125, 0.1875};
    double y[3] = {0.5, 0.25, 0.125};
    CHECK(qp_set_iers_bulletin_a(mem, 57000, 57002, dut1, x, y) == 0,
          "qp_set_iers_bulletin_a failed");
  }
  memcpy_ = qp_copy_memory(mem);
  CHECK(memcpy_ != NULL && memcpy_->init, "qp_copy_memory failed");
  CHECK(memcpy_->bulletinA.entries != mem->bulletinA.entries,
        "qp_copy_memory should copy the bulletin, not alias it");
  CHECK(memcpy_->bulletinA.mjd_min == mem->bulletinA.mjd_min,
        "the copied bulletin lost its date range");
  {
    double d1 = 0, x1 = 0, y1 = 0;
    qp_get_iers_bulletin_a(memcpy_, 57001., &d1, &x1, &y1);
    CHECK(d1 == 0.25, "the copied bulletin lost its dut1");
    CHECK(x1 == 0.125 && y1 == 0.25, "the copied bulletin lost its polar motion");
  }

done:
  failed = (msglen && msg[0]) ? step : 0;
  if (det) qp_free_det(det);
  if (ddet) qp_free_det(ddet);
  if (dets) qp_free_detarr(dets);
  if (pnt) qp_free_point(pnt);
  if (pnt2) qp_free_point(pnt2);
  if (map3) qp_free_map(map3);
  if (map2) qp_free_map(map2);
  if (pmap) qp_free_map(pmap);
  if (map) qp_free_map(map);
  if (hash2) qp_free_pixhash(hash2);
  if (hash) qp_free_pixhash(hash);
  if (memcpy_) qp_free_memory(memcpy_);
  if (mem) qp_free_memory(mem);
  return failed;
}
