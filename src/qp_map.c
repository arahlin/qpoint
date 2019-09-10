#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include "qpoint.h"
#include "fast_math.h"
#include "vec3.h"
#include "quaternion.h"
#include "chealpix.h"
#ifdef _OPENMP
#include <omp.h>
#endif

qp_det_t * qp_init_det(quat_t q_off, double weight, double gain, mueller_t mueller) {
  qp_det_t *det = malloc(sizeof(*det));

  memcpy(det->q_off, q_off, sizeof(quat_t));
  det->weight = weight;
  det->gain = gain;
  memcpy(det->mueller, mueller, sizeof(mueller_t));

  det->n = 0;

  det->tod_init = 0;
  det->tod = NULL;

  det->flag_init = 0;
  det->flag = NULL;

  det->weights_init = 0;
  det->weights = NULL;

  det->init = QP_STRUCT_INIT | QP_STRUCT_MALLOC;
  return det;
}

qp_det_t * qp_default_det(void) {
  quat_t q = {1, 0, 0, 0};
  mueller_t m = {1, 1, 0, 1};
  return qp_init_det(q, 1.0, 1.0, m);
}

void qp_init_det_tod(qp_det_t *det, size_t n) {
  det->n = n;
  det->tod = calloc(n, sizeof(double));
  det->tod_init = QP_ARR_MALLOC_1D;
}

void qp_init_det_tod_from_array(qp_det_t *det, double *tod, size_t n, int copy) {
  if (copy) {
    qp_init_det_tod(det, n);
    memcpy(det->tod, tod, n * sizeof(double));
    return;
  }

  det->n = n;
  det->tod = tod;
  det->tod_init = QP_ARR_INIT_PTR;
}

void qp_init_det_flag(qp_det_t *det, size_t n) {
  det->n = n;
  det->flag = calloc(n, sizeof(uint8_t));
  det->flag_init = QP_ARR_MALLOC_1D;
}

void qp_init_det_flag_from_array(qp_det_t *det, uint8_t *flag, size_t n, int copy) {
  if (copy) {
    qp_init_det_flag(det, n);
    memcpy(det->flag, flag, n * sizeof(uint8_t));
    return;
  }

  det->n = n;
  det->flag = flag;
  det->flag_init = QP_ARR_INIT_PTR;
}

void qp_init_det_weights(qp_det_t *det, size_t n) {
  det->n = n;
  det->weights = calloc(n, sizeof(double));
  det->weights_init = QP_ARR_MALLOC_1D;
}

void qp_init_det_weights_from_array(qp_det_t *det, double *weights, size_t n,
                                    int copy) {
  if (copy) {
    qp_init_det_weights(det, n);
    memcpy(det->weights, weights, n * sizeof(double));
    return;
  }

  det->n = n;
  det->weights = weights;
  det->weights_init = QP_ARR_INIT_PTR;
}

void qp_free_det(qp_det_t *det) {
  if (det->tod_init & QP_ARR_MALLOC_1D)
    free(det->tod);
  if (det->flag_init & QP_ARR_MALLOC_1D)
    free(det->flag);
  if (det->weights_init & QP_ARR_MALLOC_1D)
    free(det->weights);
  if (det->init & QP_STRUCT_MALLOC)
    free(det);
}

qp_detarr_t * qp_init_detarr(quat_t *q_off, double *weight, double *gain,
                             mueller_t *mueller, size_t n) {
  qp_detarr_t *dets = malloc(sizeof(*dets));
  qp_det_t *det;
  dets->n = n;
  dets->init = QP_STRUCT_INIT | QP_STRUCT_MALLOC;
  dets->arr = malloc(n * sizeof(*det));
  dets->arr_init = QP_ARR_MALLOC_1D;
  dets->diff = 0;

  for (size_t ii = 0; ii < n; ii++) {
    det = dets->arr + ii;
    memcpy(det->q_off, q_off[ii], sizeof(quat_t));
    det->weight = weight[ii];
    det->gain = gain[ii];
    memcpy(det->mueller, mueller[ii], sizeof(mueller_t));
    det->n = 0;
    det->tod_init = 0;
    det->tod = NULL;
    det->flag_init = 0;
    det->flag = NULL;
    det->weights_init = 0;
    det->weights = NULL;
    det->init = QP_STRUCT_INIT;
  }

  return dets;
}

void qp_init_detarr_tod(qp_detarr_t *dets, size_t n) {
  for (size_t ii = 0; ii < dets->n; ii++) {
    qp_init_det_tod(dets->arr + ii, n);
  }
}

void qp_init_detarr_tod_from_array(qp_detarr_t *dets, double **tod,
                                   size_t n, int copy) {
  for (size_t ii = 0; ii < dets->n; ii++) {
    qp_init_det_tod_from_array(dets->arr + ii, tod[ii], n, copy);
  }
}

void qp_init_detarr_tod_from_array_1d(qp_detarr_t *dets, double *tod,
                                      size_t n_chunk, int copy) {
  for (size_t ii = 0; ii < dets->n; ii++) {
    qp_init_det_tod_from_array(dets->arr + ii, tod + ii * n_chunk,
                               n_chunk, copy);
  }
}

void qp_init_detarr_flag(qp_detarr_t *dets, size_t n) {
  for (size_t ii = 0; ii < dets->n; ii++) {
    qp_init_det_flag(dets->arr + ii, n);
  }
}

void qp_init_detarr_flag_from_array(qp_detarr_t *dets, uint8_t **flag,
                                    size_t n, int copy) {
  for (size_t ii = 0; ii < dets->n; ii++) {
    qp_init_det_flag_from_array(dets->arr + ii, flag[ii], n, copy);
  }
}

void qp_init_detarr_flag_from_array_1d(qp_detarr_t *dets, uint8_t *flag,
                                       size_t n_chunk, int copy) {
  for (size_t ii = 0; ii < dets->n; ii++) {
    qp_init_det_flag_from_array(dets->arr + ii, flag + ii * n_chunk,
                                n_chunk, copy);
  }
}

void qp_init_detarr_weights(qp_detarr_t *dets, size_t n) {
  for (size_t ii = 0; ii < dets->n; ii++) {
    qp_init_det_weights(dets->arr + ii, n);
  }
}

void qp_init_detarr_weights_from_array(qp_detarr_t *dets, double **weights,
                                       size_t n, int copy) {
  for (size_t ii = 0; ii < dets->n; ii++) {
    qp_init_det_weights_from_array(dets->arr + ii, weights[ii], n, copy);
  }
}

void qp_init_detarr_weights_from_array_1d(qp_detarr_t *dets, double *weights,
                                          size_t n_chunk, int copy) {
  for (size_t ii = 0; ii < dets->n; ii++) {
    qp_init_det_weights_from_array(dets->arr + ii, weights + ii * n_chunk,
                                   n_chunk, copy);
  }
}

void qp_free_detarr(qp_detarr_t *dets) {
  for (size_t ii = 0; ii < dets->n; ii++) {
    qp_free_det(dets->arr + ii);
  }
  if (dets->arr_init & QP_ARR_MALLOC_1D)
    free(dets->arr);
  if (dets->init & QP_STRUCT_MALLOC)
    free(dets);
  else
    memset(dets, 0, sizeof(*dets));
}

qp_point_t * qp_init_point(size_t n, int time, int pol) {
  qp_point_t *pnt = malloc(sizeof(*pnt));

  pnt->n = n;

  pnt->ctime_init = QP_ARR_MALLOC_1D;
  if (time)
    pnt->ctime = malloc(n * sizeof(double));

  pnt->q_hwp_init = QP_ARR_MALLOC_1D;
  if (pol)
    pnt->q_hwp = malloc(n * sizeof(quat_t));

  pnt->q_bore_init = QP_ARR_MALLOC_1D;
  pnt->q_bore = malloc(n * sizeof(quat_t));

  pnt->init = QP_STRUCT_INIT | QP_STRUCT_MALLOC;
  return pnt;
}

qp_point_t *qp_init_point_from_arrays(quat_t *q_bore, double *ctime, quat_t *q_hwp,
                                      size_t n, int copy) {
  if (copy) {
    qp_point_t *pnt = qp_init_point(n, ctime ? 1 : 0, q_hwp ? 1 : 0);

    memcpy(pnt->q_bore, q_bore, n * sizeof(quat_t));
    memcpy(pnt->ctime, ctime, n * sizeof(double));
    memcpy(pnt->q_hwp, q_hwp, n * sizeof(quat_t));

    return pnt;
  }

  qp_point_t *pnt = malloc(sizeof(*pnt));

  pnt->n = n;
  pnt->q_bore_init = QP_ARR_INIT_PTR;
  pnt->q_bore = q_bore;

  if (ctime) {
    pnt->ctime_init = QP_ARR_INIT_PTR;
    pnt->ctime = ctime;
  } else {
    pnt->ctime_init = 0;
    pnt->ctime = NULL;
  }

  if (q_hwp) {
    pnt->q_hwp_init = QP_ARR_INIT_PTR;
    pnt->q_hwp = q_hwp;
  } else {
    pnt->q_hwp_init = 0;
    pnt->q_hwp = NULL;
  }

  pnt->init = QP_STRUCT_INIT | QP_STRUCT_MALLOC;
  return pnt;
}

void qp_free_point(qp_point_t *pnt) {
  if (pnt->q_bore_init & QP_ARR_MALLOC_1D)
    free(pnt->q_bore);
  if (pnt->q_hwp_init & QP_ARR_MALLOC_1D)
    free(pnt->q_hwp);
  if (pnt->ctime_init & QP_ARR_MALLOC_1D)
    free(pnt->ctime);
  if (pnt->init & QP_STRUCT_MALLOC)
    free(pnt);
  else
    memset(pnt, 0, sizeof(*pnt));
}

void qp_num_maps(qp_vec_mode vec_mode, qp_proj_mode proj_mode,
                 size_t *num_vec, size_t *num_proj) {
  size_t nm = 0;
  switch (vec_mode) {
    case QP_VEC_TEMP:
      nm = 1;
      break;
    case QP_VEC_D1:
    case QP_VEC_POL:
      nm = 3;
      break;
    case QP_VEC_VPOL:
      nm = 4;
      break;
    case QP_VEC_D1_POL:
      nm = 6;
      break;
    case QP_VEC_D2:
      nm = 9;
      break;
    case QP_VEC_D2_POL:
      nm = 18;
      break;
    default:
      nm = 0;
  }
  *num_vec = nm;

  size_t np = 0;
  switch (proj_mode) {
    case QP_PROJ_TEMP:
      np = 1;
      break;
    case QP_PROJ_POL:
      np = 6;
      break;
    case QP_PROJ_VPOL:
      np = 10;
      break;
    default:
      np = 0;
  }
  *num_proj = np;
}

// if npix != 0 then partial map
qp_map_t * qp_init_map(size_t nside, size_t npix, qp_vec_mode vec_mode,
                       qp_proj_mode proj_mode) {
  qp_map_t *map = malloc(sizeof(*map));

  map->nside = nside;
  map->npix = (npix == 0) ? nside2npix(nside) : (long) npix;
  map->partial = (npix > 0);
  map->pixinfo_init = 0;
  map->pixinfo = NULL;
  map->pixhash_init = 0;
  map->pixhash = NULL;

  qp_num_maps(vec_mode, proj_mode, &map->num_vec, &map->num_proj);

  map->vec_mode = vec_mode;
  if (map->num_vec) {
    map->vec = malloc(map->num_vec * sizeof(double *));
    for (size_t ii = 0; ii < map->num_vec; ii++)
      map->vec[ii] = calloc(map->npix, sizeof(double));
    map->vec_init = QP_ARR_MALLOC_1D | QP_ARR_MALLOC_2D;
  } else {
    map->vec_init = 0;
  }
  map->vec1d_init = 0;
  map->vec1d = NULL;

  map->proj_mode = proj_mode;
  if (map->num_proj) {
    map->proj = malloc(map->num_proj * sizeof(double *));
    for (size_t ii = 0; ii < map->num_proj; ii++)
      map->proj[ii] = calloc(map->npix, sizeof(double));
    map->proj_init = QP_ARR_MALLOC_1D | QP_ARR_MALLOC_2D;
  } else {
    map->proj_init = 0;
  }
  map->proj1d_init = 0;
  map->proj1d = NULL;

  map->init = QP_STRUCT_INIT | QP_STRUCT_MALLOC;
  return map;
}

qp_map_t * qp_init_map_from_arrays(double **vec, double **proj, size_t nside,
                                   size_t npix, qp_vec_mode vec_mode,
                                   qp_proj_mode proj_mode, int copy) {
  if (copy) {
    qp_map_t *map = qp_init_map(nside, npix, vec_mode, proj_mode);

    if (map->num_vec)
      for (size_t ii = 0; ii < map->num_vec; ii++)
        memcpy(map->vec[ii], vec[ii], map->npix * sizeof(double));
    if (map->num_proj)
      for (size_t ii = 0; ii < map->num_proj; ii++)
        memcpy(map->proj[ii], proj[ii], map->npix * sizeof(double));

    return map;
  }

  qp_map_t *map = malloc(sizeof(*map));

  map->nside = nside;
  map->npix = (npix == 0) ? nside2npix(nside) : (long) npix;
  map->partial = (npix > 0);
  map->pixinfo_init = 0;
  map->pixinfo = NULL;
  map->pixhash_init = 0;
  map->pixhash = NULL;

  qp_num_maps(vec_mode, proj_mode, &map->num_vec, &map->num_proj);
  map->vec_mode = vec_mode;
  map->proj_mode = proj_mode;

  if (map->num_vec) {
    map->vec = vec;
    map->vec_init = QP_ARR_INIT_PTR;
  } else {
    map->vec_init = 0;
  }
  map->vec1d_init = 0;
  map->vec1d = NULL;

  if (map->num_proj) {
    map->proj = proj;
    map->proj_init = QP_ARR_INIT_PTR;
  } else {
    map->proj_init = 0;
  }
  map->proj1d_init = 0;
  map->proj1d = NULL;

  map->init = QP_STRUCT_INIT | QP_STRUCT_MALLOC;
  return map;
}

qp_map_t *
qp_init_map_from_arrays_1d(double *vec, double *proj, size_t nside, size_t npix,
                           qp_vec_mode vec_mode, qp_proj_mode proj_mode, int copy) {
  if (copy) {
    qp_map_t *map = qp_init_map(nside, npix, vec_mode, proj_mode);

    if (map->num_vec)
      for (size_t ii = 0; ii < map->num_vec; ii++)
        memcpy(map->vec[ii], vec + ii * map->npix, map->npix * sizeof(double));
    if (map->num_proj)
      for (size_t ii = 0; ii < map->num_proj; ii++)
        memcpy(map->proj[ii], proj + ii * map->npix, map->npix * sizeof(double));

    return map;
  }

  qp_map_t *map = malloc(sizeof(*map));

  map->nside = nside;
  map->npix = (npix == 0) ? nside2npix(nside) : (long) npix;
  map->partial = (npix > 0);
  map->pixinfo_init = 0;
  map->pixinfo = NULL;
  map->pixhash_init = 0;
  map->pixhash = NULL;

  qp_num_maps(vec_mode, proj_mode, &map->num_vec, &map->num_proj);
  map->vec_mode = vec_mode;
  map->proj_mode = proj_mode;

  if (map->num_vec) {
    map->vec = malloc(map->num_vec * sizeof(double *));
    for (size_t ii = 0; ii < map->num_vec; ii++)
      map->vec[ii] = vec + ii * map->npix;
    map->vec_init = QP_ARR_MALLOC_1D;
  } else {
    map->vec_init = 0;
  }
  map->vec1d_init = 0;
  map->vec1d = NULL;

  if (map->num_proj) {
    map->proj = malloc(map->num_proj * sizeof(double *));
    for (size_t ii = 0; ii < map->num_proj; ii++)
      map->proj[ii] = proj + ii * map->npix;
    map->proj_init = QP_ARR_MALLOC_1D;
  } else {
    map->proj_init = 0;
  }
  map->proj1d_init = 0;
  map->proj1d = NULL;

  map->init = QP_STRUCT_INIT | QP_STRUCT_MALLOC;
  return map;
}

qp_map_t * qp_init_map_1d(size_t nside, size_t npix, qp_vec_mode vec_mode,
                          qp_proj_mode proj_mode) {
  qp_map_t *map = malloc(sizeof(*map));

  map->nside = nside;
  map->npix = (npix == 0) ? nside2npix(nside) : (long) npix;
  map->partial = (npix > 0);
  map->pixinfo_init = 0;
  map->pixinfo = NULL;
  map->pixhash_init = 0;
  map->pixhash = NULL;

  qp_num_maps(vec_mode, proj_mode, &map->num_vec, &map->num_proj);

  map->vec_mode = vec_mode;
  if (map->num_vec) {
    map->vec1d = calloc(map->num_vec * map->npix, sizeof(double));
    map->vec1d_init = QP_ARR_MALLOC_1D;
    map->vec = malloc(map->num_vec * sizeof(double *));
    for (size_t ii = 0; ii < map->num_vec; ii++)
      map->vec[ii] = map->vec1d + ii * map->npix;
    map->vec_init = QP_ARR_MALLOC_1D;
  } else {
    map->vec1d_init = 0;
    map->vec_init = 0;
  }

  map->proj_mode = proj_mode;
  if (map->num_proj) {
    map->proj1d = calloc(map->num_proj * map->npix, sizeof(double));
    map->proj1d_init = QP_ARR_MALLOC_1D;
    map->proj = malloc(map->num_proj * sizeof(double *));
    for (size_t ii = 0; ii < map->num_proj; ii++)
      map->proj[ii] = map->proj1d + ii * map->npix;
    map->proj_init = QP_ARR_MALLOC_1D;
  } else {
    map->proj1d_init = 0;
    map->proj_init = 0;
  }

  map->init = QP_STRUCT_INIT | QP_STRUCT_MALLOC;
  return map;
}

// if blank, malloc fresh arrays
// otherwise, if copy, copy arrays
// otherwise, point to arrays
// pixhash is copied if it exists
qp_map_t * qp_init_map_from_map(qp_map_t *map, int blank, int copy) {
  size_t npix = map->partial ? map->npix : 0;
  qp_map_t *new_map;

  if (blank)
    new_map = qp_init_map(map->nside, npix, map->vec_mode, map->proj_mode);
  else
    new_map = qp_init_map_from_arrays(map->vec, map->proj, map->nside, npix,
                                      map->vec_mode, map->proj_mode, copy);

  if (map->pixhash_init) {
    new_map->pixhash = qp_copy_pixhash(map->pixhash);
    new_map->pixhash_init = new_map->pixhash->init;
  }

  return new_map;
}

// convert 1d map to 2d
int qp_reshape_map(qp_map_t *map) {
  if (map->vec1d_init) {
    if (map->vec_init & QP_ARR_MALLOC_2D) {
      for (size_t ii = 0; ii < map->num_vec; ii++)
        free(map->vec[ii]);
      map->vec_init &= ~QP_ARR_MALLOC_2D;
    }
    if (!(map->vec_init & QP_ARR_MALLOC_1D)) {
      map->vec = malloc(map->num_vec * sizeof(double *));
      map->vec_init |= QP_ARR_MALLOC_1D;
    }
    for (size_t ii = 0; ii < map->num_vec; ii++)
      map->vec[ii] = map->vec1d + ii * map->npix;
  }

  if (map->proj1d_init) {
    if (map->proj_init & QP_ARR_MALLOC_2D) {
      for (size_t ii = 0; ii < map->num_proj; ii++)
        free(map->proj[ii]);
      map->proj_init &= ~QP_ARR_MALLOC_2D;
    }
    if (!(map->proj_init & QP_ARR_MALLOC_1D)) {
      map->proj = malloc(map->num_proj * sizeof(double *));
      map->proj_init |= QP_ARR_MALLOC_1D;
    }
    for (size_t ii = 0; ii < map->num_proj; ii++)
      map->proj[ii] = map->proj1d + ii * map->npix;
  }

  return 0;
}

int qp_init_map_pixhash(qp_map_t *map, long *pix, size_t npix) {
  if (!map->init)
    return QP_ERROR_INIT;
  if (npix != map->npix)
    return QP_ERROR_INIT;
  map->pixhash = qp_init_pixhash(pix, npix);
  map->pixhash_init = map->pixhash->init;
  return 0;
}

int qp_init_map_pixinfo(qp_map_t *map) {
  if (!map->init)
    return QP_ERROR_INIT;
  map->pixinfo = qp_init_pixinfo(map->nside, 0);
  map->pixinfo_init = map->pixinfo->init;
  return 0;
}

void qp_free_map(qp_map_t *map) {
  if (map->vec1d_init & QP_ARR_MALLOC_1D)
    free(map->vec1d);
  if (map->vec_init & QP_ARR_MALLOC_2D)
    for (size_t ii = 0; ii < map->num_vec; ii++)
      free(map->vec[ii]);
  if (map->vec_init & QP_ARR_MALLOC_1D)
    free(map->vec);

  if (map->proj1d_init & QP_ARR_MALLOC_1D)
    free(map->proj1d);
  if (map->proj_init & QP_ARR_MALLOC_2D)
    for (size_t ii = 0; ii < map->num_proj; ii++)
      free(map->proj[ii]);
  if (map->proj_init & QP_ARR_MALLOC_1D)
    free(map->proj);

  if (map->pixinfo_init)
    qp_free_pixinfo(map->pixinfo);

  if (map->pixhash_init)
    qp_free_pixhash(map->pixhash);

  if (map->init & QP_STRUCT_MALLOC)
    free(map);
  else
    memset(map, 0, sizeof(*map));
}

int qp_add_map(qp_memory_t *mem, qp_map_t *map, qp_map_t *maploc) {

  if (qp_check_error(mem, !map->init, QP_ERROR_INIT,
                     "qp_add_map: map not initialized."))
    return mem->error_code;
  if (qp_check_error(mem, !maploc->init, QP_ERROR_INIT,
                     "qp_add_map: maploc not initialized."))
    return mem->error_code;
  if (qp_check_error(mem, map->vec_mode != maploc->vec_mode, QP_ERROR_MAP,
                     "qp_add_map: vec_modes differ."))
    return mem->error_code;
  if (qp_check_error(mem, map->proj_mode != maploc->proj_mode, QP_ERROR_MAP,
                     "qp_add_map: proj_modes differ."))
    return mem->error_code;
  if (qp_check_error(mem, map->nside != maploc->nside, QP_ERROR_MAP,
                     "qp_add_map: nsides differ."))
    return mem->error_code;
  if (qp_check_error(mem, map->npix != maploc->npix, QP_ERROR_MAP,
                     "qp_add_map: npixs differ."))
    return mem->error_code;

  if (map->vec_init && maploc->vec_init && map->vec_mode) {
    for (size_t ii = 0; ii < map->num_vec; ii++)
      for (size_t ipix = 0; ipix < map->npix; ipix++)
        if (maploc->vec[ii][ipix] != 0)
          map->vec[ii][ipix] += maploc->vec[ii][ipix];
  }

  if (map->proj_init && maploc->proj_init && map->proj_mode) {
    for (size_t ii = 0; ii < map->num_proj; ii++)
      for (size_t ipix = 0; ipix < map->npix; ipix++)
        if (maploc->proj[ii][ipix] != 0)
          map->proj[ii][ipix] += maploc->proj[ii][ipix];
  }

  return 0;
}

int qp_tod2map1_diff(qp_memory_t *mem, qp_det_t *det, qp_det_t *det_pair,
                     qp_point_t *pnt, qp_map_t *map) {

  double spp, cpp, spp_p, cpp_p, ctime, delta;
  double alpha = 0, beta = 0, gamma = 0;
  double walpha = 0, wbeta = 0, wgamma = 0;
  long ipix, ipix_p;
  quat_t q,q_p;
  double w0 = det->weight;
  double g = det->gain;
  double *m = det->mueller;
  double w = w0;

  double w0_p = det_pair->weight;
  double g_p = det_pair->gain;
  double *m_p = det_pair->mueller;
  double w_p = w0_p;
  double wd = 0.5 * (w + w_p);
  double mtd = 0.5 * (m[0] + m_p[0]);

  if (qp_check_error(mem, !mem->init, QP_ERROR_INIT,
                     "qp_tod2map1_diff: mem not initialized."))
    return mem->error_code;
  if (qp_check_error(mem, !det->init, QP_ERROR_INIT,
                     "qp_tod2map1_diff: det not initialized."))
    return mem->error_code;
  if (qp_check_error(mem, !det_pair->init, QP_ERROR_INIT,
                     "qp_tod2map1_diff: det not initialized."))
    return mem->error_code;
  if (qp_check_error(mem, !pnt->init, QP_ERROR_INIT,
                     "qp_tod2map1_diff: pnt not initialized."))
    return mem->error_code;
  if (qp_check_error(mem, !map->init, QP_ERROR_INIT,
                     "qp_tod2map1_diff: map not initialized."))
    return mem->error_code;
  if (qp_check_error(mem, map->partial && !map->pixhash_init, QP_ERROR_INIT,
                     "qp_tod2map1_diff: map pixhash not initialized."))
    return mem->error_code;
  if (qp_check_error(mem, !mem->mean_aber && !pnt->ctime_init, QP_ERROR_POINT,
                     "qp_tod2map1_diff: ctime required if not mean_aber"))
    return mem->error_code;

  if (map->vec1d_init && !map->vec_init)
    if (qp_check_error(mem, qp_reshape_map(map), QP_ERROR_INIT,
                       "qp_tod2map1_diff: reshape error"))
      return mem->error_code;

  for (size_t ii = 0; ii < pnt->n; ii++) {
    /* if either samples are flagged then skip */
    if (det->flag_init || det_pair->flag_init){
      if(det->flag[ii] || det_pair->flag[ii]){
	continue;
      }
    }
    ctime = pnt->ctime_init ? pnt->ctime[ii] : 0;
    if (pnt->q_hwp_init){
      qp_bore2det_hwp(mem, det->q_off, ctime, pnt->q_bore[ii],
                      pnt->q_hwp[ii], q);
      qp_bore2det_hwp(mem, det_pair->q_off, ctime, pnt->q_bore[ii],
                      pnt->q_hwp[ii], q_p);
    }else{
      qp_bore2det(mem, det->q_off, ctime, pnt->q_bore[ii], q);
      qp_bore2det(mem, det_pair->q_off, ctime, pnt->q_bore[ii], q_p);
    }
    qp_quat2pix(mem, q, map->nside, &ipix, &spp, &cpp);
    qp_quat2pix(mem, q_p, map->nside, &ipix_p, &spp_p, &cpp_p);

    if (map->partial) {
      ipix = qp_repixelize(map->pixhash, ipix);
      if (ipix < 0) {
        if (mem->error_missing) {
          qp_set_error(mem, QP_ERROR_MAP,
                       "qp_tod2map1_diff: pixel out of bounds");
          return mem->error_code;
        }
        continue;
      }
      ipix_p = qp_repixelize(map->pixhash, ipix_p);
      if (ipix_p < 0) {
        if (mem->error_missing) {
          qp_set_error(mem, QP_ERROR_MAP,
                       "qp_tod2map1_diff: pair pixel out of bounds");
          return mem->error_code;
        }
        continue;
      }
    }

    if (det->weights_init)
      w = w0 * det->weights[ii];
    if (det_pair->weights_init)
      w_p = w0_p * det_pair->weights[ii];
    if (det->weights_init | det_pair->weights_init)
      wd = 0.5 * (w + w_p);

    if ((map->vec_mode >= QP_VEC_POL) || (map->proj_mode >= QP_PROJ_POL)) {
      alpha = m[1] * cpp - m[2] * spp - (m_p[1] * cpp_p - m_p[2] * spp_p);
      beta = m[2] * cpp + m[1] * spp - (m_p[2] * cpp_p + m_p[1] * spp_p);
      if (!mem->polconv)
        beta *= -1;
      walpha = 0.5 * wd * alpha;
      wbeta = 0.5 * wd * beta;
    }

    if ((map->vec_mode == QP_VEC_VPOL) || (map->proj_mode == QP_PROJ_VPOL)) {
      gamma = m[3] * cpp - m_p[3] * cpp_p;
      wgamma = 0.5 * wd * gamma;
    }

    if (det->tod_init && det_pair->tod_init && map->vec_init) {
      delta = g * det->tod[ii] - g_p * det_pair->tod[ii];

      switch (map->vec_mode) {
      case QP_VEC_VPOL:
        map->vec[3][ipix] += wgamma * delta;
        /* fall through */
      case QP_VEC_POL:
        map->vec[1][ipix] += walpha * delta;
        map->vec[2][ipix] += wbeta * delta;
	/* fall through */
      case QP_VEC_TEMP:
        map->vec[0][ipix] += 0.5 * wd *
          (g * m[0] * det->tod[ii] + g_p * m_p[0] * det_pair->tod[ii]);
	break;
      default:
	break;
      }
    }

    if (map->proj_init) {
      switch(map->proj_mode) {
        case QP_PROJ_VPOL:
          map->proj[0][ipix] += wd * mtd;
          map->proj[1][ipix] += 0.;
          map->proj[2][ipix] += 0.;
          map->proj[3][ipix] += 0.;
          map->proj[4][ipix] += walpha * alpha;
          map->proj[5][ipix] += walpha * beta;
          map->proj[6][ipix] += walpha * gamma;
          map->proj[7][ipix] += wbeta * beta;
          map->proj[8][ipix] += wbeta * gamma;
          map->proj[9][ipix] += wgamma * gamma;
          break;
        case QP_PROJ_POL:
          map->proj[1][ipix] += 0.;
          map->proj[2][ipix] += 0.;
          map->proj[3][ipix] += walpha * alpha;
          map->proj[4][ipix] += walpha * beta;
          map->proj[5][ipix] += wbeta * beta;
	  /* fall through */
        case QP_PROJ_TEMP:
          map->proj[0][ipix] += wd * mtd;
          break;
        default:
          break;
      }
    }
  }

  return 0;
}

int qp_tod2map1(qp_memory_t *mem, qp_det_t *det, qp_point_t *pnt, qp_map_t *map) {

  double spp, cpp, ctime;
  long ipix;
  quat_t q;
  double w0 = det->weight;
  double g = det->gain, gd;
  double *m = det->mueller;
  double w1, mt = m[0], mq = 0, mu = 0, mv = m[3];
  double wmt = w0 * m[0], wmq = 0, wmu = 0, wmv = w0 * m[3];

  if (qp_check_error(mem, !mem->init, QP_ERROR_INIT,
                     "qp_tod2map1: mem not initialized."))
    return mem->error_code;
  if (qp_check_error(mem, !det->init, QP_ERROR_INIT,
                     "qp_tod2map1: det not initialized."))
    return mem->error_code;
  if (qp_check_error(mem, !pnt->init, QP_ERROR_INIT,
                     "qp_tod2map1: pnt not initialized."))
    return mem->error_code;
  if (qp_check_error(mem, !map->init, QP_ERROR_INIT,
                     "qp_tod2map1: map not initialized."))
    return mem->error_code;
  if (qp_check_error(mem, map->partial && !map->pixhash_init, QP_ERROR_INIT,
                     "qp_tod2map1: map pixhash not initialized."))
    return mem->error_code;
  if (qp_check_error(mem, !mem->mean_aber && !pnt->ctime_init, QP_ERROR_POINT,
                     "qp_tod2map1: ctime required if not mean_aber"))
    return mem->error_code;

  if (map->vec1d_init && !map->vec_init)
    if (qp_check_error(mem, qp_reshape_map(map), QP_ERROR_INIT,
                       "qp_tod2map1: reshape error"))
      return mem->error_code;

  for (size_t ii = 0; ii < pnt->n; ii++) {
    if (det->flag_init && det->flag[ii])
      continue;
    ctime = pnt->ctime_init ? pnt->ctime[ii] : 0;
    
    if (pnt->q_hwp_init)
      qp_bore2det_hwp(mem, det->q_off, ctime, pnt->q_bore[ii],
                      pnt->q_hwp[ii], q);
    else
      qp_bore2det(mem, det->q_off, ctime, pnt->q_bore[ii], q);

    qp_quat2pix(mem, q, map->nside, &ipix, &spp, &cpp);

    if (map->partial) {
      ipix = qp_repixelize(map->pixhash, ipix);
      if (ipix < 0) {
        if (mem->error_missing) {
          qp_set_error(mem, QP_ERROR_MAP,
                       "qp_tod2map1: pixel out of bounds");
          return mem->error_code;
        }
        continue;
      }
    }

    if (det->weights_init) {
      w1 = w0 * det->weights[ii];
      wmt = w1 * m[0];
      if ((map->vec_mode == QP_VEC_VPOL) || (map->proj_mode == QP_PROJ_VPOL)) {
        wmv = w1 * m[3];
      }
    } else {
      w1 = w0;
    }

    if ((map->vec_mode >= QP_VEC_POL) || (map->proj_mode >= QP_PROJ_POL)) {
      mq = m[1] * cpp - m[2] * spp;
      mu = m[2] * cpp + m[1] * spp;
      if (!mem->polconv)
        mu *= -1;
      wmq = w1 * mq;
      wmu = w1 * mu;
    }

    if (det->tod_init && map->vec_init) {
      gd = g * det->tod[ii];

      switch (map->vec_mode) {
        case QP_VEC_VPOL:
          map->vec[3][ipix] += wmv * gd;
          /* fall through */
        case QP_VEC_POL:
          map->vec[1][ipix] += wmq * gd;
          map->vec[2][ipix] += wmu * gd;
          /* fall through */
        case QP_VEC_TEMP:
          map->vec[0][ipix] += wmt * gd;
          break;
        default:
          break;
      }
    }

    if (map->proj_init) {
      switch(map->proj_mode) {
        case QP_PROJ_VPOL:
          map->proj[0][ipix] += wmt * mt;
          map->proj[1][ipix] += wmt * mq;
          map->proj[2][ipix] += wmt * mu;
          map->proj[3][ipix] += wmt * mv;
          map->proj[4][ipix] += wmq * mq;
          map->proj[5][ipix] += wmq * mu;
          map->proj[6][ipix] += wmq * mv;
          map->proj[7][ipix] += wmu * mu;
          map->proj[8][ipix] += wmu * mv;
          map->proj[9][ipix] += wmv * mv;
          break;
        case QP_PROJ_POL:
          map->proj[1][ipix] += wmt * mq;
          map->proj[2][ipix] += wmt * mu;
          map->proj[3][ipix] += wmq * mq;
          map->proj[4][ipix] += wmq * mu;
          map->proj[5][ipix] += wmu * mu;
          /* fall through */
        case QP_PROJ_TEMP:
          map->proj[0][ipix] += wmt * mt;
          break;
        default:
          break;
      }
    }
  }

  return 0;
}


int qp_tod2map(qp_memory_t *mem, qp_detarr_t *dets, qp_point_t *pnt,
               qp_map_t *map) {

  if (qp_check_error(mem, !mem->init, QP_ERROR_INIT,
                     "qp_tod2map: mem not initialized."))
    return mem->error_code;
  if (qp_check_error(mem, !dets->init, QP_ERROR_INIT,
                     "qp_tod2map: dets not initialized."))
    return mem->error_code;
  if (qp_check_error(mem, !pnt->init, QP_ERROR_INIT,
                     "qp_tod2map: pnt not initialized."))
    return mem->error_code;
  if (qp_check_error(mem, !map->init, QP_ERROR_INIT,
                     "qp_tod2map: map not initialized."))
    return mem->error_code;
  if (qp_check_error(mem, map->partial && !map->pixhash_init, QP_ERROR_INIT,
                     "qp_tod2map: map pixhash not initialized."))
    return mem->error_code;
  if (qp_check_error(mem, !mem->mean_aber && !pnt->ctime_init, QP_ERROR_POINT,
                     "qp_tod2map: ctime required if not mean_aber"))
    return mem->error_code;

  if (dets->diff == 1){
    /* reset ndet to half its value*/
    dets->n = dets->n/2;
  }

#ifdef _OPENMP
  int num_threads = (int) dets->n < mem->num_threads ? (int) dets->n : mem->num_threads;
  omp_set_num_threads(num_threads);
#endif

  int err = 0;

  if (map->vec1d_init && !map->vec_init)
    if (qp_check_error(mem, qp_reshape_map(map), QP_ERROR_INIT,
                       "qp_tod2map: reshape error"))
      return mem->error_code;

#ifdef DEBUG
  qp_print_memory(mem);
#endif

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    qp_memory_t *memloc = qp_copy_memory(mem);
    const int nthreads = qp_get_opt_num_threads(memloc);

#ifdef DEBUG
    qp_print_memory(memloc);
#endif

    qp_map_t *maploc;
    int errloc = 0;
    if (nthreads > 1)
      maploc = qp_init_map_from_map(map, 1, 0);
    else
      maploc = map;

#ifdef _OPENMP
#pragma omp for
#endif
    for (size_t idet = 0; idet < dets->n; idet++) {
      if (!errloc && !err){
        if(dets->diff == 0){
	  errloc = qp_tod2map1(memloc, dets->arr + idet, pnt, maploc);
	}else{
	  errloc = qp_tod2map1_diff(memloc, dets->arr + idet, dets->arr + idet + dets->n, pnt, maploc);
	}
      }
    }

    if (nthreads > 1) {
      if (!errloc && !err) {
#ifdef _OPENMP
#pragma omp critical
#endif
        errloc = qp_add_map(memloc, map, maploc);
        if (errloc)
#ifdef _OPENMP
#pragma omp atomic
#endif
          err += errloc;
      }
      qp_free_map(maploc);
    }

    if (errloc) {
#ifdef _OPENMP
#pragma omp atomic
#endif
      err += errloc;
#ifdef _OPENMP
#pragma omp critical
#endif
      {
        mem->error_code = memloc->error_code;
        mem->error_string = memloc->error_string;
      }
    }

    qp_free_memory(memloc);
  }

  return err;
}

#define _DATUM(n) (map->vec[n][ipix])
#define DATUM(n) (mt * _DATUM(n))
#define POLDATUM(n) \
  (mt * _DATUM(n) + mq * _DATUM(n+1) + mu * _DATUM(n+2))
#define VPOLDATUM(n) \
  (POLDATUM(n) + mv * _DATUM(n+3))
#define _IDATUM(n) \
  (map->vec[n][pix[0]] * weight[0] + \
   map->vec[n][pix[1]] * weight[1] + \
   map->vec[n][pix[2]] * weight[2] + \
   map->vec[n][pix[3]] * weight[3])
#define IDATUM(n) (mt * _IDATUM(n))
#define IPOLDATUM(n) \
  (mt * _IDATUM(n) + mq * _IDATUM(n+1) + mu * _IDATUM(n+2))
#define IVPOLDATUM(n) \
  (IPOLDATUM(n) + mv * _IDATUM(n+3))

int qp_map2tod1(qp_memory_t *mem, qp_det_t *det, qp_point_t *pnt,
                qp_map_t *map) {

  if (qp_check_error(mem, !mem->init, QP_ERROR_INIT,
                     "qp_map2tod1: mem not initialized."))
    return mem->error_code;
  if (qp_check_error(mem, !det->init, QP_ERROR_INIT,
                     "qp_map2tod1: det not initialized."))
    return mem->error_code;
  if (qp_check_error(mem, !det->tod_init, QP_ERROR_INIT,
                     "qp_map2tod1: det.tod not initialized."))
    return mem->error_code;
  if (qp_check_error(mem, !pnt->init, QP_ERROR_INIT,
                     "qp_map2tod1: pnt not initialized."))
    return mem->error_code;
  if (qp_check_error(mem, !map->init, QP_ERROR_INIT,
                     "qp_map2tod1: map not initialized."))
    return mem->error_code;
  if (qp_check_error(mem, map->partial && !map->pixhash_init, QP_ERROR_INIT,
                     "qp_map2tod1: map pixhash not initialized."))
    return mem->error_code;
  if (qp_check_error(mem, !mem->mean_aber && !pnt->ctime_init, QP_ERROR_POINT,
                     "qp_map2tod1: ctime required if not mean_aber"))
    return mem->error_code;

  double ra, dec, spp, cpp, ctime, dtheta, dphi;
  long ipix;
  quat_t q;
  long pix[4];
  double weight[4];
  double g = det->gain;
  double *m = det->mueller;
  double mt = m[0], mq = 0, mu = 0, mv = m[3];
  int do_interp = (mem->interp_pix &&               \
                   (map->vec_mode == QP_VEC_TEMP || \
                    map->vec_mode == QP_VEC_POL ||  \
                    map->vec_mode == QP_VEC_VPOL));
  int jj, kk, bad_pix = 0;
  double norm1, norm2;

  if (map->vec1d_init && !map->vec_init)
    if (qp_check_error(mem, qp_reshape_map(map), QP_ERROR_INIT,
                       "qp_map2tod1: reshape error"))
      return mem->error_code;

  if (do_interp && !map->pixinfo_init)
    if (qp_check_error(mem, qp_init_map_pixinfo(map), QP_ERROR_INIT,
                       "qp_map2tod1: pixinfo init error"))
      return mem->error_code;

  for (size_t ii = 0; ii < pnt->n; ii++) {
    if (det->flag_init && det->flag[ii])
      continue;
    ctime = pnt->ctime_init ? pnt->ctime[ii] : 0;
    if (pnt->q_hwp_init)
      qp_bore2det_hwp(mem, det->q_off, ctime, pnt->q_bore[ii],
                      pnt->q_hwp[ii], q);
    else
      qp_bore2det(mem, det->q_off, ctime, pnt->q_bore[ii], q);

    if ((map->vec_mode >= QP_VEC_D1) || do_interp) {
      qp_quat2radec(mem, q, &ra, &dec, &spp, &cpp);
      ipix = qp_radec2pix(mem, ra, dec, map->nside);
      qp_pixel_offset(mem, map->nside, ipix, ra, dec, &dtheta, &dphi);
      if (do_interp)
        qp_get_interpol(mem, map->pixinfo, ra, dec, pix, weight);
    } else {
      qp_quat2pix(mem, q, map->nside, &ipix, &spp, &cpp);
    }

    if (map->partial) {
      ipix = qp_repixelize(map->pixhash, ipix);
      if (ipix < 0) {
        if (mem->error_missing) {
          qp_set_error(mem, QP_ERROR_MAP,
                       "qp_map2tod1: pixel out of bounds");
          return mem->error_code;
        } else if (mem->nan_missing) {
          det->tod[ii] = 0.0 / 0.0;
        }
        continue;
      }
      if (do_interp) {
          bad_pix = 0;
          for (jj = 0; jj < 4; jj++) {
            pix[jj] = qp_repixelize(map->pixhash, pix[jj]);
            if (pix[jj] < 0) {
              if (mem->error_missing) {
                qp_set_error(mem, QP_ERROR_MAP,
                             "qp_map2tod1: neighbor pixel out of bounds");
                return mem->error_code;
              } else if (mem->interp_missing) {
                /* compute normalization with/without bad neighbors */
                norm1 = 0.0;
                norm2 = 0.0;
                for (kk = 0; kk < 4; kk++) {
                  norm1 += weight[kk];
                  if (kk != jj)
                    norm2 += weight[kk];
                }
                /* clear pix/weight for bad neighbor */
                pix[jj] = 0;
                weight[jj] = 0;
                /* renormalize */
                for (kk = 0; kk < 4; kk++) {
                  if (kk != jj)
                    weight[kk] *= norm1 / norm2;
                }
                /* count bad neighbors */
                bad_pix += 1;
              } else {
                bad_pix = 1;
                break;
              }
            }
          }
          /* at least one good neighbor remains, so this sample is ok */
          if (mem->interp_missing && bad_pix < 4)
            bad_pix = 0;
          /* fill bad sample with nan or just skip it */
          if (bad_pix) {
            if (mem->nan_missing)
              det->tod[ii] = 0.0 / 0.0;
            continue;
          }
      }
    }

    if ((map->vec_mode >= QP_VEC_POL) || (map->proj_mode >= QP_PROJ_POL)) {
      mq = m[1] * cpp - m[2] * spp;
      mu = m[2] * cpp + m[1] * spp;
      if (!mem->polconv)
        mu *= -1;
    }

    switch (map->vec_mode) {
      case QP_VEC_VPOL:
        if (do_interp)
          det->tod[ii] += g * IVPOLDATUM(0);
        else
          det->tod[ii] += g * VPOLDATUM(0);
        break;
      case QP_VEC_D2_POL:
        det->tod[ii] += g * (dphi * dphi * POLDATUM(15)
                             + dtheta * dphi * POLDATUM(12)
                             + dtheta * dtheta * POLDATUM(9));
        /* fall through */
      case QP_VEC_D1_POL:
        det->tod[ii] += g * (dphi * POLDATUM(6) + dtheta * POLDATUM(3));
        /* fall through */
      case QP_VEC_POL:
        if (do_interp)
          det->tod[ii] += g * IPOLDATUM(0);
        else
          det->tod[ii] += g * POLDATUM(0);
        break;
      case QP_VEC_D2:
        det->tod[ii] += g * (dphi * dphi * DATUM(5) + dtheta * dphi * DATUM(4)
                             + dtheta * dtheta * DATUM(3));
        /* fall through */
      case QP_VEC_D1:
        det->tod[ii] += g * (dphi * DATUM(2) + dtheta * DATUM(1));
        /* fall through */
      case QP_VEC_TEMP:
        if (do_interp)
          det->tod[ii] += g * IDATUM(0);
        else
          det->tod[ii] += g * DATUM(0);
        break;
      default:
        break;
    }
  }

  return 0;
}

int qp_map2tod(qp_memory_t *mem, qp_detarr_t *dets, qp_point_t *pnt,
               qp_map_t *map) {

  if (qp_check_error(mem, !mem->init, QP_ERROR_INIT,
                     "qp_map2tod: mem not initialized."))
    return mem->error_code;
  if (qp_check_error(mem, !dets->init, QP_ERROR_INIT,
                     "qp_map2tod: det not initialized."))
    return mem->error_code;
  if (qp_check_error(mem, !pnt->init, QP_ERROR_INIT,
                     "qp_map2tod: pnt not initialized."))
    return mem->error_code;
  if (qp_check_error(mem, !map->init, QP_ERROR_INIT,
                     "qp_map2tod: map not initialized."))
    return mem->error_code;
  if (qp_check_error(mem, map->partial && !map->pixhash_init, QP_ERROR_INIT,
                     "qp_map2tod: map pixhash not initialized."))
    return mem->error_code;
  if (qp_check_error(mem, !mem->mean_aber && !pnt->ctime_init, QP_ERROR_POINT,
                     "qp_map2tod: ctime required if not mean_aber"))
    return mem->error_code;

#ifdef _OPENMP
  int num_threads = (int) dets->n < mem->num_threads ? (int) dets->n : mem->num_threads;
  omp_set_num_threads(num_threads);
#endif

  int err = 0;

  if (map->vec1d_init && !map->vec_init)
    if (qp_check_error(mem, qp_reshape_map(map), QP_ERROR_INIT,
                       "qp_map2tod: reshape error"))
      return mem->error_code;

#ifdef DEBUG
  qp_print_memory(mem);
#endif

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    qp_memory_t *memloc = qp_copy_memory(mem);
    int errloc = 0;

#ifdef DEBUG
    qp_print_memory(memloc);
#endif

#ifdef _OPENMP
#pragma omp for nowait
#endif
    for (size_t idet = 0; idet < dets->n; idet++) {
      if (!errloc && !err)
        errloc = qp_map2tod1(memloc, dets->arr + idet, pnt, map);
    }

    if (errloc) {
#ifdef _OPENMP
#pragma omp atomic
#endif
      err += errloc;
#ifdef _OPENMP
#pragma omp critical
#endif
      {
        mem->error_code = memloc->error_code;
        mem->error_string = memloc->error_string;
      }
    }

    qp_free_memory(memloc);
  }

  return err;
}

void qp_set_opt_num_threads(qp_memory_t *mem, int num_threads) {
  if (num_threads == 0) {
#ifdef _OPENMP
#pragma omp parallel
    {
      num_threads = omp_get_num_threads();
    }
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
