
/*
 IERS Bulletin A interpolation
 2012 Mike Nolta <mike@nolta.net>
 2015 ASR <arahlin@princeton.edu>
 2017 ASR <arahlin@fnal.gov>
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "qpoint.h"

static qp_bulletina_entry_t bulletinA_factory[] = {
  {0., 0., 0.}, // NULL
};

static const int mjd_min_factory = 0;
static const int mjd_max_factory = 0;

int qp_get_iers_bulletin_a( qp_memory_t *mem, double mjd,
                           double *dut1, double *x, double *y )
{
  qp_bulletina_t *B = &mem->bulletinA;

  if (B->entries == NULL) {
    B->entries = bulletinA_factory;
    B->mjd_min = mjd_min_factory;
    B->mjd_max = mjd_max_factory;
  }
  if ( (B->mjd_min <= mjd) && (mjd < B->mjd_max) )
  {
    double mjd_floor;
    double r = modf( mjd, &mjd_floor );
    int k = (int) mjd_floor - B->mjd_min;
    qp_bulletina_entry_t a = B->entries[k];
    qp_bulletina_entry_t b = B->entries[k+1];
    // Detect leap seconds
    double leap = b.dut1 - a.dut1;
    if (leap > 0.5)
      leap = 1.;
    else if (leap < -0.5)
      leap = -1.;
    else
      leap = 0;

    *dut1 = (1-r)*a.dut1 + r*(b.dut1-leap);
    *x = (1-r)*a.x + r*b.x;
    *y = (1-r)*a.y + r*b.y;
  }
  else
  {
    // silent failure!
    *dut1 = 0;
    *x = 0;
    *y = 0;
    return 1;
  }

  return 0;
}

int qp_set_iers_bulletin_a( qp_memory_t *mem, int mjd_min_, int mjd_max_,
                           double *dut1, double *x, double *y)
{
  qp_bulletina_t *B = &mem->bulletinA;

  if (dut1 == NULL) {
    // Destroy
    if (B->entries != bulletinA_factory && B->entries != NULL) {
      free(B->entries);
      B->entries = NULL;
    }
    return 0;
  }
  // Create
  B->mjd_min = mjd_min_;
  B->mjd_max = mjd_max_;

  int n_mjd = B->mjd_max - B->mjd_min + 1;
  B->entries = malloc(n_mjd*sizeof(*(B->entries)));
  if (B->entries == NULL)
    return 1;

  for (int k=0; k<n_mjd; k++) {
    B->entries[k].x = x[k];
    B->entries[k].y = y[k];
    B->entries[k].dut1 = dut1[k];
  }
  return 0;
}

int qp_copy_iers_bulletin_a(qp_memory_t *memdest, qp_memory_t *memsrc) {
  qp_bulletina_t *bdest = &memdest->bulletinA;
  qp_bulletina_t *bsrc = &memsrc->bulletinA;


  // nothing to do, since both source and dest point to factory
  if (bdest->entries == bsrc->entries && bsrc->entries == bulletinA_factory)
    return 0;

  if (bsrc->entries != bulletinA_factory && bsrc->entries != NULL) {
    int n_mjd = bsrc->mjd_max - bsrc->mjd_min + 1;
    bdest->mjd_max = bsrc->mjd_max;
    bdest->mjd_min = bsrc->mjd_min;
    if (bdest->entries != bulletinA_factory && bdest->entries != NULL) {
      if (bdest->entries != bsrc->entries)
        free(bdest->entries);
      bdest->entries = NULL;
    }
    bdest->entries = malloc(n_mjd*sizeof(*(bdest->entries)));
    if (bdest->entries == NULL)
      return 1;
    memcpy(bdest->entries, bsrc->entries, n_mjd*sizeof(*(bsrc->entries)));
  }

  return 0;
}
