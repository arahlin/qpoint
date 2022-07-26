#!/usr/bin/env python
#
# The following IERS data are stored:
#
#   * Bull. A PM-x (sec. of arc)
#   * Bull. A PM-y (sec. of arc)
#   * Bull. A UT1-UTC (sec. of time)
#
# 2012 Mike Nolta <mike@nolta.net>
# 2015 ASR <arahlin@princeton.edu>
# 2017 ASR <arahlin@fnal.gov> -- use astropy utilities for auto update

header = """
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

"""

footer = r"""
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
"""

if __name__ == '__main__':

    f = open( "qp_iers_bulletin_a.c", "w" )
    f.write( header )

    f.write( "static qp_bulletina_entry_t bulletinA_factory[] = {\n" )

    try:
        from astropy.utils.iers import IERS_Auto

    except ImportError:
        from warnings import warn, simplefilter
        simplefilter('always')
        warn('Compatible Astropy not found, creating an empty IERS-A table. '
             'Install astropy v1.2 or newer for accurate polar motion '
             'and UT1 corrections',
             ImportWarning)

        f.write( " {0., 0., 0.}, // NULL\n")
        mjd_min = 0
        mjd_max = 0

    else:
        import numpy as np

        columns = ['year', 'month', 'day', 'MJD', 'PM_x', 'PM_y', 'UT1_UTC']
        iers_table = IERS_Auto.open()[columns].as_array()

        # check year
        year = iers_table['year'] + 1900
        wraps, = np.where(np.ediff1d(year) < 0)
        for idx in wraps:
            year[idx + 1:] += 100
        iers_table['year'] = year
        iers_table = iers_table[year >= 2000]

        # check MJD
        mjds = iers_table['MJD']
        mjd_min = int(mjds.min())
        mjd_max = int(mjds.max())

        # populate factory table
        for row in iers_table:
            f.write( "  {%.6f, %.6f, %.7f},  // %d-%02d-%02d\n" % (
                row['PM_x'], row['PM_y'], row['UT1_UTC'],
                row['year'], row['month'], row['day']) )

    f.write( "};\n\n" )
    
    f.write( "static const int mjd_min_factory = %d;\n" % mjd_min )
    f.write( "static const int mjd_max_factory = %d;\n" % mjd_max )
    
    f.write( footer )
    f.close()
