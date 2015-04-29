#!/usr/bin/env python
#
# The file:
#
#   iers_bulletin_a.dat
#
# is from:
#
#   http://maia.usno.navy.mil/search/search.html
#
# with options:
#
#   * Bull. A PM-x (sec. of arc)
#   * Bull. A PM-y (sec. of arc)
#   * Bull. A UT1-UTC (sec. of time)
#
# 2012 Mike Nolta <mike@nolta.net>

header = """
/*
 IERS Bulletin A interpolation
 2012 Mike Nolta <mike@nolta.net>
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "qpoint.h"

"""

footer = r"""
int get_iers_bulletin_a( qp_memory_t *mem, double mjd,
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

int set_iers_bulletin_a( qp_memory_t *mem, int mjd_min_, int mjd_max_,
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

int copy_iers_bulletin_a(qp_memory_t *memdest, qp_memory_t *memsrc) {
  qp_bulletina_t *bdest = &memdest->bulletinA;
  qp_bulletina_t *bsrc = &memsrc->bulletinA;

  if (bsrc->entries != bulletinA_factory && bsrc->entries != NULL) {
    int n_mjd = bsrc->mjd_max - bsrc->mjd_min + 1;
    bdest->mjd_max = bsrc->mjd_max;
    bdest->mjd_min = bsrc->mjd_min;
    if (bdest->entries != bulletinA_factory && bdest->entries != NULL) {
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

    f = open( "iers_bulletin_a.c", "w" )
    f.write( header )

    mjds = []
    f.write( "static qp_bulletina_entry_t bulletinA_factory[] = {\n" )

    for line in open('iers_bulletin_a.dat'):
        stripped_line = line.strip()
        if len(stripped_line) == 0 or stripped_line[0] == "#":
            continue
        s = line.rsplit( None, 4 )
        mjd = int(float(s[1]))
        mjds.append( mjd )
        f.write( "  {%sf, %sf, %sf},  // %s\n" % (s[2],s[3],s[4],s[0].strip()) )

    f.write( "};\n\n" )
    
    f.write( "static const int mjd_min_factory = %d;\n" % min(mjds) )
    f.write( "static const int mjd_max_factory = %d;\n" % max(mjds) )
    
    f.write( footer )
    f.close()
