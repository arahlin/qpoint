#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include "qpoint.h"
#include "fast_math.h"
#include "vec3.h"
#include "quaternion.h"
#include <chealpix.h>
#ifdef _OPENMP
#include <omp.h>
#endif

/* Compute healpix pixel number for given nside and ra/dec */
long qp_radec2pix(qp_memory_t *mem, double nside, double ra, double dec) {
  long pix;
  if (mem->pix_order == QP_ORDER_NEST)
    ang2pix_nest(nside, M_PI_2 - deg2rad(dec), deg2rad(ra), &pix);
  ang2pix_ring(nside, M_PI_2 - deg2rad(dec), deg2rad(ra), &pix);
  return pix;
}

/* Compute pointing matrix map for given boresight timestream and detector
   offset. pmap is a 6-x-npix array containing (hits, p01, p02, p11, p12, p22) */
void qp_bore2map_single(qp_memory_t *mem,
			double delta_az, double delta_el, double delta_psi,
			double *ctime, quat_t *q_bore, int n,
			pixel_t *pmap, int nside) {
  double ra, dec, sin2psi, cos2psi;
  long ipix;
  
  quat_t q_off, q;
  qp_det_offset(delta_az, delta_el, delta_psi, q_off);
  
  for (int ii=0; ii<n; ii++) {
    qp_bore2det(mem, q_off, ctime[ii], q_bore[ii], q);
    qp_quat2radec(mem, q, &ra, &dec, &sin2psi, &cos2psi);
    ipix = qp_radec2pix(mem, nside, ra, dec);
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
void qp_bore2map_pair(qp_memory_t *mem,
		      double delta_az, double delta_el, double delta_psi,
		      double *ctime, quat_t *q_bore, int n,
		      pixel_t *pmap, int nside) {
  double ra, dec, sin2psi, cos2psi;
  long ipix;
  
  quat_t q_off, q;
  qp_det_offset(delta_az, delta_el, delta_psi, q_off);
  
  for (int ii=0; ii<n; ii++) {
    qp_bore2det(mem, q_off, ctime[ii], q_bore[ii], q);
    qp_quat2radec(mem, q, &ra, &dec, &sin2psi, &cos2psi);
    ipix = qp_radec2pix(mem, nside, ra, dec);
    pmap[ipix][0] += 2;
    pmap[ipix][3] += 2*cos2psi*cos2psi;
    pmap[ipix][4] += 2*cos2psi*sin2psi;
    pmap[ipix][5] += 2*sin2psi*sin2psi;
  }
}

/* Compute pointing matrix map for given boresight timestream and many detector
   offsets. pmap a 6-x-npix array containing (hits, p01, p02, p11, p12, p22).
   openMP-parallelized. */
void qp_bore2map(qp_memory_t *mem,
		 double *delta_az, double *delta_el, double *delta_psi, int ndet,
		 double *ctime, quat_t *q_bore, int n,
		 pixel_t *pmap, int nside) {
  
  long npix = nside2npix(nside);

#pragma omp parallel
  {
    
    // local map array
    pixel_t *pmaploc;
    if (mem->num_threads > 1) {
      pmaploc = calloc(npix,sizeof(pixel_t));
    } else {
      pmaploc = pmap;
    }
    
    // local copy of memory
    qp_memory_t memloc = *mem;
    
#pragma omp for nowait
    for (int idet=0; idet<ndet; idet++) {
#ifdef DEBUG
      printf("thread %d, det %d\n", omp_get_thread_num(), idet);
      printf("offset %f %f %f\n", delta_az[idet], delta_el[idet], delta_psi[idet]);
#endif
      if (mem->pair_dets)
	qp_bore2map_pair(&memloc, delta_az[idet], delta_el[idet], delta_psi[idet],
			 ctime, q_bore, n, pmaploc, nside);
      else
	qp_bore2map_single(&memloc, delta_az[idet], delta_el[idet], delta_psi[idet],
			   ctime, q_bore, n, pmaploc, nside);
    }
    
    if (mem->num_threads > 1) {
      // reduce
#pragma omp critical
      {
	for (int ipix=0; ipix<npix; ipix++)
	  for (int ii=0; ii<6; ii++)
	    pmap[ipix][ii] += pmaploc[ipix][ii];
      }
      
      free(pmaploc);
    }
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
  return mem->num_threads;
}
