#include <stdlib.h>
#include "qpoint.h"

int main(int argc, char *argv[]) {
  
  /* parameters */
  
  int n = 10000;     // number of samples
  int accuracy = 1;  // "low" accuracy
  int mean_aber = 0; // apply aberration correction per detector
  int fast = 1;      // fast math
  int polconv = 0;   // healpix polarization convention
  
  double rate_daber = QP_DO_ALWAYS;  // update diurnal aberration rotation
  double rate_lonlat = QP_DO_ALWAYS; // update lon/lat rotation
  double rate_wobble = QP_DO_NEVER;  // update wobble correction
  double rate_dut1 = QP_DO_NEVER;    // update ut1-utc time correction
  double rate_erot = QP_DO_ALWAYS;   // update ERA angle
  double rate_npb = 10;              // update nutation/precession/bias correction
  double rate_aaber = 100;           // update annual aberration correction
  double rate_ref = QP_DO_NEVER;     // update refraction correction
  
  // detector offsets in degrees
  double delta_az = 1.0;
  double delta_el = -1.0;
  double delta_psi = 22.5;
  
  /* allocate all the things */
  
  double *az = malloc(sizeof(double)*n);
  double *el = malloc(sizeof(double)*n);
  double *pitch = malloc(sizeof(double)*n);
  double *roll = malloc(sizeof(double)*n);
  double *lon = malloc(sizeof(double)*n);
  double *lat = malloc(sizeof(double)*n);
  double *ctime = malloc(sizeof(double)*n);
  quat_t *q_bore = malloc(sizeof(double)*n);
  
  double *ra = malloc(sizeof(double)*n);
  double *dec = malloc(sizeof(double)*n);
  double *sin2psi = malloc(sizeof(double)*n);
  double *cos2psi = malloc(sizeof(double)*n);
  
  /* initialize memory */
  
  qp_memory_t *mem = qp_init_memory();
  
  // set all parameters at once
  /*
  qp_set_rates(mem, rate_daber, rate_lonlat, rate_wobble, rate_dut1,
	       rate_erot, rate_npb, rate_aaber, rate_ref);
  qp_set_options(mem, accuracy, mean_aber, fast_math, polconv);
  */
  
  // or only a few
  qp_set_accuracy(mem, accuracy);
  qp_set_fast_math(mem, fast_math);
  
  /* Calculate boresight pointing */
  
  // TODO
  // read in az/el/etc data
  
  qp_azel2bore(mem, az, el, pitch, roll, lon, lat, ctime, q_bore, n);
  
  /* Calculate detector pointing */
  
  qp_bore2radec(mem, delta_az, delta_el, delta_psi, ctime, q_bore,
		ra, dec, sin2psi, cos2psi, n);
  
  /* free all the things */
  
  free(az);
  free(el);
  free(pitch);
  free(roll);
  free(lon);
  free(lat);
  free(ctime);
  
  free(ra);
  free(dec);
  free(sin2psi);
  free(cos2psi);
  
  qp_free_memory(mem)
  
  return 0;
}
