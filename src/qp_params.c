#include <assert.h>
#include <stdio.h>
#include "qpoint.h"

qp_params_t qp_params;

#define PARAMFUNCD(param)	 \
  void qp_set_##param(int val) { \
    qp_params.param = val;	 \
  }				 \
  int qp_get_##param(void) {	 \
    return qp_params.param;	 \
  }
PARAMFUNCD(accuracy)
PARAMFUNCD(mean_aber)
PARAMFUNCD(fast_math)
PARAMFUNCD(polconv)

#define STEPFUNCD(step)				  \
  void qp_set_##step(double rate) {		  \
    if (rate != qp_params.s_##step.update_rate) { \
      qp_params.s_##step.update_rate = rate;	  \
      qp_params.s_##step.ctime_last = -1;	  \
    }						  \
  }						  \
  void qp_reset_##step(void) {			  \
    qp_params.s_##step.ctime_last = -1;		  \
  }						  \
  double qp_get_##step(void) {			  \
    return qp_params.s_##step.update_rate;	  \
  }
STEPFUNCD(daber)
STEPFUNCD(lonlat)
STEPFUNCD(wobble)
STEPFUNCD(dut1)
STEPFUNCD(erot)
STEPFUNCD(npb)
STEPFUNCD(aaber)

void qp_print_debug(const char *tag, quat_t q) {
  printf("%s quat: [%.6g, %.6g, %.6g, %.6g]\n",
	 tag, q[0], q[1], q[2], q[3]);
}

void qp_init_params(void) {
  if (qp_params.initialized) return;
  qp_set_daber (QP_DO_ALWAYS);
  qp_set_lonlat(QP_DO_ALWAYS);
  qp_set_wobble(QP_DO_NEVER);
  qp_set_dut1  (QP_DO_NEVER);
  qp_set_erot  (QP_DO_ALWAYS);
  qp_set_npb   (10);
  qp_set_aaber (100);
  qp_params.accuracy = 0;
  qp_params.mean_aber = 0;
  qp_params.fast_math = 0;
  qp_params.polconv = 0;
  qp_params.initialized = 1;
}

void qp_set_params(double daber_rate,
		   double lonlat_rate,
		   double wobble_rate,
		   double dut1_rate,
		   double erot_rate,
		   double npb_rate,
		   double aaber_rate,
		   int accuracy,
		   int mean_aber,
		   int fast_math,
		   int polconv) {
  qp_init_params();
  qp_set_daber (daber_rate);
  qp_set_lonlat(lonlat_rate);
  qp_set_wobble(wobble_rate);
  qp_set_dut1  (dut1_rate);
  qp_set_erot  (erot_rate);
  qp_set_npb   (npb_rate);
  qp_set_aaber (aaber_rate);
  qp_params.accuracy = accuracy;
  qp_params.mean_aber = mean_aber;
  qp_params.fast_math = fast_math;
  qp_params.polconv = 0;
}

void qp_reset_all(void) {
  qp_reset_daber();
  qp_reset_lonlat();
  qp_reset_wobble();
  qp_reset_dut1();
  qp_reset_erot();
  qp_reset_npb();
  qp_reset_aaber();
}

// return 0 to skip, 1 to apply
int qp_check_update(qp_step_t *step, double ctime) {
  // don't update if set to never
  if (step->update_rate == QP_DO_NEVER) return 0;
  // don't update if set to once and already done
  if ( (step->update_rate == QP_DO_ONCE) &&
       (step->ctime_last > 0) ) return 0;
  // update if hasn't been checked yet (likely first time)
  if (step->ctime_last <= 0) {
    step->ctime_last = ctime;
    return 1;
  }
  // update if time is discontinuous
  if (ctime < step->ctime_last) {
    step->ctime_last = ctime;
    return 1;
  }
  // update if enough ticks have passed
  if ( (ctime - step->ctime_last) >= step->update_rate ) {
    step->ctime_last = ctime;
    return 1;
  }
  // otherwise don't update
  return 0;
}

int qp_check_apply(qp_step_t *step) {
  if (step->update_rate == QP_DO_NEVER) return 0;
  return 1;
}
