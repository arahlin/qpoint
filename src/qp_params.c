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
  void qp_set_rate_##step(double rate) {	  \
    if (rate != qp_params.s_##step.update_rate) { \
      qp_params.s_##step.update_rate = rate;	  \
      qp_params.s_##step.ctime_last = -1;	  \
    }						  \
  }						  \
  void qp_reset_rate_ ##step(void) {		  \
    qp_params.s_##step.ctime_last = -1;		  \
  }						  \
  double qp_get_rate_##step(void) {		  \
    return qp_params.s_##step.update_rate;	  \
  }
STEPFUNCD(daber)
STEPFUNCD(lonlat)
STEPFUNCD(wobble)
STEPFUNCD(dut1)
STEPFUNCD(erot)
STEPFUNCD(npb)
STEPFUNCD(aaber)
STEPFUNCD(refro)

void qp_print_debug(const char *tag, quat_t q) {
  printf("%s quat: [%.6g, %.6g, %.6g, %.6g]\n",
	 tag, q[0], q[1], q[2], q[3]);
}

void qp_init_params(void) {
  if (qp_params.initialized) return;
  qp_set_rate_daber (QP_DO_ALWAYS);
  qp_set_rate_lonlat(QP_DO_ALWAYS);
  qp_set_rate_wobble(QP_DO_NEVER);
  qp_set_rate_dut1  (QP_DO_NEVER);
  qp_set_rate_erot  (QP_DO_ALWAYS);
  qp_set_rate_npb   (10);
  qp_set_rate_aaber (100);
  qp_set_rate_refro (QP_DO_NEVER);
  qp_params.accuracy = 0;
  qp_params.mean_aber = 0;
  qp_params.fast_math = 0;
  qp_params.polconv = 0;
  qp_params.refro_data.height = 30000.;
  qp_params.refro_data.temperature = 0.;
  qp_params.refro_data.pressure = 10.;
  qp_params.refro_data.humidity = 0.;
  qp_params.refro_data.frequency = 150.;
  qp_params.refro_data.lapse_rate = 0.0065;
  qp_params.refro_data.tolerance = 1.0e-8;
  qp_params.refro_data.corr = 0.;
  qp_params.initialized = 1;
}

void qp_set_params(double daber_rate,
		   double lonlat_rate,
		   double wobble_rate,
		   double dut1_rate,
		   double erot_rate,
		   double npb_rate,
		   double aaber_rate,
		   double refro_rate,
		   int accuracy,
		   int mean_aber,
		   int fast_math,
		   int polconv) {
  qp_init_params();
  qp_set_rate_daber (daber_rate);
  qp_set_rate_lonlat(lonlat_rate);
  qp_set_rate_wobble(wobble_rate);
  qp_set_rate_dut1  (dut1_rate);
  qp_set_rate_erot  (erot_rate);
  qp_set_rate_npb   (npb_rate);
  qp_set_rate_aaber (aaber_rate);
  qp_set_rate_refro (refro_rate);
  if (accuracy != qp_params.accuracy) {
    qp_params.accuracy = accuracy;
    qp_reset_rate_npb();
  }
  if (mean_aber != qp_params.mean_aber) {
    qp_params.mean_aber = mean_aber;
    qp_reset_rate_aaber();
  }
  qp_params.fast_math = fast_math;
  qp_params.polconv = polconv;
}

void qp_reset_rate_all(void) {
  qp_reset_rate_daber();
  qp_reset_rate_lonlat();
  qp_reset_rate_wobble();
  qp_reset_rate_dut1();
  qp_reset_rate_erot();
  qp_reset_rate_npb();
  qp_reset_rate_aaber();
  qp_reset_rate_refro();
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

// set/get refraction parameters
#define REFROFUNCD(param)		     \
  void qp_set_refro_##param(double val) {    \
    if (val != qp_params.refro_data.param) { \
      qp_params.refro_data.param = val;	     \
      qp_reset_rate_refro();		     \
    }					     \
  }					     \
  double qp_get_refro_##param() {	     \
    return qp_params.refro_data.param;	     \
  }
REFROFUNCD(height)
REFROFUNCD(temperature)
REFROFUNCD(pressure)
REFROFUNCD(humidity)
REFROFUNCD(frequency)
REFROFUNCD(lapse_rate)
REFROFUNCD(tolerance)

// use with caution
void qp_set_refro_corr(double val) {
  qp_params.refro_data.corr = val;
}
double qp_get_refro_corr() {
  return qp_params.refro_data.corr;
}

// update all refro_data parameters
void qp_set_refro_data(double height, double temperature, double pressure,
		       double humidity, double frequency, double lapse_rate,
		       double tolerance) {
  qp_set_refro_height(height);
  qp_set_refro_temperature(temperature);
  qp_set_refro_pressure(pressure);
  qp_set_refro_humidity(humidity);
  qp_set_refro_frequency(frequency);
  qp_set_refro_lapse_rate(lapse_rate);
  qp_set_refro_tolerance(tolerance);
}
