#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "qpoint.h"
#include <omp.h>

void qp_init_state(qp_state_t *state, double rate) {
  state->update_rate = rate;
  state->ctime_last = -1;
}

qp_memory_t * qp_init_memory(void) {
  qp_memory_t *mem = malloc(sizeof(*mem));
  qp_init_state(&mem->state_daber , QP_DO_ALWAYS);
  qp_init_state(&mem->state_lonlat, QP_DO_ALWAYS);
  qp_init_state(&mem->state_wobble, QP_DO_NEVER);
  qp_init_state(&mem->state_dut1  , QP_DO_NEVER);
  qp_init_state(&mem->state_erot  , QP_DO_ALWAYS);
  qp_init_state(&mem->state_npb   , 10);
  qp_init_state(&mem->state_aaber , 100);
  qp_init_state(&mem->state_ref , QP_DO_NEVER);
  mem->accuracy = 0;
  mem->mean_aber = 0;
  mem->fast_math = 0;
  mem->polconv = 0;
  mem->pair_dets = 0;
  mem->pix_order = 0;
  qp_set_opt_num_threads(mem, 0);
  mem->weather.height = 30000.;
  mem->weather.temperature = 0.;
  mem->weather.pressure = 10.;
  mem->weather.humidity = 0.;
  mem->weather.frequency = 150.;
  mem->weather.lapse_rate = 0.0065;
  mem->ref_tol = 1.0e-8;
  mem->ref_delta = 0.;
  mem->dut1 = 0.;
  memset(mem->q_lonlat,   0, 4);
  memset(mem->q_wobble,   0, 4);
  memset(mem->q_npb,      0, 4);
  memset(mem->q_erot,     0, 4);
  memset(mem->beta_earth, 0, 3);
  mem->bulletinA.entries = NULL;
  mem->initialized = 1;
  return mem;
}

void qp_free_memory(qp_memory_t *mem) {
  set_iers_bulletin_a(mem, 0, 0, NULL, NULL, NULL);
  free(mem);
}

#define RATEFUNCD(state)				    \
  void qp_set_rate_##state(qp_memory_t *mem, double rate) { \
    if (rate != mem->state_##state.update_rate) {	    \
      mem->state_##state.update_rate = rate;		    \
      mem->state_##state.ctime_last = -1;		    \
    }							    \
  }							    \
  void qp_reset_rate_##state(qp_memory_t *mem) {	    \
    mem->state_##state.ctime_last = -1;			    \
  }							    \
  double qp_get_rate_##state(qp_memory_t *mem) {	    \
    return mem->state_##state.update_rate;		    \
  }
RATEFUNCD(daber)
RATEFUNCD(lonlat)
RATEFUNCD(wobble)
RATEFUNCD(dut1)
RATEFUNCD(erot)
RATEFUNCD(npb)
RATEFUNCD(aaber)
RATEFUNCD(ref)

void qp_print_debug(const char *tag, quat_t q) {
  printf("%s quat: [%.6g, %.6g, %.6g, %.6g]\n",
	 tag, q[0], q[1], q[2], q[3]);
}

void qp_set_rates(qp_memory_t *mem,
		  double daber_rate,
		  double lonlat_rate,
		  double wobble_rate,
		  double dut1_rate,
		  double erot_rate,
		  double npb_rate,
		  double aaber_rate,
		  double ref_rate) {
  qp_set_rate_daber (mem, daber_rate);
  qp_set_rate_lonlat(mem, lonlat_rate);
  qp_set_rate_wobble(mem, wobble_rate);
  qp_set_rate_dut1  (mem, dut1_rate);
  qp_set_rate_erot  (mem, erot_rate);
  qp_set_rate_npb   (mem, npb_rate);
  qp_set_rate_aaber (mem, aaber_rate);
  qp_set_rate_ref   (mem, ref_rate);
}

void qp_reset_rates(qp_memory_t *mem) {
  qp_reset_rate_daber (mem);
  qp_reset_rate_lonlat(mem);
  qp_reset_rate_wobble(mem);
  qp_reset_rate_dut1  (mem);
  qp_reset_rate_erot  (mem);
  qp_reset_rate_npb   (mem);
  qp_reset_rate_aaber (mem);
  qp_reset_rate_ref   (mem);
}

// return 0 to skip, 1 to apply
int qp_check_update(qp_state_t *state, double ctime) {
  // don't update if set to never
  if (state->update_rate == QP_DO_NEVER) return 0;
  // don't update if set to once and already done
  if ( (state->update_rate == QP_DO_ONCE) &&
       (state->ctime_last > 0) ) return 0;
  // update if hasn't been checked yet (likely first time)
  if (state->ctime_last <= 0) {
    state->ctime_last = ctime;
    return 1;
  }
  // update if time is discontinuous
  if (ctime < state->ctime_last) {
    state->ctime_last = ctime;
    return 1;
  }
  // update if enough ticks have passed
  if ( (ctime - state->ctime_last) >= state->update_rate ) {
    state->ctime_last = ctime;
    return 1;
  }
  // otherwise don't update
  return 0;
}

int qp_check_apply(qp_state_t *state) {
  if (state->update_rate == QP_DO_NEVER) return 0;
  return 1;
}

#define OPTIONFUNCS(opt)			     \
  void qp_set_opt_##opt(qp_memory_t *mem, int val) { \
    mem->opt = val;				     \
  }
#define OPTIONFUNCSR(opt, rate)			     \
  void qp_set_opt_##opt(qp_memory_t *mem, int val) { \
    if (val != mem->opt) {			     \
      mem->opt = val;				     \
      qp_reset_rate_##rate(mem);		     \
    }						     \
  }
#define OPTIONFUNCG(opt)		   \
  int qp_get_opt_##opt(qp_memory_t *mem) { \
    return mem->opt;			   \
  }
#define OPTIONFUNCD(opt) \
  OPTIONFUNCS(opt)	 \
  OPTIONFUNCG(opt)
#define OPTIONFUNCDR(opt, rate) \
  OPTIONFUNCSR(opt, rate)	\
  OPTIONFUNCG(opt)

OPTIONFUNCDR(accuracy, npb)
OPTIONFUNCDR(mean_aber, aaber)
OPTIONFUNCD(fast_math)
OPTIONFUNCD(polconv)
OPTIONFUNCD(pair_dets)
OPTIONFUNCD(pix_order)

void qp_set_options(qp_memory_t *mem,
		    int accuracy,
		    int mean_aber,
		    int fast_math,
		    int polconv,
		    int pair_dets,
		    int pix_order,
		    int num_threads) {
  qp_set_opt_accuracy   (mem, accuracy);
  qp_set_opt_mean_aber  (mem, mean_aber);
  qp_set_opt_fast_math  (mem, fast_math);
  qp_set_opt_polconv    (mem, polconv);
  qp_set_opt_pair_dets  (mem, pair_dets);
  qp_set_opt_pix_order  (mem, pix_order);
  qp_set_opt_num_threads(mem, num_threads);
}

// update all ref_data parameters
void qp_set_weather(qp_memory_t *mem,
		    double height, double temperature, double pressure,
		    double humidity, double frequency, double lapse_rate) {
  qp_set_weather_height     (mem, height);
  qp_set_weather_temperature(mem, temperature);
  qp_set_weather_pressure   (mem, pressure);
  qp_set_weather_humidity   (mem, humidity);
  qp_set_weather_frequency  (mem, frequency);
  qp_set_weather_lapse_rate (mem, lapse_rate);
}

// set/get refraction parameters
#define WEATHFUNCD(param)				      \
  void qp_set_weather_##param(qp_memory_t *mem, double val) { \
    if (val != mem->weather.param) {			      \
      mem->weather.param = val;				      \
      qp_reset_rate_ref(mem);				      \
    }							      \
  }							      \
  double qp_get_weather_##param(qp_memory_t *mem) {	      \
    return mem->weather.param;				      \
  }
WEATHFUNCD(height)
WEATHFUNCD(temperature)
WEATHFUNCD(pressure)
WEATHFUNCD(humidity)
WEATHFUNCD(frequency)
WEATHFUNCD(lapse_rate)

#define DOUBLEFUNCS(param)			      \
  void qp_set_##param(qp_memory_t *mem, double val) { \
    mem->param = val;				      \
  }
#define DOUBLEFUNCSR(param, rate)		      \
  void qp_set_##param(qp_memory_t *mem, double val) { \
    if (val != mem->param) {			      \
      mem->param = val;				      \
      qp_reset_rate_##rate(mem);		      \
    }						      \
  }
#define DOUBLEFUNCSRR(param, rate1, rate2)	      \
  void qp_set_##param(qp_memory_t *mem, double val) { \
    if (val != mem->param) {			      \
      mem->param = val;				      \
      qp_reset_rate_##rate1(mem);		      \
      qp_reset_rate_##rate2(mem);		      \
    }						      \
  }
#define DOUBLEFUNCG(param)		    \
  double qp_get_##param(qp_memory_t *mem) { \
    return mem->param;			    \
  }
#define DOUBLEFUNCD(param) \
  DOUBLEFUNCS(param)	   \
  DOUBLEFUNCG(param)
#define DOUBLEFUNCDR(param, rate) \
  DOUBLEFUNCSR(param, rate)	  \
  DOUBLEFUNCG(param)
#define DOUBLEFUNCDRR(param, rate1, rate2) \
  DOUBLEFUNCSRR(param, rate1, rate2)	   \
  DOUBLEFUNCG(param)

DOUBLEFUNCDR(ref_tol, ref)
DOUBLEFUNCD(ref_delta)
DOUBLEFUNCDRR(dut1, erot, wobble)
