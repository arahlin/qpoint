#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "qpoint.h"
#include <omp.h>

const int QP_DO_ALWAYS = 0;
const int QP_DO_ONCE = -1;
const int QP_DO_NEVER = -999;

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
  qp_init_state(&mem->state_ref   , QP_DO_NEVER);
  mem->accuracy = 0;
  mem->mean_aber = 0;
  mem->fast_math = 0;
  mem->polconv = 0;
  mem->pair_dets = 0;
  mem->pix_order = 0;
  mem->fast_pix = 0;
  mem->gal_init = 0;
  mem->thread_num = 0;
  qp_set_opt_num_threads(mem, 0);
  mem->weather.height = 35000.;
  mem->weather.temperature = 0.;
  mem->weather.pressure = 10.;
  mem->weather.humidity = 0.;
  mem->weather.frequency = 150.;
  mem->weather.lapse_rate = 0.0065;
  mem->ref_tol = 1.0e-8;
  mem->ref_delta = 0.;
  mem->dut1 = 0.;
  memset(mem->q_lonlat,   0, sizeof(quat_t));
  memset(mem->q_wobble,   0, sizeof(quat_t));
  memset(mem->q_npb,      0, sizeof(quat_t));
  memset(mem->q_erot,     0, sizeof(quat_t));
  memset(mem->q_ref,      0, sizeof(quat_t));
  memset(mem->q_gal,      0, sizeof(quat_t));
  memset(mem->q_gal_inv,  0, sizeof(quat_t));
  memset(mem->beta_rot,   0, sizeof(vec3_t));
  memset(mem->beta_earth, 0, sizeof(vec3_t));
  mem->bulletinA.entries = NULL;
  mem->initialized = 1;
  return mem;
}

qp_memory_t * qp_copy_memory(qp_memory_t *memsrc) {
  qp_memory_t *memdest = malloc(sizeof(*memdest));
  *memdest = *memsrc;
  copy_iers_bulletin_a(memdest, memsrc);
  return memdest;
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

void qp_print_quat(const char *tag, quat_t q) {
  printf("%s quat: [%.6g, %.6g, %.6g, %.6g]\n",
	 tag, q[0], q[1], q[2], q[3]);
}

void qp_print_vec3(const char *tag, vec3_t v) {
  printf("%s vec3: [%.6g, %.6g, %.6g]\n",
	 tag, v[0], v[1], v[2]);
}

void qp_print_state(const char *tag, qp_state_t *state) {
  printf("%s state: rate %.6g, ctime %20.6f\n",
         tag, state->update_rate, state->ctime_last);
}

void qp_print_weather(qp_weather_t *w) {
  printf("Height [m]: %.6g\n", w->height);
  printf("Temperature [C]: %.6g\n", w->temperature);
  printf("Pressure [mbar]: %.6g\n", w->pressure);
  printf("Humidity [%%]: %.6g\n", w->humidity * 100);
  printf("Frequency [GHz]: %.6g\n", w->frequency);
  printf("Lapse rate [K/m]: %.6g\n", w->lapse_rate);
}

void qp_print_quat_mp(int th, const char *tag, quat_t q) {
  printf("[%d]  %s quat: [%.6g, %.6g, %.6g, %.6g]\n",
	 th, tag, q[0], q[1], q[2], q[3]);
}

void qp_print_vec3_mp(int th, const char *tag, vec3_t v) {
  printf("[%d]  %s vec3: [%.6g, %.6g, %.6g]\n",
	 th, tag, v[0], v[1], v[2]);
}

void qp_print_state_mp(int th, const char *tag, qp_state_t *state) {
  printf("[%d]  %s state: rate %.6g, ctime %20.6f\n",
         th, tag, state->update_rate, state->ctime_last);
}

void qp_print_weather_mp(int th, qp_weather_t *w) {
  printf("[%d]  Height [m]: %.6g\n", th, w->height);
  printf("[%d]  Temperature [C]: %.6g\n", th, w->temperature);
  printf("[%d]  Pressure [mbar]: %.6g\n", th, w->pressure);
  printf("[%d]  Humidity [%%]: %.6g\n", th, w->humidity * 100);
  printf("[%d]  Frequency [GHz]: %.6g\n", th, w->frequency);
  printf("[%d]  Lapse rate [K/m]: %.6g\n", th, w->lapse_rate);
}

void qp_print_memory(qp_memory_t *mem) {
  int thread = qp_get_opt_thread_num(mem);

  printf("[%d]  ========== QPOINT MEMORY ==========\n", thread);
  qp_print_state_mp(thread, "ref", &mem->state_ref);
  qp_print_state_mp(thread, "daber", &mem->state_daber);
  qp_print_state_mp(thread, "lonlat", &mem->state_lonlat);
  qp_print_state_mp(thread, "wobble", &mem->state_wobble);
  qp_print_state_mp(thread, "dut1", &mem->state_dut1);
  qp_print_state_mp(thread, "erot", &mem->state_erot);
  qp_print_state_mp(thread, "npb", &mem->state_npb);
  qp_print_state_mp(thread, "aaber", &mem->state_aaber);

  qp_print_weather_mp(thread, &mem->weather);

  qp_print_quat_mp(thread, "ref", mem->q_ref);
  qp_print_quat_mp(thread, "lonlat", mem->q_lonlat);
  qp_print_quat_mp(thread, "wobble", mem->q_wobble);
  qp_print_quat_mp(thread, "npb", mem->q_npb);
  qp_print_quat_mp(thread, "erot", mem->q_erot);
  qp_print_quat_mp(thread, "gal", mem->q_gal);
  qp_print_quat_mp(thread, "gal_inv", mem->q_gal_inv);
  qp_print_vec3_mp(thread, "beta earth", mem->beta_earth);
  qp_print_vec3_mp(thread, "beta rot", mem->beta_rot);

  printf("[%d]  ref_tol %.6g\n", thread, mem->ref_tol);
  printf("[%d]  ref_delta %.6g\n", thread, mem->ref_delta);
  printf("[%d]  dut1 %.6g\n", thread, mem->dut1);

  printf("[%d]  accuracy: %s\n", thread, mem->accuracy ? "low" : "full");
  printf("[%d]  mean aber: %s\n", thread, mem->mean_aber ? "yes" : "no");
  printf("[%d]  fast math: %s\n", thread, mem->fast_math ? "yes" : "no");
  printf("[%d]  polconv: %s\n", thread, mem->polconv ? "IAU" : "healpix");
  printf("[%d]  pair dets: %s\n", thread, mem->pair_dets ? "yes" : "no");
  printf("[%d]  fast pix: %s\n", thread, mem->fast_pix ? "yes" : "no");
  printf("[%d]  gal init: %s\n", thread, mem->gal_init ? "yes" : "no");

  printf("[%d]  num threads: %d\n", thread, qp_get_opt_num_threads(mem));
  printf("[%d]  thread num: %d\n", thread, qp_get_opt_thread_num(mem));
  printf("[%d]  initialized: %s\n", thread, mem->initialized ? "yes" : "no");

  printf("[%d]  ===================================\n", thread);
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
OPTIONFUNCD(fast_pix)

void qp_set_options(qp_memory_t *mem,
		    int accuracy,
		    int mean_aber,
		    int fast_math,
		    int polconv,
		    int pair_dets,
		    int pix_order,
                    int fast_pix,
		    int num_threads) {
  qp_set_opt_accuracy   (mem, accuracy);
  qp_set_opt_mean_aber  (mem, mean_aber);
  qp_set_opt_fast_math  (mem, fast_math);
  qp_set_opt_polconv    (mem, polconv);
  qp_set_opt_pair_dets  (mem, pair_dets);
  qp_set_opt_pix_order  (mem, pix_order);
  qp_set_opt_fast_pix   (mem, fast_pix);
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
