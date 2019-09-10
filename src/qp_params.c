#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "qpoint.h"

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
  qp_init_state(&mem->state_daber_inv , QP_DO_ALWAYS);
  qp_init_state(&mem->state_lonlat_inv, QP_DO_ALWAYS);
  qp_init_state(&mem->state_wobble_inv, QP_DO_NEVER);
  qp_init_state(&mem->state_dut1_inv  , QP_DO_NEVER);
  qp_init_state(&mem->state_erot_inv  , QP_DO_ALWAYS);
  qp_init_state(&mem->state_npb_inv   , 10);
  qp_init_state(&mem->state_aaber_inv , 100);
  qp_init_state(&mem->state_ref_inv   , QP_DO_NEVER);
  mem->accuracy = 0;
  mem->mean_aber = 1;
  mem->fast_math = 0;
  mem->polconv = 0;
  mem->pix_order = 0;
  mem->interp_pix = 0;
  mem->fast_pix = 0;
  mem->error_missing = 1;
  mem->nan_missing = 0;
  mem->interp_missing = 0;
  mem->gal_init = 0;
  mem->dipole_init = 0;
  mem->thread_num = 0;
#ifndef ENABLE_LITE
  qp_set_opt_num_threads(mem, 0);
#endif
  mem->weather.temperature = 0.;
  mem->weather.pressure = 10.;
  mem->weather.humidity = 0.;
  mem->weather.frequency = 150.;
  mem->ref_delta = 0.;
  mem->dut1 = 0.;
  memset(mem->q_lonlat,   0, sizeof(quat_t));
  memset(mem->q_wobble,   0, sizeof(quat_t));
  memset(mem->q_npb,      0, sizeof(quat_t));
  memset(mem->q_erot,     0, sizeof(quat_t));
  memset(mem->q_ref,      0, sizeof(quat_t));
  memset(mem->q_lonlat_inv,   0, sizeof(quat_t));
  memset(mem->q_wobble_inv,   0, sizeof(quat_t));
  memset(mem->q_npb_inv,      0, sizeof(quat_t));
  memset(mem->q_erot_inv,     0, sizeof(quat_t));
  memset(mem->q_ref_inv,      0, sizeof(quat_t));
  memset(mem->q_gal,      0, sizeof(quat_t));
  memset(mem->q_gal_inv,  0, sizeof(quat_t));
  memset(mem->v_dipole,   0, sizeof(vec3_t));
  memset(mem->beta_rot,   0, sizeof(vec3_t));
  memset(mem->beta_earth, 0, sizeof(vec3_t));
  mem->bulletinA.entries = NULL;
  mem->error_code = 0;
  mem->error_string = NULL;
  mem->init = 1;
  return mem;
}

qp_memory_t * qp_copy_memory(qp_memory_t *memsrc) {
  qp_memory_t *memdest = malloc(sizeof(*memdest));
  *memdest = *memsrc;
  qp_copy_iers_bulletin_a(memdest, memsrc);
  return memdest;
}

void qp_free_memory(qp_memory_t *mem) {
  qp_set_iers_bulletin_a(mem, 0, 0, NULL, NULL, NULL);
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
RATEFUNCD(daber_inv)
RATEFUNCD(lonlat_inv)
RATEFUNCD(wobble_inv)
RATEFUNCD(dut1_inv)
RATEFUNCD(erot_inv)
RATEFUNCD(npb_inv)
RATEFUNCD(aaber_inv)
RATEFUNCD(ref_inv)

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
  printf("Temperature [C]: %.6g\n", w->temperature);
  printf("Pressure [mbar]: %.6g\n", w->pressure);
  printf("Humidity [%%]: %.6g\n", w->humidity * 100);
  printf("Frequency [GHz]: %.6g\n", w->frequency);
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
  printf("[%d]  %s state: rate %.6g, ctime %-20.6f\n",
         th, tag, state->update_rate, state->ctime_last);
}

void qp_print_weather_mp(int th, qp_weather_t *w) {
  printf("[%d]  weather: temperature [C]: %.6g\n", th, w->temperature);
  printf("[%d]  weather: pressure [mbar]: %.6g\n", th, w->pressure);
  printf("[%d]  weather: humidity [%%]: %.6g\n", th, w->humidity * 100);
  printf("[%d]  weather: frequency [GHz]: %.6g\n", th, w->frequency);
}

void qp_print_memory(qp_memory_t *mem) {
#ifndef ENABLE_LITE
  int thread = qp_get_opt_thread_num(mem);
#else
  int thread = mem->thread_num;
#endif

  printf("[%d]  ========== QPOINT MEMORY ==========\n", thread);
  qp_print_state_mp(thread, "ref", &mem->state_ref);
  qp_print_quat_mp(thread, "ref", mem->q_ref);
  printf("[%d]  ref delta %.6g\n", thread, mem->ref_delta);

  qp_print_state_mp(thread, "daber", &mem->state_daber);
  qp_print_state_mp(thread, "daber inv", &mem->state_daber_inv);
  qp_print_vec3_mp(thread, "daber beta rot", mem->beta_rot);

  qp_print_state_mp(thread, "lonlat", &mem->state_lonlat);
  qp_print_quat_mp(thread, "lonlat", mem->q_lonlat);
  qp_print_state_mp(thread, "lonlat inv", &mem->state_lonlat_inv);
  qp_print_quat_mp(thread, "lonlat inv", mem->q_lonlat_inv);

  qp_print_state_mp(thread, "wobble", &mem->state_wobble);
  qp_print_quat_mp(thread, "wobble", mem->q_wobble);
  qp_print_state_mp(thread, "wobble inv", &mem->state_wobble_inv);
  qp_print_quat_mp(thread, "wobble inv", mem->q_wobble_inv);

  qp_print_state_mp(thread, "dut1", &mem->state_dut1);
  qp_print_state_mp(thread, "dut1 inv", &mem->state_dut1_inv);
  printf("[%d]  dut1 %.6g\n", thread, mem->dut1);

  qp_print_state_mp(thread, "erot", &mem->state_erot);
  qp_print_quat_mp(thread, "erot", mem->q_erot);
  qp_print_state_mp(thread, "erot inv", &mem->state_erot_inv);
  qp_print_quat_mp(thread, "erot inv", mem->q_erot_inv);

  qp_print_state_mp(thread, "npb", &mem->state_npb);
  qp_print_quat_mp(thread, "npb", mem->q_npb);
  qp_print_state_mp(thread, "npb inv", &mem->state_npb_inv);
  qp_print_quat_mp(thread, "npb inv", mem->q_npb_inv);

  qp_print_state_mp(thread, "aaber", &mem->state_aaber);
  qp_print_state_mp(thread, "aaber inv", &mem->state_aaber_inv);
  qp_print_vec3_mp(thread, "aaber beta earth", mem->beta_earth);

  printf("[%d]  gal init: %s\n", thread, mem->gal_init ? "yes" : "no");
  qp_print_quat_mp(thread, "gal", mem->q_gal);
  qp_print_quat_mp(thread, "gal inv", mem->q_gal_inv);

  // qp_print_weather_mp(thread, &mem->weather);

  printf("[%d]  opt: accuracy: %s\n", thread, mem->accuracy ? "low" : "full");
  printf("[%d]  opt: mean aber: %s\n", thread, mem->mean_aber ? "yes" : "no");
  printf("[%d]  opt: fast math: %s\n", thread, mem->fast_math ? "yes" : "no");
  printf("[%d]  opt: polconv: %s\n", thread, mem->polconv ? "IAU" : "healpix");
  printf("[%d]  opt: interp pix: %s\n", thread, mem->interp_pix ? "yes" : "no");
  printf("[%d]  opt: fast pix: %s\n", thread, mem->fast_pix ? "yes" : "no");
  printf("[%d]  opt: error missing: %s\n", thread, mem->error_missing ? "yes" : "no");
  printf("[%d]  opt: nan missing: %s\n", thread, mem->nan_missing ? "yes" : "no");
  printf("[%d]  opt: interp missing: %s\n", thread, mem->interp_missing ? "yes" : "no");

#ifndef ENABLE_LITE
  printf("[%d]  opt: num threads: %d\n", thread, qp_get_opt_num_threads(mem));
  printf("[%d]  thread num: %d\n", thread, qp_get_opt_thread_num(mem));
#endif
  printf("[%d]  initialized: %s\n", thread, mem->init ? "yes" : "no");

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

void qp_set_inv_rates(qp_memory_t *mem,
		      double daber_rate,
		      double lonlat_rate,
		      double wobble_rate,
		      double dut1_rate,
		      double erot_rate,
		      double npb_rate,
		      double aaber_rate,
		      double ref_rate) {
  qp_set_rate_daber_inv (mem, daber_rate);
  qp_set_rate_lonlat_inv(mem, lonlat_rate);
  qp_set_rate_wobble_inv(mem, wobble_rate);
  qp_set_rate_dut1_inv  (mem, dut1_rate);
  qp_set_rate_erot_inv  (mem, erot_rate);
  qp_set_rate_npb_inv   (mem, npb_rate);
  qp_set_rate_aaber_inv (mem, aaber_rate);
  qp_set_rate_ref_inv   (mem, ref_rate);
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

void qp_reset_inv_rates(qp_memory_t *mem) {
  qp_reset_rate_daber_inv (mem);
  qp_reset_rate_lonlat_inv(mem);
  qp_reset_rate_wobble_inv(mem);
  qp_reset_rate_dut1_inv  (mem);
  qp_reset_rate_erot_inv  (mem);
  qp_reset_rate_npb_inv   (mem);
  qp_reset_rate_aaber_inv (mem);
  qp_reset_rate_ref_inv   (mem);
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
OPTIONFUNCD(pix_order)
OPTIONFUNCD(interp_pix)
OPTIONFUNCD(fast_pix)
OPTIONFUNCD(error_missing)
OPTIONFUNCD(nan_missing)
OPTIONFUNCD(interp_missing)

void qp_set_options(qp_memory_t *mem,
		    int accuracy,
		    int mean_aber,
		    int fast_math,
		    int polconv,
		    int pix_order,
                    int interp_pix,
                    int fast_pix,
                    int error_missing,
                    int nan_missing,
                    int interp_missing,
		    int num_threads) {
  qp_set_opt_accuracy      (mem, accuracy);
  qp_set_opt_mean_aber     (mem, mean_aber);
  qp_set_opt_fast_math     (mem, fast_math);
  qp_set_opt_polconv       (mem, polconv);
  qp_set_opt_pix_order     (mem, pix_order);
  qp_set_opt_interp_pix    (mem, interp_pix);
  qp_set_opt_fast_pix      (mem, fast_pix);
  qp_set_opt_error_missing (mem, error_missing);
  qp_set_opt_nan_missing   (mem, nan_missing);
  qp_set_opt_interp_missing(mem, interp_missing);
#ifndef ENABLE_LITE
  qp_set_opt_num_threads   (mem, num_threads);
#endif
}

// update all ref_data parameters
void qp_set_weather(qp_memory_t *mem, double temperature, double pressure,
		    double humidity, double frequency) {
  qp_set_weather_temperature(mem, temperature);
  qp_set_weather_pressure   (mem, pressure);
  qp_set_weather_humidity   (mem, humidity);
  qp_set_weather_frequency  (mem, frequency);
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
WEATHFUNCD(temperature)
WEATHFUNCD(pressure)
WEATHFUNCD(humidity)
WEATHFUNCD(frequency)

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

DOUBLEFUNCD(ref_delta)
DOUBLEFUNCDRR(dut1, erot, wobble)
