#pragma once

#include <array>
#include <cstddef>
#include <string_view>

namespace qp {

// Sentinel update rates, matching QP_DO_* in qp_params.c.
inline constexpr double kDoAlways = 0.;
inline constexpr double kDoOnce = -1.;
inline constexpr double kDoNever = -999.;

// One correction's update-rate cache. check() mutates ctime_last, so the
// per-sample loop over these is inherently sequential -- which is why
// vectorization stays serial in the bindings.
class UpdateState {
 public:
  constexpr explicit UpdateState(double rate = kDoAlways) : rate_(rate) {}

  constexpr double rate() const { return rate_; }

  void set_rate(double rate) {
    if (rate != rate_) {
      rate_ = rate;
      ctime_last_ = -1;
    }
  }

  void reset() { ctime_last_ = -1; }

  // Returns true if the correction should be recomputed at this sample.
  bool check(double ctime) {
    if (rate_ == kDoNever) return false;
    if (rate_ == kDoOnce && ctime_last_ > 0) return false;
    if (ctime_last_ <= 0) {
      ctime_last_ = ctime;
      return true;
    }
    // recompute on a backwards time jump
    if (ctime < ctime_last_) {
      ctime_last_ = ctime;
      return true;
    }
    if ((ctime - ctime_last_) >= rate_) {
      ctime_last_ = ctime;
      return true;
    }
    return false;
  }

  constexpr bool should_apply() const { return rate_ != kDoNever; }

 private:
  double rate_;
  double ctime_last_ = -1;
};

enum class Rate {
  daber = 0,
  lonlat,
  wobble,
  dut1,
  erot,
  npb,
  aaber,
  ref,
  COUNT,
};

inline constexpr std::size_t kNumRates = static_cast<std::size_t>(Rate::COUNT);

// Name and default rate per correction, in Rate order -- one table rather
// than two to keep in alignment. Defaults are from qp_init_memory; the
// forward and inverse sets are identical.
struct RateDesc {
  std::string_view name;
  double def;
};

inline constexpr std::array<RateDesc, kNumRates> kRates = {{
    {"daber", kDoAlways},
    {"lonlat", kDoAlways},
    {"wobble", kDoNever},
    {"dut1", kDoNever},
    {"erot", kDoAlways},
    {"npb", 10.},
    {"aaber", 100.},
    {"ref", kDoNever},
}};

// The initial cache state, so Pointing needs no constructor of its own.
inline constexpr std::array<UpdateState, kNumRates> kInitialRateStates = {{
    UpdateState(kRates[0].def), UpdateState(kRates[1].def),
    UpdateState(kRates[2].def), UpdateState(kRates[3].def),
    UpdateState(kRates[4].def), UpdateState(kRates[5].def),
    UpdateState(kRates[6].def), UpdateState(kRates[7].def),
}};

struct Weather {
  double temperature = 0.;
  double pressure = 10.;
  double humidity = 0.;
  double frequency = 150.;
};

struct Options {
  int accuracy = 0;
  int mean_aber = 1;
  int fast_aber = 1;
  int fast_math = 0;
  int polconv = 0;
  int pix_order = 0;
  int interp_pix = 0;
  // On by default here, unlike the C. It takes the pixel from the pointing
  // vector instead of going round through ra/dec, which is 26-34% of
  // tod2map, and in qpoint2 it costs nothing: the polarization angles are
  // bit-identical to the slow path, and the pixel differs only inside
  // 2.1e-8 rad of a pole, where the fast path is the correct one.
  int fast_pix = 1;
  int error_missing = 1;
  int nan_missing = 0;
  int interp_missing = 0;
};

// Descriptor tables replacing the RATEFUNC / OPTIONFUNC / WEATHFUNC /
// DOUBLEFUNC macro families in qp_params.c, which between them generate
// around a hundred near-identical functions. Setting a value resets the
// associated forward rate state, and only when the value actually changed;
// Rate::COUNT means there is nothing to reset.
struct OptionDesc {
  std::string_view name;
  int Options::*mem;
  Rate reset = Rate::COUNT;
};

inline constexpr OptionDesc kOptionParams[] = {
    {"accuracy", &Options::accuracy, Rate::npb},
    {"mean_aber", &Options::mean_aber, Rate::aaber},
    {"fast_aber", &Options::fast_aber},
    {"fast_math", &Options::fast_math},
    {"polconv", &Options::polconv},
    {"pix_order", &Options::pix_order},
    {"interp_pix", &Options::interp_pix},
    {"fast_pix", &Options::fast_pix},
    {"error_missing", &Options::error_missing},
    {"nan_missing", &Options::nan_missing},
    {"interp_missing", &Options::interp_missing},
};

// Every weather parameter resets the refraction rate, so there is nothing
// to tabulate but the name and the member.
struct WeatherDesc {
  std::string_view name;
  double Weather::*mem;
};

inline constexpr WeatherDesc kWeatherParams[] = {
    {"temperature", &Weather::temperature},
    {"pressure", &Weather::pressure},
    {"humidity", &Weather::humidity},
    {"frequency", &Weather::frequency},
};

}  // namespace qp
