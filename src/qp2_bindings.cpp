#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <cctype>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <memory>
#include <mutex>
#include <new>
#include <optional>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "cpp/cshims.hpp"
#include "cpp/error.hpp"
#include "cpp/map.hpp"
#include "cpp/pixinfo.hpp"
#include "cpp/pointing.hpp"
#include "cpp/quat.hpp"
#include "cpp/span.hpp"
#include "cpp/state.hpp"

namespace py = pybind11;
using namespace qp;

namespace {

// ---------------------------------------------------------------------------
// Argument handling
//
// Every per-sample argument is either a scalar or a full-length 1-D array;
// there is no true N-dimensional broadcasting anywhere in this API. That makes
// the whole thing resolvable here, with no Python-side preparation: a scalar
// becomes a stride-0 column, so it costs one read per sample and no memory.
// The alternative -- np.broadcast_arrays plus a contiguity copy -- allocates a
// full-length array for every scalar argument on every call.
//
// An existing array is never converted. pybind11's usual options all copy
// silently: forcecast copies on any dtype or layout mismatch, and even a
// plain py::array_t<double> routes through PyArray_FromAny, so a "safe" cast
// like float32 -> float64 also copies. Declaring the parameters as
// py::object and checking here gives .noconvert() semantics -- raise, never
// convert -- with errors that name the argument and say what was wrong.
//
// Two things are exempt, because neither is a silent copy of a buffer the
// caller already owns:
//   - scalars, where widening one value is exact and allocates nothing;
//   - lists and tuples, which have no buffer to share, so there is nothing
//     to alias and materializing one is the only option. Convenient in
//     tests, and the cost is bounded by what was literally written out.
// ---------------------------------------------------------------------------

std::string shape_str(const py::array &a) {
  std::string s = "(";
  for (py::ssize_t i = 0; i < a.ndim(); ++i) {
    if (i) s += ", ";
    s += std::to_string(a.shape(i));
  }
  return s + ")";
}

void require_double_array(const py::array &a, const char *name) {
  if (!a.dtype().is(py::dtype::of<double>()))
    throw py::type_error(std::string(name) + " must be a float64 array, got " +
                         py::str(a.dtype()).cast<std::string>() +
                         " (convert it yourself; this API never copies an "
                         "array you already own)");
  if (!(a.flags() & py::array::c_style))
    throw py::value_error(std::string(name) + " must be C-contiguous");
  if (!(a.flags() & py::detail::npy_api::NPY_ARRAY_ALIGNED_))
    throw py::value_error(std::string(name) + " must be aligned");
}

bool is_convertible_sequence(const py::handle &o) {
  return py::isinstance<py::sequence>(o) && !py::isinstance<py::str>(o) &&
         !py::isinstance<py::bytes>(o);
}

// Array to read from, for something that is array-like. A numpy array is
// used in place after validation; a list or tuple is materialized.
py::array as_array(const py::handle &o, const char *name) {
  if (py::isinstance<py::array>(o)) {
    auto a = py::reinterpret_borrow<py::array>(o);
    if (a.ndim() != 0) require_double_array(a, name);
    return a;
  }
  auto conv =
      py::array_t<double, py::array::c_style | py::array::forcecast>::ensure(o);
  if (!conv)
    throw py::type_error(std::string(name) +
                         " could not be converted to a float64 array");
  return conv;
}

// A parsed argument, before the sample count is known.
template <class T>
struct RawCol {
  const T *ptr = nullptr;
  py::ssize_t len = 0;  // 0 when there is no array behind it
  T value{};
  bool present = false;
  // Whether the argument arrived with a sample axis -- a 1-d array of
  // doubles, or an (n, 4) of quaternions. A scalar, a 0-d array and a bare
  // (4,) quaternion are all one sample without an axis, and the outputs
  // follow: those degrade, an array does not, even at length 1.
  bool vector = false;
  // Owns the array `ptr` points into. Always set, so the buffer's lifetime
  // does not depend on the caller keeping the argument alive -- which it
  // would not, for an array we materialized from a list.
  py::object keep;
};

RawCol<double> parse_col(const py::handle &o, const char *name) {
  RawCol<double> r;
  if (o.is_none()) return r;
  r.present = true;

  if (!py::isinstance<py::array>(o) && !is_convertible_sequence(o)) {
    try {
      r.value = o.cast<double>();
    } catch (const py::cast_error &) {
      throw py::type_error(std::string(name) +
                           " must be a float64 array, a sequence, or a scalar");
    }
    return r;
  }

  py::array a = as_array(o, name);
  if (a.ndim() == 0) {
    // 0-d array: a scalar, so dtype is not load-bearing
    r.value = a.attr("item")().cast<double>();
    return r;
  }
  if (a.ndim() != 1)
    throw py::value_error(std::string(name) +
                          " must be 1-dimensional or scalar, got " +
                          shape_str(a));
  r.ptr = static_cast<const double *>(a.data());
  r.len = a.shape(0);
  r.vector = true;
  r.keep = std::move(a);
  return r;
}

RawCol<Quat> parse_quat_col(const py::handle &o, const char *name) {
  RawCol<Quat> r;
  if (o.is_none()) return r;
  r.present = true;

  if (!py::isinstance<py::array>(o) && !is_convertible_sequence(o))
    throw py::type_error(std::string(name) +
                         " must be a float64 array or sequence of shape "
                         "(4,) or (n, 4)");
  py::array a = as_array(o, name);

  // Quat is standard-layout and four contiguous doubles, and alignment was
  // just checked, so an (n,4) buffer can be viewed as Quats directly.
  if ((a.ndim() == 1 && a.shape(0) == 4) ||
      (a.ndim() == 2 && a.shape(1) == 4)) {
    r.ptr = reinterpret_cast<const Quat *>(a.data());
    r.len = (a.ndim() == 1) ? 1 : a.shape(0);
    r.vector = a.ndim() == 2;
    r.keep = std::move(a);
    return r;
  }
  throw py::value_error(std::string(name) +
                        " must have shape (4,) or (n, 4), got " +
                        shape_str(a));
}

// Resolve against the sample count: full length reads normally, a single
// element repeats via stride 0.
template <class T>
Col<T> fix(const RawCol<T> &r, py::ssize_t n, const char *name) {
  Col<T> c;
  c.present = r.present;
  c.value = r.value;
  if (!r.ptr) return c;
  if (r.len == n)
    c.stride = 1;
  else if (r.len == 1)
    c.stride = 0;
  else
    throw py::value_error(std::string(name) + " has length " +
                          std::to_string(r.len) + ", expected " +
                          std::to_string(n) + " or 1");
  c.ptr = r.ptr;
  return c;
}

Quat load_quat(const py::handle &o, const char *name) {
  if (!py::isinstance<py::array>(o) && !is_convertible_sequence(o))
    throw py::type_error(std::string(name) +
                         " must be a float64 array or sequence of 4 elements");
  py::array a = as_array(o, name);
  if (a.size() != 4)
    throw py::value_error(std::string(name) + " must have 4 elements, got " +
                          shape_str(a));
  return Quat::load(static_cast<const double *>(a.data()));
}

// An array written in place. No stride-0 broadcasting here: a repeated
// column cannot be an output, and the length must match exactly.
//
// The n-less overload takes the length from the array itself, for the
// in-place entry points where that array is the only thing that sets the
// sample count -- parsing it a second time just to count it was pure
// duplicate validation.
template <class T>
Span<T> out_span(const py::handle &o, const char *name) {
  if (!py::isinstance<py::array>(o))
    throw py::type_error(std::string(name) +
                         " must be a float64 array; it is written in place");
  auto a = py::reinterpret_borrow<py::array>(o);
  require_double_array(a, name);
  if (!a.writeable())
    throw py::value_error(std::string(name) +
                          " must be writeable; it is written in place");
  const py::ssize_t want = static_cast<py::ssize_t>(sizeof(T) / sizeof(double));
  // For a quaternion span, the same shapes parse_quat_col accepts -- a
  // count that merely divides evenly is not enough.
  const bool shape_ok =
      want == 1 || (a.ndim() == 1 && a.shape(0) == want) ||
      (a.ndim() == 2 && a.shape(1) == want);
  if (!shape_ok)
    throw py::value_error(std::string(name) + " must have shape (" +
                          std::to_string(want) + ",) or (n, " +
                          std::to_string(want) + "), got " + shape_str(a));
  return {reinterpret_cast<T *>(a.mutable_data()),
          static_cast<size_t>(a.size() / want)};
}

template <class T>
Span<T> out_span(const py::handle &o, const char *name, py::ssize_t n) {
  if (!py::isinstance<py::array>(o))
    throw py::type_error(std::string(name) +
                         " must be a float64 array; it is written in place");
  auto a = py::reinterpret_borrow<py::array>(o);
  require_double_array(a, name);
  if (!a.writeable())
    throw py::value_error(std::string(name) +
                          " must be writeable; it is written in place");
  const py::ssize_t want = static_cast<py::ssize_t>(sizeof(T) / sizeof(double));
  const py::ssize_t have = a.size() / want;
  if (have != n || a.size() % want != 0)
    throw py::value_error(std::string(name) + " has length " +
                          std::to_string(have) + ", expected " +
                          std::to_string(n));
  return {reinterpret_cast<T *>(a.mutable_data()), static_cast<size_t>(n)};
}

struct MapArray {
  double *ptr = nullptr;
  std::size_t nrow = 0;
  std::size_t npix = 0;
  py::object keep;
};

MapArray parse_map(const py::handle &o, const char *name, bool need_write) {
  MapArray r;
  if (o.is_none() || (py::isinstance<py::bool_>(o) && !o.cast<bool>()))
    return r;
  py::array a = as_array(o, name);
  if (a.ndim() != 2)
    throw py::value_error(std::string(name) +
                          " must have shape (nrow, npix), got " + shape_str(a));
  if (need_write && !(a.flags() & py::detail::npy_api::NPY_ARRAY_WRITEABLE_))
    throw py::value_error(std::string(name) +
                          " must be writeable; it is accumulated into");
  r.nrow = static_cast<std::size_t>(a.shape(0));
  r.npix = static_cast<std::size_t>(a.shape(1));
  // tod2map only reads the source map, so a read-only buffer is fine there.
  r.ptr = static_cast<double *>(
      const_cast<void *>(static_cast<const void *>(a.data())));
  r.keep = std::move(a);
  return r;
}

// Output conventions: a call made entirely of scalars degrades to a scalar,
// and to a bare (4,) for a quaternion. A call made with arrays keeps its
// axis, even when the array holds one sample.
//
// The distinction cannot be made from the output, which is length 1 either
// way -- it has to come from how the arguments arrived, which is what
// SampleCount carries. qpoint decides it from the output size and so cannot
// tell `f(1.0)` from `f(np.array([1.0]))`; this is a deliberate divergence.
//
// py::cast takes the Python type from T, so a double still arrives as a
// float and a long as an int.
template <class T>
py::object maybe_scalar(py::array_t<T> a, bool degrade) {
  if (degrade && a.size() == 1) return py::cast(*a.data());
  return std::move(a);
}

// Not the same operation: this degrades a lone (1, 4) to a bare (4,)
// rather than to a scalar.
py::object maybe_scalar_quat(py::array_t<double> q, bool degrade) {
  if (degrade && q.shape(0) == 1) return q[py::int_(0)];
  return std::move(q);
}

// ---------------------------------------------------------------------------
// The two halves of every vectorized entry point
//
// Arguments: a ColSet collects the per-sample columns, then resolves them
// against the sample count in one step. Columns are handed out as
// references and filled by resolve(), so add them all before the loop.
//
// Outputs: vec_out/long_out/quat_out allocate, take the raw pointers while
// the GIL is still held, drop it, run the loop, and apply the single-sample
// scalar degradation on the way out. Those four steps have to happen in
// that order, which is exactly why they belong in one place.
// ---------------------------------------------------------------------------

// What resolve() reports: the sample count, and whether every argument
// arrived as a scalar. The second is what decides whether the outputs keep
// their axis, and only the arguments can answer it -- see maybe_scalar.
//
// Converts to ssize_t so the loops and length checks that only want the
// count read naturally.
struct SampleCount {
  py::ssize_t n = 1;
  bool all_scalar = true;

  constexpr operator py::ssize_t() const { return n; }
};

class ColSet {
 public:
  Col<double> &add(const py::handle &o, const char *name) {
    dbl_.push_back({parse_col(o, name), {}, name});
    return dbl_.back().col;
  }

  // required rejects None: some arguments have no sensible zero.
  Col<Quat> &add_quat(const py::handle &o, const char *name,
                      bool required = false) {
    quat_.push_back({parse_quat_col(o, name), {}, name});
    if (required && !quat_.back().raw.present)
      throw py::value_error(std::string(name) + " is required");
    return quat_.back().col;
  }

  // The sample count is the longest array argument; all-scalar means one
  // sample. Every column then reads at stride 1 or, if it is a lone
  // element, stride 0.
  SampleCount resolve() {
    py::ssize_t n = 1;
    bool all_scalar = true;
    for (const auto &e : dbl_) {
      if (e.raw.len > n) n = e.raw.len;
      if (e.raw.vector) all_scalar = false;
    }
    for (const auto &e : quat_) {
      if (e.raw.len > n) n = e.raw.len;
      if (e.raw.vector) all_scalar = false;
    }
    for (auto &e : dbl_) e.col = fix(e.raw, n, e.name);
    for (auto &e : quat_) e.col = fix(e.raw, n, e.name);
    return {n, all_scalar};
  }

 private:
  template <class T>
  struct Entry {
    RawCol<T> raw;
    Col<T> col;
    const char *name;
  };
  // deque, not vector: add() hands out references that must survive later
  // additions.
  std::deque<Entry<double>> dbl_;
  std::deque<Entry<Quat>> quat_;
};

template <std::size_t K, class Fn>
py::object vec_out(SampleCount count, Fn &&body) {
  const py::ssize_t n = count.n;
  const bool degrade = count.all_scalar;
  std::array<py::array_t<double>, K> arr;
  std::array<double *, K> op{};
  for (std::size_t k = 0; k < K; ++k) {
    arr[k] = py::array_t<double>(n);
    op[k] = arr[k].mutable_data();
  }
  {
    py::gil_scoped_release nogil;
    for (py::ssize_t i = 0; i < n; ++i) body(i, op);
  }
  if constexpr (K == 1) return maybe_scalar(std::move(arr[0]), degrade);
  py::tuple out(K);
  for (std::size_t k = 0; k < K; ++k)
    out[k] = maybe_scalar(std::move(arr[k]), degrade);
  return std::move(out);
}

template <class Fn>
py::object long_out(SampleCount count, Fn &&body) {
  const py::ssize_t n = count.n;
  py::array_t<long> arr(n);
  long *op = arr.mutable_data();
  {
    py::gil_scoped_release nogil;
    for (py::ssize_t i = 0; i < n; ++i) op[i] = body(i);
  }
  return maybe_scalar(std::move(arr), count.all_scalar);
}

// The sky-coordinate output shapes. Both end in the polarization pair --
// p1 always, p2 only when sin2psi/cos2psi were asked for -- and differ
// only in what they put in front of it. p2 reaches the body as nullptr in
// PA mode, which is what quat2radec and quat2pix already expect.
template <class Fn>
py::object radec_out(SampleCount count, bool return_pa, Fn &&body) {
  const py::ssize_t n = count.n;
  const bool degrade = count.all_scalar;
  py::array_t<double> ra(n), dec(n), p1(n), p2(return_pa ? 0 : n);
  double *ora = ra.mutable_data(), *odec = dec.mutable_data();
  double *op1 = p1.mutable_data();
  double *op2 = return_pa ? nullptr : p2.mutable_data();
  {
    py::gil_scoped_release nogil;
    for (py::ssize_t i = 0; i < n; ++i) body(i, ora, odec, op1, op2);
  }
  if (return_pa)
    return py::make_tuple(maybe_scalar(std::move(ra), degrade),
                          maybe_scalar(std::move(dec), degrade),
                          maybe_scalar(std::move(p1), degrade));
  return py::make_tuple(maybe_scalar(std::move(ra), degrade),
                        maybe_scalar(std::move(dec), degrade),
                        maybe_scalar(std::move(p1), degrade),
                        maybe_scalar(std::move(p2), degrade));
}

template <class Fn>
py::object pix_out(SampleCount count, bool return_pa, Fn &&body) {
  const py::ssize_t n = count.n;
  const bool degrade = count.all_scalar;
  py::array_t<long> pix(n);
  py::array_t<double> p1(n), p2(return_pa ? 0 : n);
  long *opix = pix.mutable_data();
  double *op1 = p1.mutable_data();
  double *op2 = return_pa ? nullptr : p2.mutable_data();
  {
    py::gil_scoped_release nogil;
    for (py::ssize_t i = 0; i < n; ++i) body(i, opix, op1, op2);
  }
  if (return_pa)
    return py::make_tuple(maybe_scalar(std::move(pix), degrade),
                          maybe_scalar(std::move(p1), degrade));
  return py::make_tuple(maybe_scalar(std::move(pix), degrade),
                        maybe_scalar(std::move(p1), degrade),
                        maybe_scalar(std::move(p2), degrade));
}

// degrade: whether a single sample comes back as a bare (4,). Boresight
// quaternions keep their leading axis; det_offset and friends do not.
template <class Fn>
py::object quat_out(SampleCount count, bool degrade, Fn &&body) {
  const py::ssize_t n = count.n;
  py::array_t<double> arr({n, py::ssize_t(4)});
  double *op = arr.mutable_data();
  {
    py::gil_scoped_release nogil;
    for (py::ssize_t i = 0; i < n; ++i) body(i).store(op + 4 * i);
  }
  // degrade is the entry point's own choice -- azel2bore always returns
  // (n, 4) -- and applies only to a call made of scalars.
  return maybe_scalar_quat(std::move(arr), degrade && count.all_scalar);
}

// ---------------------------------------------------------------------------
// Pointing wrapper
//
// num_threads is configuration for the (Phase 2) OpenMP driver in this layer,
// so it lives here rather than on the core -- unlike qp_memory_t, which
// carries both num_threads and a per-thread thread_num.
// ---------------------------------------------------------------------------

// 0 means "use every available thread". qp_set_opt_num_threads resolves
// this by opening a parallel region purely to count threads; querying the
// max directly is equivalent and does not spawn one.
int resolve_num_threads(int n) {
  if (n != 0) return n;
#ifdef _OPENMP
  return omp_get_max_threads();
#else
  return 1;
#endif
}

class PointingWrap {
 public:
  Pointing core;
  int num_threads = resolve_num_threads(0);
};

std::string lower(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(),
                 [](unsigned char c) { return std::tolower(c); });
  return s;
}

template <class T, std::size_t N>
const T *find_desc(const T (&table)[N], const std::string &key) {
  for (const auto &d : table)
    if (d.name == key) return &d;
  return nullptr;
}

// The three options the API presents as strings, in one table so that
// set_param and get_param cannot disagree about the spellings. An
// unrecognized string reads as the off value, which is what the C does.
struct EnumOpt {
  std::string_view name, off, on, alias;
};

constexpr EnumOpt kEnumOpts[] = {
    {"accuracy", "high", "low", ""},
    {"polconv", "cosmo", "iau", ""},
    {"pix_order", "ring", "nest", "nested"},
};

const EnumOpt *find_enum(const std::string &key) {
  for (const auto &d : kEnumOpts)
    if (d.name == key) return &d;
  return nullptr;
}

void set_param(PointingWrap &w, const std::string &key, py::object val) {
  if (key.size() > 5 && key.compare(0, 5, "rate_") == 0) {
    double rate;
    if (py::isinstance<py::str>(val)) {
      const auto s = val.cast<std::string>();
      if (s == "always")
        rate = kDoAlways;
      else if (s == "once")
        rate = kDoOnce;
      else if (s == "never")
        rate = kDoNever;
      else
        throw py::key_error("Unknown rate state: " + s);
    } else {
      rate = val.cast<double>();
    }
    w.core.set_rate(key.substr(5), rate);
    return;
  }

  // num_threads drives the OpenMP driver in this layer, not the core.
  if (key == "num_threads") {
    int ival = 0;
    if (!val.is_none() && !py::isinstance<py::str>(val)) ival = val.cast<int>();
    w.num_threads = resolve_num_threads(ival);
    return;
  }

  if (find_desc(kOptionParams, key)) {
    int ival = 0;
    if (!val.is_none()) {
      if (py::isinstance<py::str>(val)) {
        const auto s = lower(val.cast<std::string>());
        const EnumOpt *e = find_enum(key);
        if (e && (s == e->on || (!e->alias.empty() && s == e->alias))) ival = 1;
      } else {
        ival = val.cast<int>();
      }
    }
    w.core.set_opt(key, ival);
    return;
  }

  if (find_desc(kWeatherParams, key)) {
    w.core.set_weather_param(key, val.cast<double>());
    return;
  }
  if (key == "ref_delta" || key == "dut1") {
    w.core.set_double(key, val.cast<double>());
    return;
  }
  throw py::key_error("Unknown parameter: " + key);
}

py::object get_param(const PointingWrap &w, const std::string &key) {
  if (key.size() > 5 && key.compare(0, 5, "rate_") == 0) {
    const double rate = w.core.get_rate(key.substr(5));
    if (rate == kDoAlways) return py::cast(std::string("always"));
    if (rate == kDoOnce) return py::cast(std::string("once"));
    if (rate == kDoNever) return py::cast(std::string("never"));
    return py::cast(rate);
  }
  if (key == "num_threads") return py::cast(w.num_threads);
  if (find_desc(kOptionParams, key)) {
    const int ival = w.core.get_opt(key);
    if (const EnumOpt *e = find_enum(key))
      return py::cast(std::string(ival ? e->on : e->off));
    return py::cast(bool(ival));
  }
  if (find_desc(kWeatherParams, key))
    return py::cast(w.core.get_weather_param(key));
  if (key == "ref_delta" || key == "dut1")
    return py::cast(w.core.get_double(key));
  throw py::key_error("Unknown parameter: " + key);
}

// ---------------------------------------------------------------------------
// Vectorized entry points
//
// The n-loop lives here, not in the core: core methods are scalar. The loop
// is serial because the correction rate-cache is sequential across samples,
// and it runs with the GIL released.
// ---------------------------------------------------------------------------

py::object py_det_offset(py::object daz, py::object del_, py::object dpsi) {
  ColSet cs;
  auto &az = cs.add(daz, "delta_az");
  auto &el = cs.add(del_, "delta_el");
  auto &psi = cs.add(dpsi, "delta_psi");
  const auto n = cs.resolve();

  return quat_out(n, true, [&](py::ssize_t i) {
    return det_offset(az[i], el[i], psi[i]);
  });
}

py::object py_refraction(py::object el, py::object temp, py::object press,
                         py::object hum, py::object freq) {
  ColSet cs;
  auto &vel = cs.add(el, "el");
  auto &vtemp = cs.add(temp, "temp");
  auto &vpress = cs.add(press, "press");
  auto &vhum = cs.add(hum, "hum");
  auto &vfreq = cs.add(freq, "freq");
  const auto n = cs.resolve();

  return vec_out<1>(n, [&](py::ssize_t i, auto &o) {
    o[0][i] = refraction(vel[i], vtemp[i], vpress[i], vhum[i], vfreq[i]);
  });
}

py::object py_hwp_quat(py::object theta) {
  ColSet cs;
  auto &th = cs.add(theta, "theta");
  const auto n = cs.resolve();

  return quat_out(n, true, [&](py::ssize_t i) { return hwp_quat(th[i]); });
}

py::object py_gmst(PointingWrap &w, py::object ctime) {
  ColSet cs;
  auto &ct = cs.add(ctime, "ctime");
  const auto n = cs.resolve();

  return vec_out<1>(n, [&](py::ssize_t i, auto &o) {
    o[0][i] = w.core.gmst(ct[i]);
  });
}

py::object py_lmst(PointingWrap &w, py::object ctime, py::object lon) {
  ColSet cs;
  auto &ct = cs.add(ctime, "ctime");
  auto &lo = cs.add(lon, "lon");
  const auto n = cs.resolve();

  return vec_out<1>(n, [&](py::ssize_t i, auto &o) {
    o[0][i] = w.core.lmst(ct[i], lo[i]);
  });
}

// One entry point covers both azel2bore and azelpsi2bore: psi is just an
// optional argument. The C needs two functions because the loop lives there.
py::object py_azel2bore(PointingWrap &w, py::object az, py::object el,
                        py::object psi, py::object pitch, py::object roll,
                        py::object lon, py::object lat, py::object ctime) {
  ColSet cs;
  auto &vaz = cs.add(az, "az");
  auto &vel = cs.add(el, "el");
  auto &vpsi = cs.add(psi, "psi");
  auto &vpitch = cs.add(pitch, "pitch");
  auto &vroll = cs.add(roll, "roll");
  auto &vlon = cs.add(lon, "lon");
  auto &vlat = cs.add(lat, "lat");
  auto &vct = cs.add(ctime, "ctime");
  const auto n = cs.resolve();

  // Boresight quaternions keep their leading axis even for one sample,
  // unlike det_offset/hwp_quat. Both reference implementations do this.
  return quat_out(n, false, [&](py::ssize_t i) {
    // azelpsi2quat accumulates, so seed with identity
    Quat q = Quat::identity();
    w.core.azelpsi2quat(vaz[i], vel[i], vpsi[i], vpitch[i], vroll[i], vlon[i],
                        vlat[i], vct[i], q);
    return q;
  });
}

// The azel*2radec* family forces mean aberration on for the duration of the
// call. Going through set_opt rather than the field keeps the side effect
// the C has: changing the option also resets the annual-aberration rate
// state, so it is reset both on entry and on exit.
class MeanAberGuard {
 public:
  explicit MeanAberGuard(Pointing &p) : p_(p), saved_(p.get_opt("mean_aber")) {
    p_.set_opt("mean_aber", 1);
  }
  ~MeanAberGuard() { p_.set_opt("mean_aber", saved_); }

 private:
  Pointing &p_;
  int saved_;
};

// Covers the whole azel2radec / azelpsi2radec family, including the hwp and
// sindec and pa variants -- sixteen C entry points behind one loop.
py::object py_azel2radec(PointingWrap &w, double daz, double del, double dpsi,
                         py::object az, py::object el, py::object psi,
                         py::object pitch, py::object roll, py::object lon,
                         py::object lat, py::object ctime, py::object hwp,
                         bool sindec, bool return_pa) {
  ColSet cs;
  auto &vaz = cs.add(az, "az");
  auto &vel = cs.add(el, "el");
  auto &vpsi = cs.add(psi, "psi");
  auto &vpitch = cs.add(pitch, "pitch");
  auto &vroll = cs.add(roll, "roll");
  auto &vlon = cs.add(lon, "lon");
  auto &vlat = cs.add(lat, "lat");
  auto &vct = cs.add(ctime, "ctime");
  auto &vhwp = cs.add(hwp, "hwp");
  const auto n = cs.resolve();

  const DecOut dmode = sindec ? DecOut::SinDec : DecOut::Dec;
  const PolOut pmode = return_pa ? PolOut::PA : PolOut::SinCos;
  const Quat q_off = det_offset(daz, del, dpsi);

  // Brackets the loop, as it did when it sat inside the nogil block: it
  // touches only the core, never Python.
  MeanAberGuard guard(w.core);
  return radec_out(n, return_pa,
                   [&](py::ssize_t i, double *ra_, double *dec_, double *p1_,
                       double *p2_) {
                     // seeded with the offset, because azelpsi2quat
                     // accumulates
                     Quat q = q_off;
                     if (vhwp) mul_right(q, hwp_quat(vhwp.ref(i)));
                     w.core.azelpsi2quat(vaz[i], vel[i], vpsi[i], vpitch[i],
                                         vroll[i], vlon[i], vlat[i], vct[i], q);
                     w.core.quat2radec(q, dmode, pmode, ra_[i], dec_[i], p1_[i],
                                       p2_ ? &p2_[i] : nullptr);
                   });
}

// Collapses the six qp_bore2radec* entry points. The dec and pol axes are
// independent, so (sindec, pa) works too -- the C has no such function.
py::object py_bore2radec(PointingWrap &w, py::object q_off, py::object ctime,
                         py::object q_bore, py::object q_hwp, bool sindec,
                         bool return_pa) {
  const Quat off = load_quat(q_off, "q_off");
  ColSet cs;
  auto &vct = cs.add(ctime, "ctime");
  auto &vbore = cs.add_quat(q_bore, "q_bore", true);
  auto &vhwp = cs.add_quat(q_hwp, "q_hwp");
  const auto n = cs.resolve();

  const DecOut dmode = sindec ? DecOut::SinDec : DecOut::Dec;
  const PolOut pmode = return_pa ? PolOut::PA : PolOut::SinCos;
  return radec_out(
      n, return_pa,
      [&](py::ssize_t i, double *ra, double *dec, double *p1, double *p2) {
        const Quat q = vhwp
                           ? w.core.bore2det_hwp(off, vct[i], vbore.ref(i), vhwp.ref(i))
                           : w.core.bore2det(off, vct[i], vbore.ref(i));
        w.core.quat2radec(q, dmode, pmode, ra[i], dec[i], p1[i],
                          p2 ? &p2[i] : nullptr);
      });
}

py::object py_bore2azel(PointingWrap &w, py::object q_bore, py::object lon,
                        py::object lat, py::object ctime) {
  ColSet cs;
  auto &vbore = cs.add_quat(q_bore, "q_bore", true);
  auto &vlon = cs.add(lon, "lon");
  auto &vlat = cs.add(lat, "lat");
  auto &vct = cs.add(ctime, "ctime");
  const auto n = cs.resolve();

  return vec_out<3>(n, [&](py::ssize_t i, auto &o) {
    w.core.quat2azel(vbore.ref(i), vlon[i], vlat[i], vct[i], o[0][i], o[1][i],
                     o[2][i]);
  });
}

py::object py_radec2azel(PointingWrap &w, py::object ra, py::object dec,
                         py::object pa, py::object lon, py::object lat,
                         py::object ctime) {
  ColSet cs;
  auto &vra = cs.add(ra, "ra");
  auto &vdec = cs.add(dec, "dec");
  auto &vpa = cs.add(pa, "pa");
  auto &vlon = cs.add(lon, "lon");
  auto &vlat = cs.add(lat, "lat");
  auto &vct = cs.add(ctime, "ctime");
  const auto n = cs.resolve();

  return vec_out<3>(n, [&](py::ssize_t i, auto &o) {
    const Quat q = w.core.radecpa2quat(vra[i], vdec[i], vpa[i]);
    w.core.quat2azel(q, vlon[i], vlat[i], vct[i], o[0][i], o[1][i], o[2][i]);
  });
}

py::object py_radecpa2quat(PointingWrap &w, py::object ra, py::object dec,
                           py::object pa) {
  ColSet cs;
  auto &vra = cs.add(ra, "ra");
  auto &vdec = cs.add(dec, "dec");
  auto &vpa = cs.add(pa, "pa");
  const auto n = cs.resolve();

  return quat_out(n, true, [&](py::ssize_t i) {
    return w.core.radecpa2quat(vra[i], vdec[i], vpa[i]);
  });
}

py::object py_quat2radecpa(PointingWrap &w, py::object quat) {
  ColSet cs;
  auto &vq = cs.add_quat(quat, "quat", true);
  const auto n = cs.resolve();

  return vec_out<3>(n, [&](py::ssize_t i, auto &o) {
    w.core.quat2radec(vq.ref(i), DecOut::Dec, PolOut::PA, o[0][i], o[1][i],
                      o[2][i], nullptr);
  });
}

py::object py_update_ref(PointingWrap &w, py::object q) {
  ColSet cs;
  auto &vq = cs.add_quat(q, "q", true);
  const auto n = cs.resolve();

  return vec_out<1>(n, [&](py::ssize_t i, auto &o) {
    o[0][i] = w.core.update_ref(vq.ref(i));
  });
}

void py_set_bulletin_a(PointingWrap &w, int mjd_min, int mjd_max,
                            py::object dut1, py::object x, py::object y) {
  const py::ssize_t n = mjd_max - mjd_min + 1;
  auto rd = parse_col(dut1, "dut1");
  auto rx = parse_col(x, "x");
  auto ry = parse_col(y, "y");
  for (auto *r : {&rd, &rx, &ry})
    if (r->len != n)
      throw py::value_error("bulletin A arrays must have length " +
                            std::to_string(n) + " to span mjd_min..mjd_max");
  w.core.bulletin().set(mjd_min, mjd_max, rd.ptr, rx.ptr, ry.ptr);
}

py::object py_get_bulletin_a(PointingWrap &w, py::object mjd) {
  ColSet cs;
  auto &vm = cs.add(mjd, "mjd");
  const auto n = cs.resolve();

  return vec_out<3>(n, [&](py::ssize_t i, auto &o) {
    const IersValues v = w.core.bulletin().interp(vm[i]);
    o[0][i] = v.dut1;
    o[1][i] = v.x;
    o[2][i] = v.y;
  });
}

// ---------------------------------------------------------------------------
// Dipole
// ---------------------------------------------------------------------------

py::object py_dipole(PointingWrap &w, py::object ctime, py::object ra,
                     py::object dec) {
  ColSet cs;
  auto &vct = cs.add(ctime, "ctime");
  auto &vra = cs.add(ra, "ra");
  auto &vdec = cs.add(dec, "dec");
  const auto n = cs.resolve();

  return vec_out<1>(n, [&](py::ssize_t i, auto &o) {
    o[0][i] = w.core.dipole(vct[i], vra[i], vdec[i]);
  });
}

py::object py_bore2dipole(PointingWrap &w, py::object q_off, py::object ctime,
                          py::object q_bore) {
  const Quat off = load_quat(q_off, "q_off");
  ColSet cs;
  auto &vct = cs.add(ctime, "ctime");
  auto &vbore = cs.add_quat(q_bore, "q_bore", true);
  const auto n = cs.resolve();

  return vec_out<1>(n, [&](py::ssize_t i, auto &o) {
    const Quat q = w.core.bore2det(off, vct[i], vbore.ref(i));
    o[0][i] = w.core.quat2dipole(vct[i], q);
  });
}

// ---------------------------------------------------------------------------
// Boresight offset and gyro integration
// ---------------------------------------------------------------------------

// Adjusts q_bore in place by a per-sample offset. post selects whether the
// offset is applied in the detector frame or the sky frame.
void py_bore_offset(PointingWrap &w, py::object q_bore, py::object ang1,
                    py::object ang2, py::object ang3, bool post) {
  auto vbore = out_span<Quat>(q_bore, "q_bore");
  const py::ssize_t n = static_cast<py::ssize_t>(vbore.len);

  const auto a1 = fix(parse_col(ang1, "ang1"), n, "ang1");
  const auto a2 = fix(parse_col(ang2, "ang2"), n, "ang2");
  const auto a3 = fix(parse_col(ang3, "ang3"), n, "ang3");

  {
    py::gil_scoped_release nogil;
    for (py::ssize_t i = 0; i < n; ++i) {
      if (post) {
        const Quat q = w.core.radecpa2quat(a1[i], a2[i], a3[i]);
        mul_left(q, vbore[i]);
      } else {
        mul_right(vbore[i], det_offset(a1[i], a2[i], a3[i]));
      }
    }
  }
}

// Integrates gyro rates into an attitude timestream. The loop carries the
// attitude quaternion forward, so it is sequential by construction -- but
// that state is a local here, not something the core holds.
// Takes no Pointing state at all: the attitude is a local. Bound as a
// lambda that drops the wrapper argument.
py::object py_omega2azelpsi(double init_az, double init_el, double init_psi,
                            py::object omega_x, py::object omega_y,
                            py::object omega_z, double dt) {
  ColSet cs;
  auto &ox = cs.add(omega_x, "omega_x");
  auto &oy = cs.add(omega_y, "omega_y");
  auto &oz = cs.add(omega_z, "omega_z");
  const auto n = cs.resolve();

  // The attitude carries across samples, so it is declared outside the
  // body the loop calls.
  Quat attitude = azelpsi_quat(init_az, init_el, init_psi, 0.0, 0.0);
  return vec_out<3>(n, [&](py::ssize_t i, auto &o) {
    const Vec3 omega{{ox[i], oy[i], oz[i]}};
    const double mag = omega.norm();
    if (std::fabs(mag) > 1e-12) {
      const Vec3 axis{{omega[0] / mag, omega[1] / mag, omega[2] / mag}};
      mul_right(attitude, Quat::rot(mag * dt, axis));
      attitude.unit();
    }
    quat_azelpsi(attitude, o[0][i], o[1][i], o[2][i]);
  });
}

// ---------------------------------------------------------------------------
// Pixelization
// ---------------------------------------------------------------------------

py::object py_radec2pix(PointingWrap &w, py::object ra, py::object dec,
                        int nside) {
  ColSet cs;
  auto &vra = cs.add(ra, "ra");
  auto &vdec = cs.add(dec, "dec");
  const auto n = cs.resolve();

  return long_out(n, [&](py::ssize_t i) {
    return w.core.radec2pix(vra[i], vdec[i], nside);
  });
}

py::object py_quat2pix(PointingWrap &w, py::object quat, int nside,
                       bool return_pa) {
  ColSet cs;
  auto &vq = cs.add_quat(quat, "quat", true);
  const auto n = cs.resolve();
  const PolOut pmode = return_pa ? PolOut::PA : PolOut::SinCos;

  return pix_out(n, return_pa,
                 [&](py::ssize_t i, long *pix, double *p1, double *p2) {
                   w.core.quat2pix(vq.ref(i), nside, pmode, pix[i], p1[i],
                                   p2 ? &p2[i] : nullptr);
                 });
}

py::object py_bore2pix(PointingWrap &w, py::object q_off, py::object ctime,
                       py::object q_bore, py::object q_hwp, int nside,
                       bool return_pa) {
  const Quat off = load_quat(q_off, "q_off");
  ColSet cs;
  auto &vct = cs.add(ctime, "ctime");
  auto &vbore = cs.add_quat(q_bore, "q_bore", true);
  auto &vhwp = cs.add_quat(q_hwp, "q_hwp");
  const auto n = cs.resolve();

  const PolOut pmode = return_pa ? PolOut::PA : PolOut::SinCos;
  return pix_out(n, return_pa,
                 [&](py::ssize_t i, long *pix, double *p1, double *p2) {
                   const Quat q =
                       vhwp ? w.core.bore2det_hwp(off, vct[i], vbore.ref(i),
                                                  vhwp.ref(i))
                            : w.core.bore2det(off, vct[i], vbore.ref(i));
                   w.core.quat2pix(q, nside, pmode, pix[i], p1[i],
                                   p2 ? &p2[i] : nullptr);
                 });
}

// The ring table is built once per call and shared across the loop, rather
// than rebuilt per sample as qp_get_interp_valn does.
py::object py_get_interp_val(PointingWrap &w, int nside, py::object map,
                             py::object ra, py::object dec) {
  auto rmap = parse_col(map, "map");
  if (!rmap.ptr) throw py::value_error("map must be an array");
  const long npix = ::nside2npix(nside);
  if (rmap.len != npix)
    throw py::value_error("map has length " + std::to_string(rmap.len) +
                          ", expected " + std::to_string(npix) +
                          " for nside " + std::to_string(nside));

  ColSet cs;
  auto &vra = cs.add(ra, "ra");
  auto &vdec = cs.add(dec, "dec");
  const auto n = cs.resolve();
  const bool nest = w.core.opt().pix_order == 1;

  // One ring table for the whole call, rather than one per sample as
  // qp_get_interp_valn does.
  const PixInfo pixinfo(nside);
  return vec_out<1>(n, [&](py::ssize_t i, auto &o) {
    o[0][i] = pixinfo.interp_val(rmap.ptr, vra[i], vdec[i], nest);
  });
}

// Resample a polarized map into another coordinate frame.
//
// For each output pixel: take its direction, rotate backwards to find where
// it came from, read T/Q/U there, then rotate the polarization basis
// forwards. Matches healpy's Rotator.rotate_map_pixel to ~5e-7 with
// interpolation on.
//
// Unlike the C, this checks the row count. qp_rotate_map dereferences
// map_in[1] and map_in[2] unconditionally, so a T-only map reads past the
// end of the row table and segfaults.
py::object py_rotate_map(PointingWrap &w, int nside, py::object map_in,
                         bool to_gal) {
  MapArray in = parse_map(map_in, "map_in", false);
  if (!in.ptr) throw py::value_error("map_in is required");
  if (in.nrow != 3)
    throw py::value_error(
        "rotate_map needs a polarized map with exactly 3 rows (T, Q, U), got " +
        std::to_string(in.nrow));

  const long npix = ::nside2npix(nside);
  if (static_cast<long>(in.npix) != npix)
    throw py::value_error("map_in has " + std::to_string(in.npix) +
                          " pixels, expected " + std::to_string(npix) +
                          " for nside " + std::to_string(nside));

  const double *ti = in.ptr, *qi = in.ptr + npix, *ui = in.ptr + 2 * npix;
  py::array_t<double> out({py::ssize_t(3), py::ssize_t(npix)});
  double *to = out.mutable_data();
  std::memset(to, 0, static_cast<std::size_t>(3 * npix) * sizeof(double));
  double *qo = to + npix, *uo = to + 2 * npix;

  const bool interp = w.core.opt().interp_pix;
  const bool nest = w.core.opt().pix_order == 1;
  {
    py::gil_scoped_release nogil;
    std::optional<PixInfo> pixinfo;
    if (interp) pixinfo.emplace(nside);

    for (long ii = 0; ii < npix; ++ii) {
      double ra, dec;
      w.core.pix2radec(nside, ii, ra, dec);
      double sin2psi = 0., cos2psi = 1.;

      // backwards, to the direction this pixel came from
      w.core.rotate_coord(ra, dec, sin2psi, cos2psi, !to_gal);

      double t, q, u;
      if (interp) {
        long pix[4];
        double wt[4];
        pixinfo->interpol(ra, dec, nest, pix, wt);
        t = q = u = 0;
        for (int jj = 0; jj < 4; ++jj) {
          t += ti[pix[jj]] * wt[jj];
          q += qi[pix[jj]] * wt[jj];
          u += ui[pix[jj]] * wt[jj];
        }
      } else {
        const long pix = w.core.radec2pix(ra, dec, nside);
        t = ti[pix];
        q = qi[pix];
        u = ui[pix];
      }

      if (t == 0 && q == 0 && u == 0) continue;

      // forwards, for the polarization basis at the source position
      sin2psi = 0.;
      cos2psi = 1.;
      w.core.rotate_coord(ra, dec, sin2psi, cos2psi, to_gal);

      to[ii] = t;
      qo[ii] = q * cos2psi + u * sin2psi;
      uo[ii] = u * cos2psi - q * sin2psi;
    }
  }
  return std::move(out);
}

// ---------------------------------------------------------------------------
// Galactic rotation, all in place
// ---------------------------------------------------------------------------

void py_rotate_quat(PointingWrap &w, py::object quat, bool to_gal) {
  auto q = out_span<Quat>(quat, "quat");
  const py::ssize_t n = static_cast<py::ssize_t>(q.len);
  {
    py::gil_scoped_release nogil;
    for (py::ssize_t i = 0; i < n; ++i) {
      if (to_gal)
        w.core.radec2gal_quat(q[i]);
      else
        w.core.gal2radec_quat(q[i]);
    }
  }
}

void py_rotate_coord(PointingWrap &w, py::object ra, py::object dec,
                     py::object pa, py::object sin2psi, py::object cos2psi,
                     bool to_gal) {
  const bool do_pa = !pa.is_none();
  auto vra = out_span<double>(ra, "ra");
  const py::ssize_t n = static_cast<py::ssize_t>(vra.len);

  auto vdec = out_span<double>(dec, "dec", n);
  auto vp1 = out_span<double>(do_pa ? pa : sin2psi, do_pa ? "pa" : "sin2psi", n);
  Span<double> vp2;
  if (!do_pa) vp2 = out_span<double>(cos2psi, "cos2psi", n);

  {
    py::gil_scoped_release nogil;
    for (py::ssize_t i = 0; i < n; ++i) {
      Quat q = do_pa ? w.core.radecpa2quat(vra[i], vdec[i], vp1[i])
                     : w.core.radec2quat(vra[i], vdec[i], vp1[i], vp2[i]);
      if (to_gal)
        w.core.radec2gal_quat(q);
      else
        w.core.gal2radec_quat(q);
      w.core.quat2radec(q, DecOut::Dec, do_pa ? PolOut::PA : PolOut::SinCos,
                        vra[i], vdec[i], vp1[i], do_pa ? nullptr : &vp2[i]);
    }
  }
}

// ---------------------------------------------------------------------------
// Mapmaking structures
//
// Each wrapper owns the py::object for every array its Map/Det/PointData
// views into, so the buffer outlives the call regardless of what the caller
// keeps. The core structs hold only non-owning spans.
// ---------------------------------------------------------------------------

bool is_false(const py::handle &o) {
  return py::isinstance<py::bool_>(o) && !o.cast<bool>();
}

// A flag the caller may have no opinion on: None takes the default.
// py::object rather than std::optional<bool>, since the latter needs
// pybind11/stl.h, whose container casters would add exactly the silent
// copies the array contract exists to prevent.
bool flag_or(const py::handle &o, bool deflt) {
  return o.is_none() ? deflt : o.cast<bool>();
}

// chealpix's npix2nside uses an integer sqrt and returns -1 when npix is
// not 12*nside^2, so this is exact where a float sqrt would not be.
long checked_npix2nside(py::ssize_t npix) {
  const long n = ::npix2nside(static_cast<long>(npix));
  if (n < 0) throw py::value_error("Invalid npix, must be 12 * nside**2");
  return n;
}

// Mirrors qmap_class.check_map: a C-contiguous float64 (nrow, npix)
// array. A 1-D map becomes one row.
//
// Row count of whichever mode, so one template can serve both components.
std::size_t mode_rows(VecMode m) { return num_vec(m); }
std::size_t mode_rows(ProjMode m) { return num_proj(m); }

// map.hpp owns the mode <-> row-count table; these two add the Python-facing
// error, which the core has no business raising.
VecMode infer_vec_mode(std::size_t nrow, bool pol, bool vpol) {
  const VecMode m = vec_mode_for(nrow, pol, vpol);
  if (m == VecMode::None)
    throw py::value_error("Cannot infer vec mode from " + std::to_string(nrow) +
                          " map rows");
  return m;
}

ProjMode infer_proj_mode(std::size_t nrow) {
  const ProjMode m = proj_mode_for(nrow);
  if (m == ProjMode::None)
    throw py::value_error("Cannot infer proj mode from " +
                          std::to_string(nrow) + " projection rows");
  return m;
}

class MapWrap {
 public:
  Map core;
  bool inited = false;

  // Each of vec and proj takes a map, None to allocate one of zeros, or
  // False for a map with no such component -- a source map has no
  // projection, and a dest map can accumulate either alone.
  //
  // A supplied map fixes npix, and its row count fixes the mode. A default
  // is sized from the other component when there is one, so the two always
  // describe the same number of map components, and from pol/vpol when
  // there isn't.
  //
  // pol and vpol are optional because that is the truth about them: they
  // decide a mode only where a row count leaves one open -- 3 rows being
  // T,Q,U or T with a first derivative, 4 rows being T,Q,U,V or nothing
  // valid -- and are overridden outright wherever a supplied map settles
  // it. None asks for that reading; the fallbacks here are the documented
  // API defaults, applied at the point of use rather than in the
  // signature, so no caller has to pass a flag it has no opinion on.
  void setup(py::object vec, py::object proj, py::object pixels,
             py::object nside, py::object pol_arg, py::object vpol_arg) {
    reset();
    bool pol = flag_or(pol_arg, true);
    bool vpol = flag_or(vpol_arg, false);
    pol_ = pol;
    vpol_ = vpol;

    const bool want_vec = !is_false(vec);
    const bool want_proj = !is_false(proj);
    if (!want_vec && !want_proj)
      throw py::value_error("one of vec or proj must not be False");

    const bool partial = !pixels.is_none();
    if (partial && nside.is_none())
      throw py::value_error("nside required for partial maps");

    const py::ssize_t expect_npix = partial ? py::len(pixels) : 0;

    MapArray v, p;
    const bool vec_given = want_vec && !vec.is_none();
    const bool proj_given = want_proj && !proj.is_none();
    if (vec_given) v = parse_map(vec, "vec", false);
    if (proj_given) p = parse_map(proj, "proj", false);
    if (vec_given && proj_given && v.npix != p.npix)
      throw py::value_error("vec and proj must have the same npix");

    // npix and nside, from whichever of a supplied map, the pixel list or
    // nside is available. The 256 default is the documented API default.
    if (vec_given || proj_given) {
      core.npix = vec_given ? v.npix : p.npix;
      core.nside =
          nside.is_none() ? static_cast<std::size_t>(checked_npix2nside(
                                static_cast<py::ssize_t>(core.npix)))
                          : nside.cast<std::size_t>();
    } else {
      core.nside = nside.is_none() ? 256 : nside.cast<std::size_t>();
      core.npix = partial ? static_cast<std::size_t>(expect_npix)
                          : static_cast<std::size_t>(
                                nside2npix(static_cast<long>(core.nside)));
    }

    // A partial map is the one case a shape cannot be judged on its own --
    // a (3, 2) map of two pixels is taller than it is wide and still right
    // -- so the pixel list is what catches a transposed one.
    if (partial && core.npix != static_cast<std::size_t>(expect_npix))
      throw py::value_error(
          "partial map has " + std::to_string(core.npix) + " columns for " +
          std::to_string(expect_npix) + " pixels");

    if (want_vec) {
      // p.nrow is 0 unless a proj was supplied, in which case the default
      // vec follows it rather than pol/vpol -- and so does the mode, since
      // a 10-row proj describes a TQUV map whatever vpol said. Without
      // that, a defaulted vec was sized from the proj and then rejected
      // for not matching the flags.
      if (!vec_given) {
        v = parse_map(zeros_map(vec_rows(p.nrow)), "vec", false);
        if (proj_given) {
          pol = v.nrow >= 3;
          vpol = v.nrow == 4;
          pol_ = pol;
          vpol_ = vpol;
        }
      }
      core.vec_mode = infer_vec_mode(v.nrow, pol, vpol);
      core.vec = v.ptr;
      vec_keep_ = std::move(v.keep);
    }
    if (want_proj) {
      if (!proj_given)
        p = parse_map(zeros_map(proj_rows(v.nrow)), "proj", false);
      core.proj_mode = infer_proj_mode(p.nrow);
      core.proj = p.ptr;
      proj_keep_ = std::move(p.keep);
    }

    // An N-component map has N vec rows and N*(N+1)/2 proj rows. Supplying
    // both and getting that wrong hands the kernels two maps describing
    // different things.
    if (want_vec && want_proj && num_proj_for_vec(v.nrow) != p.nrow)
      throw py::value_error("proj shape incompatible with vec: " +
                            std::to_string(p.nrow) + " rows for a " +
                            std::to_string(v.nrow) + "-row vec");

    if (partial) {
      py::array a = py::cast<py::array>(pixels);
      if (!a.dtype().is(py::dtype::of<long>()))
        throw py::type_error("pixels must be an int64 array");
      if (!(a.flags() & py::array::c_style))
        throw py::value_error("pixels must be C-contiguous");
      pixhash_.emplace(Span<const long>{static_cast<const long *>(a.data()),
                                        static_cast<std::size_t>(a.size())});
      core.pixhash = &pixhash_.value();
      pix_keep_ = std::move(a);
    }

    inited = true;
  }

  // The same three argument forms as setup, against an initialized map:
  // a map replaces the installed buffer and must match its shape, since
  // re-deciding the mode would change what the kernels read out from under
  // a half-filled accumulator; None installs a fresh map of zeros; False
  // disables the component.
  void update_vec(py::object arr) {
    update_component(arr, "vec", core.vec, core.vec_mode, vec_keep_,
                     core.has_proj(), vec_rows(num_proj(core.proj_mode)),
                     [&](std::size_t nrow) {
                       return infer_vec_mode(nrow, pol_, vpol_);
                     });
  }

  void update_proj(py::object arr) {
    update_component(arr, "proj", core.proj, core.proj_mode, proj_keep_,
                     core.has_vec(), proj_rows(num_vec(core.vec_mode)),
                     infer_proj_mode);
  }

  void reset() {
    core = Map();
    pol_ = true;
    vpol_ = false;
    vec_keep_ = py::none();
    proj_keep_ = py::none();
    pix_keep_ = py::none();
    pixhash_.reset();
    pixinfo_.reset();
    inited = false;
  }

  // Called before entering a parallel region. PixInfo fills its ring table
  // lazily, which mutates the cache on read, so it must be populated up
  // front before threads share it.
  void ensure_pixinfo() {
    if (!pixinfo_) {
      pixinfo_ = std::make_unique<PixInfo>(static_cast<long>(core.nside));
      core.pixinfo = pixinfo_.get();
    }
    pixinfo_->populate();
  }

  py::object get_vec() const { return vec_keep_; }
  py::object get_proj() const { return proj_keep_; }
  bool is_init() const { return inited; }
  bool is_partial() const { return core.partial(); }
  bool has_vec() const { return core.has_vec(); }
  bool has_proj() const { return core.has_proj(); }
  // What the map contains, which is a different question from the one the
  // kernels ask. at_least(vec_mode, Pol) is an ordering test over an enum
  // that interleaves the derivative modes, so it calls D1 and D2
  // polarized; and a dest can be proj-only, with no vec mode to read at
  // all, so the proj has to answer when the vec cannot.
  bool is_pol() const {
    return is_polarized(core.vec_mode) ||
           at_least(core.proj_mode, ProjMode::Pol);
  }
  bool is_vpol() const {
    return core.vec_mode == VecMode::VPol || core.proj_mode == ProjMode::VPol;
  }

 private:
  // update_vec and update_proj are the same three branches over different
  // members. other_present is whether the map would still have a component
  // after this one is switched off; deflt is the shape to allocate when it
  // is switched back on.
  template <class Mode, class Infer>
  void update_component(py::object arr, const char *name, double *&ptr,
                        Mode &mode, py::object &keep, bool other_present,
                        std::size_t deflt, Infer &&infer) {
    if (is_false(arr)) {
      if (!other_present)
        throw py::value_error("one of vec or proj must not be False");
      ptr = nullptr;
      mode = Mode::None;
      keep = py::none();
      return;
    }

    const bool present = ptr != nullptr && mode != Mode::None;
    if (arr.is_none()) {
      // At the installed shape, or the one the mode implies when the
      // component is switched off. Replacing the buffer rather than
      // zeroing it leaves whatever the caller still holds alone.
      arr = zeros_map(present ? mode_rows(mode) : deflt);
    }

    MapArray a = parse_map(arr, name, false);
    if (a.npix != core.npix || (present && a.nrow != mode_rows(mode)))
      throw py::value_error(std::string(name) +
                            " shape does not match the existing map");
    mode = infer(a.nrow);
    ptr = a.ptr;
    keep = std::move(a.keep);
  }

  // Zeros at the map's pixel count, ready for parse_map: numpy hands back
  // a C-contiguous float64 block, and zeros() is what the Python layer
  // used to allocate here.
  py::object zeros_map(std::size_t nrow) const {
    py::module_ np = py::module_::import("numpy");
    return np.attr("zeros")(py::make_tuple(nrow, core.npix));
  }

  // Rows of a default component. Passing the other component's row count
  // matches the two up -- an N-component map has N vec rows and
  // N*(N+1)/2 proj rows -- and 0 means there is no other component to
  // follow, so pol/vpol as given to setup decide.
  std::size_t vec_rows(std::size_t nproj) const {
    if (const std::size_t n = num_vec_for_proj(nproj)) return n;
    return vpol_ ? 4 : pol_ ? 3 : 1;
  }

  std::size_t proj_rows(std::size_t nvec) const {
    if (const std::size_t n = num_proj_for_vec(nvec)) return n;
    return vpol_ ? 10 : pol_ ? 6 : 1;
  }

  // As passed to setup, so a component re-enabled later is sized the way
  // the caller asked for originally. The installed modes cannot stand in:
  // a 3-row source map with pol off is D1, not T,Q,U.
  bool pol_ = true;
  bool vpol_ = false;
  py::object vec_keep_ = py::none();
  py::object proj_keep_ = py::none();
  py::object pix_keep_ = py::none();
  std::optional<PixHash> pixhash_;
  std::unique_ptr<PixInfo> pixinfo_;
};

class PointWrap {
 public:
  PointData core;

  void set_bore(py::object q_bore) {
    auto r = parse_quat_col(q_bore, "q_bore");
    if (!r.present) throw py::value_error("q_bore is required");
    core.q_bore = {r.ptr, static_cast<std::size_t>(r.len)};
    bore_keep_ = r.keep;
  }

  void set_ctime(py::object ctime) {
    auto r = parse_col(ctime, "ctime");
    if (!r.ptr) throw py::value_error("ctime must be an array");
    if (static_cast<std::size_t>(r.len) != core.n())
      throw py::value_error("ctime length does not match q_bore");
    core.ctime = {r.ptr, static_cast<std::size_t>(r.len)};
    ctime_keep_ = r.keep;
  }

  void clear_ctime() {
    core.ctime = {};
    ctime_keep_ = py::none();
  }

  void set_hwp(py::object q_hwp) {
    auto r = parse_quat_col(q_hwp, "q_hwp");
    if (!r.present) throw py::value_error("q_hwp is required");
    if (static_cast<std::size_t>(r.len) != core.n())
      throw py::value_error("q_hwp length does not match q_bore");
    core.q_hwp = {r.ptr, static_cast<std::size_t>(r.len)};
    hwp_keep_ = r.keep;
  }

  void clear_hwp() {
    core.q_hwp = {};
    hwp_keep_ = py::none();
  }

  void reset() {
    core = PointData();
    bore_keep_ = py::none();
    ctime_keep_ = py::none();
    hwp_keep_ = py::none();
  }

  bool is_init() const { return core.q_bore.size() > 0; }
  std::size_t n_samples() const { return core.n(); }

 private:
  py::object bore_keep_ = py::none();
  py::object ctime_keep_ = py::none();
  py::object hwp_keep_ = py::none();
};

class DetArrWrap {
 public:
  std::vector<Det> dets;
  bool diff = false;
  bool inited = false;

  // Detector pairs are folded together in diff mode, so the loop runs over
  // half the array. The C halves dets->n in place instead, which means
  // calling it twice on the same detarr halves it again.
  std::size_t n_loop() const { return diff ? dets.size() / 2 : dets.size(); }

  void setup(py::object q_off, py::object weight, py::object gain,
             py::object mueller, py::object tod, py::object flag,
             py::object weights, py::ssize_t n_samp, bool do_diff,
             bool write) {
    reset();

    auto roff = parse_quat_col(q_off, "q_off");
    if (!roff.present) throw py::value_error("q_off is required");
    const py::ssize_t ndet = roff.len;
    qoff_keep_ = roff.keep;

    auto rw = parse_col(weight, "weight");
    auto rg = parse_col(gain, "gain");
    weight_keep_ = rw.keep;
    gain_keep_ = rg.keep;

    const double *mu = nullptr;
    if (!mueller.is_none()) {
      py::array a = as_array(mueller, "mueller");
      if (a.size() != ndet * 4)
        throw py::value_error("mueller shape does not match q_off");
      mu = static_cast<const double *>(a.data());
      mueller_keep_ = std::move(a);
    }

    // tod2map only reads the tod; map2tod writes it, so that path needs a
    // writable buffer and allocates one when none was supplied.
    double *tod_ptr = nullptr;
    if (!tod.is_none()) {
      MapArray t = parse_map(tod, "tod", write);
      if (!t.ptr) throw py::value_error("tod must be a 2-D array");
      if (static_cast<py::ssize_t>(t.nrow) != ndet ||
          static_cast<py::ssize_t>(t.npix) != n_samp)
        throw py::value_error("tod shape does not match q_off and n_samp");
      tod_ptr = t.ptr;
      tod_keep_ = std::move(t.keep);
    } else if (write) {
      py::array_t<double> a({ndet, n_samp});
      std::memset(a.mutable_data(), 0,
                  static_cast<std::size_t>(ndet * n_samp) * sizeof(double));
      tod_ptr = a.mutable_data();
      tod_keep_ = std::move(a);
    }

    const std::uint8_t *flag_ptr = nullptr;
    if (!flag.is_none()) {
      py::array a = py::cast<py::array>(flag);
      if (!a.dtype().is(py::dtype::of<std::uint8_t>()))
        throw py::type_error("flag must be a uint8 array");
      if (!(a.flags() & py::array::c_style))
        throw py::value_error("flag must be C-contiguous");
      if (a.size() != ndet * n_samp)
        throw py::value_error("flag shape does not match q_off and n_samp");
      flag_ptr = static_cast<const std::uint8_t *>(a.data());
      flag_keep_ = std::move(a);
    }

    const double *wts_ptr = nullptr;
    if (!weights.is_none()) {
      MapArray a = parse_map(weights, "weights", false);
      if (static_cast<py::ssize_t>(a.nrow) != ndet ||
          static_cast<py::ssize_t>(a.npix) != n_samp)
        throw py::value_error("weights shape does not match q_off and n_samp");
      wts_ptr = a.ptr;
      weights_keep_ = std::move(a.keep);
    }

    const auto ns = static_cast<std::size_t>(n_samp);
    dets.resize(static_cast<std::size_t>(ndet));
    for (py::ssize_t i = 0; i < ndet; ++i) {
      Det &d = dets[static_cast<std::size_t>(i)];
      d.q_off = roff.ptr[i];
      d.weight = rw.present ? (rw.ptr ? rw.ptr[i] : rw.value) : 1.0;
      d.gain = rg.present ? (rg.ptr ? rg.ptr[i] : rg.value) : 1.0;
      if (mu)
        d.mueller = {{mu[i * 4], mu[i * 4 + 1], mu[i * 4 + 2], mu[i * 4 + 3]}};
      else
        d.mueller = {{1., 1., 0., 1.}};
      if (tod_ptr) d.tod = {tod_ptr + i * ns, ns};
      if (flag_ptr) d.flag = {flag_ptr + i * ns, ns};
      if (wts_ptr) d.weights = {wts_ptr + i * ns, ns};
    }

    diff = do_diff;
    inited = true;
  }

  void reset() {
    dets.clear();
    diff = false;
    inited = false;
    qoff_keep_ = py::none();
    weight_keep_ = py::none();
    gain_keep_ = py::none();
    mueller_keep_ = py::none();
    tod_keep_ = py::none();
    flag_keep_ = py::none();
    weights_keep_ = py::none();
  }

  py::object get_tod() const { return tod_keep_; }
  bool is_init() const { return inited; }

 private:
  py::object qoff_keep_ = py::none();
  py::object weight_keep_ = py::none();
  py::object gain_keep_ = py::none();
  py::object mueller_keep_ = py::none();
  py::object tod_keep_ = py::none();
  py::object flag_keep_ = py::none();
  py::object weights_keep_ = py::none();
};

// ---------------------------------------------------------------------------
// Mapmaking drivers
//
// OpenMP lives here rather than in the core: the core stays scalar and
// single-threaded, and this layer owns the thread-local copies.
// ---------------------------------------------------------------------------

// At most one thread per detector, and never fewer than one.
int clamp_threads(int nthr, std::ptrdiff_t ndet) {
  if (nthr > ndet) nthr = static_cast<int>(ndet);
  return nthr < 1 ? 1 : nthr;
}

// A zeroed double buffer, from calloc rather than a vector.
//
// The thread-local accumulators are as long as the whole map while a scan
// touches a fraction of a percent of it. vector::assign writes every byte
// -- 906 MB per thread at nside 1024 -- where calloc hands back pages the
// kernel zeroes lazily, so the untouched rest is never faulted in. That
// only pays off because the merge visits recorded pixels rather than
// scanning for nonzero ones; scanning would fault the whole thing in.
class ZeroBuf {
 public:
  ZeroBuf() = default;
  explicit ZeroBuf(std::size_t n)
      : ptr_(n ? static_cast<double *>(std::calloc(n, sizeof(double)))
                : nullptr) {
    if (n && !ptr_) throw std::bad_alloc();
  }
  ~ZeroBuf() { std::free(ptr_); }
  ZeroBuf(ZeroBuf &&o) noexcept : ptr_(o.ptr_) { o.ptr_ = nullptr; }
  ZeroBuf &operator=(ZeroBuf &&o) noexcept {
    std::swap(ptr_, o.ptr_);
    return *this;
  }
  ZeroBuf(const ZeroBuf &) = delete;
  ZeroBuf &operator=(const ZeroBuf &) = delete;

  double *get() const { return ptr_; }

 private:
  double *ptr_ = nullptr;
};

// A zeroed map with the same geometry, backed by thread-private buffers.
// This is the one copy the design cannot avoid -- each thread needs private
// storage to accumulate into -- but it is internal scratch, not data shared
// with Python.
//
// nslot is the row stride: g.npix for a directly-indexed accumulator, or
// the accumulator's capacity for a compact one, whose rows are indexed by
// slot. The pixhash is inherited either way, so a partial destination map
// still translates sky pixel to map index before any slot lookup.
Map blank_like(const Map &g, std::size_t nslot, ZeroBuf &vbuf, ZeroBuf &pbuf) {
  Map m = g;
  m.npix = nslot;
  if (g.has_vec()) {
    vbuf = ZeroBuf(num_vec(g.vec_mode) * nslot);
    m.vec = vbuf.get();
  }
  if (g.has_proj()) {
    pbuf = ZeroBuf(num_proj(g.proj_mode) * nslot);
    m.proj = pbuf.get();
  }
  return m;
}

// Whether a thread's private accumulator should be compacted.
//
// reach bounds the distinct pixels one thread can touch: one per sample
// per detector it handles, doubled for the diff kernel, which remaps a
// pair. Compacting costs a lookup per sample and saves a private copy of
// everything the thread did not touch, so it only pays when the map is
// much bigger than that -- a full-sky destination for a scan that covers
// a sliver of it. For a map that is already a compact footprint, reach is
// the same size as the map and this stays off.
bool want_compact(std::size_t npix, std::size_t reach) {
  return npix > 4 * reach;
}

// The detector loop, with all of the OpenMP plumbing in one place.
//
// per_thread() runs once inside each thread and returns that thread's
// state; body(state, i) handles one detector; finish(state) merges it.
// The state is returned by value, which is safe even if that moves: a
// moved std::vector keeps its buffer, so a Map pointing into one stays
// valid.
//
// Three things this encodes. An exception must not cross the structured
// block -- that is undefined behaviour -- so the first one is captured and
// rethrown afterwards, with the GIL reacquired, because raising through
// pybind's translator needs it. You cannot break out of an omp for, so the
// remaining iterations no-op instead. And num_threads is a clause rather
// than the C's process-wide omp_set_num_threads.
template <class Setup, class Body, class Finish>
void parallel_dets(int nthr, std::ptrdiff_t ndet, const char *what, bool balance,
                   Setup &&per_thread, Body &&body, Finish &&finish) {
  std::atomic<bool> failed{false};
  std::mutex errmtx;
  std::string errmsg;

#ifdef _OPENMP
  // Detectors cost the same as each other, but cores do not: split them
  // evenly across a machine with both performance and efficiency cores and
  // the fast cores sit waiting for the slow ones. Dynamic scheduling fixes
  // that -- 22% at ten threads on a 4+6 machine -- but it also lets one
  // thread take any number of the detectors, which a compact accumulator
  // cannot be sized for without assuming the worst. So the caller says
  // which of the two it needs.
  omp_set_schedule(balance ? omp_sched_dynamic : omp_sched_static,
                   balance ? 1 : 0);
#endif

  {
    py::gil_scoped_release nogil;
#ifdef _OPENMP
#pragma omp parallel num_threads(nthr)
#endif
    {
      auto state = per_thread();
#ifdef _OPENMP
#pragma omp for schedule(runtime)
#endif
      for (std::ptrdiff_t i = 0; i < ndet; ++i) {
        if (failed.load(std::memory_order_relaxed)) continue;
        try {
          body(state, i);
        } catch (const std::exception &e) {
          if (!failed.exchange(true)) {
            std::lock_guard<std::mutex> lk(errmtx);
            errmsg = e.what();
          }
        } catch (...) {
          if (!failed.exchange(true)) {
            std::lock_guard<std::mutex> lk(errmtx);
            errmsg = std::string("unknown error in ") + what;
          }
        }
      }
      if (!failed.load()) finish(state);
    }
  }

  if (failed) throw QpMapError(errmsg);
}

// One thread's share of a tod2map. The accumulator is private only when
// there is more than one thread; with one, it writes the real map directly
// and there is nothing to merge.
struct Tod2MapLocal {
  Pointing core;
  ZeroBuf vbuf, pbuf;
  // Heap-allocated so the pointers lmap holds into them survive this
  // struct being moved out of the setup callable.
  std::unique_ptr<PixTouch> touched;
  std::unique_ptr<PixAccum> accum;
  Map lmap;
  bool reduce = false;
};

void py_tod2map(PointingWrap &w, DetArrWrap &dets, PointWrap &pnt,
                MapWrap &dest) {
  if (!dets.is_init()) throw QpInitError("tod2map: detectors not initialized");
  if (!pnt.is_init()) throw QpInitError("tod2map: pointing not initialized");
  if (!dest.is_init()) throw QpInitError("tod2map: map not initialized");

  Map &global = dest.core;
  const PointData pd = pnt.core;
  const std::ptrdiff_t ndet = static_cast<std::ptrdiff_t>(dets.n_loop());
  const std::ptrdiff_t nhalf = static_cast<std::ptrdiff_t>(dets.dets.size() / 2);

  // Whether to compact has to be settled before the parallel region,
  // because it chooses the schedule. It uses the even split, which is
  // what a thread gets under the static schedule that compacting forces.
  const int nthr = clamp_threads(w.num_threads, ndet);
  const std::size_t mult = dets.diff ? 2 : 1;
  const std::size_t even_split =
      (static_cast<std::size_t>(ndet) + nthr - 1) / nthr;
  const bool compacting =
      nthr > 1 && want_compact(global.npix, even_split * pd.n() * mult);

  parallel_dets(
      nthr, ndet, "tod2map", !compacting,
      [&] {
        Tod2MapLocal st;
        st.core = w.core;  // O(1): the IERS table is shared
        int nthreads = 1;
#ifdef _OPENMP
        nthreads = omp_get_num_threads();
#endif
        st.reduce = nthreads > 1;
        if (!st.reduce) {
          // One thread writes the real map, so there is nothing to record
          // and nothing to merge.
          st.lmap = global;
          return st;
        }

        if (compacting) {
          // Sized from the thread count the runtime actually gave us,
          // which is what the static schedule will divide the detectors by.
          const std::size_t cap =
              (static_cast<std::size_t>(ndet) + nthreads - 1) / nthreads *
              pd.n() * mult;
          st.accum = std::make_unique<PixAccum>(cap);
          st.lmap = blank_like(global, cap, st.vbuf, st.pbuf);
          st.lmap.accum = st.accum.get();
        } else {
          st.touched = std::make_unique<PixTouch>(global.npix);
          st.lmap = blank_like(global, global.npix, st.vbuf, st.pbuf);
          st.lmap.touched = st.touched.get();
        }
        return st;
      },
      [&](Tod2MapLocal &st, std::ptrdiff_t i) {
        if (dets.diff)
          tod2map1_diff(st.core, dets.dets[i], dets.dets[i + nhalf], pd,
                        st.lmap);
        else
          tod2map1(st.core, dets.dets[i], pd, st.lmap);
      },
      [&](Tod2MapLocal &st) {
        if (!st.reduce) return;
#ifdef _OPENMP
#pragma omp critical(qp_map_accum)
#endif
        if (st.accum)
          add_map(global, st.lmap, *st.accum);
        else
          add_map(global, st.lmap, st.touched.get());
      });
}

void py_map2tod(PointingWrap &w, DetArrWrap &dets, PointWrap &pnt,
                MapWrap &source) {
  if (!dets.is_init()) throw QpInitError("map2tod: detectors not initialized");
  if (!pnt.is_init()) throw QpInitError("map2tod: pointing not initialized");
  if (!source.is_init()) throw QpInitError("map2tod: map not initialized");

  // built here, before the parallel region, so threads only read it
  if (w.core.opt().interp_pix) source.ensure_pixinfo();

  const Map &global = source.core;
  const PointData pd = pnt.core;
  const std::ptrdiff_t ndet = static_cast<std::ptrdiff_t>(dets.dets.size());

  // No reduction -- detectors write disjoint rows of the tod -- so there is
  // no accumulator to size and nothing stopping the balanced schedule.
  parallel_dets(
      clamp_threads(w.num_threads, ndet), ndet, "map2tod", true,
      [&] { return w.core; },
      [&](Pointing &local, std::ptrdiff_t i) {
        map2tod1(local, dets.dets[i], pd, global);
      },
      [](Pointing &) {});
}

}  // namespace

PYBIND11_MODULE(_libqpoint2, m) {
  m.doc() = "qpoint2: C++ core with vectorizing, zero-copy pybind11 bindings";

  // Whether this build can actually run the detector loops in parallel.
  // False on macOS, where linking an OpenMP runtime collides with the one
  // healpy bundles -- see the revert of the Apple-clang fallback.
#ifdef _OPENMP
  m.attr("HAS_OPENMP") = true;
#else
  m.attr("HAS_OPENMP") = false;
#endif

  // Derived from RuntimeError so that `except RuntimeError` keeps working:
  // that is what the C-backed packages raise, and what qmap_class already
  // raises from Python.
  auto base = py::register_exception<QpError>(m, "QpError", PyExc_RuntimeError);
  py::register_exception<QpInitError>(m, "QpInitError", base);
  py::register_exception<QpPointError>(m, "QpPointError", base);
  py::register_exception<QpMapError>(m, "QpMapError", base);

  m.def("det_offset", &py_det_offset, py::arg("delta_az"), py::arg("delta_el"),
        py::arg("delta_psi"),
           R"doc(
Quaternion for the requested detector centroid offset from boresight.

Arguments
---------
delta_az : array_like
    Azimuthal centroid offset of the detector in degrees
delta_el : array_like
    Elevation centroid offset of the detector in degrees
delta_psi : array_like
    Polarization offset of the detector from vertical in degrees

Returns
-------
q : array_like
    Detector centroid offset quaternion for each detector
)doc");
  m.def("refraction", &py_refraction, py::arg("el"), py::arg("temp"),
        py::arg("press"), py::arg("hum"), py::arg("freq"),
           R"doc(
Atmospheric refraction correction, in degrees.

Stateless: it reads none of the stored weather and updates none of it,
which is what distinguishes it from the QPoint.refraction method.

Arguments
---------
el : array_like
    Observer elevation angle in degrees
temp : array_like
    Ambient temperature in Celcius
press : array_like
    Ambient pressure in mbar
hum : array_like
    Relative humidity, as a fraction
freq : array_like
    Observing frequency in GHz

Returns
-------
delta : array_like
    Refraction correction in degrees, for each elevation
)doc");
  m.def("hwp_quat", &py_hwp_quat, py::arg("theta"),
           R"doc(
Quaternion for rotation by 2*theta (physical HWP angle).

Arguments
---------
theta : array_like
    HWP physical angle in degrees

Returns
-------
q : array_like
    Quaternion for each hwp angle
)doc");

  py::class_<PointingWrap>(m, "Pointing")
      .def(py::init<>())
      .def("set_param", &set_param, py::arg("key"), py::arg("val"),
           R"doc(
Set one parameter by name.

Arguments
---------
key : str
    Name of a rate, option, weather or double-valued parameter.
val : object
    The value. Rates take 'never', 'once', 'always' or an interval in
    seconds; the string-valued options take their own spellings.
)doc")
      .def("get_param", &get_param, py::arg("key"),
           R"doc(
Return one parameter by name.

Arguments
---------
key : str
    Name of the parameter.

Returns
-------
val : object
    Its value, in the same spelling set accepts.
)doc")
      .def("reset_rates", [](PointingWrap &w) { w.core.reset_rates(); },
           R"doc(
Reset every forward correction cache, so each is recomputed.
)doc")
      .def("reset_inv_rates", [](PointingWrap &w) { w.core.reset_inv_rates(); },
           R"doc(
Reset every inverse correction cache.
)doc")
      .def("reset_rate",
           [](PointingWrap &w, const std::string &n) { w.core.reset_rate(n); },
           R"doc(
Reset one correction cache, leaving its rate alone.

Arguments
---------
name : str
    Correction name, without the rate_ prefix.
)doc")
      .def("print_memory",
           [](PointingWrap &w) {
             py::print("qpoint2 Pointing:");
             for (const auto &d : kRates) {
               const std::string name(d.name);
               py::print("  rate", name, w.core.get_rate(name), "| inv",
                         w.core.get_rate(name + "_inv"));
             }
             for (const auto &d : kOptionParams)
               py::print("  opt", std::string(d.name), w.core.get_opt(d.name));
             py::print("  num_threads", w.num_threads);
             for (const auto &d : kWeatherParams)
               py::print("  weather", std::string(d.name),
                         w.core.get_weather_param(d.name));
             py::print("  ref_delta", w.core.get_double("ref_delta"));
             py::print("  dut1", w.core.get_double("dut1"));
           },
           R"doc(
Print the correction state and options to stdout.
)doc")
      .def("azel2radec", &py_azel2radec, py::arg("delta_az"),
           py::arg("delta_el"), py::arg("delta_psi"), py::arg("az"),
           py::arg("el"), py::arg("psi"), py::arg("pitch"), py::arg("roll"),
           py::arg("lon"), py::arg("lat"), py::arg("ctime"), py::arg("hwp"),
           py::arg("sindec"), py::arg("return_pa"),
           R"doc(
Sky coordinates from az/el boresight and a detector offset.

Arguments
---------
delta_az : float
    Azimuthal offset of the detector in degrees
delta_el : float
    Elevation offset of the detector in degrees
delta_psi : float
    Polarization offset of the detector in degrees
az : array_like
    Boresight azimuth in degrees
el : array_like
    Boresight elevation in degrees
pitch : array_like
    Boresight pitch in degrees.  If None, this term is ignored.
roll : array_like
    Boresight roll in degrees.  If None, this term is ignored.
lon : array_like
    Observer longitude in degrees.
lat : array_like
    Observer latitude in degrees.
ctime : array_like
    Unix time in seconds UTC
hwp : array_like, optional
    HWP angles in degrees
sindec : bool, optional
    If `True`, return sin(dec) instead of dec in degrees (default False)
return_pa : bool, optional
    If `True`, return pa instead of sin2psi/cos2psi
psi : array_like
    See :meth:`qpoint2.QPoint.azel2radec`.

Returns
-------
ra : array_like
    Detector right ascension in degrees
dec/sindec : array_like
    Detector declination in degrees
pa : array_like
    Detector position angle, if `return_pa` is True
sin2psi : array_like
    Detector polarization orientation, if `return_pa` is False
cos2psi : array_like
    Detector polarization orientation, if `return_pa` is False
)doc")
      .def("gmst", &py_gmst, py::arg("ctime"),
           R"doc(
Greenwich mean sidereal time, in hours.

Arguments
---------
ctime : array_like
    Unix time in seconds UTC

Returns
-------
gmst : array_like
    Greenwich mean sidereal time of the observer
)doc")
      .def("lmst", &py_lmst, py::arg("ctime"), py::arg("lon"),
           R"doc(
Local mean sidereal time, in hours.

Arguments
---------
ctime : array_like
    Unix time in seconds UTC
lon : array_like
    Observer longitude (degrees)

Returns
-------
lmst : array_like
    Local mean sidereal time of the observer
)doc")
      .def("azel2bore", &py_azel2bore, py::arg("az"), py::arg("el"),
           py::arg("psi"), py::arg("pitch"), py::arg("roll"), py::arg("lon"),
           py::arg("lat"), py::arg("ctime"),
           R"doc(
Boresight quaternion from az/el/pitch/roll/lon/lat/ctime.

Arguments
---------
az : array_like
    Boresight azimuth in degrees
el : array_like
    Boresight elevation in degrees
pitch : array_like
    Boresight pitch in degrees.  If `None`, this term is ignored.
roll : array_like
    Boresight roll in degrees.  If `None`, this term is ignored.
lon : array_like
    Observer longitude in degrees
lat : array_like
    Observer latitude in degrees
ctime : array_like
    Unix time in seconds UTC
psi : array_like
    See :meth:`qpoint2.QPoint.azel2bore`.

Returns
-------
q : array_like
    Nx4 numpy array of quaternions for each supplied timestamp.
)doc")
      .def("bore2radec", &py_bore2radec, py::arg("q_off"), py::arg("ctime"),
           py::arg("q_bore"), py::arg("q_hwp"), py::arg("sindec"),
           py::arg("return_pa"),
           R"doc(
Sky coordinates for a detector offset and a boresight timestream.

Arguments
---------
q_off : quaternion
    Detector offset quaternion for a single detector, calculated using
    :meth:`det_offset`.
ctime : array_like
    Unix time in seconds UTC, broadcastable to shape (N,),
    the long dimension of `q_bore`.
q_bore : quaternion or array of quaternions
    Nx4 array of quaternions encoding the boresight orientation on the
    sky (as output by :meth:`azel2radec` or equivalent)
q_hwp : quaternion or array of quaternions, optional
    HWP angle quaternions calculated using :meth:`hwp_quat`.  Must be
    broadcastable to the same shape as `q_bore`.
sindec : bool, optional
    If `True`, return sin(dec) instead of dec in degrees
    (default False).
return_pa : bool, optional
    If `True`, return pa instead of sin2psi / cos2psi

Returns
-------
ra : array_like
    Detector right ascension in degrees
dec/sindec : array_like
    Detector declination in degrees or sin(dec) if `sindec` is `True`.
pa/sin2psi : array_like
    Detector polarization orientation if `return_pa` is `True`, or
    sin(2*pa) if `return_pa` is `False`.
cos2psi : array_like
    detector polarization orientation cos(2*pa), if `return_pa` is `False`.
)doc")
      .def("bore2azel", &py_bore2azel, py::arg("q_bore"), py::arg("lon"),
           py::arg("lat"), py::arg("ctime"),
           R"doc(
Horizon coordinates from a boresight quaternion timestream.

Arguments
---------
q_bore : array_like
    Nx4 array of boresight quaternions (as output by :meth:`azel2bore`).
lon : array_like
    Observer longitude in degrees.
lat : array_like
    Observer latitude in degrees.
ctime : array_like
    Unix time in seconds UTC

Returns
-------
az : array_like
    Azimuth in degrees
el : array_like
    Elevation in degrees
pa : array_like
    Position angle in horizon coordinates
)doc")
      .def("radec2azel", &py_radec2azel, py::arg("ra"), py::arg("dec"),
           py::arg("pa"), py::arg("lon"), py::arg("lat"), py::arg("ctime"),
           R"doc(
Horizon coordinates from equatorial coordinates.

Arguments
---------
ra : array_like
    Right ascension angle
dec : array_like
    Declination angle
pa : array_like
    Position angle in equatorial coordinates
lon : array_like
    Observer longitude in degrees.
lat : array_like
    Observer latitude in degrees.
ctime : array_like
    Unix time in seconds UTC

Returns
-------
az : array_like
    Azimuth in degrees
el : array_like
    Elevation in degrees
hpa : array_like
    Position angle in horizon coordinates
)doc")
      .def("radecpa2quat", &py_radecpa2quat, py::arg("ra"), py::arg("dec"),
           py::arg("pa"),
           R"doc(
Quaternion from ra/dec/position angle, in degrees.

Arguments
---------
ra : array_like
    Right ascension angle
dec : array_like
    Declination angle
pa : array_like
    Position angle

Returns
-------
q : array_like
    Quaternion constructed from the input angles.
)doc")
      .def("quat2radecpa", &py_quat2radecpa, py::arg("quat"),
           R"doc(
ra/dec/position angle from a quaternion, in degrees.

Arguments
---------
quat : quaternion or array of quaternions
    Orientation quaternions, of shape (N, 4).

Returns
-------
ra : array_like
    Right ascension in degrees.
dec : array_like
    Declination in degrees.
pa : array_like
    Position angle in degrees.
)doc")
      .def("update_ref", &py_update_ref, py::arg("q"),
           R"doc(
Compute and store the refraction correction for an orientation.

Arguments
---------
q : array_like
    Boresight quaternion, of shape (4,).

Returns
-------
delta : float
    The correction in degrees, also stored as ref_delta.
)doc")
      .def("set_bulletin_a", &py_set_bulletin_a, py::arg("mjd_min"),
           py::arg("mjd_max"), py::arg("dut1"), py::arg("x"), py::arg("y"),
           R"doc(
Load an IERS Bulletin A table.

Arguments
---------
mjd_min : int
mjd_max : int
    Inclusive range of dates the table covers.
dut1 : array_like
    UT1 - UTC in seconds, one per day.
x : array_like
y : array_like
    Polar motion in arcseconds, one per day.
)doc")
      .def("get_bulletin_a", &py_get_bulletin_a, py::arg("mjd"),
           R"doc(
Interpolated (dut1, x, y) from the loaded Bulletin A table.

Arguments
---------
mjd : array_like
    Modified Julian date, of any shape.

Returns
-------
dut1 : array_like
    UT1 - UTC in seconds, interpolated to `mjd`.
x : array_like
y : array_like
    Polar motion in arcseconds, interpolated to `mjd`.

Notes
-----
A date outside the loaded table returns zeros rather than raising,
which is what the C does and what callers who never loaded a bulletin
rely on. Zeros are also what the parameters mean when no correction is
applied.
)doc")
      .def("radec2pix", &py_radec2pix, py::arg("ra"), py::arg("dec"),
           py::arg("nside"),
           R"doc(
HEALPix pixel numbers for the given sky coordinates.

Arguments
---------
ra : array_like
    Right ascension angle
dec : array_like
    Declination angle
nside : int
    HEALpix resolution parameter

Returns
-------
pix : array_like
    Pixel number(s) corresponding to the input positions(s).
)doc")
      .def("quat2pix", &py_quat2pix, py::arg("quat"), py::arg("nside"),
           py::arg("return_pa"),
           R"doc(
Pixel number and polarization angle for a quaternion.

Arguments
---------
quat : quaternion or array of quaternions
    Pointing orientation(s)
nside : int, optional
    HEALpix resolution parameter
return_pa : array_like
    See :meth:`qpoint2.QPoint.quat2pix`.

Returns
-------
pix : array_like
    Pixel number(s) for the given input quaternion(s)
sin2psi : array_like
cos2psi : array_like
    Polarization coefficients, if `pol` is `True`.
)doc")
      .def("bore2pix", &py_bore2pix, py::arg("q_off"), py::arg("ctime"),
           py::arg("q_bore"), py::arg("q_hwp"), py::arg("nside"),
           py::arg("return_pa"),
           R"doc(
Pixel and polarization timestreams for a detector offset.

Arguments
---------
q_off : quaternion
    Detector offset quaternion for a single detector,
    calculated using :meth:`det_offset`.
ctime : array_like
    Unix times in seconds UTC, broadcastable to shape (N,),
    the long dimenions of `q_bore`.
q_bore : quaternion or array of quaternions
    Nx4 array of quaternions encoding the boresight orientation on the
    sky (as output by :meth:`azel2radec` or equivalent)
q_hwp : quaternion or array of quaternions, optional
    HWP angle quaternions calculated using :meth:`hwp_quat`.  Must be
    broadcastable to the same shape as `q_bore`.
nside : int, optional
    HEALpix map dimension.  Default: 256.
return_pa : bool, optional
    If `True`, return pa instead of sin2psi / cos2psi

Returns
-------
pix : array_like
    Detector pixel number
pa/sin2psi : array_like
    Detector polarization orientation if `return_pa` is `True`, or
    sin(2*pa) if `return_pa` is `False`.
cos2psi : array_like
    detector polarization orientation cos(2*pa), if `return_pa` is `False`.
)doc")
      .def("dipole", &py_dipole, py::arg("ctime"), py::arg("ra"),
           py::arg("dec"),
           R"doc(
CMB dipole amplitude in the given equatorial direction, in K.

Arguments
---------
ctime : array_like
    Unix time in seconds UTC
ra : array_like
    Right ascension on the sky, in degrees.
dec : array_like
    Declination on the sky, in degrees

Returns
-------
dipole : array_like
    Dipole amplitude in K
)doc")
      .def("bore2dipole", &py_bore2dipole, py::arg("q_off"), py::arg("ctime"),
           py::arg("q_bore"),
           R"doc(
CMB dipole timestream for a detector offset and boresight.

Arguments
---------
q_off : quaternion
    Detector offset quaternion for a single detector, calculated using
    :meth:`det_offset`
ctime : array_like
    Array of unix times in seconds UTC
q_bore : quaternion or array of quaternions
    Array of quaternions encoding the boresight orientation on the sky
    (as output by :meth:`azel2radec` or similar).  Broadcastable to the
    same length as `ctime`.

Returns
-------
dipole : array_like
    Dipole amplitude in K
)doc")
      .def("bore_offset", &py_bore_offset, py::arg("q_bore"), py::arg("ang1"),
           py::arg("ang2"), py::arg("ang3"), py::arg("post"),
           R"doc(
Apply a fixed or per-sample offset to a boresight quaternion.

Arguments
---------
q_bore : array_like
    boresight pointing quaternion
ang1 : array_like, optional
    Azimuthal or ra offset in degrees
ang2 : array_like, optional
    Elevation or dec offset in degrees
ang3 : array_like, optional
    Position angle offset in degrees
post : bool, optional
    If False, apply offset as an az/el/pa pre-rotation
    If True, apply offset as an ra/dec/pa post-rotation

Returns
-------
q_bore : array_like
    Offset boresight quaternion
)doc")
      // a method on Pointing for API parity, though it uses no state
      .def(
          "omega2azelpsi",
          [](PointingWrap &, double init_az, double init_el, double init_psi,
             py::object ox, py::object oy, py::object oz, double dt) {
            return py_omega2azelpsi(init_az, init_el, init_psi, ox, oy, oz, dt);
          },
          py::arg("init_az"), py::arg("init_el"), py::arg("init_psi"),
          py::arg("omega_x"), py::arg("omega_y"), py::arg("omega_z"),
          py::arg("dt"),
           R"doc(
Integrate gyro rates into an az/el/psi attitude timestream.

Arguments
---------
init_az : float
    Initial azimuth in degrees.
init_el : float
    Initial elevation in degrees.
init_psi : float
    Initial rotation about the boresight in degrees.
omega_x : array_like
omega_y : array_like
omega_z : array_like
    Body-frame angular rates in degrees per second, of shape (N,).
dt : array_like
    See :meth:`qpoint2.QPoint.omega2azelpsi`.

Returns
-------
az : array_like
el : array_like
psi : array_like
    The integrated attitude in degrees, of shape (N,).
)doc")
      .def("get_interp_val", &py_get_interp_val, py::arg("nside"),
           py::arg("map"), py::arg("ra"), py::arg("dec"),
           R"doc(
Bilinearly interpolate a map at the given sky coordinates.

Arguments
---------
ra : array_like
    Right ascension in degrees, of shape (N,).
dec : array_like
    Declination in degrees, of shape (N,).
nside : array_like
    See :meth:`qpoint2.QPoint.get_interp_val`.
map : array_like
    See :meth:`qpoint2.QPoint.get_interp_val`.

Returns
-------
val : array_like
    Bilinearly interpolated map values, of shape (N,) for one map or
    (nmap, N) for several. A single sample degrades to a scalar.
)doc")
      .def("rotate_map", &py_rotate_map, py::arg("nside"),
           py::arg("map_in"), py::arg("to_gal"),
           R"doc(
Resample a polarized (3, npix) map into another coordinate frame.

Arguments
---------
map_in : array_like
    Input map, of shape (3, N)
nside : array_like
    See :meth:`qpoint2.QPoint.rotate_map`.
to_gal : array_like
    See :meth:`qpoint2.QPoint.rotate_map`.

Returns
-------
map_out : array_like
    Rotated output map.
)doc")
      .def("rotate_quat", &py_rotate_quat, py::arg("quat"), py::arg("to_gal"),
           R"doc(
Rotate a quaternion between celestial and galactic coordinates.

Arguments
---------
quat : array_like
    array of quaternions, of shape (n, 4)
to_gal : array_like
    See :meth:`qpoint2.QPoint.rotate_quat`.

Returns
-------
quat : array_like
    rotated quaternion array
)doc")
      .def("rotate_coord", &py_rotate_coord, py::arg("ra"), py::arg("dec"),
           py::arg("pa"), py::arg("sin2psi"), py::arg("cos2psi"),
           py::arg("to_gal"),
           R"doc(
Rotate sky coordinates between celestial and galactic.

Arguments
---------
ra : array_like
    Right ascension in degrees, of shape (N,). Rotated in place unless
    `inplace` is False.
dec : array_like
    Declination in degrees, of shape (N,). Rotated in place unless
    `inplace` is False.
pa : array_like, optional
    Position angle in degrees. Supply this or the `sin2psi`/`cos2psi`
    pair, not both.
sin2psi : array_like, optional
    sin(2*pa), if the polarization angle is carried as a pair.
cos2psi : array_like, optional
    cos(2*pa), paired with `sin2psi`.
to_gal : array_like
    See :meth:`qpoint2.QPoint.rotate_coord`.

Returns
-------
ra : array_like
    Rotated right ascension in degrees.
dec : array_like
    Rotated declination in degrees.
pa/sin2psi : array_like
    Rotated position angle, or sin(2*pa) if the pair was supplied.
cos2psi : array_like
    Rotated cos(2*pa), if the pair was supplied.
)doc");

  // Mapmaking takes a Pointing rather than belonging to one, so these are
  // module functions. Binding them onto Pointing put tod2map and map2tod
  // on QPoint, which has no maps to make.
  m.def("tod2map", &py_tod2map, py::arg("pointing"), py::arg("detarr"),
        py::arg("point"), py::arg("map_out"),
        R"doc(
Accumulate timestreams into a destination map.

Arguments
---------
pointing : Pointing
    Correction state to project with.
detarr : QpDetArr
    Detectors, their timestreams and their per-sample data.
point : QpPoint
    Boresight pointing for the same samples.
map_out : QpMap
    Destination, accumulated into rather than overwritten.
)doc");
  m.def("map2tod", &py_map2tod, py::arg("pointing"), py::arg("detarr"),
        py::arg("point"), py::arg("map_in"),
        R"doc(
Scan a source map into timestreams.

Arguments
---------
pointing : Pointing
    Correction state to project with.
detarr : QpDetArr
    Detectors, whose timestreams are accumulated into with +=.
point : QpPoint
    Boresight pointing for the same samples.
map_in : QpMap
    Source map.
)doc");

  py::class_<MapWrap>(m, "QpMap")
      .def(py::init<>())
      .def("setup", &MapWrap::setup, py::arg("vec"), py::arg("proj"),
           py::arg("pixels"), py::arg("nside"), py::arg("pol") = py::none(),
           py::arg("vpol") = py::none(),
           R"doc(
Install the map components and fix their shapes.

Arguments
---------
vec : array_like, None or False
    The signal map, a fresh map of zeros, or no such component.
proj : array_like, None or False
    The projection matrix, likewise.
pixels : array_like, optional
    Pixel numbers for a partial map, of shape (npix,).
nside : int
    HEALPix resolution.
pol : bool, optional
    Whether the map carries polarization. Consulted only where a row
    count leaves the mode open, and overridden by a supplied map that
    settles it. None takes the documented default, T,Q,U.
vpol : bool, optional
    Whether it carries the V component as well. Same reading; None is
    False, so a 4-row map needs this set unless a proj settles it.
)doc")
      .def("update_vec", &MapWrap::update_vec, py::arg("arr") = py::none(),
           R"doc(
Replace the signal map, or switch it off.

Arguments
---------
arr : array_like, None or False
    A map of the installed row count, None for a fresh map of zeros,
    or False to remove the component.
)doc")
      .def("update_proj", &MapWrap::update_proj, py::arg("arr") = py::none(),
           R"doc(
Replace the projection matrix, or switch it off.

Arguments
---------
arr : array_like, None or False
    As for update_vec.
)doc")
      .def("get_vec", &MapWrap::get_vec,
           R"doc(
The installed signal map.

Returns
-------
vec : array_like or None
    The array itself, not a copy.
)doc")
      .def("get_proj", &MapWrap::get_proj,
           R"doc(
The installed projection matrix.

Returns
-------
proj : array_like or None
    The array itself, not a copy.
)doc")
      .def("reset", &MapWrap::reset,
           R"doc(
Release both components.
)doc")
      .def("is_init", &MapWrap::is_init,
           R"doc(
Whether the map has been set up.

Returns
-------
bool
)doc")
      .def("is_partial", &MapWrap::is_partial,
           R"doc(
Whether the map covers a pixel list rather than the sphere.

Returns
-------
bool
)doc")
      .def("has_vec", &MapWrap::has_vec,
           R"doc(
Whether a signal map is installed.

Returns
-------
bool
)doc")
      .def("has_proj", &MapWrap::has_proj,
           R"doc(
Whether a projection matrix is installed.

Returns
-------
bool
)doc")
      .def("is_pol", &MapWrap::is_pol,
           R"doc(
Whether the map is polarized, from either component. The unpolarized
derivative modes D1 and D2 are not, whatever their row count.

Returns
-------
bool
)doc")
      .def("is_vpol", &MapWrap::is_vpol,
           R"doc(
Whether the map carries the V component, from either component.

Returns
-------
bool
)doc");

  py::class_<PointWrap>(m, "QpPoint")
      .def(py::init<>())
      .def("set_bore", &PointWrap::set_bore, py::arg("q_bore"),
           R"doc(
Install the boresight quaternions, which fixes the sample count.

Arguments
---------
q_bore : array_like
    Quaternions, of shape (N, 4).
)doc")
      .def("set_ctime", &PointWrap::set_ctime, py::arg("ctime"),
           R"doc(
Install the timestamps.

Arguments
---------
ctime : array_like
    Unix time in seconds UTC, of shape (N,).
)doc")
      .def("clear_ctime", &PointWrap::clear_ctime,
           R"doc(
Drop the timestamps, leaving the pointing without time.
)doc")
      .def("set_hwp", &PointWrap::set_hwp, py::arg("q_hwp"),
           R"doc(
Install the half-wave plate angles.

Arguments
---------
q_hwp : array_like
    HWP quaternions, of shape (N, 4).
)doc")
      .def("clear_hwp", &PointWrap::clear_hwp,
           R"doc(
Drop the half-wave plate angles.
)doc")
      .def("reset", &PointWrap::reset,
           R"doc(
Release the pointing.
)doc")
      .def("is_init", &PointWrap::is_init,
           R"doc(
Whether the boresight has been installed.

Returns
-------
bool
)doc")
      .def("n_samples", &PointWrap::n_samples,
           R"doc(
Number of samples the installed pointing covers.

Returns
-------
int
)doc");

  py::class_<DetArrWrap>(m, "QpDetArr")
      .def(py::init<>())
      .def("setup", &DetArrWrap::setup, py::arg("q_off"), py::arg("weight"),
           py::arg("gain"), py::arg("mueller"), py::arg("tod"), py::arg("flag"),
           py::arg("weights"), py::arg("n_samp"), py::arg("do_diff"),
           py::arg("write"),
           R"doc(
Install the detectors and the arrays they carry.

Arguments
---------
q_off : array_like
    Offset quaternions, of shape (ndet, 4).
weight : array_like
    Per-detector mapmaking weight, of shape (ndet,).
gain : array_like
    Per-detector gain, of shape (ndet,).
mueller : array_like
    Per-detector polarization efficiency, of shape (ndet, 4).
tod : array_like, optional
    Timestreams, of shape (ndet, n_samp).
flag : array_like, optional
    Per-sample flags, of shape (ndet, n_samp). A flagged sample is
    skipped.
weights : array_like, optional
    Per-sample weights, of shape (ndet, n_samp).
n_samp : int
    Samples per detector, which must match the pointing.
do_diff : bool
    Pair the first half of the detectors with the second half and
    difference them.
write : bool
    Whether the timestreams are written to rather than read.
)doc")
      .def("get_tod", &DetArrWrap::get_tod,
           R"doc(
The installed timestreams.

Returns
-------
tod : array_like or None
    The array itself, not a copy.
)doc")
      .def("reset", &DetArrWrap::reset,
           R"doc(
Release the detectors.
)doc")
      .def("is_init", &DetArrWrap::is_init,
           R"doc(
Whether the detectors have been set up.

Returns
-------
bool
)doc");
}
