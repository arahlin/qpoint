#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <string>

#include "cpp/error.hpp"
#include "cpp/pointing.hpp"
#include "cpp/quat.hpp"
#include "cpp/span.hpp"

namespace py = pybind11;
using namespace qp;

namespace {

// ---------------------------------------------------------------------------
// Checked array accessors
//
// Array arguments are declared as untyped py::array and validated here rather
// than as py::array_t<double, forcecast>. Both of the convenient alternatives
// copy silently: forcecast copies on any dtype or layout mismatch, and plain
// py::array_t<double> still routes through PyArray_FromAny, so a safe cast
// such as float32 -> float64 also copies. A copy is a perf trap for a plain
// input and a correctness trap for anything persistent or in-place, where the
// caller would never see the buffer that was actually written.
//
// Coercion stays in Python, where qpoint already does it.
// ---------------------------------------------------------------------------

std::string shape_str(const py::array &a) {
  std::string s = "(";
  for (py::ssize_t i = 0; i < a.ndim(); ++i) {
    if (i) s += ", ";
    s += std::to_string(a.shape(i));
  }
  return s + ")";
}

void require_double(const py::array &a, const char *name) {
  if (!a.dtype().is(py::dtype::of<double>()))
    throw py::type_error(std::string(name) + " must be float64, got " +
                         py::str(a.dtype()).cast<std::string>());
  if (!(a.flags() & py::array::c_style))
    throw py::value_error(std::string(name) + " must be C-contiguous");
  if (!(a.flags() & py::detail::npy_api::NPY_ARRAY_ALIGNED_))
    throw py::value_error(std::string(name) + " must be aligned");
}

py::ssize_t length_of(const py::array &a, const char *name) {
  require_double(a, name);
  if (a.ndim() != 1)
    throw py::value_error(std::string(name) + " must be 1-dimensional, got " +
                          shape_str(a));
  return a.shape(0);
}

Span<const double> in_dbl(const py::array &a, const char *name,
                          py::ssize_t n = -1) {
  require_double(a, name);
  const py::ssize_t sz = a.size();
  if (n >= 0 && sz != n)
    throw py::value_error(std::string(name) + " has length " +
                          std::to_string(sz) + ", expected " +
                          std::to_string(n));
  return {static_cast<const double *>(a.data()), static_cast<size_t>(sz)};
}

// Absent optional array -> empty span, which the core reads as "not supplied".
Span<const double> in_dbl_opt(const py::object &o, const char *name,
                              py::ssize_t n) {
  if (o.is_none()) return {};
  return in_dbl(py::cast<py::array>(o), name, n);
}

// View an (n,4) buffer as Quats. Sound because Quat is standard-layout, four
// contiguous doubles, and numpy guarantees the alignment we checked above.
Span<const Quat> in_quat(const py::array &a, const char *name,
                         py::ssize_t n = -1) {
  require_double(a, name);
  if (a.ndim() != 2 || a.shape(1) != 4)
    throw py::value_error(std::string(name) +
                          " must have shape (n, 4), got " + shape_str(a));
  if (n >= 0 && a.shape(0) != n)
    throw py::value_error(std::string(name) + " has length " +
                          std::to_string(a.shape(0)) + ", expected " +
                          std::to_string(n));
  return {reinterpret_cast<const Quat *>(a.data()),
          static_cast<size_t>(a.shape(0))};
}

// A single (4,) quaternion, by value.
Quat load_quat(const py::array &a, const char *name) {
  require_double(a, name);
  if (a.size() != 4)
    throw py::value_error(std::string(name) + " must have 4 elements, got " +
                          shape_str(a));
  return Quat::load(static_cast<const double *>(a.data()));
}

// Match qpoint's output convention: a single sample degrades to a bare (4,).
py::object maybe_scalar_quat(py::array_t<double> q) {
  if (q.shape(0) == 1) return q[py::int_(0)];
  return std::move(q);
}

// ---------------------------------------------------------------------------
// Vectorized entry points
//
// The n-loop lives here, not in the core: core methods are scalar. The loop
// is serial because the correction rate-cache is sequential across samples,
// and it runs with the GIL released.
// ---------------------------------------------------------------------------

py::object py_det_offset(py::array daz, py::array del, py::array dpsi) {
  const py::ssize_t n = length_of(daz, "delta_az");
  auto vaz = in_dbl(daz, "delta_az", n);
  auto vel = in_dbl(del, "delta_el", n);
  auto vpsi = in_dbl(dpsi, "delta_psi", n);

  py::array_t<double> out({n, py::ssize_t(4)});
  double *op = out.mutable_data();
  {
    py::gil_scoped_release nogil;
    for (py::ssize_t i = 0; i < n; ++i)
      det_offset(vaz[i], vel[i], vpsi[i]).store(op + 4 * i);
  }
  return maybe_scalar_quat(std::move(out));
}

py::object py_hwp_quat(py::array theta) {
  const py::ssize_t n = length_of(theta, "theta");
  auto vt = in_dbl(theta, "theta", n);

  py::array_t<double> out({n, py::ssize_t(4)});
  double *op = out.mutable_data();
  {
    py::gil_scoped_release nogil;
    for (py::ssize_t i = 0; i < n; ++i) hwp_quat(vt[i]).store(op + 4 * i);
  }
  return maybe_scalar_quat(std::move(out));
}

}  // namespace

PYBIND11_MODULE(_libqpoint2, m) {
  m.doc() = "qpoint2: C++ core with vectorizing, zero-copy pybind11 bindings";

  py::register_exception<QpError>(m, "QpError");

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
}
