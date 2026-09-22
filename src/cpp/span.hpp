#pragma once

#include <cstddef>

namespace qp {

// Non-owning view over memory owned elsewhere -- in practice always a numpy
// buffer held alive by a binding wrapper. An empty span is the "argument not
// supplied" sentinel, replacing the NULL-pointer tests in the C library.
template <class T>
struct Span {
  T *ptr = nullptr;
  std::size_t len = 0;

  constexpr Span() = default;
  constexpr Span(T *p, std::size_t n) : ptr(p), len(n) {}

  constexpr std::size_t size() const { return len; }
  constexpr explicit operator bool() const { return ptr != nullptr; }

  constexpr T &operator[](std::size_t i) const { return ptr[i]; }
};

// A per-sample input column: either a real array (stride 1), a single array
// element repeated across every sample (stride 0), or a stored scalar with no
// array behind it at all.
//
// Stride 0 is what makes scalar arguments free. The Python API lets any
// argument be a scalar standing in for a full-length column, and materializing
// those -- np.broadcast_arrays followed by a contiguity copy -- would allocate
// a full-length array per scalar argument on every call. Reading the same
// element n times costs nothing and allocates nothing.
//
template <class T>
struct Col {
  const T *ptr = nullptr;
  std::size_t stride = 0;
  T value{};
  bool present = false;

  // By value: T is only ever double or Quat, and returning a reference to
  // the stored `value` made a temporary Col a dangling read waiting to
  // happen.
  constexpr T operator[](std::size_t i) const {
    return ptr ? ptr[i * stride] : value;
  }

  // The same read without the copy, for a Quat, which is four doubles and
  // costs about 4% of bore2radec to copy once per sample. Safe on a named
  // Col -- which is what ColSet hands out, and what every loop holds --
  // and not on a temporary, that being why operator[] returns by value.
  constexpr const T &ref(std::size_t i) const {
    return ptr ? ptr[i * stride] : value;
  }
  // Whether the caller supplied anything. Tracked separately from ptr,
  // because a supplied scalar also has no array behind it.
  constexpr explicit operator bool() const { return present; }
};

}  // namespace qp
