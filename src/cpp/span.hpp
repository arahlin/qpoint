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

}  // namespace qp
