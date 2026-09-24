#pragma once

#include <stdexcept>
#include <string>

namespace qp {

// Replaces the error_code / error_string members of qp_memory_t. The four
// types mirror the qp_error_codes enum so the bindings can map them onto
// distinct Python exception types.
class QpError : public std::runtime_error {
 public:
  explicit QpError(const std::string &msg) : std::runtime_error(msg) {}
};

class QpInitError : public QpError {
 public:
  explicit QpInitError(const std::string &msg) : QpError(msg) {}
};

class QpPointError : public QpError {
 public:
  explicit QpPointError(const std::string &msg) : QpError(msg) {}
};

class QpMapError : public QpError {
 public:
  explicit QpMapError(const std::string &msg) : QpError(msg) {}
};

}  // namespace qp
