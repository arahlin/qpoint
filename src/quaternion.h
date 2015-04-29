/*
  Quaternion library
  2012 Mike Nolta <mike@nolta.net>
*/


#pragma once

#ifdef __cplusplus
extern "C" {
#endif

  typedef double Quaternion[4];

  // q = a + b
  static inline void
  Quaternion_add(Quaternion q, const Quaternion a, const Quaternion b)
  {
    q[0] = a[0] + b[0];
    q[1] = a[1] + b[1];
    q[2] = a[2] + b[2];
    q[3] = a[3] + b[3];
  }

  // q = q*
  static inline void
  Quaternion_conj(Quaternion q)
  {
    q[1] = -q[1];
    q[2] = -q[2];
    q[3] = -q[3];
  }

  // q = a
  static inline void
  Quaternion_copy(Quaternion q, const Quaternion a)
  {
    q[0] = a[0];
    q[1] = a[1];
    q[2] = a[2];
    q[3] = a[3];
  }

  // q = (1,0,0,0)
  static inline void
  Quaternion_identity(Quaternion q)
  {
    q[0] = 1.;
    q[1] = 0.;
    q[2] = 0.;
    q[3] = 0.;
  }

  // q = 1/q = q*/|q|^2
  void Quaternion_inv(Quaternion q);

  // q = a * b
  void Quaternion_mul(Quaternion q, const Quaternion a, const Quaternion b);

  // q = a * q
  void Quaternion_mul_left(const Quaternion a, Quaternion q);

  // q = q * a
  void Quaternion_mul_right(Quaternion q, const Quaternion a);

  // q = {w,x,y,z}
  static inline void
  Quaternion_new(Quaternion q, double w, double x, double y, double z)
  {
    q[0] = w;
    q[1] = x;
    q[2] = y;
    q[3] = z;
  }

  // |q|
  double Quaternion_norm(const Quaternion q);

  // |q|^2
  static inline double
  Quaternion_norm2(const Quaternion q)
  {
    return q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3];
  }

  // q = (rotate by angle around arbitrary vector v)
  int Quaternion_rot(Quaternion q, double angle, const double v[3]);

  // q = R_i(angle)
  void Quaternion_r1(Quaternion q, double angle);
  void Quaternion_r2(Quaternion q, double angle);
  void Quaternion_r3(Quaternion q, double angle);

  // q = R_i(angle) * q
  void Quaternion_r1_mul(double angle, Quaternion q);
  void Quaternion_r2_mul(double angle, Quaternion q);
  void Quaternion_r3_mul(double angle, Quaternion q);

  // q = q * scale
  static inline void
  Quaternion_scale(Quaternion q, double scale)
  {
    q[0] *= scale;
    q[1] *= scale;
    q[2] *= scale;
    q[3] *= scale;
  }

  // q = a - b
  static inline void
  Quaternion_sub(Quaternion q, const Quaternion a, const Quaternion b)
  {
    q[0] = a[0] - b[0];
    q[1] = a[1] - b[1];
    q[2] = a[2] - b[2];
    q[3] = a[3] - b[3];
  }

  void Quaternion_to_matrix(const Quaternion q, double mat[3][3]);

  // same as Quaternion_to_matrix, but only returns single column
  // note: does not normalize quaternion before converting

  static inline void
  Quaternion_to_matrix_col1(const Quaternion u, double col1[3])
  {
    // make sure quaternion is normalized before calling
    double a2 = u[0]*u[0], b2 = u[1]*u[1], c2 = u[2]*u[2], d2 = u[3]*u[3];
    col1[0] = a2 + b2 - c2 - d2;
    col1[1] = 2.*(u[1]*u[2] + u[0]*u[3]);
    col1[2] = 2.*(u[1]*u[3] - u[0]*u[2]);
  }

  static inline void
  Quaternion_to_matrix_col2(const Quaternion u, double col2[3])
  {
    // make sure quaternion is normalized before calling
    double a2 = u[0]*u[0], b2 = u[1]*u[1], c2 = u[2]*u[2], d2 = u[3]*u[3];
    col2[0] = 2.*(u[1]*u[2] - u[0]*u[3]);
    col2[1] = a2 - b2 + c2 - d2;
    col2[2] = 2.*(u[2]*u[3] + u[0]*u[1]);
  }

  static inline void
  Quaternion_to_matrix_col3(const Quaternion u, double col3[3])
  {
    // make sure quaternion is normalized before calling
    double a2 = u[0]*u[0], b2 = u[1]*u[1], c2 = u[2]*u[2], d2 = u[3]*u[3];
    col3[0] = 2.*(u[1]*u[3] + u[0]*u[2]);
    col3[1] = 2.*(u[2]*u[3] - u[0]*u[1]);
    col3[2] = a2 - b2 - c2 + d2;
  }

  void Quaternion_unit(Quaternion q);

  //

  typedef struct
  {
    Quaternion q0, q1;
    double alpha, sin_alpha;
  }
    QuaternionSlerp;

  void
  QuaternionSlerp_init(QuaternionSlerp *slerp, const Quaternion a, const Quaternion b);

  void
  QuaternionSlerp_interpolate(const QuaternionSlerp *slerp, double t, Quaternion z);

#ifdef __cplusplus
}
#endif

