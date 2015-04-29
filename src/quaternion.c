/*
  Quaternion library
  2012 Mike Nolta <mike@nolta.net>
 */

#include <math.h>
#include "quaternion.h"

#ifndef M_PI
#define M_PI		3.14159265358979323846	// pi
#define M_PI_2		1.57079632679489661923	// pi/2
#endif

#ifndef invsqrt
static double inline invsqrt( double x ) { return 1./sqrt(x); };
#endif

void
Quaternion_inv(Quaternion q)
{
  double norm2 = Quaternion_norm2(q);
  q[0] = q[0]/norm2;
  q[1] = -q[1]/norm2;
  q[2] = -q[2]/norm2;
  q[3] = -q[3]/norm2;
}

void
Quaternion_mul(Quaternion q, const Quaternion a, const Quaternion b)
{
  q[0] = a[0]*b[0] - a[1]*b[1] - a[2]*b[2] - a[3]*b[3];
  q[1] = a[0]*b[1] + a[1]*b[0] + a[2]*b[3] - a[3]*b[2];
  q[2] = a[0]*b[2] - a[1]*b[3] + a[2]*b[0] + a[3]*b[1];
  q[3] = a[0]*b[3] + a[1]*b[2] - a[2]*b[1] + a[3]*b[0];
}

void
Quaternion_mul_left(const Quaternion a, Quaternion q)
{
  Quaternion b;
  Quaternion_copy(b, q);
  Quaternion_mul(q, a, b);
}

void
Quaternion_mul_right(Quaternion q, const Quaternion a)
{
  Quaternion b;
  Quaternion_copy(b, q);
  Quaternion_mul(q, b, a);
}

double
Quaternion_norm(const Quaternion q)
{
  return sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
}

int
Quaternion_rot(Quaternion q, double angle, const double v[3])
{
  double angle_2 = 0.5*angle;
  double s = sin(angle_2);
  double norm = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
  if (norm <= 0.)
    return 1;
  q[0] = cos(angle_2);
  q[1] = s*v[0]/norm;
  q[2] = s*v[1]/norm;
  q[3] = s*v[2]/norm;
  return 0;
}

void
Quaternion_r1(Quaternion q, double angle)
{
  double angle_2 = 0.5*angle;
  q[0] = cos(angle_2);
  q[1] = sin(angle_2);
  q[2] = 0.;
  q[3] = 0.;
}

void
Quaternion_r2(Quaternion q, double angle)
{
  double angle_2 = 0.5*angle;
  q[0] = cos(angle_2);
  q[1] = 0.;
  q[2] = sin(angle_2);
  q[3] = 0.;
}

void
Quaternion_r3(Quaternion q, double angle)
{
  double angle_2 = 0.5*angle;
  q[0] = cos(angle_2);
  q[1] = 0.;
  q[2] = 0.;
  q[3] = sin(angle_2);
}

void
Quaternion_r1_mul(double angle, Quaternion q)
{
  Quaternion a, b;
  Quaternion_r1(a, angle);
  Quaternion_copy(b, q);

  q[0] = a[0]*b[0] - a[1]*b[1]; // - a[2]*b[2] - a[3]*b[3];
  q[1] = a[0]*b[1] + a[1]*b[0]; // + a[2]*b[3] - a[3]*b[2];
  q[2] = a[0]*b[2] - a[1]*b[3]; // + a[2]*b[0] + a[3]*b[1];
  q[3] = a[0]*b[3] + a[1]*b[2]; // - a[2]*b[1] + a[3]*b[0];
}

void
Quaternion_r2_mul(double angle, Quaternion q)
{
  Quaternion b;
  double angle_2 = 0.5*angle;
  double c = cos(angle_2), s = sin(angle_2);
  Quaternion_copy(b, q);

  q[0] = c*b[0] - s*b[2];
  q[1] = c*b[1] + s*b[3];
  q[2] = c*b[2] + s*b[0];
  q[3] = c*b[3] - s*b[1];

  /*
    Quaternion_r2(angle, a);

    q[0] = a[0]*b[0] - a[2]*b[2]; //- a[3]*b[3];
    q[1] = a[0]*b[1] + a[2]*b[3]; //- a[3]*b[2];
    q[2] = a[0]*b[2] + a[2]*b[0]; //+ a[3]*b[1];
    q[3] = a[0]*b[3] - a[2]*b[1]; //+ a[3]*b[0];
  */
}

void
Quaternion_r3_mul(double angle, Quaternion q)
{
  Quaternion a, b;
  Quaternion_r3(a, angle);
  Quaternion_copy(b, q);

  q[0] = a[0]*b[0] /*- a[1]*b[1] - a[2]*b[2]*/ - a[3]*b[3];
  q[1] = a[0]*b[1] /*+ a[1]*b[0] + a[2]*b[3]*/ - a[3]*b[2];
  q[2] = a[0]*b[2] /*- a[1]*b[3] + a[2]*b[0]*/ + a[3]*b[1];
  q[3] = a[0]*b[3] /*+ a[1]*b[2] - a[2]*b[1]*/ + a[3]*b[0];
}

void
Quaternion_to_matrix(const Quaternion q, double mat[3][3])
{
  Quaternion u;
  Quaternion_copy(u, q);
  Quaternion_unit(u);

  double a2 = u[0]*u[0], b2 = u[1]*u[1], c2 = u[2]*u[2], d2 = u[3]*u[3];
  mat[0][0] = a2 + b2 - c2 - d2;
  mat[1][1] = a2 - b2 + c2 - d2;
  mat[2][2] = a2 - b2 - c2 + d2;
  mat[0][1] = 2.*(u[1]*u[2] - u[0]*u[3]);
  mat[0][2] = 2.*(u[1]*u[3] + u[0]*u[2]);
  mat[1][2] = 2.*(u[2]*u[3] - u[0]*u[1]);
  mat[1][0] = 2.*(u[1]*u[2] + u[0]*u[3]);
  mat[2][0] = 2.*(u[1]*u[3] - u[0]*u[2]);
  mat[2][1] = 2.*(u[2]*u[3] + u[0]*u[1]);
}

void
Quaternion_unit(Quaternion q)
{
  double norm2 = Quaternion_norm2(q);
  double invnorm = invsqrt(norm2);
  Quaternion_scale(q, invnorm);
}

void
QuaternionSlerp_init(QuaternionSlerp *slerp, const Quaternion a, const Quaternion b)
{
  double cos_alpha = a[0]*b[0] + a[1]*b[1] + a[2]*b[2] + a[3]*b[3];
  slerp->sin_alpha = sqrt(1. - cos_alpha*cos_alpha);
  Quaternion_copy(slerp->q0, a);
  Quaternion_copy(slerp->q1, b);

  if (cos_alpha < 0.) {
    slerp->alpha = acos(-cos_alpha);
    slerp->q1[0] = -slerp->q1[0];
    slerp->q1[1] = -slerp->q1[1];
    slerp->q1[2] = -slerp->q1[2];
    slerp->q1[3] = -slerp->q1[3];
  } else {
    slerp->alpha = acos(cos_alpha);
  }
}

void
QuaternionSlerp_interpolate(const QuaternionSlerp *slerp, double t, Quaternion q)
{
  double s0 = sin((1.-t)*slerp->alpha)/slerp->sin_alpha;
  double s1 = sin(t*slerp->alpha)/slerp->sin_alpha;
  for (int i = 0; i != 4; ++i)
    q[i] = s0*slerp->q0[i] + s1*slerp->q1[i];
}

