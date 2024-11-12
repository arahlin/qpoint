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


// Standard vector rotation using quaternions. 32 multiplications.
void 
Quaternion_rot_vector(Quaternion q, const double v[3], double v_rot[3]) 
{
    Quaternion vector_quat = {0, v[0], v[1], v[2]};
    Quaternion temp, rotated_quat;

    Quaternion_mul(temp, q, vector_quat);
    Quaternion_conj(q);
    Quaternion_mul(rotated_quat, temp, q);
    Quaternion_conj(q); // Restore original quaternion

    v_rot[0] = rotated_quat[1];
    v_rot[1] = rotated_quat[2];
    v_rot[2] = rotated_quat[3];
}

/* Vector rotation using quaternions with fewer operations. 15 multiplcations.
this computes v' = v + q_0 x t + (q_v  x t)
where t = 2 x (q_v x v_v) */
void 
Quaternion_rot_vector_fast(Quaternion q, const double v[3], double v_rot[3])
{
    double t[3];
    t[0] = 2.0 * (q[1] * v[2] - q[3] * v[1]);
    t[1] = 2.0 * (q[2] * v[0] - q[1] * v[2]);
    t[2] = 2.0 * (q[3] * v[0] - q[2] * v[1]);
    v_rot[0] = v[0] + q[0] * t[0] + (q[2] * t[2] - q[3] * t[1]);
    v_rot[1] = v[1] + q[0] * t[1] + (q[3] * t[0] - q[1] * t[2]);
    v_rot[2] = v[2] + q[0] * t[2] + (q[1] * t[1] - q[2] * t[0]);
}


void 
apply_angular_velocity(Quaternion attitude, double omega_x, double omega_y, double omega_z, double delta_t) 
{
    // Compute the magnitude of the angular velocity vector
    double omega_mag = sqrt(omega_x * omega_x + omega_y * omega_y + omega_z * omega_z);
    double angle = omega_mag * delta_t;
    
    if (omega_mag == 0.0 || angle == 0.0) {
        return; // No change to attitude, so return early
    }

    // Normalize the omega vector to get the rotation axis
    double axis[3] = {omega_x / omega_mag, omega_y / omega_mag, omega_z / omega_mag};
    // printf("Axis: [%f, %f, %f]\n", axis[0], axis[1], axis[2]);
    Quaternion q_delta;
    Quaternion_rot(q_delta, angle, axis);
    // printf("Quaternion delta: [%f, %f, %f, %f]\n", q_delta[0], q_delta[1], q_delta[2], q_delta[3]);

    // Corrected: Right-multiply q_delta to attitude (attitude = attitude * q_delta)
    Quaternion temp;
    Quaternion_mul(temp, attitude, q_delta);
    attitude[0] = temp[0];
    attitude[1] = temp[1];
    attitude[2] = temp[2];
    attitude[3] = temp[3];
    // Normalize the updated attitude quaternion to prevent drift
    double norm = Quaternion_norm(attitude);
    attitude[0] /= norm;
    attitude[1] /= norm;
    attitude[2] /= norm;
    attitude[3] /= norm;
    // printf("Updated quaternion: [%f, %f, %f, %f]\n", attitude[0], attitude[1], attitude[2], attitude[3]);
}


