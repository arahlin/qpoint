/*
  Vector math
  2012 Mike Nolta <mike@nolta.net>
 */

#pragma once

#include <math.h> 

#ifdef __cplusplus
extern "C" {
#endif

#ifndef M_PI
#define M_PI		3.14159265358979323846	// pi
#define M_PI_2		1.57079632679489661923	// pi/2
#endif

#ifndef invsqrt
static double inline invsqrt( double x ) { return 1./sqrt(x); };
#endif

static inline double
vec3_norm(const double v[3])
{
    return sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
}

static inline void
matrix_times_vec3(double a[3], double m[3][3], const double v[3])
{
    a[0] = m[0][0]*v[0] + m[0][1]*v[1] + m[0][2]*v[2];
    a[1] = m[1][0]*v[0] + m[1][1]*v[1] + m[1][2]*v[2];
    a[2] = m[2][0]*v[0] + m[2][1]*v[1] + m[2][2]*v[2];
}

static inline void
vec3_cross_product(double axb[3], const double a[3], const double b[3])
{
    axb[0] = a[1]*b[2] - a[2]*b[1];
    axb[1] = a[2]*b[0] - a[0]*b[2];
    axb[2] = a[0]*b[1] - a[1]*b[0];
}

static inline double
vec3_dot_product(const double a[3], const double b[3])
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

static inline double
vec3_invnorm(const double v[3])
{
    return invsqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
}

static inline void
vec3_unit(double a[3])
{
    double invnorm = vec3_invnorm(a);
    a[0] *= invnorm;
    a[1] *= invnorm;
    a[2] *= invnorm;
}

#ifdef __cplusplus
}
#endif

