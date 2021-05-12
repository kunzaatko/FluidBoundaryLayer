#ifndef UTIL_H
#define UTIL_H

#include <math.h>
#include <stddef.h>
#include <stdlib.h>

// #define DEBUG = 1
// #define DEBUG_RK4 = 1
// #define DEBUG_USINT = 1
// #define DEBUG_DFDZ = 1
// #define DEBUG_NRRF = 1
// #define DEBUG_IBIS = 1

#define UNUSED(x) (void)(x)

/* allocates a float vector with subscript range v[nl..nh] */
float *vector(long nl, long nh);

/* free a float vector allocated with vector() */
void free_vector(float *v, long nl, long nh);

/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
float **matrix(long nrl, long nrh, long ncl, long nch);

/* free a float matrix allocated by matrix() */
void free_matrix(float **m, long nrl, long nrh, long ncl, long nch);

/* Evaluation of the Z(x + h) using the fourth-order Rundge-Kutta method.
 *  z[] - value of Z(x)
 *  dzstart[] - derivative of z at point x (it is inserted as an argument
 *              because it is used in determining the step size so as not
 *              to allocate more times than needed)
 *  n - size of vector z (number of first-order equations that are calculated)
 *  x - point we are advancing from h - step size
 *  derivs - function for derivative of Z derivs(x,z,dzdx)
 *  zout[] - [output] value of Z(x + h)
 */
void rk4(float z[], float dzstart[], int n, float x, float h, float zout[],
         void (*derivs)(float, float[], float[]));

/* Evaluation of the function from a to b using uniform steps and the fourth
 * order Rundge-Kutta method.
 *  z[] - value of Z(a)
 *  n - size of vector z (number of first-order equations that are calculated)
 *  s - number of uniform steps to take between a and b
 *  a - initial point
 *  b - end point
 *  derivs - function for derivative of Z
 *  zout[] - [output] value of Z(b)
 */
void usint(float z[], int n, int s, float a, float b, float zout[],
           void (*derivs)(float, float[], float[]));

/* Same as usint but prints the output to a csvfile
 * z[] - value of Z(a)
 *  n - size of vector z (number of first-order equations that are calculated)
 *  s - number of uniform steps to take between a and b
 *  a - initial point
 *  b - end point
 *  derivs - function for derivative of Z
 *  zout[] - [output] value of Z(b)
 *  char* csvfile - file name for the output file
 */
void usint_print(float z[], int n, int s, float a, float b, float zout[],
                 char *csvfile, void (*derivs)(float, float[], float[]));

/* Approximate the derivative of function f by z using delta dz.
 *  z - value where the derivative should be evaluated
 *  dz - delta to use for the evaluation
 *  *dfdz - [output] derivative at z
 *  *fz - value of f at z (for use in the Newton Raphson)
 *  f - function of z with arguments (z, *fout)
 */
void dfdz(float z, float dz, float *dfdz, float *fz,
          void (*f)(float x, float *fout));

/* Use the Newton-Rhapson method for finding the root
 *  zstart - initial value to start algorithm from
 *  acc - accuracy to be accieved
 *  zout - [output] root
 *  dfdz - function of the derivative z with arguments (z, *fz, *dfdz)
 */
void nrrf(float zstart, float acc, float *zout,
          void (*dfdz)(float z, float *fz, float *dfdz));

/* Use the interval bisection method for finding the root
 *  zstart_1,zstart_2 - initial values to start algorithm (it is necessary that
 * sign(f(zstart_1)) = -sign(f(zstart_2)), that means the root must lie in the
 * interval)
 *  acc - accuracy to be accieved zout - [output] root func - function
 */
void ibis(float zstart_1, float zstart_2, float acc, float *zout,
          void (*func)(float x, float *fout));

#endif  // UTIL_H
