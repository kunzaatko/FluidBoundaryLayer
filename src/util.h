#ifndef UTIL_H
#define UTIL_H
#include <math.h>
#include <stddef.h>
#include <stdlib.h>

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
 *  derivs - function for derivative of Z [derivs(x,y,dydz)]
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

/* Evaluate the derivative of function f by z using delta dz.
 *  z - value where the derivative should be evaluated
 *  dz - delta to use for the evaluation
 *  dfdz - [output] derivative at z
 *  fz - value of f at z (for use in the Newton Raphson)
 *  f - function of z with arguments (z, *fout)
 */
void dfdz(float z, float dz, float *dfdz, float *fz, void (*f)(float, float *));

/* Use the Newton-Rhapson method for finding the root
 *  zstart - initial value to start algorithm from
 *  acc - accuracy to be accieved
 *  zout - [output] root
 *  dfdz - function of the derivative z with arguments (z, *fz, *dfdz)
 */
void nrrf(float zstart, float acc, float *zout,
          void (*dfdz)(float z, float *, float *));

#endif // UTIL_H
