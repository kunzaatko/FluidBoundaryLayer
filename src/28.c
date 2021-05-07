#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

/* allocates a float vector with subscript range v[nl..nh] */
float *vector(long nl, long nh) {
  // {{{
  float *v;
  v = (float *)malloc((size_t)((nh - nl + 1) * sizeof(float)));
  return v - nl;
} // }}}

/* free a float vector allocated with vector() */
void free_vector(float *v, long nl, long nh) {
  // {{{
  free((char *)(v + nl));
} // }}}

/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
float **matrix(long nrl, long nrh, long ncl, long nch) {
  // {{{
  long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
  float **m;
  /* allocate pointers to rows */
  m = (float **)malloc((size_t)((nrow) * sizeof(float *)));
  m -= nrl;
  /* allocate rows and set pointers to them */
  m[nrl] = (float *)malloc((size_t)((nrow * ncol) * sizeof(float)));
  m[nrl] -= ncl;
  for (i = nrl + 1; i <= nrh; i++)
    m[i] = m[i - 1] + ncol;
  /* return pointer to array of pointers to rows */
  return m;
} // }}}

/* free a float matrix allocated by matrix() */
void free_matrix(float **m, long nrl, long nrh, long ncl, long nch) {
  // {{{
  free((char *)(m[nrl] + ncl));
  free((char *)(m + nrl));
} // }}}

// clang-format off
/* Evaluation of the Z(x + h) using the fourth-order Rundge-Kutta method.
 *  z[] - value of Z(x)
 *  dzstart[] - derivative of z at point x (it is inserted as an argument
 *              because it is used in determining the step size so as not
 *              to allocate more times than needed) 
 *  n - size of vector z (number of first-order equations that are calculated) 
 *  x - point we are advancing from h - step size 
 *  derivs - function for derivative of Z
 *  zout[] - [output] value of Z(x + h) 
 */
// clang-format on
void rk4(float z[], float dzstart[], int n, float x, float h, float zout[],
         void (*derivs)(float, float[], float[])) {
  // {{{
  /* TODO: join these two lines?! <07-05-21, kunzaatko> */
  float *zt;
  zt = vector(1, n); // function value at midpoints (only storage)

  float hh = h * 0.5; // half step
  for (int i = 1; i <= n; i++) {
    // function at midpoint by startpoint derivative (midpoint1)
    zt[i] = z[i] + hh * dzstart[i];
  }
  float *dz1; // derivative at midpoint1
  (*derivs)(x + hh, zt, dz1);

  for (int i = 1; i <= n; i++) {
    // function at midpoint by midpoint1 derivative (midpoint2)
    zt[i] = z[i] + hh * dz1[i];
  }
  float *dz2; // derivative at midpoint2
  (*derivs)(x + hh, zt, dz2);

  for (int i = 1; i <= n; i++) {
    zt[i] = z[i] + h * dz2[i];
    dz2[i] += dz1[i];
  }
  float *dzmid = dz2; // accumulated midpoint derivative
  float *dzend = dz1; // end point derivative
  (*derivs)(x + h, zt, dzend);

  // evaluation of `yout` by the Rundge-Kutta method
  for (int i = 1; i <= n; i++) {
    zout[i] = z[i] + h / 6 * (dzstart[i] + dzend[i] + 2.0 * dzmid[i]);
  }

  // free memory
  free_vector(zt, 1, n);
  free_vector(dzend, 1, n);
  free_vector(dzmid, 1, n);
} // }}}

// clang-format off
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
// clang-format on
void usint(float z[], int n, int s, float a, float b, float zout[],
           void (*derivs)(float, float[], float[])) {
  // {{{
  float *dzstart;
  dzstart = vector(1, n); // derivative at the start of the interval
  float ss = (b - a) / s; // step size

  // initializing moving variables
  float x = a;   // moving value of x
  float *zstart; // moving value of z at the beginning of the interval
  zstart = vector(1, n);
  for (int i = 1; i <= n; i++) {
    zstart[i] = z[i];
  }
  float *zend; // moving value of z at the end of the interval
  zend = vector(1, n);

  // evaluation from a to b
  while (x < b) {
    derivs(x, zstart, dzstart); // derivative at the start of the interval
    rk4(zstart, dzstart, n, x, ss, zend, derivs); // evaluating the next zend
    for (int i = 1; i <= n; i++) {
      zstart[i] = zend[i]; // next zend is the start of the next iteration
    }
    x += ss;
  }
  for (int i = 1; i <= n; i++) {
    zout[i] = zend[i];
  }

  // free memory
  free_vector(zend, 1, n);
  free_vector(zstart, 1, n);
  free_vector(dzstart, 1, n);
} // }}}
