#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

float *vector(long nl, long nh) {/*{{{*/
  /* allocate a float vector with subscript range v[nl..nh] */
  float *v;
  v = (float *)malloc((size_t)((nh - nl + 1) * sizeof(float)));
  return v - nl;
}/*}}}*/

void free_vector(float *v, long nl, long nh)/*{{{*/
/* free a float vector allocated with vector() */
{
  free((char *)(v + nl));
}/*}}}*/

float **matrix(long nrl, long nrh, long ncl, long nch)/*{{{*/
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
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
}/*}}}*/

void free_matrix(float **m, long nrl, long nrh, long ncl, long nch)/*{{{*/
/* free a float matrix allocated by matrix() */
{
  free((char *)(m[nrl] + ncl));
  free((char *)(m + nrl));
}/*}}}*/

void rk4(float y[], float dydx[], int n, float x, float h, float yout[],/*{{{*/
         void (*derivs)(float, float[], float[]))
/* Given values for the variables y[1..n] and their derivatives dydx[1..n]
   known at x, use the fourth-order Runge-Kutta method to advance the
   solution over an interval h and return the incremented variables as
   yout[1..n], which need not be a distinct array from y. The user */
{
  int i;
  float xh, hh, h6, *dym, *dyt, *yt;
  yt = vector(1, n);
  hh = h * 0.5;
  h6 = h / 6.0;
  xh = x + hh;
  for (i = 1; i <= n; i++)
    yt[i] = y[i] + hh * dydx[i];
  (*derivs)(xh, yt, dyt);
  for (i = 1; i <= n; i++)
    yt[i] = y[i] + hh * dyt[i];
  (*derivs)(xh, yt, dym);
  for (i = 1; i <= n; i++) {
    yt[i] = y[i] + h * dym[i];
    dym[i] += dyt[i];
  }
  (*derivs)(x + h, yt, dyt);
  for (i = 1; i <= n; i++)
    yout[i] = y[i] + h6 * (dydx[i] + dyt[i] + 2.0 * dym[i]);
  free_vector(yt, 1, n);
  free_vector(dyt, 1, n);
  free_vector(dym, 1, n);
}/*}}}*/

