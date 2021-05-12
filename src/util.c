#include "util.h"

#include <stdio.h>

float *vector(long nl, long nh) {
  // {{{
  float *v;
  v = (float *)malloc((size_t)((nh - nl + 1) * sizeof(float)));
  return v - nl;
}  // }}}

void free_vector(float *v, long nl, long nh) {
  // {{{
  UNUSED(nh);
  free((char *)(v + nl));
}  // }}}

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
  for (i = nrl + 1; i <= nrh; i++) m[i] = m[i - 1] + ncol;
  /* return pointer to array of pointers to rows */
  return m;
}  // }}}

void free_matrix(float **m, long nrl, long nrh, long ncl, long nch) {
  // {{{
  UNUSED(nrh);
  UNUSED(nch);
  free((char *)(m[nrl] + ncl));
  free((char *)(m + nrl));
}  // }}}

void rk4(float z[], float dzstart[], int n, float x, float h, float zout[],
         void (*derivs)(float, float[], float[])) {
  // {{{
#ifdef DEBUG_RK4
  printf("\t--- rk4 ---\n");
#endif
  float *zt;
  zt = vector(1, n);  // function value at midpoints (only storage)

  float hh = h * 0.5;  // half step
  for (int i = 1; i <= n; i++) {
    // function at midpoint by startpoint derivative (midpoint1)
    zt[i] = z[i] + hh * dzstart[i];
  }
  float *dz1;  // derivative at midpoint1
  dz1 = vector(1, 2);
  (*derivs)(x + hh, zt, dz1);

  for (int i = 1; i <= n; i++) {
    // function at midpoint by midpoint1 derivative (midpoint2)
    zt[i] = z[i] + hh * dz1[i];
  }
  float *dz2;  // derivative at midpoint2
  dz2 = vector(1, 2);
  (*derivs)(x + hh, zt, dz2);

  for (int i = 1; i <= n; i++) {
    zt[i] = z[i] + h * dz2[i];
    dz2[i] += dz1[i];
  }
  float *dzmid = dz2;  // accumulated midpoint derivative
  float *dzend = dz1;  // end point derivative
  (*derivs)(x + h, zt, dzend);

  // evaluation of `yout` by the Rundge-Kutta method
  for (int i = 1; i <= n; i++) {
    zout[i] = z[i] + h / 6 * (dzstart[i] + dzend[i] + 2.0 * dzmid[i]);
  }

#ifdef DEBUG_RK4
  printf(
      "\tRundge-Kutta 4 step from Z(%.6f) = [%.2f,%.2f,%.2f] to Z(%.6f) = "
      "[%.2f,%.2f,%.2f]\n",
      x, z[1], z[2], z[3], x + h, zout[1], zout[2], zout[3]);
#endif

  // free memory
  free_vector(zt, 1, n);
  free_vector(dzend, 1, n);
  free_vector(dzmid, 1, n);
}  // }}}

void usint(float z[], int n, int s, float a, float b, float zout[],
           void (*derivs)(float, float[], float[])) {
  // {{{
#ifdef DEBUG_USINT
  printf("\t--- usint ---\n");
#endif

  float *dzstart;
  dzstart = vector(1, n);  // derivative at the start of the interval
  float ss = (b - a) / s;  // step size

  // initializing moving variables
  float x = a;    // moving value of x
  float *zstart;  // moving value of z at the beginning of the interval
  zstart = vector(1, n);
  for (int i = 1; i <= n; i++) {
    zstart[i] = z[i];
  }
  float *zend;  // moving value of z at the end of the interval
  zend = vector(1, n);

  // evaluation from a to b
  while (x < b) {
    derivs(x, zstart, dzstart);  // derivative at the start of the interval
    rk4(zstart, dzstart, n, x, ss, zend, derivs);  // evaluating the next zend

    for (int i = 1; i <= n; i++) {
      zstart[i] = zend[i];  // next zend is the start for the next iteration
    }
    x += ss;
  }

  for (int i = 1; i <= n; i++) {
    zout[i] = zend[i];
  }

#ifdef DEBUG_USINT
  printf(
      "\tInterval evaluation from Z(%f) = [%.2f,%.2f,%.2f] to Z(%f) = "
      "[%.2f,%.2f,%.2f]\n",
      a, z[1], z[2], z[3], b, zout[1], zout[2], zout[3]);
#endif

  // free memory
  free_vector(zend, 1, n);
  free_vector(zstart, 1, n);
  free_vector(dzstart, 1, n);
}  // }}}

void usint_print(float z[], int n, int s, float a, float b, float zout[],
                 char *csvfile, void (*derivs)(float, float[], float[])) {
  // {{{
#ifdef DEBUG_USINT
  printf("\t--- usint_print ---\n");
#endif

  FILE *ftp;
  ftp = fopen(csvfile, "w+");
  fprintf(ftp, "x, z_0, z_1, z_2\n");

  float *dzstart;
  dzstart = vector(1, n);  // derivative at the start of the interval
  float ss = (b - a) / s;  // step size

  // initializing moving variables
  float x = a;    // moving value of x
  float *zstart;  // moving value of z at the beginning of the interval
  zstart = vector(1, n);
  for (int i = 1; i <= n; i++) {
    zstart[i] = z[i];
  }

  fprintf(ftp, "%f, %f, %f, %f\n", x, zstart[1], zstart[2], zstart[3]);

  float *zend;  // moving value of z at the end of the interval
  zend = vector(1, n);

  // evaluation from a to b
  while (x < b) {
    derivs(x, zstart, dzstart);  // derivative at the start of the interval
    rk4(zstart, dzstart, n, x, ss, zend, derivs);  // evaluating the next zend

    for (int i = 1; i <= n; i++) {
      zstart[i] = zend[i];  // next zend is the start for the next iteration
    }
    x += ss;

    fprintf(ftp, "%f, %f, %f, %f\n", x, zstart[1], zstart[2], zstart[3]);
  }

  for (int i = 1; i <= n; i++) {
    zout[i] = zend[i];
  }

  fprintf(ftp, "%f, %f, %f, %f\n", b, zout[1], zout[2], zout[3]);

#ifdef DEBUG_USINT
  printf(
      "\tInterval evaluation from Z(%f) = [%.2f,%.2f,%.2f] to Z(%f) = "
      "[%.2f,%.2f,%.2f]\n",
      a, z[1], z[2], z[3], b, zout[1], zout[2], zout[3]);
#endif

  // free memory
  free_vector(zend, 1, n);
  free_vector(zstart, 1, n);
  free_vector(dzstart, 1, n);
}  // }}}

void dfdz(float z, float dz, float *dfdz, float *fz,
          void (*f)(float, float *)) {
  // {{{
#ifdef DEBUG_DFDZ
  printf("\t--- dfdz ---\n");
#endif
  float fzpdz;  // f(z+dz), f(z)
  f(z + dz, &fzpdz);
  f(z, fz);
  /* TODO: Tady by bylo asi vhodnější vracet hodnotu, protože to nic nestojí.
   * Nastavení proměnné je stejně náročné jako jí vytvořit protože float žije na
   * stacku. (Copy x Clone) <08-05-21, kunzaatko> */
  *dfdz = (fzpdz - *fz) / dz;

#ifdef DEBUG_DFDZ
  printf(
      "\tdfdz: derivative dfdz(%.2f) = (%.6f - %.6f)/%.6f = %.6f/%.6f = %.2f\n",
      z, fzpdz, *fz, dz, fzpdz - *fz, dz, *dfdz);
#endif
}  // }}}

void nrrf(float zstart, float acc, float *zout,
          void (*dfdz)(float, float *, float *)) {
  // {{{
#ifdef DEBUG_NRRF
  printf("\t--- nrrf ---\n");
#endif
  float fz, ddfdz;   // function value, derivative at z
  float z = zstart;  // moving point z
  for (;;) {
    dfdz(z, &fz, &ddfdz);
#ifdef DEBUG_NRRF
    printf("\tNewton-Raphson method dfdz = %.2f: z1 = %.2f, f(z1) = %.2f",
           ddfdz, z, fz);
#endif
    if (fabsf(fz) <= acc) {
#ifdef DEBUG
      printf("\n");
#endif
      break;
    } else {
      z -= fz / ddfdz;  // Newton-Raphson step for one dimensional root finding
#ifdef DEBUG_NRRF
      printf("\tnrrf: Moving variable to x = %.2f\n", z);
#endif
    }
  }
#ifdef DEBUG
  printf("\tNewton-Raphson method returning root %f\n", z);
#endif
  *zout = z;
}  // }}}

void ibis(float zstart_1, float zstart_2, float acc, float *zout,
          void (*func)(float x, float *fout)) {
  // {{{
#ifdef DEBUG_IBIS
  printf("\t--- ibis ---\n");
  float val_at_root;
#endif

  float fout_1, fout_2, fout_new, znew;
  func(zstart_1, &fout_1);
  func(zstart_2, &fout_2);

  for (;;) {
    if (fabs(zstart_1 - zstart_2) < acc) {
      break;
    }

    znew = (zstart_1 + zstart_2) / 2;
    func(znew, &fout_new);

    if (fout_new * fout_1 > 0) {
      zstart_1 = znew;
      continue;
    } else if (fout_new == 0)  // we should check for the neiche but possible
      break;  // outcome that the new point landed exactly at the root
    else {
      zstart_2 = znew;
    }

#ifdef DEBUG_IBIS
    printf("\tBisection at interval (%.2f,%.2f)\n", zstart_1, zstart_2);
#endif
  }

  *zout = (zstart_1 + zstart_2) / 2;

#ifdef DEBUG_IBIS
  func(*zout, &val_at_root);
  printf(
      "\tBisection root finding exiting with interval (%.2f,%.2f) returning "
      "value %.3f where f(%.3f) = %f\n",
      zstart_1, zstart_2, *zout, *zout, val_at_root);
#endif
}  // }}}
