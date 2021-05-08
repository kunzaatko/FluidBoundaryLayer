#include "blais.h"

void blaisderivs(float x, float z[], float dzdx[]) {
  // {{{
  UNUSED(x);
  dzdx[1] = z[1];
  dzdx[2] = z[2];
  dzdx[3] = -LAMBDA * (1 - pow(z[1], 2)) - z[2] * z[0];
} // }}}

void blaisB2(float zA2start, float *blaisB2out) {
  // {{{
  float *zstart;
  zstart = vector(1, 3);
  zstart[1] = ZA0;
  zstart[2] = ZA1;
  zstart[3] = zA2start;

  float *zout;
  zout = vector(1, 3);
  usint(zstart, 3, STEPS, A, B, zout, blaisderivs);

  *blaisB2out = zout[3] - 1;

  free_vector(zout, 1, 3);
} // }}}

void blaisdB2dz(float zA2start, float *fz, float *dB2dz) {
  // {{{
  dfdz(zA2start, EPS, dB2dz, fz, *blaisB2);
} // }}}

void blais(float zA2start, float acc, char *csvfilename, float zout[]) {
  // {{{
  float zA2out;
  nrrf(zA2start, acc, &zA2out, *blaisdB2dz);

  float *zAout;
  zAout = vector(1, 3);
  zAout[1] = ZA0;
  zAout[2] = ZA1;
  zAout[3] = zA2out;

  UNUSED(csvfilename);

  usint(zAout, 3, STEPS, A, B, zout, *blaisderivs);
  free_vector(zAout, 1, 3);
} // }}}
