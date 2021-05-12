#include "blais.h"

void blaisnrrf(float *zA2start, float acc, float zAout[], float zBout[]) {
  // {{{
  float zA2out;  // output root of bliasB2dz (the end condition)
  nrrf(*zA2start, acc, &zA2out, *blaisdB2dz);

#ifdef DEBUG
  printf("\tblaisnrrf: Root of blais end boundary condition ZA[2] = %f\n",
         zA2out);
#endif

  zAout[1] = ZA0;
  zAout[2] = ZA1;
  zAout[3] = zA2out;

  *zA2start = zA2out;

  usint(zAout, 3, STEPS, A, B, zBout, *blaisderivs);
}  // }}}

void blaisibis(float zA2start_1, float zA2start_2, float acc, char *csvfile,
               float zAout[], float zBout[]) {
  // {{{
  float zA2out;  // output root of bliasB2dz (the end condition)
  ibis(zA2start_1, zA2start_2, acc, &zA2out, *blaisB2);

  zAout[1] = ZA0;
  zAout[2] = ZA1;
  zAout[3] = zA2out;

#ifdef DEBUG
  printf("\tblaisibis: Root of blais end boundary condition ZA[2] = %f\n",
         zA2out);
#endif

  if (strcmp(csvfile, "")) {
    usint_print(zAout, 3, STEPS, A, B, zBout, csvfile, *blaisderivs);
  } else {
    usint(zAout, 3, STEPS, A, B, zBout, *blaisderivs);
  }
}  // }}}

void blaisderivs(float x, float z[], float dzdx[]) {
  // {{{
  UNUSED(x);
  dzdx[1] = z[2];
  dzdx[2] = z[3];
  dzdx[3] = (-LAMBDA) * (1 - pow(z[2], 2)) - z[3] * z[1];
}  // }}}

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

  *blaisB2out = zout[2] - ZB1;

  free_vector(zout, 1, 3);
}  // }}}

void blaisdB2dz(float zA2start, float *fz, float *dB2dz) {
  // {{{
  dfdz(zA2start, EPS, dB2dz, fz, *blaisB2);
}  // }}}
