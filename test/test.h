#ifndef TEST_H
#define TEST_H

#include <math.h>

#include "../src/blais.h"
#include "../src/util.h"
#include "vendor/unity.h"

void derivs_euler(float x, float z[], float dydx[]);
void dydx_poly(float x, float *fx, float *dydx);

void poly(float x, float *y);
void dydx_poly_dfdz(float x, float *y, float *dydx);

void monom_3ord(float x, float *fout);

#endif  // TEST_H
