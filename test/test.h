#ifndef TEST_H
#define TEST_H

#include "../src/blais.h"
#include "../src/util.h"
#include "vendor/unity.h"
#include <math.h>

void derivs_euler(float x, float z[], float dydx[]);
void dydx_poly(float x, float *fx, float *dydx);

void poly(float x, float *y);
void dydx_poly_dfdz(float x, float *y, float *dydx);

#endif // TEST_H
