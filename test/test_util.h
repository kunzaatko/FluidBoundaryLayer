#ifndef TEST_UTIL_H
#define TEST_UTIL_H

const float EPS = 1e-6;

void derivs_euler(float x, float z[], float dydx[]);
void dydx_poly(float x, float *fx, float *dydx);

void poly(float x, float *y);
void dydx_poly_dfdz(float x, float *y, float *dydx);

#endif // TEST_UTIL_H
