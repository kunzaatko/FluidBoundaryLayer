#include "test_util.h"
#include "../src/util.h"
#include "vendor/unity.h"
#include <math.h>

void setUp(void) {}

void tearDown(void) {}

// {{{ TEST_EULER_USINT
// solution of diffeq analytically is y(x) = e^{-3x} + 2e^{2x}
void derivs_euler(float x, float z[], float dydx[]) {
  UNUSED(x);
  dydx[1] = z[2];
  dydx[2] = 6 * z[1] - z[2];
}
static void test_euler_usint(void) {
  float *z;
  z = vector(1, 2);
  z[1] = 3;
  z[2] = 1;
  float *zout;
  zout = vector(1, 2);
  int nsteps = 10;
  // true value is ~ 14.828
  usint(z, 2, nsteps, 0, 1, zout, *derivs_euler);
  printf("\tValue of y(1) = e^{-3} + 2e^{2} estimated from %d steps is %.3f "
         "(true value 14.828)\n",
         nsteps, zout[1]);
  TEST_ASSERT_TRUE(zout[1] < 16);
  TEST_ASSERT_TRUE(zout[1] > 14);
  free_vector(z, 1, 2);
  free_vector(zout, 1, 2);
}
// }}}

// {{{ TEST_POLY_NEWTON_RAPHSON
void dydx_poly(float x, float *fx, float *dydx) {
  *fx = pow(x, 4) - pow(x, 3) + 5 * pow(x, 2) - 6;
  *dydx = 4 * pow(x, 3) - 3 * pow(x, 2) + 10 * pow(x, 1);
}
static void test_poly_nrrf(void) {
  float x1 = 0.1;
  float root_right;
  nrrf(x1, EPS, &root_right, (*dydx_poly));
  printf("\tRight root of polynom x^{4} - x^{3} + 5x^{2} - 6 is found as %.5f"
         "(true value 1.0854)\n",
         root_right);
  TEST_ASSERT_TRUE(root_right < 1.086);
  TEST_ASSERT_TRUE(root_right > 1.085);

  float x2 = -90;
  float root_left;
  nrrf(x2, EPS, &root_left, (*dydx_poly));
  printf("\tLeft root of polynom x^{4} - x^{3} + 5x^{2} - 6 is found as %.5f"
         "(true value -0.93809)\n",
         root_left);
  TEST_ASSERT_TRUE(root_left < -0.9380);
  TEST_ASSERT_TRUE(root_left > -0.9381);
}
// }}}

// {{{ TEST_POLY_DFDZ_NRRF
void poly(float x, float *y) { *y = pow(x, 4) - pow(x, 3) + 5 * pow(x, 2) - 6; }
void dydx_poly_dfdz(float x, float *y, float *dydx) {
  dfdz(x, EPS, dydx, y, *poly);
}
static void test_poly_dfdz_nrrf(void) {
  float x1 = 0.1;
  float root_right;
  nrrf(x1, EPS, &root_right, (*dydx_poly_dfdz));
  printf("\tRight root of polynom x^{4} - x^{3} + 5x^{2} - 6 is found as %.5f"
         "(true value 1.0854)\n",
         root_right);
  TEST_ASSERT_TRUE(root_right < 1.086);
  TEST_ASSERT_TRUE(root_right > 1.085);

  float x2 = -5;
  float root_left;
  nrrf(x2, EPS, &root_left, (*dydx_poly_dfdz));
  printf("\tLeft root of polynom x^{4} - x^{3} + 5x^{2} - 6 is found as %.5f"
         "(true value -0.93809)\n",
         root_left);
  TEST_ASSERT_TRUE(root_left < -0.9380);
  TEST_ASSERT_TRUE(root_left > -0.9381);
}
// }}}

int main(void) {
  UnityBegin("test/test_util.c");

  RUN_TEST(test_euler_usint);
  printf("--------------\n");
  RUN_TEST(test_poly_nrrf);
  printf("--------------\n");
  RUN_TEST(test_poly_dfdz_nrrf);

  return UnityEnd();
}
