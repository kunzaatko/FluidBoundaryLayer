#include "../src/28.h"
#include "vendor/unity.h"
#include "test_28.h"

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
  // true value is ~ 14.828
  usint(z, 2, 500, 0, 1, zout, *derivs_euler);
  TEST_ASSERT_TRUE(zout[1] < 16);
  TEST_ASSERT_TRUE(zout[1] > 14);
  free_vector(z, 1, 2);
  free_vector(zout, 1, 2);
}
// }}}

int main(void) {
  UnityBegin("test/test_28.c");

  RUN_TEST(test_euler_usint);

  return UnityEnd();
}
