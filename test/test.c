#include "test.h"

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
#ifdef DEBUG
  printf("\t--- test_euler_usint ---\n");
#endif
  float *z;
  z = vector(1, 2);
  z[1] = 3;
  z[2] = 1;
  float *zout;
  zout = vector(1, 2);
  int nsteps = 10;
  // true value is ~ 14.828
  usint(z, 2, nsteps, 0, 1, zout, *derivs_euler);
#ifdef DEBUG
  printf(
      "\tValue of y(1) = e^{-3} + 2e^{2} estimated from %d steps is %.3f "
      "(true value 14.828)\n",
      nsteps, zout[1]);
#endif
  TEST_ASSERT_TRUE(zout[1] < 16);
  TEST_ASSERT_TRUE(zout[1] > 14);

  free_vector(z, 1, 2);
  free_vector(zout, 1, 2);
}
// }}}

// {{{ TEST_POLY_DFDZ_NRRF
void poly(float x, float *y) { *y = pow(x, 4) - pow(x, 3) + 5 * pow(x, 2) - 6; }
void dydx_poly_dfdz(float x, float *y, float *dydx) {
  float acc = 1e-6;
  dfdz(x, acc, dydx, y, *poly);
}
static void test_poly_dfdz_nrrf(void) {
#ifdef DEBUG
  printf("\t--- test_poly_dfdz_nrrf ---\n");
#endif
  float acc = 1e-6;

  float x1 = 0.1;
  float root_right;
  nrrf(x1, acc, &root_right, (*dydx_poly_dfdz));
#ifdef DEBUG
  printf(
      "\tRight root of polynom x^{4} - x^{3} + 5x^{2} - 6 is found as %.5f"
      "(true value 1.0854)\n",
      root_right);
#endif
  TEST_ASSERT_TRUE(root_right < 1.086);
  TEST_ASSERT_TRUE(root_right > 1.085);

  float x2 = -5;
  float root_left;
  nrrf(x2, acc, &root_left, (*dydx_poly_dfdz));
#ifdef DEBUG
  printf(
      "\tLeft root of polynom x^{4} - x^{3} + 5x^{2} - 6 is found as %.5f"
      "(true value -0.93809)\n",
      root_left);
#endif
  TEST_ASSERT_TRUE(root_left > -0.9381);
}
// }}}

// {{{ TEST_POLY_NEWTON_RAPHSON
void dydx_poly(float x, float *fx, float *dydx) {
  *fx = pow(x, 4) - pow(x, 3) + 5 * pow(x, 2) - 6;
  *dydx = 4 * pow(x, 3) - 3 * pow(x, 2) + 10 * pow(x, 1);
}
static void test_poly_nrrf(void) {
#ifdef DEBUG
  printf("\t--- test_poly_nrrf ---\n");
#endif
  float acc = 1e-6;

  float x1 = 0.1;
  float root_right;
  nrrf(x1, acc, &root_right, (*dydx_poly));
#ifdef DEBUG
  printf(
      "\tRight root of polynom x^{4} - x^{3} + 5x^{2} - 6 is found as %.5f"
      "(true value 1.0854)\n",
      root_right);
#endif
  TEST_ASSERT_TRUE(root_right < 1.086);
  TEST_ASSERT_TRUE(root_right > 1.085);

  float x2 = -90;
  float root_left;
  nrrf(x2, acc, &root_left, (*dydx_poly));
#ifdef DEBUG
  printf(
      "\tLeft root of polynom x^{4} - x^{3} + 5x^{2} - 6 is found as %.5f"
      "(true value -0.93809)\n",
      root_left);
#endif
  TEST_ASSERT_TRUE(root_left < -0.9380);
  TEST_ASSERT_TRUE(root_left > -0.9381);
}
// }}}

// {{{ TEST_POLY_IBIS
static void test_poly_ibis(void) {
#ifdef DEBUG
  printf("\t--- test_poly_ibis ---\n");
#endif
  float acc = 1e-6;

  float xr1 = 0.1;
  float xr2 = 2;
  float root_right;
  ibis(xr1, xr2, acc, &root_right, *poly);
#ifdef DEBUG
  printf(
      "\tRight root of polynom x^{4} - x^{3} + 5x^{2} - 6 is found as %.5f"
      "(true value 1.0854)\n",
      root_right);
#endif
  TEST_ASSERT(root_right < 1.086);
  TEST_ASSERT(root_right > 1.085);

  float xl1 = -34;
  float xl2 = -0.1;
  float root_left;
  ibis(xl1, xl2, acc, &root_left, *poly);
#ifdef DEBUG
  printf(
      "\tLeft root of polynom x^{4} - x^{3} + 5x^{2} - 6 is found as %.5f"
      "(true value -0.93809)\n",
      root_left);
#endif
  TEST_ASSERT(root_left < -0.9380);
  TEST_ASSERT(root_left > -0.9381);
}
// }}}

// {{{ TEST_MONOM_3ORD_IBIS
void monom_3ord(float x, float *fout) { *fout = powf(x - 2, 3); }
static void test_monom_3ord_ibis(void) {
#ifdef DEBUG
  printf("\t--- test_monom_3ord_ibis ---\n");
#endif
  float acc = 1e-6;

  float x1 = 5;
  float x2 = -5;
  float root;
  ibis(x1, x2, acc, &root, *monom_3ord);
#ifdef DEBUG
  printf(
      "\tRoot of polynom (x-2)^{3} is found as %.5f"
      "(true value 2)\n",
      root);
#endif
  TEST_ASSERT(root < 2.01);
  TEST_ASSERT(root > 1.99);
}
// }}}

// {{{ TEST_BLAIS_NRRF
static void test_blais_nrrf(void) {
#ifdef DEBUG
  printf("\t--- test_blais_nrrf ---\n");
#endif
  float zA2start = 1.6;
  float acc = 0.1;
  float *zAout;
  zAout = vector(1, 3);
  float *zBout;
  zBout = vector(1, 3);
  blaisnrrf(&zA2start, acc, zAout, zBout);
#ifdef DEBUG
  printf(
      "\tWith starting parameter Z2(%.1f) = %.5f and accuracy %.2f, "
      "Z(%.1f)=[%.2f,%.2f,%.2f]\n",
      A, zA2start, acc, B, zBout[1], zBout[2], zBout[3]);
#endif
  TEST_ASSERT(zA2start > 0.4)
  TEST_ASSERT(zA2start < 0.5)
  free_vector(zAout, 1, 3);
  free_vector(zBout, 1, 3);
}
// }}}

// {{{ TEST_BLAIS_IBIS
static void test_blais_ibis(void) {
#ifdef DEBUG
  printf("\t--- test_blais_ibis ---\n");
#endif
  float zA2start_1 = 0;
  float zA2start_2 = 3;
  float acc = 0.0001;
  float *zAout;
  zAout = vector(1, 3);
  float *zBout;
  zBout = vector(1, 3);
  blaisibis(zA2start_1, zA2start_2, acc, "", zAout, zBout);
#ifdef DEBUG
  printf(
      "\tWith starting interval at x = %.1f (%.2f,%.2f) and accuracy %.2f, "
      "Z(%.1f) = [%.2f,%.2f,%.2f], Z(%.1f)=[%.2f,%.2f,%.2f]\n",
      A, zA2start_1, zA2start_2, acc, A, zAout[1], zAout[2], zAout[3], B,
      zBout[1], zBout[2], zBout[3]);
#endif
  TEST_ASSERT(zAout[3] > 0.4)
  TEST_ASSERT(zAout[3] < 0.5)
  free_vector(zAout, 1, 3);
  free_vector(zBout, 1, 3);
}
// }}}

int main(void) {
  UnityBegin("test/test.c");

  printf("\n");
  printf("   UTIL TEST\n");
  printf("==============\n\n");
  RUN_TEST(test_euler_usint);
  printf("--------------\n");
  RUN_TEST(test_poly_nrrf);
  printf("--------------\n");
  RUN_TEST(test_poly_dfdz_nrrf);
  printf("--------------\n");
  RUN_TEST(test_poly_ibis);
  printf("--------------\n");
  RUN_TEST(test_monom_3ord_ibis);
  printf("==============\n\n");

  printf("   BLAIS TEST\n");
  printf("==============\n");
  RUN_TEST(test_blais_nrrf);
  printf("==============\n");
  RUN_TEST(test_blais_ibis);
  printf("==============\n");

  return UnityEnd();
}
