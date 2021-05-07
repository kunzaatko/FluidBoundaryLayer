#include "28.h"
#include "stdio.h"

void derivs_euler(float x, float z[], float dydx[]) {
  dydx[1] = z[2];
  dydx[2] = 6 * z[1] - z[2];
}

int main() {
  float *z;
  z = vector(1, 2);
  z[1] = 3;
  z[2] = 1;
  float *zout;
  zout = vector(1, 2);
  usint(z, 2, 50, 0, 1, zout, derivs_euler);
  printf("The result at x = 1 is %f", zout[2]);
}
