#include "blais.h"
#include "util.h"
#include <stdio.h>

int main() {
  float zA2start = 0;
  float acc = 0.1;
  float *zout;
  zout = vector(1, 3);
  blais(zA2start, acc, "", zout);
  printf("\tWith starting parameter Z2(%.0f) = %.0f and accuracy %.2f, "
         "Z(%.0f)=[%.2f,%.2f,%.2f]",
         A, zA2start, acc, B, zout[1], zout[2], zout[3]);
  free_vector(zout, 1, 3);
}
