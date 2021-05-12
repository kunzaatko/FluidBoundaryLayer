#include "blais.h"
#include "util.h"
#include <stdio.h>

int main() {
  float zA2start_1 = 0.1;
  float zA2start_2 = 3;
  float acc = 0.001;
  float *zAout;
  zAout = vector(1, 3);
  float *zBout;
  zBout = vector(1, 3);
  blaisibis(zA2start_1, zA2start_2, acc, "./lam_0.1.csv", zAout, zBout);
  free_vector(zAout, 1, 3);
  free_vector(zBout, 1, 3);
}
