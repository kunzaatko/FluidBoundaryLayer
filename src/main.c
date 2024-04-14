#include "blais.h"
#include "util.h"
#include <stdio.h>

int main() {
  float zA2start_1 = 0.1; // lower bound to y''(A)
  float zA2start_2 = 3; // upper bound to y''(A)
  float acc = 0.001; // accuracy of the end condition `ZB1`
  float *zAout;
  zAout = vector(1, 3);
  float *zBout;
  zBout = vector(1, 3);
  char* filename = "./lam_0_0.csv"; // output file of the data
  blaisibis(zA2start_1, zA2start_2, acc, filename, zAout, zBout);
  free_vector(zAout, 1, 3);
  free_vector(zBout, 1, 3);
}
