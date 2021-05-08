#ifndef BLAIS_H
#define BLAIS_H

#include "util.h"
#include <math.h>
#include <stdio.h>

// initial conditions
static const float ZA0 = 0;
static const float ZA1 = 0;

// final condition
static const float ZB1 = 1;

// interval
static const float A = 0;
static const float B = 10;

static const float LAMBDA = 0.5;

static const int STEPS = 1e4;  // number of steps to make in shot
static const float EPS = 1e-4; // delta to evalutate the derivative

// file to write the output to
// static FILE *ftp;

/* Evaluate the blais function from A to B with conditions ZA0 and ZA1 at A and
 * condition ZB1 at B
 *  zA2start - initial condition to start with
 *  acc - accuracy on the final condition
 *  csv - the file name we want to log to
 *  zout[] - [output] the final z
 */
void blais(float zA2start, float acc, char *csvfilename, float zout[]);

/* Derivative of blais (Z) at x and z
 *  x - x for evaluation
 *  z[] - z for evaluation
 *  dzdx[] - [output] derivatives dzdx at z and x
 */
void blaisderivs(float x, float z[], float dzdx[]);

/* Cost function for blais at end point
 *  zA2start - starting parameter for the shot
 *  *blaisB2out - [output] cost function value for zAstart
 */
void blaisB2(float zA2start, float *blaisB2out);

/* Derivative of cost function for blais at end point
 *  zA2start - starting parameter for the shot
 *  *fz - [output] function value at zA2start
 *  *dB2dz - [output] dervative value at zA2start
 */
void blaisdB2dz(float zA2start, float *fz, float *dB2dz);

#endif // BLAIS_H
