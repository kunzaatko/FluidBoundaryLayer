#ifndef BLAIS_H
#define BLAIS_H

#include <math.h>
#include <stdio.h>
#include <string.h>

#include "util.h"

// initial conditions
static const float ZA0 = 0; // y(A)
static const float ZA1 = 0; // y'(A)

// final condition
static const float ZB1 = 1; // y'(B)

// interval
static const float A = 0;
static const float B = 10;

static const float LAMBDA = 0;

static const int STEPS = 1e3;   // number of steps to make in shot
static const float EPS = 1e-4;  // delta to evaluate the derivative

/* Evaluate the blais function from A to B with conditions ZA0 and ZA1 at A and
 * condition ZB1 at B using Newton-Raphson root finding algorithm
 *  zA2start - initial condition to start with
 *  acc - accuracy on the final condition
 *  zAout[] - [output] z at A
 *  zBout[] - [output] z at B
 */
void blaisnrrf(float *zA2start, float acc, float zAout[], float zBout[]);

/* Evaluate the blais function from A to B with conditions ZA0 and ZA1 at A and
 * condition ZB1 at B using interval bisection root finding algorithm
 *  zA2start_1 - initial condition interval boundary 1
 *  zA2start_2 - initial condition interval boundary 2
 *  acc - accuracy on the final condition
 *  zAout[] - [output] z at A
 *  zBout[] - [output] z at B
 */
void blaisibis(float zA2start_1, float zA2start_2, float acc, char *csvfile,
               float zAout[], float zBout[]);

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

#endif  // BLAIS_H
