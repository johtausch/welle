/*
 *  kernelWelle.c 
 *  define routines related to the kernel in the wave equation
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "globals.h"
#include "kernel.h"

/*
 * Single layer potential kernel. Each routine corresponds to one term in the telescoping sum
 * representation of the convolution of the ansatz function in time.
 *   kernelSingleC0 .. kernelSingleC2   piecewise constant functions in time
 *   kernelSingleL0 .. kernelSingleL4   piecewise linear functions in time
 *
 * Parameters
 *   x, y        evaluation point
 *   nrmX, nrmY  normals at evaluation points (not referenced for the SLP).
 *   param       parameters to the kernel. See resetKernel()
 *   r           the denominator |x-y|  (return value)
 *   fcn         the numerator (return value)
 */
void kernelSingle(double *x, double *y, double *nrmX, double *nrmY, double *param, double *r, double *fcn){
  double rr[3];

  rr[0] = x[0] - y[0];
  rr[1] = x[1] - y[1];
  rr[2] = x[2] - y[2];
  double r2 = rr[0]*rr[0] + rr[1]*rr[1] + rr[2]*rr[2];
  r[0] = sqrt(r2); 
  fcn[0] = FOURPII;
} /* kernelSingle */

void kernelSingleC0(double *x, double *y, double *nrmX, double *nrmY, double *param, double *r, double *fcn){
  double rr[3];

  rr[0] = x[0] - y[0];
  rr[1] = x[1] - y[1];
  rr[2] = x[2] - y[2];
  double r2 = rr[0]*rr[0] + rr[1]*rr[1] + rr[2]*rr[2];
  r[0] = sqrt(r2); 
  fcn[0] = -FOURPII;
} /* kernelSingleC0 */

void kernelSingleC1(double *x, double *y, double *nrmX, double *nrmY, double *param, double *r, double *fcn){
  double rr[3];

  rr[0] = x[0] - y[0];
  rr[1] = x[1] - y[1];
  rr[2] = x[2] - y[2];
  double r2 = rr[0]*rr[0] + rr[1]*rr[1] + rr[2]*rr[2];
  r[0] = sqrt(r2); 
  fcn[0] = TWOPII;
} /* kernelSingleC1 */

void kernelSingleC2(double *x, double *y, double *nrmX, double *nrmY, double *param, double *r, double *fcn){
  double rr[3];

  rr[0] = x[0] - y[0];
  rr[1] = x[1] - y[1];
  rr[2] = x[2] - y[2];
  double r2 = rr[0]*rr[0] + rr[1]*rr[1] + rr[2]*rr[2];
  r[0] = sqrt(r2); 
  fcn[0] = -FOURPII;
} /* kernelSingleC2 */



void kernelSingleL0(double *x, double *y, double *nrmX, double *nrmY, double *param, double *r, double *fcn){
  double rr[3];

  rr[0] = x[0] - y[0];
  rr[1] = x[1] - y[1];
  rr[2] = x[2] - y[2];
  double r2 = rr[0]*rr[0] + rr[1]*rr[1] + rr[2]*rr[2];
  r[0] = sqrt(r2); 
  double rho = r[0]/param[0] - param[1];
  double rho2 = rho + 2;
  double p0 = -0.5*rho2*rho2;

  fcn[0] = FOURPII*p0;
} /* kernelSingleL0 */


void kernelSingleL1(double *x, double *y, double *nrmX, double *nrmY, double *param, double *r, double *fcn){
  double rr[3];

  rr[0] = x[0] - y[0];
  rr[1] = x[1] - y[1];
  rr[2] = x[2] - y[2];
  double r2 = rr[0]*rr[0] + rr[1]*rr[1] + rr[2]*rr[2];
  r[0] = sqrt(r2); 
  double rho = r[0]/param[0] - param[1];
  double rho1 = rho + 1;
  double p1 = 2*rho1*rho1;
 
  fcn[0] = FOURPII*p1;
} /* kernelSingleL1 */



/*
 * This evaluates the smooth part of the single layer potential 
 */
void kernelSingleL2(double *x, double *y, double *nrmX, double *nrmY, double *param, double *r, double *fcn){
  double rr[3];

  rr[0] = x[0] - y[0];
  rr[1] = x[1] - y[1];
  rr[2] = x[2] - y[2];
  double r2 = rr[0]*rr[0] + rr[1]*rr[1] + rr[2]*rr[2];
  r[0] = sqrt(r2); 
  double rho = r[0]/param[0] - param[1];
  double p2 = -3*rho*rho;
  
  fcn[0] = FOURPII*p2;
} /* kernelSingleL2 */



/*
 * This evaluates the smooth part of the single layer potential 
 */
void kernelSingleL3(double *x, double *y, double *nrmX, double *nrmY, double *param, double *r, double *fcn){
  double rr[3];

  rr[0] = x[0] - y[0];
  rr[1] = x[1] - y[1];
  rr[2] = x[2] - y[2];
  double r2 = rr[0]*rr[0] + rr[1]*rr[1] + rr[2]*rr[2];
  r[0] = sqrt(r2); 
  double rho = r[0]/param[0] - param[1];
  double rho1 = rho - 1;
  double p3 = 2*rho1*rho1;

  fcn[0] = FOURPII*p3;
} /* kernelSingleL3 */


/*
 * This evaluates the smooth part of the single layer potential 
 */
void kernelSingleL4(double *x, double *y, double *nrmX, double *nrmY, double *param, double *r, double *fcn){
  double rr[3];

  rr[0] = x[0] - y[0];
  rr[1] = x[1] - y[1];
  rr[2] = x[2] - y[2];
  double r2 = rr[0]*rr[0] + rr[1]*rr[1] + rr[2]*rr[2];
  r[0] = sqrt(r2); 
  double rho = r[0]/param[0] - param[1];
  double rho2 = rho - 2;
  double p4 = -0.5*rho2*rho2;

  fcn[0] = FOURPII*p4;
} /* kernelSingleL4 */



/*
 * Reset the kernel. 
 */ 
void resetKernel(kernel *krnl, int type, double h, int d){
  if ( krnl->ansatzT == 0 ) { //piecewise constant functions in time
    switch(type) {
    case -1:
      krnl->fcn = kernelSingleC0;
      break; 
    case 0:
      krnl->fcn = kernelSingleC1;
      break; 
    case 1:
      krnl->fcn = kernelSingleC2;
      break; 
    default:
      krnl->fcn = NULL;
    }
  }
  else if ( krnl->ansatzT == 1 ) { //piecewise linear functions in time
    switch(type) {
    case -2:
      krnl->fcn = kernelSingleL0;
      break; 
    case -1:
      krnl->fcn = kernelSingleL1;
      break; 
    case 0:
      krnl->fcn = kernelSingleL2;
      break; 
    case 1:
      krnl->fcn = kernelSingleL3;
      break; 
    case 2:
      krnl->fcn = kernelSingleL4;
      break; 
    default:
      krnl->fcn = NULL;
    }
  }
  else {
    krnl->fcn = NULL;
  }
  krnl->param[0] = h;
  krnl->param[1] = d;
} /* resetKernel */



/*
 * Initialize the kernel.
 * ansatzT == 0: use piecewise constant functions in time
 * ansatzT == 1: use piecewise linear functions in time
 */ 
kernel *initKernel(int ansatzT){
  kernel *krnl;

  ASSERT(ansatzT==0); //so far, only this has been implemented
  
  krnl = malloc( sizeof(kernel) );
  krnl->nParam = 2;
  krnl->param = calloc( 2, sizeof(double) );
  krnl->ansatzT = ansatzT;
  krnl->fcn = kernelSingle;
  return krnl;
} /* initKernel */
