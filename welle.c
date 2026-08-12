/*
 * welle:  main driver 
 *
 *  usage:
 *    welle  [options] panelfile [options]
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include "globals.h"
#include "kernel.h"
#include "mesh.h"
#include "quadRule.h"
#include "polytope.h"
#include "solver.h"

int verbose = 0;
bool doVtxApx;
kernel *krnl;

quadRule **qrJacobi, **qrJ1, **qrJ2, **qrSimplex;
poly *fX, *fY;
double **rMin;
double **rMinV;
double **rMax;
double ***apexX;
double ***apexY;
double ***apexXV;
double ***apexYV;
double **apices;
double *mBarycentX, *mBarycentY;
int nIntegrals = 0;
uint64_t *nIntegralsD;


void initQuadJacobiOLD(int order);
void initQuadSimplexOLD(int order);
void initQuadJacobi(int order, char *path);
void initQuadSimplex(int order, char *path);
void intFceS0(void (*fcn)(double *x,double t, double *f), face *fX, double t, double *intgr);
void fcefceIA(poly *fX, poly *fY, double R2, double *intgr);


double fcn( double t) {
  return 0.5*pow(t,4)*exp(-2*t);
} /* fcn */

double fcn1( double t) {
  double t3 = pow(t,3);
  double t4 = t*t3;
  return (4*t3 - 2*t4)*exp(-2*t);
} /* fcn1 */



/*
 * right hand side in the integral equation
 */
void fcnRHS(double *x, double t, double *f){
  f[0] = fcn(t);
} /* fcnRHS */

/*
 * right hand side in the integral equation
 */
void fcnSol(double *x, double t, double *f){
  double sum =  fcn1(t);

  for ( double t2 = t-2; t2>0; t2 -= 2 ) {
    sum += fcn1(t2);
  }
  f[0] = sum;
} /* fcnSol */





/*
 * prints out all computed integrals, to setup the first column of all A[d]'s.
 * used to check their convergence
 */
void testIntegrals(int nVtxs, int nFces, face *fces, int Nt, double tMax) {
  verbose = 2;
  double h = tMax/Nt; 
  face *fceX = &fces[0];
  setVerticesTriangle(fX, fceX->vtx[0]->x, fceX->vtx[1]->x, fceX->vtx[2]->x);
  mBarycentX = fces[0].mBarycent;
    
  for ( int jF=0; jF<nFces; jF++ ) {
    face *fceY = &fces[jF];
    setVerticesTriangle(fY, fceY->vtx[0]->x, fceY->vtx[1]->x, fceY->vtx[2]->x);
    mBarycentY = fces[jF].mBarycent;
      
    setDistancesCartesianProduct(7, fX, fY);

    int dMin = (int) (sqrt(rMin[0][0])/h);
    int dMax = (int) (sqrt(rMax[0][0])/h) + 1;
    if ( dMax >= Nt ) dMax = Nt - 1;
    int dMin1 = MAX(dMin-1,1);
      
    for ( int d=dMin1; d<=dMax+1; d++ ) {
      double tmp;
      double R = d*h;
      fcefceIA(fX, fY, R*R, &tmp);
    }      
  }
  
} /* testIntegrals */



/*
 * setup and return all matrices A[d]
 */
double **setupMatrix(int nVtxs, int nFces, face *fces, int Nt, double tMax) {
  double h = tMax/Nt;
  struct timespec time1;
  timespec_get(&time1, TIME_UTC);
  
  nIntegralsD = calloc(100,sizeof(uint64_t));

  double **A = calloc( Nt+2, sizeof(double*) );
  int nA = nFces*nFces;

  //determine the largest d to allocate the matrix 
  int dMaxMax = 0;
  for ( int iF=0; iF<nFces; iF++ ) { 
    face *fceX = &fces[iF];
    for ( int jF=0; jF<nFces; jF++ ) {
      face *fceY = &fces[jF];

      double rMax2 = 0;
      for ( int i=0; i<3; i++ ) {
        for ( int j=0; j<3; j++ ) {
          double r2 = distPointPoint(fceX->vtx[i]->x, fceY->vtx[j]->x);
          if ( r2 > rMax2 ) rMax2 = r2;
        }
      }

      int dMax = (int) (sqrt(rMax2)/h) + 1;
      if ( dMax >= Nt ) dMax = Nt - 1;
      if ( dMax>dMaxMax ) dMaxMax = dMax; 
    }
  }
  for (int k=0; k<=dMaxMax; k++ ) {
    A[k] = calloc( nA, sizeof(double) );
  }
  double *B = calloc(dMaxMax+2,sizeof(double*));
 
      
  //actual quadrature
  uint64_t nNonZeroA = 0;
  uint64_t nFcefceIA = 0;
  for ( int iF=0; iF<nFces; iF++ ) { 
    face *fceX = &fces[iF];
    setVerticesTriangle(fX, fceX->vtx[0]->x, fceX->vtx[1]->x, fceX->vtx[2]->x);
    mBarycentX = fces[iF].mBarycent;
    for ( int jF=0; jF<nFces; jF++ ) {
      face *fceY = &fces[jF];
      double areaFac = fceX->area2*fceY->area2;
      setVerticesTriangle(fY, fceY->vtx[0]->x, fceY->vtx[1]->x, fceY->vtx[2]->x);
      mBarycentY = fces[jF].mBarycent;
      
      setDistancesCartesianProduct(7, fX, fY);

      int dMin = (int) (sqrt(rMin[0][0])/h);
      int dMax = (int) (sqrt(rMax[0][0])/h) + 1;
      if ( dMax >= Nt ) dMax = Nt - 1;
      int dMin1 = MAX(dMin-1,1);
      
      for ( int d=dMin1; d<=dMax+1; d++ ) {
        double R = d*h;
        fcefceIA(fX, fY, R*R, &B[d]);
        B[d] *= areaFac;
        nIntegralsD[nIntegrals]++;
        nFcefceIA++;
      }
      
      for ( int d=dMin; d<=dMax; d++ ) {
        if ( d>=1 ) {
          double aij = -B[d-1] + 2*B[d] - B[d+1];
          A[d][iF+nFces*jF] += aij;
        }
        else { 
          double aij = -B[1];
          A[0][iF+nFces*jF] += aij;
        }
        nNonZeroA++;
      }   
    }
  }

  free(B);
  struct timespec time2;
  timespec_get(&time2, TIME_UTC);
  double setupTime = time2.tv_sec - time1.tv_sec + 1e-9*(time2.tv_nsec - time1.tv_nsec);

  if ( verbose > 0 ) {
    uint64_t nIntegralsTotal = 0;
    for (int i=0; i<50; i++) {
      //printf("nIntegrals(%d)=%le;\n", i, (double)nIntegralsD[i]);
      nIntegralsTotal += i*nIntegralsD[i];
    }
    printf("nNonZeroA=%le  nFcefceIA=%le nIntegralsTotal=%le setupTime=%le\n",
           (double)nNonZeroA, (double)nFcefceIA,(double)nIntegralsTotal, setupTime );
  }
  return A;
} /* setupMatrix */





/*
 * setup the right hand side
 */
double **setupRHS(void (*fcn)(double *x,double t, double *f), 
                  int nVtxs, int nFces, face *fces, int Nt, double h) {
  double intgr;

  double **b = calloc( Nt, sizeof(double*) );
  for (int k=0; k<Nt; k++ ) {
    b[k] = calloc( nFces, sizeof(double) );
  }

  //k=0
  for ( int iF=0; iF<nFces; iF++ ) {
    face *fX = &fces[iF];
    intFceS0(fcn, fX, 0, &intgr);
    b[0][iF] += intgr;
  } 
  
  for (int k=1; k<Nt; k++ ) {
    for ( int iF=0; iF<nFces; iF++ ) {
      face *fX = &fces[iF];
      intFceS0(fcn, fX, k*h, &intgr);
      b[k][iF] += intgr;
      b[k-1][iF] -= intgr;
    } 
  }

  //k=Nt
  for ( int iF=0; iF<nFces; iF++ ) {
    face *fX = &fces[iF];
    intFceS0(fcn, fX, Nt*h, &intgr);
    b[Nt-1][iF] -= intgr;
  }
  
  return b;
} /* setupRHS */



/*
 * calculate the average value of a piecewise constant function on a triangulation
 */
double getAverage(double *fcn, int nFces, face *fces) {
  double avg = 0;
  double area = 0.0;
  
  for ( int iF=0; iF<nFces; iF++ ) {
    double thisArea = 0.5*fces[iF].area2;
    area += thisArea;
    avg += fcn[iF]*thisArea;
  }
  
  avg /= area;
  return avg;
} /* getAverage */



int main(int nargs, char *argv[]) {
  int job = 1;
  int order = 5;
  int nTimeSteps = 100;
  double tMax = 1.0;
  doVtxApx = true;

  char meshFile[1024] = {0};
  char quadPath[1024] = "./QuadRules";

  //parse the command line 
  for (int i = 1; i < nargs; i++)
    if (argv[i][0] == '-')
      switch (argv[i][1]) {
      case 'j':
        job = atoi(argv[i] + 3);
        break;
      case 'q': order = atoi( argv[i]+3 );
        break;
      case 'Q': strcpy( quadPath, argv[i]+3);
        break;
      case 'v': verbose = atoi(argv[i] + 3);
        break;
      case 'V': doVtxApx = false;
        break;
      case 'T': tMax = atof(argv[i] + 3);
        break;
      case 'N': nTimeSteps = atoi(argv[i] + 3);
        break;
      }
     else
      strcpy(meshFile, argv[i]);

  
  krnl = initKernel(0);

  // load the mesh  
  vertex *vtxs;
  face *fces;
  int nVtxs, nFces;
  loadVF(meshFile, &nVtxs, &vtxs, &nFces, &fces);
  fillFaces(nFces, fces);

  // initialize things to setup the matrix/rhs
  fX = initTriangle();
  fY = initTriangle();
  initDistancesCartesianProduct(7);
  initQuadJacobi(order, quadPath);
  initQuadSimplex(order, quadPath);


  double h = tMax/nTimeSteps;
  
  if ( job==1 ) { //print integrals only
    testIntegrals(nVtxs, nFces, fces, nTimeSteps, tMax) ;
  }
  if ( job==2 ) { //solve the linear system
   double **rhsol = setupRHS(fcnRHS, nVtxs, nFces, fces, nTimeSteps, h);
    double **A = setupMatrix(nVtxs, nFces, fces, nTimeSteps, tMax);
    
    blockForwardElimination(nFces, nTimeSteps, A, rhsol);
    
    for (int k=0; k<nTimeSteps; k++) {//
      double avg = getAverage(rhsol[k], nFces, fces); 
      double x[3] = {1.0,0.0,0.0};
      double y0;
      fcnSol(x,h*(k+0.5),&y0);
      printf("%lf %18.15lf %18.15lf %lg\n", h*k, y0, avg, avg/y0);
    }
  }
  return 0;
} /* main */







