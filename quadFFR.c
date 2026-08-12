/*
 * quadFFR.c
 * Routines for the integration over a pair of faces where the kernel 
 * can be singular at r=0 and is only supported in 0<=r<=R
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include "globals.h"
#include "kernel.h"
#include "quadRule.h"
#include "mesh.h"
#include "polytope.h"

extern int verbose;
extern int nIntegrals;
extern bool doVtxApx;
extern kernel *krnl;

extern double *mBarycentX, *mBarycentY;
extern quadRule ***qrJacobi, **qrSimplex;
double *xLam, *wLam, *xEta, *wEta, *xXi, *wXi; 
extern double **rMin;
extern double **rMinV;
extern double **rMax;
extern double ***apexX;
extern double ***apexY;
extern double ***apexXV;
extern double ***apexYV;
extern double **apices;



/*
 * Initialize all Gauss-Jacobi Rules that are needed later on.
 */
void initQuadJacobi(int order, char *path) {
  qrJacobi = malloc( 4*sizeof(quadRule**) );
  for ( int i=0; i<4; i++ ) {
    qrJacobi[i] = malloc( 4*sizeof(quadRule*) );
    for ( int j=0; j<4; j++ ) {
      qrJacobi[i][j] = initQuadRuleJacobi(order, i, j, path); 
    }
  }

  int nPts = qrJacobi[0][0]->nPts;
  xLam = malloc( nPts*sizeof(double) );
  wLam = malloc( nPts*sizeof(double) );
  xEta = malloc( nPts*sizeof(double) );
  wEta = malloc( nPts*sizeof(double) );
  xXi = malloc( nPts*sizeof(double) );
  wXi = malloc( nPts*sizeof(double) );
}/* initQuadJacobi */




/*
 * Initialize generalized Gauss rules for the simplex of all
 * dimensions that are needed later on.
 */
void initQuadSimplex(int order, char *path) { 
  qrSimplex = malloc( 4*sizeof(quadRule*) );
  qrSimplex[0] = initQuadRuleBarycentric(0, order, path);
  qrSimplex[1] = initQuadRuleBarycentric(1, order, path);
  qrSimplex[2] = initQuadRuleBarycentric(2, order, path);
  qrSimplex[3] = initQuadRuleBarycentric(3, order, path);
} /* initQuadSimplex */





/*
 * Maps a quadrature point on the simplex. Based on barycentric coordinates
 * Parameters
 *     qr      quadrature rule, must be initiated by initQuadRuleBarycentric()
 *     dimS    dimension of the vertices
 *     v       the vertices
 *     k       the index of the quadrature point
 *     x       the mapped point
 *     
 */
void mapSimplexBarycentric(quadRule *qr, int dimS, double **v, int k, double *x) {
  int dim = qr->dim;
  double *lam = &qr->x[(dim+1)*k];

  for (int i=0; i<dimS; i++ ) {
    x[i] = lam[0]*v[0][i];
  }
  
  for (int k=1; k<=dim; k++ ) {
    for (int i=0; i<dimS; i++ ) {
      x[i] += lam[k]*v[k][i];
    }   
  }

} /* mapSimplexBarycentric */


void integrateNonSng(poly *fX, poly *fY, double *intgr) {
  double **v = fX->vtxs;
  double **w = fY->vtxs;
  int dimX = fX->dim;
  int dimY = fY->dim;
  quadRule *qrX = qrSimplex[dimX];
  quadRule *qrY = qrSimplex[dimY];
  void (*kernelFcn)() = krnl->fcn;

  double x[3], y[3];
  for ( int kx=0; kx<qrX->nPts; kx++ ) {
    mapSimplexBarycentric(qrX, 3, v, kx, x);
    for ( int ky=0; ky<qrY->nPts; ky++ ) {
      mapSimplexBarycentric(qrY, 3, w, ky, y);
      double wXY = qrX->w[kx]*qrY->w[ky];
      double r, kFcn;
      kernelFcn(x, y, NULL, NULL, krnl->param, &r, &kFcn);
      kFcn /= r;
      kFcn *= wXY;
      intgr[0] += kFcn;
    }
  }
  nIntegrals++;
} /* integrateNonSng */




/*
 * Perform an integral over a full convex hull. 
 */
void integrateF(int nApex, int nSngApex, double **apices, poly *fX, poly *fY, double *intgr) {
  double **v = fX->vtxs;
  double **w = fY->vtxs;
  int dimX = fX->dim;
  int dimY = fY->dim;
  int dimB = dimX + dimY;
  int dimS = nSngApex - 1;
  int dimR = nApex - nSngApex - 1;
  int dimA = nApex - 1;
  quadRule *qrX = qrSimplex[dimX];
  quadRule *qrY = qrSimplex[dimY];
  double x[3], y[3], xB[3], yB[3], zA[6], zS[6], zR[6];
  void (*kernelFcn)() = krnl->fcn;

  int quadType = 0;
  double integral = 0.0;

  if ( nSngApex==0 ) { //only regular apices
    quadType = 1;
    quadRule *qrA = qrSimplex[dimA];
    quadRule *qrJ = qrJacobi[dimB][dimA];
    for ( int kx=0; kx<qrX->nPts; kx++ ) {
      mapSimplexBarycentric(qrX, 3, v, kx, xB); 
      for ( int ky=0; ky<qrY->nPts; ky++ ) {
        mapSimplexBarycentric(qrY, 3, w, ky, yB); 
        double wXY = qrX->w[kx]*qrY->w[ky];
        for ( int kz=0; kz<qrA->nPts; kz++ ) {
          mapSimplexBarycentric(qrA, 6, apices, kz, zA); 
          double wXYZ = wXY*qrA->w[kz];
          for ( int l=0; l<qrJ->nPts; l++ ) {
            double lam  = qrJ->x[l];
            double lam1 = 1.0 - lam;
            for (int i=0; i<3; i++ ){
              x[i] = lam1*zA[i] + lam*xB[i];
              y[i] = lam1*zA[i+3] + lam*yB[i];
            }
            double r, kFcn;
            kernelFcn(x, y, NULL, NULL, krnl->param, &r, &kFcn);
            kFcn /= r;
            kFcn *= qrJ->w[l]*wXYZ;
            integral += kFcn;
          }
        }
      }
    }
  }
  else {
    quadRule *qrS = qrSimplex[dimS];
    double **aS = apices;
    double **aR = &apices[nSngApex];

    if ( dimR < 0 ){ //only singular apices
     quadType = 2;
     quadRule *qrJ  = qrJacobi[dimB-1][0];

      for ( int kx=0; kx<qrX->nPts; kx++ ) {
        mapSimplexBarycentric(qrX, 3, v, kx, xB); 
        for ( int ky=0; ky<qrY->nPts; ky++ ) {
          mapSimplexBarycentric(qrY, 3, w, ky, yB); 
          double wXY = qrX->w[kx]*qrY->w[ky];
          for ( int ks=0; ks<qrS->nPts; ks++ ) {
            mapSimplexBarycentric(qrS, 6, apices, ks, zS); 
            double wSimplex = wXY*qrS->w[ks];
            for ( int l=0; l<qrJ->nPts; l++ ) {
              double lam  = qrJ->x[l];
              double lam1 = 1.0 - lam;
              double facLam = pow(lam1,dimS);
              for (int i=0; i<3; i++ ){
                x[i] = lam1*zS[i] + lam*xB[i];
                y[i] = lam1*zS[i+3] + lam*yB[i];
              }
              double r, kFcn;
              kernelFcn(x, y, NULL, NULL, krnl->param, &r, &kFcn);
              r /= lam;
              kFcn *= qrJ->w[l]*wSimplex*facLam/r;
              integral += kFcn;
            }
          }
        }
      }
    }
    else { //singular and regular apices
      quadType = 3;
      quadRule *qrR = qrSimplex[dimR];
      quadRule *qrXi  = qrJacobi[dimR+dimB][dimS];
      quadRule *qrEt  = qrJacobi[dimB][dimR];

      for ( int kx=0; kx<qrX->nPts; kx++ ) {
        mapSimplexBarycentric(qrX, 3, v, kx, xB); 
        for ( int ky=0; ky<qrY->nPts; ky++ ) {
          mapSimplexBarycentric(qrY, 3, w, ky, yB); 
          for ( int kr=0; kr<qrR->nPts; kr++ ) {
            mapSimplexBarycentric(qrR, 6, aR, kr, zR); 
            for ( int ks=0; ks<qrS->nPts; ks++ ) {
              mapSimplexBarycentric(qrS, 6, aS, ks, zS); 
              double wSimplex = qrX->w[kx]*qrY->w[ky]*qrR->w[kr]*qrS->w[ks];
              for ( int i=0; i<qrXi->nPts; i++ ) {
                double xi = qrXi->x[i];
                double wi = qrXi->w[i];
                double lam12 = 1 - xi;
                for ( int j=0; j<qrEt->nPts; j++ ) {
                  double eta = qrEt->x[j];
                  double wj  = qrEt->w[j];
                  double lam1 = (1 - eta)*xi;
                  double lam2 = xi*eta;
 
                  for (int k=0; k<3; k++ ){
                    x[k] = lam12*zS[k]   + lam1*zR[k]   + lam2*xB[k];
                    y[k] = lam12*zS[k+3] + lam1*zR[k+3] + lam2*yB[k];
                  }
                  double r, kFcn;
                  kernelFcn(x, y, NULL, NULL, krnl->param, &r, &kFcn);
                  r /= xi;
                  kFcn *= wi*wj*wSimplex/r;
                  integral += kFcn;
                }
              }
            }
          }
        }
      }
    }
  }

  intgr[0] += integral;
  nIntegrals++;
  if ( verbose > 1 ) printf("%d %18.15le\n", quadType, integral);
 

} /* integrateF */ 




/*
 * get the upper bound in the lambda-integration
 */
double getLamMax(double *zA, double *xB, double *yB, double R2) {
  double *xA = zA;
  double *yA = &zA[3];

  double r2 = 0;
  double rv = 0;
  double v2 = 0;
  for ( int i=0; i<3; i++ ) {
    double r = xA[i] - yA[i];
    double v = xB[i] - yB[i] - r;
    r2 += r*r;
    rv += r*v;
    v2 += v*v;    
  }

  double rad = sqrt( (R2 - r2)*v2 + rv*rv );
  double lam = (rad - rv)/v2;


  if (!(lam>=0) ) {
    printf("lam=%lf\n", lam); fflush(stdout);
  }
  ASSERT(lam>=0);
  ASSERT(lam<ONEP);

  return lam;
} /* getLamMax */


double getLamMaxSimple(double *xB, double *yB, double R2){
  double r2 = 0;
  for ( int i=0; i<3; i++ ) {
    double r = xB[i] - yB[i];
    r2 += r*r;
  }

  double lam = sqrt(R2/r2);

  ASSERT(lam<ONEP);
  return lam;
} /* getLamMaxSimple*/


double getXiMax(double eta, double *zR, double *xB, double *yB, double R2){
  double eta1 = 1 - eta;
  double r2 = 0.0;
  for ( int i=0; i<3; i++ ) {
    double r = eta1*(zR[i] - zR[i+3]) + eta*(xB[i] - yB[i]);
    r2 += r*r;
  }
  double xi = sqrt(R2/r2);
  ASSERT( xi <= 1 );
  return xi;
} /* getXiMax */


/*
 * Perform an integral over the part of the convex hull where |x-y|^2 > R2.
 */
void integrateB(int nApex, int nSngApex, double **apices, poly *fX, poly *fY, double R2, double *intgr) {
  double **v = fX->vtxs;
  double **w = fY->vtxs;
  int dimX = fX->dim;
  int dimY = fY->dim;
  int dimB = dimX + dimY;
  int dimS = nSngApex - 1;
  int dimR = nApex - nSngApex - 1;
  int dimA = nApex - 1;
  int nPtsJ = qrJacobi[0][0]->nPts;
  quadRule *qrX = qrSimplex[dimX];
  quadRule *qrY = qrSimplex[dimY];
  double x[3], y[3], xB[3], yB[3], zA[6], zS[6], zR[6];
  void (*kernelFcn)() = krnl->fcn;

  int quadType = 0;
  double integral = 0.0;
  
  if ( nSngApex==0 ) { //only regular apices
    quadType = 4;
    quadRule *qrA = qrSimplex[dimA];
    for ( int kx=0; kx<qrX->nPts; kx++ ) {
      mapSimplexBarycentric(qrX, 3, v, kx, xB); 
      for ( int ky=0; ky<qrY->nPts; ky++ ) {
        mapSimplexBarycentric(qrY, 3, w, ky, yB); 
        double wXY = qrX->w[kx]*qrY->w[ky];   
        for ( int kz=0; kz<qrA->nPts; kz++ ) {
          mapSimplexBarycentric(qrA, 6, apices, kz, zA); 
          double wSimplex = wXY*qrA->w[kz];
          double lamStar = getLamMax(zA, xB, yB, R2);
          rescaleGJrule0(qrJacobi[dimB][0], lamStar, dimB, dimA, xLam, wLam);
         
          for ( int l=0; l<nPtsJ; l++ ) {
            double lam = xLam[l];
            double lam1 = 1.0 - lam;
            for (int i=0; i<3; i++ ){
              x[i] = lam1*zA[i] + lam*xB[i];
              y[i] = lam1*zA[i+3] + lam*yB[i];
            }
            double r, kFcn;
            kernelFcn(x, y, NULL, NULL, krnl->param, &r, &kFcn);
            kFcn *= wLam[l]*wSimplex/r;
            integral += kFcn;
          }
        }
      }
    }
  }
  else {
    double **aS = apices;
    double **aR = &apices[nSngApex];
    quadRule *qrS = qrSimplex[dimS];
 
    if ( dimR < 0 ) { //only singular apices
     quadType = 5;
     for ( int kx=0; kx<qrX->nPts; kx++ ) {
        mapSimplexBarycentric(qrX, 3, v, kx, xB); 
        for ( int ky=0; ky<qrY->nPts; ky++ ) {
          mapSimplexBarycentric(qrY, 3, w, ky, yB); 
          double wXY = qrX->w[kx]*qrY->w[ky];       
          double lamStar = getLamMaxSimple(xB, yB, R2); 
          rescaleGJrule0(qrJacobi[dimB-1][0], lamStar, dimB-1, dimS, xLam, wLam);
          for ( int ks=0; ks<qrS->nPts; ks++ ) {
            mapSimplexBarycentric(qrS, 6, apices, ks, zS); 
            double wSimplex = wXY*qrS->w[ks];
            for ( int l=0; l<nPtsJ; l++ ) {
              double lam  = xLam[l];
              double lam1 = 1.0 - lam;
              for (int i=0; i<3; i++ ){
                x[i] = lam1*zS[i] + lam*xB[i];
                y[i] = lam1*zS[i+3] + lam*yB[i];
              }
              double r, kFcn;
              kernelFcn(x, y, NULL, NULL, krnl->param, &r, &kFcn);
              r /= lam;
              kFcn *= wLam[l]*wSimplex/r;
              integral += kFcn;
            }
          }
        }
      }
    }
    else { //singular and regular apices
      quadType = 6;
      quadRule *qrR = qrSimplex[dimR];
 
      for ( int kx=0; kx<qrX->nPts; kx++ ) {
        mapSimplexBarycentric(qrX, 3, v, kx, xB); 
        for ( int ky=0; ky<qrY->nPts; ky++ ) {
          mapSimplexBarycentric(qrY, 3, w, ky, yB); 
          for ( int kr=0; kr<qrR->nPts; kr++ ) {
            mapSimplexBarycentric(qrR, 6, aR, kr, zR); 
            double etaStar = getLamMax(zR, xB, yB, R2);

            //first domain
            rescaleGJrule0(qrJacobi[dimB][0], etaStar, dimB, dimR, xEta, wEta);
            quadRule *qrXi = qrJacobi[dimR+dimB][dimS]; 

            for ( int ks=0; ks<qrS->nPts; ks++ ) {
              mapSimplexBarycentric(qrS, 6, aS, ks, zS); 
              double wSimplex = qrX->w[kx]*qrY->w[ky]*qrR->w[kr]*qrS->w[ks];
             
              for ( int i=0; i<nPtsJ; i++ ) {
                double eta = xEta[i];
                double wi  = wEta[i];
                for ( int j=0; j<nPtsJ; j++ ) {
                  double xi = qrXi->x[j];
                  double wj = qrXi->w[j];
                  double lam1 = (1-eta)*xi;
                  double lam2 = eta*xi;
                  double lam12 = 1 - xi;
                  for (int k=0; k<3; k++ ){
                    x[k] = lam12*zS[k]   + lam1*zR[k]   + lam2*xB[k];
                    y[k] = lam12*zS[k+3] + lam1*zR[k+3] + lam2*yB[k];
                  }
 
                  double r, kFcn;
                  kernelFcn(x, y, NULL, NULL, krnl->param, &r, &kFcn);
                  r /= xi;
                  kFcn *= wi*wj*wSimplex/r;
                  integral += kFcn;
                }
              }
            }

            //second domain
            rescaleGJrule1(qrJacobi[dimR][0], etaStar, dimB, dimR, xEta, wEta);
            for ( int ks=0; ks<qrS->nPts; ks++ ) {
              mapSimplexBarycentric(qrS, 6, aS, ks, zS); 
              double wSimplex = qrX->w[kx]*qrY->w[ky]*qrR->w[kr]*qrS->w[ks];
             
              for ( int i=0; i<nPtsJ; i++ ) {
                double eta = xEta[i];
                double wi  = wEta[i];
                double xiStar = getXiMax(eta, zR, xB, yB, R2);
                rescaleGJrule0(qrJacobi[dimR+dimB][0], xiStar, dimR+dimB, dimS, xXi, wXi);
                
                for ( int j=0; j<nPtsJ; j++ ) {
                  double xi = xXi[j];
                  double wj = wXi[j];
                  double lam1 = (1-eta)*xi;
                  double lam2 = eta*xi;
                  double lam12 = 1 - xi;
                  for (int k=0; k<3; k++ ){
                    x[k] = lam12*zS[k]   + lam1*zR[k]   + lam2*xB[k];
                    y[k] = lam12*zS[k+3] + lam1*zR[k+3] + lam2*yB[k];
                  }
                  double r, kFcn;
                  kernelFcn(x, y, NULL, NULL, krnl->param, &r, &kFcn);
                  r /= xi;
                  kFcn *= wi*wj*wSimplex/r;
                  integral += kFcn;
                }
              }
            }
          }
        }
      }
    }
  }

  intgr[0] += integral;
  nIntegrals++;
  if ( verbose > 1 ) printf("%d %18.15le\n", quadType, integral);
} /* integrateB */ 



/*
 * determine the factor in the recursive integration
 */
double getBaryFactor(poly *fX, poly *kX, double *lam) {
  for ( int i=0; i<3; i++ ) {
    if ( ( fX-> bCoordinates[i] == true ) && ( kX-> bCoordinates[i] == false) ) {
      return lam[i];
    }
  }
  ABORT();
  return 0;
} /* getBaryFactor */



/*
 * Integrate over the Cartesian product of two polytopes. Return zero if the minimal 
 * distance is greater than R. Otherwise, select an apex with distance < R and recurse
 * into the faces.
 *
 * Parameters
 *   nApex    level = number of apices in the recursion.
 *   nSngApex number of singular apices (must be the first in the array aprices)
 *   apices   array of apex pointers. Must be allocated sufficiently long.
 *   fX,fY    the current Cartesian product 
 *   R2       square of distance threshold.
 *   intgr    returned integral
 *
 */
void integrateRec(int nApex, int nSngApex, double **apices, poly *fX, poly *fY, double R2, double *intgr) {
  int idxX = fX->idx;
  int idxY = fY->idx;

  if ( rMin[idxX][idxY]>0 && rMax[idxX][idxY] <= R2 ) {
    integrateF(nApex, nSngApex, apices, fX, fY, intgr);
  }
  else if ( rMin[idxX][idxY] >= R2 && nApex>0 ) {
    integrateB(nApex, nSngApex, apices, fX, fY, R2, intgr);
  }
  else {
    double *apxX, *apxY;
    if ( rMinV[idxX][idxY] <= R2 && doVtxApx ) {
      apxX = apexXV[idxX][idxY];
      apxY = apexYV[idxX][idxY];
    }
    else {
      apxX = apexX[idxX][idxY];
      apxY = apexY[idxX][idxY];
    }
    
    memcpy(  apices[nApex],    apxX, 3*sizeof(double) );
    memcpy( &apices[nApex][3], apxY, 3*sizeof(double) );
    nApex++;
    if ( rMin[idxX][idxY] == 0 ) nSngApex++;
 
    
    double lamX[3], lamY[3];
    evalBarycentric(3, apxX, mBarycentX, lamX);
    evalBarycentric(3, apxY, mBarycentY, lamY);

    double intgrLocal;
    for ( int k=0; k<fX->nFces; k++ ) {
      poly *kX = fX->fces[k];
      double lam = getBaryFactor(fX, kX, lamX);
      if ( lam > ZEROP ) {
        intgrLocal = 0.0;
        integrateRec(nApex, nSngApex, apices, kX, fY, R2, &intgrLocal);
        intgr[0] += lam*intgrLocal;
      }
    }

    for ( int k=0; k<fY->nFces; k++ ) {
      poly *kY = fY->fces[k];
      double lam = getBaryFactor(fY, kY, lamY); 
      if ( lam > ZEROP ) {
        intgrLocal = 0.0;
        integrateRec(nApex, nSngApex, apices, fX, kY, R2, &intgrLocal);
        intgr[0] += lam*intgrLocal;
      }
    }
  }
}/* integrateRec */


/*
 * Top level routine for face-face interaction
 */
void fcefceIA(poly *fX, poly *fY, double R2, double *intgr){

  nIntegrals = 0;
  intgr[0] = 0;
  
  if ( rMin[fX->idx][fY->idx] >= R2 ) {
    return;
  }
  else if ( rMin[fX->idx][fY->idx] > 0 && rMax[fX->idx][fY->idx] <= R2 ) {
    integrateNonSng(fX, fY, intgr);
  }
  else {
    integrateRec(0, 0, apices, fX, fY, R2, intgr);
  }

} /* fcefceIA */



/*
 * Integrate the function fcn(x,t) over a face. The time value is fixed.
 * Parameters
 *    fcn    function in integrand
 *    fX     face
 *    t      time
 *    intgr  on return, the integral
 * The quadrature rules are globals 
 */
void intFceS0(void (*fcn)(double *x,double t, double *f), face *fX, double t, double *intgr){
  double *v[3], fcnVal;
  v[0] = fX->vtx[0]->x;
  v[1] = fX->vtx[1]->x;
  v[2] = fX->vtx[2]->x;
  quadRule *qrX = qrSimplex[2];

  intgr[0] = 0.0;
  for ( int i=0; i<qrX->nPts; i++ ) {
    double x[3];
    mapSimplexBarycentric(qrX, 3, v, i, x);
    fcn(x, t, &fcnVal);
    fcnVal *= qrX->w[i];
    intgr[0] += fcnVal;
  }

  intgr[0] *= fX->area2;
} /* intFceS0 */

