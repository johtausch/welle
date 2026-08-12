/*
 * polytope.c:
 * routines related to convex polytopes and determining
 * minimal distances for their Cartesian products
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include "globals.h"
#include "polytope.h"
#include "mesh.h"

extern double **rMin;
extern double **rMinV;
extern double **rMax;
extern double ***apexX;
extern double ***apexY;
extern double ***apexXV;
extern double ***apexYV;
extern double **apices;



/*
 * Initiate a poly structure for a triangle.
 * This allocates space for all faces, face pointers and vertex pointers
 * and assigns the face pointers.
 * The vertex pointers for an actual triangle are assigned by the
 * routine setVerticesTriangle() below.
 */
poly *initTriangle(){
  poly *triang = calloc( 7, sizeof(poly) );
  poly** fceBuf = calloc( 9, sizeof(poly*) );
  double** vtxBuf = calloc( 12, sizeof(double*));

  triang[0].dim = 2;
  triang[0].nFces = 3;
  triang[0].fces = fceBuf;
  triang[0].nVtxs = 3;
  triang[0].vtxs = vtxBuf;

  for (int i=1; i<=3; i++ ) {
    triang[i].dim = 1;
    triang[i].nFces = 2;
    triang[i].fces = &fceBuf[2*i+1];
    triang[i].nVtxs = 2;
    triang[i].vtxs = &vtxBuf[2*i+1];
  }
  for (int i=4; i<=6; i++ ) {
    triang[i].dim = 0;
    triang[i].nFces = 0;
    triang[i].nVtxs = 1;
    triang[i].vtxs = &vtxBuf[i+5];
  }

  for (int i=0; i<3; i++ ) {
    triang[0].fces[i] = &triang[i+1];
  }
  triang[1].fces[0] = &triang[5];
  triang[1].fces[1] = &triang[6];
 
  triang[2].fces[0] = &triang[4];
  triang[2].fces[1] = &triang[6];

  triang[3].fces[0] = &triang[4];
  triang[3].fces[1] = &triang[5];

  for ( int k=0; k<7; k++ ) {
    triang[k].idx = k;
  };

  memcpy(triang[0].bCoordinates, (bool[]){true,true,true}, 3*sizeof(bool) );
  memcpy(triang[1].bCoordinates, (bool[]){false,true,true}, 3*sizeof(bool) );
  memcpy(triang[2].bCoordinates, (bool[]){true,false,true}, 3*sizeof(bool) );
  memcpy(triang[3].bCoordinates, (bool[]){true,true,false}, 3*sizeof(bool) );
  memcpy(triang[4].bCoordinates, (bool[]){true,false,false}, 3*sizeof(bool) );
  memcpy(triang[5].bCoordinates, (bool[]){false,true,false}, 3*sizeof(bool) );
  memcpy(triang[6].bCoordinates, (bool[]){false,false,true}, 3*sizeof(bool) );

  return triang;
} /* initTriangle */


/*
 * fill the vertex pointers for a triangle
 */
void setVerticesTriangle(poly *triang, double *v0, double *v1, double *v2) {
  double *vtx[3] = {v0,v1,v2};
  int nAllFaces = 7;

  for ( int i=0; i<nAllFaces; i++ ) {
    int c=0;
    for ( int k=0; k<3; k++ ) {
      if ( triang[i].bCoordinates[k] ) {
        triang[i].vtxs[c++] = vtx[k];
      }
    }
  }
} /* setVerticesTriangle */


/*
 * print all faces of a polytope
 */
void dumpPoly(int nPoly, poly *polyBuf){
  for (int i=0; i<nPoly; i++ ) {
    poly *f = &polyBuf[i];
    printf("Face %d,  dim=%d  ", f->idx, f->dim);
    printf("    faces=[ ");
    for (int j=0; j<f->nFces; j++ ) {
      printf("%d ", f->fces[j]->idx);
    }
    printf("]\n");
    for (int j=0; j<f->nVtxs; j++ ) {
      printf("[ ");
      for (int k=0; k<3; k++ ) printf("%lf ", f->vtxs[j][k] );
      printf("]"); 
    }
    printf("\n");
  }
} /* dumpPoly */



/*
 * Returns the square of the distance of two points
 */
double distPointPoint(double *x, double *y){
  double r2 = 0.0;
  for (int i=0; i<3; i++ ) {
    double r = x[i] - y[i];
    r2 += r*r;
  }
  return r2;
} /* distPointPoint */


/*
 * Determines the distance of a point to an edge.
 * If the point is closer to the edge than to the endpoints the square of
 * the distance is returned, otherwise the return value is -1.
 *
 * Parameters
 *    x    the point
 *    v0   first endpoint of the edge
 *    v1   second endpoint of the edge
 *    v    the minimizing point
 *
 * v is changed only if the minimizer is in the interior of the edge
 */
double distPointEdge(double *x, double *v0, double *v1, double *v){
  double rw = 0.0;
  double ww = 0.0;
  for (int i=0; i<3; i++ ) {
    double r = x[i] - v0[i];
    double w = v1[i] - v0[i];
    rw += r*w;
    ww += w*w;
  }
  double t = rw/ww;
  if ( t<ZEROP ) return -1;
  if ( t>ONEM  ) return -1;

  double lam[2];
  lam[0] = 1.0 - t;
  lam[1] = t;

  double r2 = 0;
  for (int i=0; i<3; i++ ) {
    v[i] = lam[0]*v0[i] + lam[1]*v1[i];
    double r = x[i] - v[i];
    r2 += r*r;
  }
  
  return r2;
} /* distPointEdge */


/*
 * Determines the distance of two edges.
 * If the edges are closer than their respective endpoints the 
 * square of the distance is returned, otherwise the return value is -1.
 *
 * Parameters
 *    v0   first endpoint of edge 1
 *    v1   second endpoint of edge 1
 *    w0   first endpoint of edge 1
 *    w1   second endpoint of edge 2
 *    v    minimizing point on edge 1 
 *    w    minimizing point on edge 2
 *
 * lam v and w are changed only if the minimizer is in the interior of both edges
 */
double distEdgeEdge(double *v0, double *v1, double *w0, double *w1, double *v, double *w){
  double aa = 0.0;
  double ab = 0.0;
  double bb = 0.0;
  double ra = 0.0;
  double rb = 0.0;
  
  for (int i=0; i<3; i++ ) {
    double r = v0[i] - w0[i];
    double a = v1[i] - v0[i];
    double b = w0[i] - w1[i];
    aa += a*a;
    ab += a*b;
    bb += b*b;
    ra += r*a;
    rb += r*b;
  }
  double det = aa*bb - ab*ab;
  double s=0, t=0;
  if ( fabs(det)>1e-12 ) {
    s = (ab*rb - bb*ra)/det;
    t  = (ab*ra - aa*rb)/det;
  }

  if ( s<ZEROP ) return -1;
  if ( t<ZEROP ) return -1;
  if ( s>ONEM ) return -1;
  if ( t>ONEM ) return -1;

  double lam[2], mu[2];
  lam[0] = 1.0 - s;
  lam[1] = s;
  
  mu[0] = 1.0 - t;
  mu[1] = t;
  
  double r2 = 0;
  for (int i=0; i<3; i++ ) {
    v[i] = lam[0]*v0[i] + lam[1]*v1[i];
    w[i] = mu[0]*w0[i]  + mu[1]*w1[i];
    double r = v[i] - w[i];
    r2 += r*r;
  }
  return r2;
    
} /* distEdgeEdge */


/*
 * Determines the distance of a point to a triangular face.
 * If the face is closer than their respective edges the 
 * square of the distance is returned, otherwise the return value is -1.
 *
 * Parameters
 *    x    the point
 *    v0   first vertex of the face
 *    v1   first vertex of the face
 *    v2   first vertex of the face
 *    v    minimizing point in the face
 *
 * v is changed only if the minimizer is in the interior of the face
 */
double distPointFace(double *x, double *v0, double *v1, double *v2, double *v){
  double aa = 0.0;
  double ab = 0.0;
  double bb = 0.0;
  double ra = 0.0;
  double rb = 0.0;
  
  for (int i=0; i<3; i++ ) {
    double r = v0[i] - x[i];
    double a = v1[i] - v0[i];
    double b = v2[i] - v0[i];
    aa += a*a;
    ab += a*b;
    bb += b*b;
    ra += r*a;
    rb += r*b;
  }
  double det = aa*bb - ab*ab;
  double s = (ab*rb - bb*ra)/det;
  double t  = (ab*ra - aa*rb)/det;

  if ( s<ZEROP ) return -1;
  if ( t<ZEROP ) return -1;
  if ( s+t>ONEM ) return -1;

  double lam[3];
  lam[0] = 1 - s - t;
  lam[1] = s;
  lam[2] = t;

  double r2 = 0;
  for (int i=0; i<3; i++ ) {
    v[i] = lam[0]*v0[i] + lam[1]*v1[i] + lam[2]*v2[i];
    double r = x[i] - v[i];
    r2 += r*r;
  }
  return r2;
} /* distPointFace */


/* 
 * Allocate memory for the distance calculations
 * of two polytopes.
 */
void initDistancesCartesianProduct(int nPoly){
  rMin   = calloc(nPoly, sizeof(double*));
  rMinV  = calloc(nPoly, sizeof(double*));
  rMax   = calloc(nPoly, sizeof(double*));
  apexX  = calloc(nPoly, sizeof(double**));
  apexY  = calloc(nPoly, sizeof(double**));
  apexXV = calloc(nPoly, sizeof(double**));
  apexYV = calloc(nPoly, sizeof(double**));

  for ( int i=0; i<nPoly; i++ ) {
    rMin[i]   = calloc(nPoly, sizeof(double));
    rMinV[i]  = calloc(nPoly, sizeof(double));
    rMax[i]   = calloc(nPoly, sizeof(double));
    apexXV[i] = calloc(nPoly, sizeof(double*));
    apexYV[i] = calloc(nPoly, sizeof(double*));
    apexX[i]  = calloc(nPoly, sizeof(double*));
    apexY[i]  = calloc(nPoly, sizeof(double*));
    for ( int j=0; j<nPoly; j++ ) {
      apexX[i][j] = calloc(3, sizeof(double));
      apexY[i][j] = calloc(3, sizeof(double));
      apexXV[i][j] = calloc(3, sizeof(double));
      apexYV[i][j] = calloc(3, sizeof(double));
   }
  }
  
  apices = calloc(4, sizeof(double*));
  for (int i=0; i<4; i++ ) {
    apices[i] = calloc(6, sizeof(double));
  }
} /* initDistancesCartesianProduct */



/*
 * Find the maximum and minimum distances based on facets
 * Update the corresponding entries in the rMin, rMinV, rMax and apex arrays
 */
void minMaxFacets( poly *fX, poly *fY){
  rMin[fX->idx][fY->idx]  = 1e100;
  rMinV[fX->idx][fY->idx] = 1e100;
  rMax[fX->idx][fY->idx]  = 0;
 
  for ( int i=0; i<fX->nFces; i++ ) {
    poly *kX = fX->fces[i];

    if ( rMin[kX->idx][fY->idx]  < rMin[fX->idx][fY->idx] ) {
      rMin[fX->idx][fY->idx] = rMin[kX->idx][fY->idx];
      memcpy(apexX[fX->idx][fY->idx], apexX[kX->idx][fY->idx], 3*sizeof(double) );
      memcpy(apexY[fX->idx][fY->idx], apexY[kX->idx][fY->idx], 3*sizeof(double) );
    }
    if ( rMinV[kX->idx][fY->idx]  < rMinV[fX->idx][fY->idx] ) {
      rMinV[fX->idx][fY->idx] = rMinV[kX->idx][fY->idx];
      memcpy(apexXV[fX->idx][fY->idx], apexXV[kX->idx][fY->idx], 3*sizeof(double) );
      memcpy(apexYV[fX->idx][fY->idx], apexYV[kX->idx][fY->idx], 3*sizeof(double) );
     }
    if ( rMax[kX->idx][fY->idx] > rMax[fX->idx][fY->idx] ) {
      rMax[fX->idx][fY->idx] = rMax[kX->idx][fY->idx];
    }    
  }
  
  for ( int j=0; j<fY->nFces; j++ ) {
    poly *kY = fY->fces[j];

    if ( rMin[fX->idx][kY->idx]  < rMin[fX->idx][fY->idx] ) {
      rMin[fX->idx][fY->idx] = rMin[fX->idx][kY->idx];
      memcpy(apexX[fX->idx][fY->idx], apexX[fX->idx][kY->idx], 3*sizeof(double) );
      memcpy(apexY[fX->idx][fY->idx], apexY[fX->idx][kY->idx], 3*sizeof(double) );
    }
    if ( rMinV[fX->idx][kY->idx]  < rMinV[fX->idx][fY->idx] ) {
      rMinV[fX->idx][fY->idx] = rMinV[fX->idx][kY->idx];
      memcpy(apexXV[fX->idx][fY->idx], apexXV[fX->idx][kY->idx], 3*sizeof(double) );
      memcpy(apexYV[fX->idx][fY->idx], apexYV[fX->idx][kY->idx], 3*sizeof(double) );
    }
    if ( rMax[fX->idx][kY->idx] > rMax[fX->idx][fY->idx] ) {
      rMax[fX->idx][fY->idx] = rMax[fX->idx][kY->idx];
    }
  }

} /* minMaxFacets */



/* 
 * Evaluate the members  rMax, rMin,  rMinV, apex, apexV of the poly2 structure.
 * This can be done after calling the routines initCartesianProduct()
 * and setVerticesTriangle().
 *
 * Parameters:
 *  nPoly2 number of faces of the cartesian product
 *  fXY    the array of poly2 stuctures. The 0-th element is the Cartesian
 *         product of two poly's. The faces are sorted degreasingly by
 *         dimension. The vertices come last.
 */
void setDistancesCartesianProduct(int nPoly, poly *fXArr, poly *fYArr){
  poly *fX, *fY;
  int n1 = 3, idx1[3] = {1,2,3};   //number and which faces are one-dimensional (=edges)
  int n0 = 3, idx0[3] = {4,5,6};   //number and which faces are zero-dimensional (=vertices)
  double r, v[3], w[3]; 

  // vertex x vertex, dim=0
  for ( int ii=0; ii<n0; ii++ ) {
    int i = idx0[ii];
    fX = &fXArr[i];

    for ( int jj=0; jj<n0; jj++ ) {
      int j = idx0[jj];
      fY = &fYArr[j];
      r = distPointPoint(fX->vtxs[0],fY->vtxs[0]);
      rMax[i][j]  = r;
      rMin[i][j]  = r;
      rMinV[i][j] = r;
      memcpy( apexX[i][j], fX->vtxs[0], 3*sizeof(double) );
      memcpy( apexY[i][j], fY->vtxs[0], 3*sizeof(double) );
      memcpy( apexXV[i][j], fX->vtxs[0], 3*sizeof(double) );
      memcpy( apexYV[i][j], fY->vtxs[0], 3*sizeof(double) );
    }
  }

  
  // vertex x edge, dim=1
  for ( int ii=0; ii<n0; ii++ ) {
    int i = idx0[ii];
    fX = &fXArr[i];

    for ( int jj=0; jj<n1; jj++ ) {
      int j = idx1[jj];
      fY = &fYArr[j];
      minMaxFacets(fX, fY);

      r = distPointEdge(fX->vtxs[0], fY->vtxs[0], fY->vtxs[1], w);
      if ( r>ZEROM ) {
        rMin[i][j] = r;
        memcpy(apexX[i][j], fX->vtxs[0], 3*sizeof(double) );
        memcpy(apexY[i][j], w, 3*sizeof(double) );
      }
    }
  }
  
  // edge x vertex, dim=1
  for ( int ii=0; ii<n1; ii++ ) {
    int i = idx1[ii];
    fX = &fXArr[i];

    for ( int jj=0; jj<n0; jj++ ) {
      int j = idx0[jj];
      fY = &fYArr[j];
      minMaxFacets(fX, fY);

      r = distPointEdge(fY->vtxs[0], fX->vtxs[0], fX->vtxs[1], v);
      if ( r>ZEROM ) {
        rMin[i][j] = r;
        memcpy(apexX[i][j], v, 3*sizeof(double) );
        memcpy(apexY[i][j], fY->vtxs[0], 3*sizeof(double) );
      }
    }
  }
  
  // vertex x face, dim=2
  fY = &fYArr[0];
  for ( int ii=0; ii<n0; ii++ ) {
    int i = idx0[ii];
    fX = &fXArr[i];
    minMaxFacets(fX, fY);
    r = distPointFace(fX->vtxs[0], fY->vtxs[0], fY->vtxs[1], fY->vtxs[2], w);
    if ( r>ZEROM ) {
      rMin[i][0] = r;
      memcpy(apexX[i][0], fX->vtxs[0], 3*sizeof(double) );
      memcpy(apexY[i][0], w, 3*sizeof(double) );
    }
  }

  // edge x edge, dim=2
  for ( int ii=0; ii<n1; ii++ ) {
    int i = idx1[ii];
    fX = &fXArr[i];

    for ( int jj=0; jj<n1; jj++ ) {
      int j = idx1[jj];
      fY = &fYArr[j];
      minMaxFacets(fX, fY);
      
      r = distEdgeEdge(fX->vtxs[0], fX->vtxs[1], fY->vtxs[0], fY->vtxs[1], v, w);
      if ( r>ZEROM ) {
        rMin[i][j] = r;
        memcpy(apexX[i][j], v, 3*sizeof(double) );
        memcpy(apexY[i][j], w, 3*sizeof(double) );
      }
    }
  }

  // face x vertex, dim=2
  fX = &fXArr[0];
  for ( int jj=0; jj<n0; jj++ ) {
    int j = idx0[jj];
    fY = &fYArr[j];
    minMaxFacets(fX, fY);
    double r = distPointFace(fY->vtxs[0], fX->vtxs[0], fX->vtxs[1], fX->vtxs[2], w);
    if ( r>ZEROM ) {
      rMin[0][j] = r;
      memcpy(apexX[0][j], v, 3*sizeof(double) );
      memcpy(apexY[0][j], fY->vtxs[0], 3*sizeof(double) );
    }
 }

  // edge x face, dim=3
  fY = &fYArr[0];
  for ( int ii=0; ii<n1; ii++ ) {
    int i = idx1[ii];
    fX = &fXArr[i];
    minMaxFacets(fX, fY);
  }

  // face x edge, dim=3
  fX = &fXArr[0];
  for ( int jj=0; jj<n1; jj++ ) {
    int j = idx1[jj];
    fY = &fYArr[j];
    minMaxFacets(fX, fY);
  }

  // face x face, dim=3
  fX = &fXArr[0];
  fY = &fYArr[0];
  minMaxFacets(fX, fY);

} /* setDistancesCartesianProduct */




/*
 * Use sampling to determine the distance of two polytopes.
 * In this implementation fX and fY must be simplices of dimension less or equal two.
 */
void distPolySlow(poly *fX, poly *fY, int nSmpls, double *rMin, double *rMax){
  double a[4][3]={{0,0,0},{0,0,0},{0,0,0},{0,0,0}};
  double lam[4], x[3], y[3], r;
  double *v0 = fX->vtxs[0];
  double *w0 = fY->vtxs[0];
  double *v1 = NULL, *v2 = NULL, *w1 = NULL, *w2 = NULL;
  int nSmplsX = 0, nSmplsY = 0;
 
  if ( fX->dim > 0 ) {
    nSmplsX = nSmpls;
    v1 = fX->vtxs[1];
    if ( fX->dim > 1 ) {
      v2 = fX->vtxs[2];
      if ( fX->dim > 2 ) ABORT();
    }
  }

  if ( fY->dim > 0 ) {
    nSmplsY = nSmpls;
    w1 = fY->vtxs[1];
    if ( fY->dim > 1 ) {
      w2 = fY->vtxs[2];
      if ( fY->dim > 2 ) ABORT();
    }
  }

  for ( int k=0; k<3; k++ ) {
    if ( fX->dim > 0 ) a[0][k] = v1[k] - v0[k];
    if ( fX->dim > 1 ) a[1][k] = v2[k] - v1[k];
    if ( fY->dim > 0 ) a[2][k] = w1[k] - w0[k];
    if ( fY->dim > 1 ) a[3][k] = w2[k] - w1[k];
  }
   
    
  rMin[0] = 1e100,  rMax[0]=0;
  for ( int i0=0; i0<=nSmplsX; i0++ ) {
    lam[0] = (double)i0/nSmpls;
    int i00 = i0;
    if ( fX->dim < 2 ) i00 = 0;
    for ( int i1=0; i1<=i00; i1++ ) {
      lam[1] = (double)i1/nSmpls;
      for ( int j0=0; j0<=nSmplsY; j0++ ) {
        lam[2] = (double)j0/nSmpls;
        int j00 = j0;
        if ( fY->dim < 2 ) j00 = 0;
        for ( int j1=0; j1<=j00; j1++ ) {
          lam[3] = (double)j1/nSmpls;
          for ( int k=0; k<3; k++ ) {
            x[k] = v0[k] + lam[0]*a[0][k] + lam[1]*a[1][k];
            y[k] = w0[k] + lam[2]*a[2][k] + lam[3]*a[3][k];
          }
          r = distPointPoint(x,y);
          if ( r>rMax[0] ) rMax[0] = r;
          if ( r<rMin[0] ) rMin[0] = r;
        }
      }
    }
  }
}/* distPolySlow */

  
void testDistancesCartesianProduct(int nPoly, poly *fXArr, poly *fYArr, int nSmpls) {
  
  for ( int i=0; i<nPoly; i++ ) {
    poly *fX = &fXArr[i];
    for ( int j=0; j<nPoly; j++ ) {
      poly *fY = &fYArr[j];

      printf("Face (%d,%d) dim=%d ", fX->idx, fY->idx, fX->dim + fY->dim);
      double rMinS, rMaxS;  
      distPolySlow(fX, fY, nSmpls, &rMinS, &rMaxS);
      double eMin = rMinS - rMin[i][j];
      double eMax = rMax[i][j] - rMaxS;
      printf("rMin=%lf eMin=%lg, rMax=%lf eMax=%lg\n", rMinS, eMin, rMaxS, eMax);
    }
  }
} /* testDistancesCartesianProduct */


/*
 * Returns true if a point is on a face and false if the point is on the face
 * Parameters
 *    f   the face
 *    lam the point in barycentric coordinates
 */
bool isNotOnFace( poly *f, double *lam ){
  for ( int k=0; k<3; k++ ) {
    if ( (f->bCoordinates[k]==false) && lam[k] > ZEROP ) return true; 
  }
  return false;
} /* isNotOnFace */
