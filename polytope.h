#ifndef POLYTOPE_H
#define POLYTOPE_H

#define ZEROP 0.0001  //for testing whether a floating point number is negative
#define ZEROM -0.0001 //for testing whether a floating point number is negative
#define ONEM  0.999   //for testing whether a floating number is greater than 1.0
#define ONEP  1.0001  //for testing whether a floating number is greater than 1.0


struct poly {           // container for vertex, edge or face 
  int dim;              // dimension: 0=vertex, 1=edge, etc
  int nFces;            // number of faces 
  struct poly **fces;   // faces
  int nVtxs;            // number of vertices
  double **vtxs;        // vertices
  bool bCoordinates[3]; // describes the face in barycentric coordinates
  int idx;
};
typedef struct poly poly;



poly *initTriangle();
void setVerticesTriangle(poly *triang, double *v0, double *v1, double *v2);
void dumpPoly(int nPoly, poly *polyBuf);
void initDistancesCartesianProduct(int nPoly);
void setDistancesCartesianProduct(int nPoly, poly *fXArr, poly *fYArr);
void testDistancesCartesianProduct(int nPoly, poly *fXArr, poly *fYArr, int nSmpls);
void distPolySlow(poly *fX, poly *fY, int nSmpls, double *rMin, double *rMax);
double distPointPoint(double *x, double *y);
#endif /* POLYTOPE_H */

