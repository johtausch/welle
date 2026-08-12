/*
 * mesh.c
 *   read and refine the spatial mesh, determine neighbors etc
 *
 *   copyright by Johannes Tausch
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "globals.h"
#include "quadRule.h"
#include "mesh.h"

/*
 * Reads vertices and faces 
 *
 * Parameters
 *   basename  base name of files, the actual name must have an extra P or V at the end
 *   nVtxs     number of vertices
 *   vtxs      array of vertices
 *   nFces     number of faces
 *   fces      array of faces
 *   
 */
void loadVF(char *basename, int *nVtxsP, vertex **vtxsP, int *nFcesP, face **fcesP) {
  FILE *fp;
  char file[256];
  double x0, x1, x2;
  int nFces=0, nVtxs=0, n0, n1, n2, i, nReturned;
  vertex *vtxs;
  face *fces;
  
  /*
   * read vertices 
   */
  sprintf(file,"%sV",basename); 
  if ( ( fp = fopen(file,"r") ) == NULL ) {
    printf("\ncould not open %s\n",basename);
    exit(1);
  }
  while ( fscanf(fp,"%lf %lf %lf\n",&x0,&x1,&x2)==3 ) nVtxs++ ;
  *nVtxsP = nVtxs;
  fclose(fp);

  vtxs = *vtxsP = (vertex*)calloc( nVtxs,sizeof(vertex) );

  fp = fopen(file,"r");
  for ( i=0; i<nVtxs; i++ ){
    nReturned = fscanf(fp,"%lf %lf %lf",&vtxs[i].x[0],&vtxs[i].x[1],&vtxs[i].x[2]);
    vtxs[i].idx = i;
    ASSERT( nReturned==3 );
  }
  fclose(fp);

  /*
   * read faces 
   */
  sprintf(file,"%sF",basename); 
  if ( ( fp = fopen(file,"r") ) == NULL ) {
    printf("\ncould not open %s\n",file);
    exit(1);
  }

  while ( fscanf(fp,"%d %d %d",&n0,&n1,&n2)==3 ) nFces++;
  *nFcesP = nFces;
  fclose(fp);

  fces = *fcesP = (face*)calloc( nFces,sizeof(face) );

  fp = fopen(file,"r");
  for ( i=0; i<nFces; i++ ){
    nReturned = fscanf(fp,"%d %d %d",&n0,&n1,&n2);
    ASSERT( nReturned==3 );
    fces[i].vtx[0] = &vtxs[n0];
    fces[i].vtx[1] = &vtxs[n1];
    fces[i].vtx[2] = &vtxs[n2];
    fces[i].idx = i;
  }
  fclose(fp);

} /* loadVF */




/*
 * compute the normal and the Jacobian of a face
 */
double jacobian2(face *f, double *nrm){
  double a1[3], a2[3], jac;
  int k;

  for (k=0; k<3; k++ ) {
    a1[k] = f->vtx[1]->x[k] - f->vtx[0]->x[k];
    a2[k] = f->vtx[2]->x[k] - f->vtx[0]->x[k];
  }

  nrm[0] = a1[1]*a2[2] - a1[2]*a2[1];
  nrm[1] = a1[2]*a2[0] - a1[0]*a2[2];
  nrm[2] = a1[0]*a2[1] - a1[1]*a2[0];
  jac = sqrt(SQR(nrm[0]) + SQR(nrm[1]) + SQR(nrm[2]));
  
  nrm[0] /= jac;
  nrm[1] /= jac;
  nrm[2] /= jac;
  
  return jac;
} /* jacobian2 */



/*
 * dump vertices or a function on the vertices
 */
void dumpVtxs(int nVtxs, vertex *vtxs, double *v) {
  vertex *vtx;
  int i;

  if ( v==NULL ) {
    for ( i=0; i<nVtxs; i++  ) {
      vtx = &vtxs[i];
      printf("vtx=%d (%lf,%lf,%lf) \n", vtx->idx, vtx->x[0], vtx->x[1], vtx->x[2]);  
    }
  }
  else {
    for ( i=0; i<nVtxs; i++  ) {
      vtx = &vtxs[i];
      printf("vtx=%d (%lf %lf %lf) %lf\n", vtx->idx, vtx->x[0], vtx->x[1], vtx->x[2], v[i]);    
    }
  }
} /* dumpVtxs */


/*
 * dump faces
 */
void dumpFces(int nFces, face *fces) {
  face *f;
  int i;

  for ( i=0; i<nFces; i++  ) {
    f = &fces[i];

    printf("fce=%d  vtx=[%d %d %d] area=%lf nrm=[%lf,%lf,%lf], \n", f->idx,
           f->vtx[0]->idx, f->vtx[1]->idx, f->vtx[2]->idx, 0.5*f->area2,
           f->normal[0], f->normal[1], f->normal[2]); 
  }
} /* dumpFces */




/*
 * dump a surface triangulation into the files geomV, geomF 
 */
void dumpVF(int nVtxs, vertex *vtxs, int nFces, face *fces){
  int i;
  vertex *vtx;
  face *fce;
  FILE *fp;

  fp = fopen("geomV","w");
  for ( i=0; i<nVtxs; i++ ) {
    vtx = &vtxs[i];
    fprintf(fp, "%lf %lf %lf\n", vtx->x[0],vtx->x[1],vtx->x[2]);
  }
  fclose(fp);

  fp = fopen("geomF","w");
  for ( i=0; i<nFces; i++ ) {
    fce = &fces[i];
    fprintf(fp, "%d %d %d\n", fce->vtx[0]->idx, fce->vtx[1]->idx, fce->vtx[2]->idx);
  }
  fclose(fp);

  //printf("Wrote the surface triangulation into the files geomV and geomF\n");
} /* dumpVF */



/* 
 * computes the crossproduct of the vectors
 */
void crossPrd(double *a, double *b, double *nrm){
  nrm[0] = a[1]*b[2] - a[2]*b[1];
  nrm[1] = a[2]*b[0] - a[0]*b[2];
  nrm[2] = a[0]*b[1] - a[1]*b[0];
} /* crossPrd */



/* 
 * computes the normalized crossproduct of the vectors and returns the length of a x b
 */
double nrmCrossPrd(double *a, double *b, double *nrm){
  double len;
  nrm[0] = a[1]*b[2] - a[2]*b[1];
  nrm[1] = a[2]*b[0] - a[0]*b[2];
  nrm[2] = a[0]*b[1] - a[1]*b[0];
  len = sqrt(SQR(nrm[0]) + SQR(nrm[1]) + SQR(nrm[2]));
  nrm[0] /= len;
  nrm[1] /= len;
  nrm[2] /= len;
  return len;
} /* nrmCrssPrd */



/*
 * Set up the affine transformation that computes the barycentric
 * coordinates (or the piecewise linear shape functions) of the triangle with
 * vertices v0, v1, v2 in either two or three dimensional space.
 * The returned vector is the parameter m of the function evalBarycentric()
 * which evaluates the affine transformation for a given point. 
 */
double *setupBarycentric(int sDim, double *v0, double *v1, double *v2){

  double a1[3], a2[3];
  double a11 = 0;
  double a12 = 0;
  double a22 = 0;;
  for ( int k=0; k<sDim; k++ ) {
    a1[k] = v1[k] - v0[k];
    a2[k] = v2[k] - v0[k];
    a11 += a1[k]*a1[k];
    a12 += a1[k]*a2[k];
    a22 += a2[k]*a2[k];
  }
  
  double det = a11*a22 - a12*a12;
  double b11 = a22/det;
  double b12 = -a12/det;
  double b22 = a11/det;

  double m1[3], m2[3];

  for ( int k=0; k<sDim; k++ ) {
    m1[k] = a1[k]*b11 + a2[k]*b12;
    m2[k] = a1[k]*b12 + a2[k]*b22;
  }

  double *m = calloc( 2*sDim+2, sizeof(double) );
      
  for ( int k=0; k<sDim; k++ ) {
    m[0] -= m1[k]*v0[k];
    m[1] -= m2[k]*v0[k];
  }

  if ( sDim==2 ) {
    m[2] = m1[0];
    m[3] = m1[1];
    m[4] = m2[0];
    m[5] = m2[1];
  }
  else if ( sDim==3 ) {
    m[2] = m1[0];
    m[3] = m1[1];
    m[4] = m1[2];
    m[5] = m2[0];
    m[6] = m2[1];
    m[7] = m2[2];
  }
  else {
    ABORT();
  }
  return m;
} /* setupBarycentric */

/*
 * Evaluate the barycentric coordinates of a triangle using an affine transformation.
 * Parameters
 *     sDim   space dimension (either 2 or 3)
 *     x      evaluation point
 *     m      components of the transformation as computed by setupBarycentric()
 *     lam    barycentric coordinates
 *            if sDim==3 then lam are the coordinates of the orthogonal projection
 *            into the triangle plane.
 */
void evalBarycentric(int sDim, double *x, double *m, double *lam){
  double *m1 = &m[2];
  double *m2 = &m[2+sDim];

  lam[1] = m[0];
  lam[2] = m[1];
  for ( int k=0; k<sDim; k++ ) {
    lam[1] += m1[k]*x[k];
    lam[2] += m2[k]*x[k];
  }

  lam[0] = 1 - lam[1] - lam[2];
} /* evalBarycentric */



/*
 * Map the standard simplex to any simplex.
 * Parameters
 *     dimS   space dimension = number of components of the vertices
 *     dimX   dimension of the simplex
 *     v      coordinates of the dimX+1 vertices
 *     t      coordinates in the standard simplex
 *     x      the mapped point
 */
void mapSimplex(int dimS, int dimX, double *v, double *t, double *x) {

  double dt = 1.0 - t[0];  
  for (int i=0; i<dimS; i++ ) {
    x[i] = dt*v[i];
  }
  
  for (int k=1; k<dimX; k++ ) {
    double *vk = &v[dimS*k];
    dt = t[k-1] - t[k];
    for (int i=0; i<dimS; i++ ) {
      x[i] += dt*vk[i];
    }   
  }

  double *vk = &v[dimS*dimX];
  dt = t[dimX-1];
  for (int i=0; i<dimS; i++ ) {
    x[i] += dt*vk[i];
  }   
 
} /* mapSimplex */

/*
 * same as mapSimplex, just that the vertices are passed as a double array
 */
void mapSimplex2(int dimS, int dimX, double **v, double *t, double *x) {

  double dt = 1.0 - t[0];  
  for (int i=0; i<dimS; i++ ) {
    x[i] = dt*v[0][i];
  }
  
  for (int k=1; k<dimX; k++ ) {
    dt = t[k-1] - t[k];
    for (int i=0; i<dimS; i++ ) {
      x[i] += dt*v[k][i];
    }   
  }

  dt = t[dimX-1];
  for (int i=0; i<dimS; i++ ) {
    x[i] += dt*v[dimX][i];
  }   
 
} /* mapSimplex2 */


/*
 * compute additional face members
 */
void fillFaces(int nFces, face *fces){

  for ( int k=0; k<nFces; k++ ) {
    double *vtx0 = fces[k].vtx[0]->x;
    double *vtx1 = fces[k].vtx[1]->x;
    double *vtx2 = fces[k].vtx[2]->x;
    double a2[3], a0[3];
   
    a2[0] = vtx1[0] - vtx0[0];
    a2[1] = vtx1[1] - vtx0[1];
    a2[2] = vtx1[2] - vtx0[2];

    a0[0] = vtx2[0] - vtx1[0];
    a0[1] = vtx2[1] - vtx1[1];
    a0[2] = vtx2[2] - vtx1[2];

    fces[k].area2 = nrmCrossPrd(a2, a0, fces[k].normal);

    fces[k].mBarycent = setupBarycentric(3, vtx0, vtx1, vtx2);    
  } /* for k */
  
} /* fillFaces */



