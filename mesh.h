#ifndef MESH_H
#define MESH_H
/*
 * panel, vertex and edge definitions
 * prototypes for reading and refining the mesh
 */

struct vertex {		/* vertex */
  double x[3];          /* coordinates */
  int idx;              /* index in array */
};
typedef struct vertex vertex;


struct face {		/* face, always triangles */
  vertex *vtx[3];       /* pointers to vertices */
  double normal[3];     /* normal, determined by right-hand rule */
  double area2;         /* Jacobian = 2*area */
  double *mBarycent;    /* parameters for calculating the barycentric coordinates */
  int idx;              /* index in array */
};
typedef struct face face;


// function prototypes 
void loadVF(char *basename, int *nVtxsP, vertex **vtxsP, int *nFcesP, face **fcesP);
void fillFaces(int nFces, face *fces);
void getSurfaceMesh(int nVtxsV, vertex *vtxsV, int *nVtxsSp, vertex **vtxsSp,
                    int nFcesV, face *fcesV,  int *nFcesSp, face **fcesSp);
void dumpVtxs(int nVtxs, vertex *vtxs, double *v);
void dumpFces(int nFces, face *fces);
void dumpVF(int nVtxs, vertex *vtxs, int nFces, face *fces);
double *setupBarycentric(int sDim, double *v0, double *v1, double *v2);
void evalBarycentric(int sDim, double *x, double *sFcnParam, double *phi);
void mapSimplex(int dimS, int dimX, double *v, double *t, double *x);
void mapSimplex2(int dimS, int dimX, double **v, double *t, double *x);
#endif /* MESH_H */




