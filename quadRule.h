#ifndef QUADRULE_H
#define QUADRULE_H

struct quadRule {	/* quadrature rule */
  int dim;              /* dimension of the hypercube */
  int nPts;             /* number of points */
  double *x;            /* coordinates */
  double *w;            /* weights */
};
typedef struct quadRule quadRule;


/* function prototypes */
void Jacobi(int, double, double, double*, double*);  /* Gauss-Jacobi for [0,1] */
quadRule *initQuadRuleBarycentric(int dim, int order, char *path);
quadRule *initQuadRuleJacobi(int order, int alf, int bet, char *path);
void freeQuadRule(quadRule *qr);
void rescaleGJrule0(quadRule *qr, double a, double alf, double bet, double *x, double *w);
void rescaleGJrule1(quadRule *qr, double a, double alf, double bet, double *x, double *w);

#endif /* QUADRULE_H */
