/*
 * quadRule.c
 *    set up various quadrature rules
 *   Copyright  Johannes Tausch
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "globals.h"
#include "quadRule.h"



/*
 *  Initialize a Gauss-Jacobi quadrature rule 
 *  
 *  Parameters
 *     order   quadrature order 
 *     alf     coefficient in weight function
 *     bet     coefficient in weight function
 *     path    path where the quadrature file is located
 *  Returns a struct that contains precomputed quadrature nodes weights
 *  for this integral
 *
 *    \int_0^1 x^alf (1-x)^bet f(x) dx =~ \sum_{k=0}^{p-1} f(x[k])w[k]
 *
 */
quadRule *initQuadRuleJacobi(int order, int alf, int bet, char *path) {
  char buffer[1024];
  int tmp, nConv;
  char *str;

  quadRule *qr = malloc( sizeof(quadRule) );
  qr->dim = 1;

  if ( order%2 == 0 ) {//there are no even order rules, take the next odd order
    sprintf(buffer, "%s/quad_GJ_a%d_b%d_deg%d.txt", path, alf, bet, order+1);
  }
  else {
    sprintf(buffer, "%s/quad_GJ_a%d_b%d_deg%d.txt", path, alf, bet, order);
  }
  
  FILE *fp = fopen(buffer, "r");
  if ( fp==NULL ) {
    printf("loadQuadRuleJacobi(): cannot open %s\n",buffer);
    exit(1);
  }

  str = fgets(buffer, sizeof(buffer), fp);  //ignore the first line
  ASSERT( str != NULL );
  nConv = fscanf(fp, "%s %d\n", buffer, &tmp); //ignore the order
  nConv = fscanf(fp, "%s %d\n", buffer, &qr->dim); //dimension
  ASSERT(nConv==2);
  ASSERT(qr->dim == 1);
  nConv = fscanf(fp, "%s %d\n", buffer, &qr->nPts); //number points
  ASSERT(nConv==2);

  double *w = qr->w = calloc( qr->nPts, sizeof(double) );
  double *x = qr->x = calloc( qr->nPts, sizeof(double) );

  for ( int i=0; i<qr->nPts; i++ ) {
    nConv = fscanf(fp, "%lf %lf\n", &w[i], &x[i]);
    ASSERT(nConv==2);
  }

  return qr;
} /* initQuadRuleJacobi */



/*
 * rescale a Gauss-Jacobi rule from [0,1] to [0,a]
 * Parameters
 *    qr           input rule, must be a Gauss-Jacobi rule with weight w(x)=x^alf in [0,1]
 *    a, alf, bet  ouput rule parameters: w(x)=x^alf (1-x)^bet in [0,a]
 *    x, w         nodes and weights of the output rule
 */
void rescaleGJrule0(quadRule *qr, double a, double alf, double bet, double *x, double *w){
  int order = qr->nPts;
  double aFac = pow(a, 1+alf);
  
  for ( int i=0; i<order; i++ ) {
    x[i] = a*qr->x[i];
    double bFac = pow(1-x[i], bet);
    w[i] = aFac*bFac*qr->w[i];
  }  
} /* rescaleGJrule0 */

/*
 * rescale a Gauss-Jacobi rule from [0,1] to [a,1]
 * Parameters
 *    qr           input rule, must be a Gauss-Jacobi rule with weight w(x)=x^bet in [0,1]
 *    a, alf, bet  ouput rule parameters: w(x)=x^alf (1-x)^bet in [0,a]
 *    x, w         nodes and weights of the output rule
 */
void rescaleGJrule1(quadRule *qr, double a, double alf, double bet, double *x, double *w){
  int order = qr->nPts;
  double a1 = 1 - a;
  double bFac = pow(a1, 1+bet);
  
  for ( int i=0; i<order; i++ ) {
    x[i] = 1 - a1*qr->x[i];
    double aFac = pow(x[i], alf);
    w[i] = aFac*bFac*qr->w[i];
  }  
} /* rescaleGJrule1 */




/*
 *  Initialize a quadrature rule for the standard simplex. The nodes
 *  are returned in barycentric coordinates.
 *
 *  Parameters
 *     dim     dimension (must be <= 3)
 *     order   quadrature order 
 *     path    path where the quadrature file is located
 *  Returns a struct that contains precomputed quadrature nodes weights
 */
quadRule *initQuadRuleBarycentric(int dim, int order, char *path) {
  char buffer[1024], *str;
  double tmp;
  int nConv;

  quadRule *qr = (quadRule *)malloc( sizeof(quadRule) );
  int dim1 = dim+1;

  if ( dim==0 ) {
    qr->dim = 0;
    qr->nPts = 1;
    qr->w = malloc( sizeof(double) );
    qr->x = malloc( sizeof(double) );
    qr->w[0] = 1.0;
    qr->x[0] = 1.0;
  }
  else if ( dim==1 ) { //load Gauss-Legendre rule
    if ( order%2 == 0 ) {//there are no even order rules, take the next odd order
      sprintf(buffer, "%s/quad_GJ_a0_b0_deg%d.txt", path, order+1);
    }
    else {
      sprintf(buffer, "%s/quad_GJ_a0_b0_deg%d.txt", path, order);
    }
    FILE *fp = fopen(buffer, "r");
    if ( fp==NULL ) {
      printf("loadQuadRuleBarycentric(): cannot open %s\n",buffer);
      exit(1);
    }
    str = fgets(buffer, sizeof(buffer), fp);  //ignore the entire first line
    ASSERT( str != NULL );
    nConv = fscanf(fp, "%s %lf\n", buffer, &tmp); //ignore the second line with order
    nConv = fscanf(fp, "%s %d\n", buffer, &qr->dim); //dimension 
    ASSERT(nConv == 2);
    nConv = fscanf(fp, "%s %d\n", buffer, &qr->nPts); //number points
    ASSERT(qr->dim == 1);
    ASSERT(nConv == 2);

    double *w = qr->w = calloc( qr->nPts, sizeof(double) );
    double *x = qr->x = calloc( dim1*qr->nPts, sizeof(double) );
   
    for ( int i=0; i<qr->nPts; i++ ) {
      double t;
      nConv = fscanf(fp, "%lf %lf\n", &w[i], &t);
      ASSERT(nConv == 2);
      double *xi = &x[2*i];
      xi[0] = 1 - t;
      xi[1] = t;
    }
  }
  else { // load the proper rule for the standard simplex genrated by Arkadijs code
    sprintf(buffer, "%s/quad_simplex_deg%d_dim%d.txt", path, order, dim);
    FILE *fp = fopen(buffer, "r");
    if ( fp==NULL ) {
      printf("loadQuadRuleBarycentric(): cannot open %s\n",buffer);
      exit(1);
    }

    nConv = fscanf(fp, "%s %lf\n", buffer, &tmp); //ignore the first line with efficiency
    nConv = fscanf(fp, "%s %lf\n", buffer, &tmp); //ignore the second line with order
    nConv = fscanf(fp, "%s %d\n", buffer, &qr->dim); //dimension 
    ASSERT(nConv == 2);
    nConv = fscanf(fp, "%s %d\n", buffer, &qr->nPts); //number points
    ASSERT(qr->dim == dim);
    ASSERT(nConv == 2);
 
    double *w = qr->w = calloc( qr->nPts, sizeof(double) );
    double *x = qr->x = calloc( dim1*qr->nPts, sizeof(double) );
    double t[3];
    
    for ( int i=0; i<qr->nPts; i++ ) {
      if ( qr->dim == 2 ){
        nConv = fscanf(fp, "%s %lf %lf %lf\n", buffer, &w[i], &t[0], &t[1]);
        ASSERT(nConv == 4);
      }
      else if ( qr->dim == 3 ){
        nConv = fscanf(fp, "%s %lf %lf %lf %lf\n", buffer, &w[i], &t[0], &t[1], &t[2]);
        ASSERT(nConv == 5);
      }
      else {
        ABORT();
      }  

      double *xi = &x[i*dim1];
      xi[0] = 1 - t[0];
      for (int d=1; d<dim; d++ ) {
        xi[d] = t[d-1] - t[d];
      }
      xi[dim] = t[dim-1];  
    }
  }
  
  return qr;
} /* initQuadRuleBarycentric */



/*
 * free a quadrule
 */
void freeQuadRule(quadRule *qr){
  free(qr->x);
  free(qr->w);
} /* freeQuadRule */



/*
 * print a quadrule
 */
void printQuadRule(quadRule *qr){
  int n=qr->nPts, d=qr->dim;
  int k,l;
  
  for ( k=0; k<n; k++ ){
    printf("%d x=[",k);
    for ( l=0; l<d; l++ ) {
      printf("%lf,", qr->x[k*d+l]);
    }
    printf("] w=%lf\n", qr->w[k]);  
  }
} /* printQuadRule */



