/*
 * LAPACK based direct linear system solvers 
 * copyright   Johannes Tausch
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "solver.h"

extern int verbose;

/*
 * Solve a block lower triangular linear system with Toeplitz structure
 * Parameters
 *    nS      size of the blocks
 *    nT      number of blocks
 *    A       array of nS x nS matrices, where A[k] is the block in the k-th subdiagonal
 *    rhssol  on entrance the right hand side,  on return the solution
 */
int blockForwardElimination(int nS, int nT, double **A, double **rhsol){
  int info;
  int nRow = nS, nCol = nS, lda = nS;
  int *ipiv = malloc( nRow*sizeof(int) );
  char trans = 'N';
  struct timespec time1;
  timespec_get(&time1, TIME_UTC);
  /* Compute the LU factorization */
  dgetrf_(&nRow, &nCol, A[0], &lda, ipiv, &info);
  
  if( info != 0 ) {
    printf( "dgetrf-error, code=%d\n", info );
    return info;
  }

  int nrhs = 1;
  int inc = 1;
  int ldb = lda;
  double one = 1.0, negOne = -1.0;
  
  double *y = calloc(nS, sizeof(double));
  for (int k=0; k<nT; k++) {

    memcpy(y, rhsol[k], nS*sizeof(double) );
    for ( int l=1; l<=k; l++ ) {
      if ( A[k-l] != NULL ) {
        dgemv_(&trans, &nRow, &nCol, &negOne, A[k-l], &nRow, rhsol[l], &inc, &one, y, &inc);
      }
    }
    dgetrs_(&trans, &nRow, &nrhs, A[0], &lda, ipiv, y, &ldb, &info);
    memcpy(rhsol[k], y, nS*sizeof(double) );
  }
  free(y);
  free(ipiv);

  struct timespec time2;
  timespec_get(&time2, TIME_UTC);
  double solveTime = time2.tv_sec - time1.tv_sec + 1e-9*(time2.tv_nsec - time1.tv_nsec);
  if ( verbose > 0 ) printf("solveTime = %le\n", solveTime);
  return info;
} /* blockForwardElimination */

  

/*  printing a matrix
 *  desc is the name of the matrix as printed, e.g. desc="A = ["
 *  format is the printf-format string for a floating point number, e.g. format=" %6.2f"
 */
void print_matrix(char* desc, char* format, int m, int n, double* a, int lda ) {
  int i, j;
  printf( "\n %s\n", desc );
  for( i = 0; i < m; i++ ) {
    for( j = 0; j < n; j++ ) printf( format, a[i+j*lda] );
    printf( "\n" );
  }
} /* print_matrix */



