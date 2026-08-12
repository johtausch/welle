/*
 * solver.h
 * copyright   Johannes Tausch
 */
#ifndef SOLVER_H
#define SOLVER_H

void dgetrf_(int *m, int *n, double *A, int *lda, int *ipiv, int *info );

void dgetrs_(char *trans, int *n, int *nrhs, double *A, int *lda, int *ipiv,
             double *b, int *ldb, int *info );
void dgemv_(char *tr, int *m, int *n, double *alpha, double *A, int *lda, 
          double *x, int *incx, double *beta, double *y, int *incy);

int blockForwardElimination(int nS, int nT, double **A, double **rhsol);
void print_matrix(char* desc, char* format, int m, int n, double* a, int lda );


#endif /* SOLVER_H */

