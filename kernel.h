#ifndef KERNEL_H
#define KERNEL_H

#define SINGLE 0
#define DOUBLE 1

/* 
 * this structure gives more flexibility for what kind of kernels can be used in combination with
 * panelIA(). When running the quadrature in parallel this must be a private variable. 
 */
struct kernel {		/* kernel */
  int nParam;           /* number of parameters */
  int ansatzT;          /* temporal ansatz function: pcw const = 0, pcw linear = 1 */
  double *param;        /* extra parameters to be passed to kernel function */
                        /* kernel function */
  void (*fcn)(double *x, double *y, double *nX, double *nY, double *param, double *r, double *fcn);
};
typedef struct kernel kernel;

kernel *initKernel(int ansatzT);
void resetKernel(kernel *krnl, int type, double h, int d);
void resetKernelSlow(kernel *krnl, double h, int d);
void freeKernel(kernel *krnl);

#endif /* KERNEL_H */
