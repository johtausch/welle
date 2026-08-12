#ifndef GLOBALS_H
#define GLOBALS_H

/* Constants */
#define M_4PI     12.566370614359172954         /* (4*pi) */
#define FOURPII   0.079577471545947667882       /* (4*pi)^(-1) */
#define TWOPII    0.15915494309189533577        /* (2*pi)^(-1) */

#ifndef MAX
#define MAX(A,B)  ( (A) > (B) ? (A) : (B) )
#endif

#ifndef MIN
#define MIN(A,B)  ( (A) > (B) ? (B) : (A) )
#endif

#define ABS(A) ( ( (A) > 0 ) ? (A) : (-(A)) )
#define SQR(x) ( (x)*(x) )
#define DIST2(V,W) (SQR((V)[0]-(W)[0]) + SQR((V)[1]-(W)[1]) + SQR((V)[2]-(W)[2]))
#define PARITY(A) ( (A)%2==0 ? 1 : -1   )


#define  ABORT() {\
  fflush(stdout);\
  fprintf(stdout, "panic in file `%s' at line %d.\n", __FILE__, __LINE__);\
  exit(1);\
}

#define ASSERT(condition) if(!(condition)) ABORT()


#endif /* GLOBALS_H */
