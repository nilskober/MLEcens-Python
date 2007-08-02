/* MLEcens.h */

#ifndef MLECENS_H
#define MLECENS_H

#include <string.h>   /* for function memcpy */

#include <Rinternals.h>

/* Define constants */
/********************/

#define INFTY 1e6


/* Define functions */
/********************/

#define MemCopy(dest,source,n)  memcpy( dest, source, (size_t)( (n) * sizeof(*(dest)) ) )
#define MemAlloc(n,t)   (t *) chk_malloc( (size_t) (n), sizeof(t) )
#define MemCalloc(n,t)  (t *) chk_calloc( (size_t) (n), sizeof(t) )
#define MemFree(p)      (free( (void *)(p) ), (p) = NULL)
#define SQR(a)          (a*a)
#define IsInRectangle(p,r) ( ((p)->x >  (r)->x1) && \
                             ((p)->x <= (r)->x2) && \
                             ((p)->y >  (r)->y1) && \
                             ((p)->y <= (r)->y2) ) 

/* Typedefs: */
/*************/

/* Type definition of a point with double precision coordinates */
typedef struct { 
   double x, y; 
} SDoublePoint;

/* Type definition of a point with integer coordinates */
typedef struct { 
   int x, y; 
} SIntPoint;

/* Type definition of a rectangle with double precision coordinates
 * and specification of open/closed boundaries */
typedef struct {
  double  x1,  x2,  y1,  y2; /* x1<=x2, y1<=y2 */
  int   cx1, cx2, cy1, cy2;  /* 0=open, 1=closed */
} SRealRect;

/* Type definition of a rectangle with integer coordinates 
 * and no boundary information */
typedef struct { 
   int x1, x2, y1, y2; 
} SCanonRect;

/* Type definition of an endpoint: 
 *   value is coordinate value of endpoint
 *   closed=1 if endpoint is closed
 *   left=1 if endpoint is left endpoint
 *   index is index of the rectangle to which endpoint belongs
 */
typedef struct { 
   double value; 
   int closed, left; 
   int index;  
} SEndPoint;


/* Function declarations */
/*************************/

/* Functions in util.c */
void VerifyInputRectangles(SEXP RR, SEXP BB);
void VerifyInputCanonicalRectangles(SEXP CanonRects, SEXP RR, SEXP BB);
int  CompareEndpoints (const SEndPoint *A, const SEndPoint *B);
int  SortEndpoints(const void *A, const void *B);
void RealToCanonical(int n, double *pRR, int *pBB, SCanonRect *CanonRects, 
                     int *rx, int *ry, int *lb, int LengthBB);
void CanonicalToReal(SCanonRect *CanonRects, int mm, int m, int *ind,  
                     int n, double *pRR, int *pBB, int BBexplicit, int *rx, int *ry, 
                     double *pAnsRects, int *pAnsBounds);
void Convex_Comb(int ndata, double p, double array1[], double array2[], 
                 double result[]);
void SolveSymmetricLinearSystem(double *A, int nrowA, double *B, int ncolB,
                                int *dummy);

/* Functions in real2canon.c */
SEXP RealToCanonicalForR(SEXP RR, SEXP BB);

/* Functions in canon2real.c */
SEXP CanonicalToRealForR(SEXP InputRects, SEXP RR, SEXP BB);

/* Functions in reduc.c */
SEXP ReductionStepForR(SEXP RR, SEXP BB, SEXP hm, SEXP cm);

/* Functions in computeMLE.c */
SEXP ComuteMLEForR(SEXP RR, SEXP BB, SEXP tol);

#endif /* MLECENS_H */

/* vim:set et ts=3 sw=3: */
