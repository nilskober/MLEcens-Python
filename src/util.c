/* Util.c: collection of functions that are used by the other routines */

#include <R.h>
#include <Rinternals.h>

#include "MLEcens.h"


/* Function to check if input has the correct format */
void VerifyInputRectangles(SEXP RR, SEXP BB)
{
   int i,n;
   double *pRR;
   int *pBB;

   n = nrows(RR);

   pRR = REAL(RR);
   pBB = INTEGER(BB);

   if (!isMatrix(RR) || ncols(RR)!=4)
      error("invalid argument 'R': it must be a matrix with 4 columns\n");

   if (isMatrix(BB))
   {
      if (!(ncols(BB)==2 || ncols(BB)==4))
         error("argument 'B' has invalid number of columns\n");
      if (!(nrows(BB)==1 || nrows(BB)==n))
         error("argument 'B' has invalid number of rows\n");
   }
   else if (!(length(BB)==2 || length(BB)==4))
      error("argument 'B' has invalid length\n");

   for (i=0; i<n; i++)
   {
      if (pRR[i] > pRR[i+n])
         error("invalid argument 'R': x1 is larger than x2 in R[%d,]\n", i+1);
      else if (pRR[i] == pRR[i+n])
      {
         if (isMatrix(BB))
         {
            if (pBB[i]!=1 || pBB[i+n]!=1)
               error("x1==x2 in R[%d,], so boundaries for these endpoints need to be specified as 1=closed\n", i+1);
         }
       else if (pBB[0]!=1 || pBB[1]!=1)
            error("x1==x2 in R[%d,], so boundaries for these endpoints need to be specified as 1=closed\n", i+1);
      }

      if (pRR[i+2*n] > pRR[i+3*n])
         error("invalid argument 'R': y1 is larger than y2 in R[%d,]\n", i+1);
      else if (pRR[i+2*n] == pRR[i+3*n])
      {
         if (isMatrix(BB))
         {
            if (pBB[i+2*n]!=1 || pBB[i+3*n]!=1)
               error("y1==y2 in R[%d,],so boundaries for these endpoints need to be specified as 1=closed\n", i+1);
         }
         else if (pBB[length(BB)-2]!=1 || pBB[length(BB)-1]!=1)
            error("y1==y2 in R[%d,], so boundaries for these endpoints need to be specified as 1=closed\n", i+1);
      }
   }

   for (i=0; i<length(BB); i++)
   {
      if (pBB[i]!=0 && pBB[i]!=1)
         error("argument 'B' may only contain 0's and 1's\n");
   }
}

void VerifyInputCanonicalRectangles(SEXP CanonRects, SEXP RR, SEXP BB)
{
   int i,m,n;
   int *pCanonRects, *pBB;
   double *pRR;

   m = nrows(CanonRects);
   n = nrows(RR);

   pCanonRects = INTEGER(CanonRects);
   pRR         = REAL(RR);
   pBB         = INTEGER(BB);

   if (!isMatrix(CanonRects) || ncols(CanonRects)!=4)
      error("invalid argument 'Rcanon': it must be a matrix with 4 columns\n");

   if (!isMatrix(RR) || ncols(RR)!=4)
      error("invalid argument 'R': it must be a matrix with 4 columns\n");

   if (isMatrix(BB))
   {
      if (!(ncols(BB)==2 || ncols(BB)==4))
         error("argument 'B' has invalid number of columns\n");
      if (!(nrows(BB)==1 || nrows(BB)==n))
         error("argument 'B' has invalid number of rows\n");
   }
   else if (!(length(BB)==2 || length(BB)==4))
      error("argument 'B' has invalid length\n");

   for (i=0; i<m; i++)
   {
      if (pCanonRects[i] > pCanonRects[i+m])
         error("invalid argument 'Rcanon': x1 is larger than x2 in Rcanon[%d,]\n", i+1);
      if (pCanonRects[i+2*m] > pCanonRects[i+3*m])
         error("invalid argument 'Rcanon': y1 is larger than y2 in Rcanon[%d,]\n", i+1);
   }

   for (i=0; i<n; i++)
   {
      if (pRR[i] > pRR[i+n])
         error("invalid argument 'R': x1 is larger than x2 in R[%d,]\n", i+1);
      else if (pRR[i] == pRR[i+n])
      {
         if (isMatrix(BB))
         {
            if (pBB[i]!=1 || pBB[i+n]!=1)
               error("x1==x2 in R[%d,], so boundaries for these endpoints need to be specified as 1=closed\n", i+1);
         }
         else if (pBB[0]!=1 || pBB[1]!=1)
            error("x1==x2 in R[%d,], so boundaries for these endpoints need to be specified as 1=closed\n", i+1);
      }

      if (pRR[i+2*n] > pRR[i+3*n])
         error("invalid argument 'R': y1 is larger than y2 in R[%d,]\n", i+1);
      else if (pRR[i+2*n] == pRR[i+3*n])
      {
         if (isMatrix(BB))
         {
            if (pBB[i+2*n]!=1 || pBB[i+3*n]!=1)
               error("y1==y2 in R[%d,] so boundaries for these endpoints need to be specified as 1=closed\n", i+1);
         }
         else if (pBB[length(BB)-2]!=1 || pBB[length(BB)-1]!=1)
            error("y1==y2 in R[%d,] so boundaries for these endpoints need to be specified as 1=closed\n", i+1);
      }
   }

   for (i=0; i<length(BB); i++)
   {
      if (pBB[i]!=0 && pBB[i]!=1)
         error("argument 'B' may only contain 0's and 1's\n");
   }
}


/* Function to compare endpoints A and B.
 * If A<B, then return value is 1, otherwise the return value is 0.
 * This function is used in SortEndPoints(), which is used in RealToCanonical().
 */
int CompareEndpoints (const SEndPoint *A, const SEndPoint *B)
{
   if (A->value != B->value)
      return (A->value < B->value);

   if (A->left == B->left && A->closed == B->closed)
      return (A->index < B->index);

   if (A->left != B->left && A->closed != B->closed)
      return (B->left);

   return (A->left == A->closed);
}


/* Function to compare endpoints A and B.
 * If elem1<elem2, then return value is -1, otherwise return value is 1.
 * This function uses CompareEndpoints() to do the actual comparison,
 * and it is used in RealToCanonical().
 */
int SortEndpoints(const void *A, const void *B)
{
   if (CompareEndpoints((SEndPoint*) A, (SEndPoint*) B))
      return -1;
   else
      return 1;
}


/* Transform real rectangles to canonical rectangles.
 * Input:
 *   n: Number of observation rectangles
 *   pRR: Vector of length 4*n with coordinates of the observation rectangles
 *        The coordinates x1,x2,y1,y2 of the ith rectangle are given by
 *        pRR[i+0*n], pRR[i+1*n], pRR[i+2*n], pRR[i+3*n]
 *   pBB: Vector of length 2, 4, or 4*n with the boundary structure of the observation
 *         rectangles (0=open, 1=closed).
 *        - if length=4*n, pBB represents a 4xn matrix, stored in the same way as R
 *        - if length=4, all observation rectangles have the same boundary structure
 *        - if length=2, all x-intervals and all y-intervals have the same boundary structure
 *    LengthBB: total length of pBB (4*n, 4, or 2). This tells us how pBB is specified.
 * Output:
 *   CanonRects: nx4 matrix of canonical rectangles, stored as *SCanonRect
 *   rx: vector of length 2*n containing the indices of the rectangles
 *        corresponding to the sorted x-coordinates
 *    ry: vector of length 2*n containing the indices of the rectangles
 *        corresponding to the sorted y-coordinates
 *    lb: vector of length 2*n indicating the left/right boundaries of the
 *        sorted x-coordinates (1=left, 0=right)
 */
void RealToCanonical(int n, double *pRR, int *pBB, SCanonRect *CanonRects,
           int *rx, int *ry, int *lb, int LengthBB)
{
   int i, BBexplicit;

   SEndPoint *XEndPoints = Calloc(2*n, SEndPoint);
   SEndPoint *YEndPoints = Calloc(2*n, SEndPoint);
   int       *BBvalue    = Calloc(4, int);
   /* is used when boundary matrix BB is not given explicitly */

   BBexplicit = (LengthBB==4*n); /* indicates wether BB is given explicitly */

   if (!BBexplicit)
   {
      BBvalue[0] = pBB[0];
      BBvalue[1] = pBB[1];
      BBvalue[2] = pBB[LengthBB-2+0];
      BBvalue[3] = pBB[LengthBB-2+1];
   }

   /* Separate observations into X and Y endpoints.
   * Note that I have used more explicit names than in the paper,
   * and instead of k=1,2 to indicate left,right endpoint, I now use left=1,0.
   */
   for (i=0; i<n; i++)
   {
      XEndPoints[i*2].value  = pRR[i+0*n];
      XEndPoints[i*2].left   = 1;
      XEndPoints[i*2].index  = i;

      XEndPoints[i*2 + 1].value  = pRR[i+1*n];
      XEndPoints[i*2 + 1].left   = 0;
      XEndPoints[i*2 + 1].index  = i;

      YEndPoints[i*2].value  = pRR[i+2*n];
      YEndPoints[i*2].left   = 1;
      YEndPoints[i*2].index  = i;

      YEndPoints[i*2 + 1].value  = pRR[i+3*n];
      YEndPoints[i*2 + 1].left   = 0;
      YEndPoints[i*2 + 1].index  = i;
   }
   if (BBexplicit)
   {
      for (i=0; i<n; i++)
      {
         XEndPoints[i*2].closed     = pBB[i+0*n];
         YEndPoints[i*2].closed     = pBB[i+2*n];
         XEndPoints[i*2 + 1].closed = pBB[i+1*n];
         YEndPoints[i*2 + 1].closed = pBB[i+3*n];
      }
   }
   else
   {
      for (i=0; i<n; i++)
      {
         XEndPoints[i*2].closed     = BBvalue[0];
         YEndPoints[i*2].closed     = BBvalue[2];
         XEndPoints[i*2 + 1].closed = BBvalue[1];
         YEndPoints[i*2 + 1].closed = BBvalue[3];
      }
   }

   /* Sort endpoint arrays using quicksort with SortEndpoints() */
   qsort(XEndPoints, n*2, sizeof(SEndPoint), SortEndpoints);
   qsort(YEndPoints, n*2, sizeof(SEndPoint), SortEndpoints);

   /* Transform the sorted arrays into canonical rectangles */
   for (i=0; i<2*n; i++)
   {
      if (XEndPoints[i].left)
         CanonRects[XEndPoints[i].index].x1 = i;
      else
         CanonRects[XEndPoints[i].index].x2 = i;

      if (YEndPoints[i].left)
         CanonRects[YEndPoints[i].index].y1 = i;
      else
         CanonRects[YEndPoints[i].index].y2 = i;
   }

   for (i=0; i<2*n; i++)
   {
      rx[i] = XEndPoints[i].index;
      ry[i] = YEndPoints[i].index;
      lb[i] = XEndPoints[i].left;
   }

    Free(XEndPoints);
    Free(YEndPoints);
    Free(BBvalue);
}

/* Transform canonical rectangles back to original coordinates.
 * Input:
 *   CanonRects: matrix of canonical rectangles, stored as SCanonRect
 *   mm: number of rows of CanonRects
 *   m: number of canonical rectangles that are to be transformed
 *   ind: if (m<mm), ind is a vector of length m indicating which
 *        canonical rectangles are to be transformed
 *   n: number of original rectangles
 *   pRR: Vector of length 4*n with coordinates of the observation rectangles
 *        The coordinates x1,x2,y1,y2 of the ith rectangle are given by
 *        pRR[i+0*n], pRR[i+1*n], pRR[i+2*n], pRR[i+3*n]
 *   pBB: Vector of length 2, 4, or 4*n with the boundary structure of the observation
 *         rectangles (0=open, 1=closed).
 *        - if length=4*n, pBB represents a 4xn matrix, stored in the same way as R
 *        - if length=4, all observation rectangles have the same boundary structure
 *        - if length=2, all x-intervals and all y-intervals have the same boundary structure
 *   BBexplicit: indicates whether pBB is specified as a matrix
 *   rx: vector of length 2*n containing the indices of the rectangles
 *        corresponding to the sorted x-coordinates
 *    ry: vector of length 2*n containing the indices of the rectangles
 *        corresponding to the sorted y-coordinates
 * Output:
 *   pAnsRects: vector of length 4*m with the real coordinates of the
 *              canonical rectangles that were to be transformed (stored as pRR)
 *   pAnsBounds: boundary structure of pAnsRects,
 *               specified in the same way as pBB
 */
void CanonicalToReal(SCanonRect *CanonRects, int mm, int m, int *ind,
                int n, double *pRR, int *pBB, int BBexplicit, int *rx, int *ry,
                double *pAnsRects, int *pAnsBounds)
{
   int i;

   if (m==mm) /* all CanonRects are to be transformed */
   {
      for (i=0; i<mm; i++)
      {
         pAnsRects[i+0*mm] = pRR[rx[CanonRects[i].x1] + 0*n];
         pAnsRects[i+1*mm] = pRR[rx[CanonRects[i].x2] + 1*n];
         pAnsRects[i+2*mm] = pRR[ry[CanonRects[i].y1] + 2*n];
         pAnsRects[i+3*mm] = pRR[ry[CanonRects[i].y2] + 3*n];
      }

      if (BBexplicit)
      {
         for (i=0; i<mm; i++)
         {
            pAnsBounds[i+0*mm] = pBB[rx[CanonRects[i].x1] + 0*n];
            pAnsBounds[i+1*mm] = pBB[rx[CanonRects[i].x2] + 1*n];
            pAnsBounds[i+2*mm] = pBB[ry[CanonRects[i].y1] + 2*n];
            pAnsBounds[i+3*mm] = pBB[ry[CanonRects[i].y2] + 3*n];
         }
      }
   }

   else /* subset of CanonRects are to be transformed */
   {
      for (i=0; i<m; i++)
      {
         pAnsRects[i + 0*m] = pRR[rx[CanonRects[ind[i]].x1] + 0*n];
         pAnsRects[i + 1*m] = pRR[rx[CanonRects[ind[i]].x2] + 1*n];
         pAnsRects[i + 2*m] = pRR[ry[CanonRects[ind[i]].y1] + 2*n];
         pAnsRects[i + 3*m] = pRR[ry[CanonRects[ind[i]].y2] + 3*n];
      }

      if (BBexplicit)
      {
         for (i=0; i<m; i++)
         {
            pAnsBounds[i + 0*m] = pBB[rx[CanonRects[ind[i]].x1] + 0*n];
            pAnsBounds[i + 1*m] = pBB[rx[CanonRects[ind[i]].x2] + 1*n];
            pAnsBounds[i + 2*m] = pBB[ry[CanonRects[ind[i]].y1] + 2*n];
            pAnsBounds[i + 3*m] = pBB[ry[CanonRects[ind[i]].y2] + 3*n];
         }
      }
   }
}


/* LAPACK function used in SolveSymmetricLinearSystem() */
void F77_NAME (dspsv)(char *, int *, int *, double *, int *, double *,
                 int *, int *);


/* Compute convex combination of two arrays of doubles of length n:
 *   result = p*array1 + (1-p)*array2
 */
void Convex_Comb(int ndata, double p, double array1[], double array2[],
               double result[])
{
   int i;

   for (i=0; i<ndata; i++)
      result[i] = p*array1[i] + (1-p)*array2[i];
}


/* Solve linear system A*X = B using the FORTRAN LAPACK routines DSPSV
 * A is a SYMMETRIC SQUARE matrix in PACKED STORAGE
 *   Its size is nrowA*(nrowA+1)/2, its upper diagonal matrix
 *     (including the diagonal) is stored column wise
 *   On exit, A contains the block diagonal matrix and multipliers
 *      from the factorization of A
 * B is an (nrowA*ncolB) matrix of doubles
 *   On exit, B contains the solution matrix
 * dummy is integer array of length nrowA, that is needed in the function dspsv
 * All matrices are stored in arrays, by column, starting at zero
 */
void SolveSymmetricLinearSystem(double *A, int nrowA, double *B, int ncolB, int* dummy)
{
   int info;
   char uplo='U';

   F77_CALL (dspsv) (&uplo, &nrowA, &ncolB, A, dummy, B,
               &nrowA, &info);

   if (info!=0)
   {
      Rf_error("Error in SolveSymmetricLinearSystem");
   }
}

/* vim:set et ts=3 sw=3: */
