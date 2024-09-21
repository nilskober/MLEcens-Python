/* real2canon.c */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <R.h>
#include <Rinternals.h>

#include "MLEcens.h"

/* Transform rectangles into canonical rectangles.
 * For R, I made the canonical rectangles have coordinates in {1,..,2n}
 */
SEXP RealToCanonicalForR(SEXP RR, SEXP BB)
{
   int i, numprotected=0;
   int n = nrows(RR);   /* retrieve n from R-object RR */
   SEXP Ans;
   int  *pAns;

   double *pRR = REAL(RR);
   int    *pBB = INTEGER(BB);
   int    BBexplicit;	/* equals 1 if boundary matrix BB is given explicitly */

   SEndPoint  *XEndPoints = R_Calloc(2*n, SEndPoint);
   SEndPoint  *YEndPoints = R_Calloc(2*n, SEndPoint);
   int        *BBvalue    = R_Calloc(4, int);   
   /* BBvalue is used when boundary matrix BB is not given explicitly */

   VerifyInputRectangles(RR,BB);

   /* allocate memory for output: canonical rectangles */
   PROTECT(Ans = allocMatrix(INTSXP, n, 4));   
   numprotected++;
   pAns = INTEGER(Ans);

   BBexplicit = (isMatrix(BB) && nrows(BB)==n);
   if (!BBexplicit)
   {
      BBvalue[0] = pBB[0];
      BBvalue[1] = pBB[1];
      BBvalue[2] = pBB[length(BB)-2+0];
      BBvalue[3] = pBB[length(BB)-2+1];
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

      YEndPoints[i*2].value  = pRR[i+2*n];
      YEndPoints[i*2].left   = 1;
      YEndPoints[i*2].index  = i;

      XEndPoints[i*2 + 1].value  = pRR[i+1*n];
      XEndPoints[i*2 + 1].left   = 0;
      XEndPoints[i*2 + 1].index  = i;

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

   /* Transform the sorted arrays into canonical rectangles. 
    * For R, I let the canonical rectangles have coordinates in {1,..,2n}. 
    * Because of this, right hand side is i+1 instead of i.
    */
   for (i=0; i<2*n; i++)
   {
      if (XEndPoints[i].left)
         pAns[XEndPoints[i].index + 0*n] = i+1;
      else
         pAns[XEndPoints[i].index + 1*n] = i+1;

      if (YEndPoints[i].left)
         pAns[YEndPoints[i].index + 2*n] = i+1;
      else
         pAns[YEndPoints[i].index + 3*n] = i+1;
   }

   R_Free(XEndPoints);
   R_Free(YEndPoints);
   R_Free(BBvalue);

   UNPROTECT(numprotected);
   return Ans;
}

/* vim:set et ts=3 sw=3: */
