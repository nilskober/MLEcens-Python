/* canon2real.c */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <R.h>
#include <Rinternals.h>

#include "MLEcens.h"

SEXP CanonicalToRealForR(SEXP InputRects, SEXP RR, SEXP BB)
{
   int    i, *pAnsBounds, BBexplicit, BBvalue[4];
   int    numprotected=0;
   double *pAnsRects;
   char   *names[2] = {"rects","bounds"};
   SEXP   AnsRects, AnsBounds, ListNames, Ans;

   int m = nrows(InputRects);
   int n = nrows(RR);
 
   int    *pInputRects = INTEGER(InputRects);     
   double *pRR         = REAL(RR);
   int    *pBB         = INTEGER(BB);
 
   VerifyInputCanonicalRectangles(InputRects, RR, BB);

   SEndPoint  *XEndPoints = Calloc(2*n, SEndPoint);
   SEndPoint  *YEndPoints = Calloc(2*n, SEndPoint);

   /* allocate space for output */
   PROTECT(AnsRects = allocMatrix(REALSXP, m, 4));  
   numprotected++;
   pAnsRects = REAL(AnsRects);

   PROTECT(AnsBounds = allocMatrix(INTSXP, m, 4));  
   numprotected++;
   pAnsBounds = INTEGER(AnsBounds);

   BBexplicit = (isMatrix(BB) && nrows(BB)==n);

   if (!BBexplicit)
   {
      BBvalue[0] = pBB[0];
      BBvalue[1] = pBB[1];
      BBvalue[2] = pBB[length(BB)-2+0];
      BBvalue[3] = pBB[length(BB)-2+1];
   }
	
   /* Separate observation rectangles of RR into X and Y endpoints.
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

   for (i=0; i<m; i++)
   { 
      pAnsRects[i+0*m]  = XEndPoints[pInputRects[i+0*m]-1].value;
      pAnsRects[i+1*m]  = XEndPoints[pInputRects[i+1*m]-1].value;
      pAnsRects[i+2*m]  = YEndPoints[pInputRects[i+2*m]-1].value;
      pAnsRects[i+3*m]  = YEndPoints[pInputRects[i+3*m]-1].value;
   }

   /* Output list containing the transformed rectangles and their boundaries */
   PROTECT(ListNames = allocVector(STRSXP, 2));  
   numprotected++;
   for(i=0; i<2; i++)
      SET_STRING_ELT(ListNames,i,mkChar(names[i]));
 
   PROTECT(Ans = allocVector(VECSXP, 2));     
   numprotected++;
   SET_VECTOR_ELT(Ans,0,AnsRects);            
  
   if(BBexplicit)
   {
      for (i=0; i<m; i++){
         pAnsBounds[i+0*m] = XEndPoints[pInputRects[i+0*m]-1].closed;
         pAnsBounds[i+1*m] = XEndPoints[pInputRects[i+1*m]-1].closed;
         pAnsBounds[i+2*m] = YEndPoints[pInputRects[i+2*m]-1].closed;
         pAnsBounds[i+3*m] = YEndPoints[pInputRects[i+3*m]-1].closed;
      }
      SET_VECTOR_ELT(Ans,1,AnsBounds);    
   }
   else
      SET_VECTOR_ELT(Ans,1,BB);
   
   setAttrib(Ans,R_NamesSymbol,ListNames);   /* attach vector names */

   Free(XEndPoints);
   Free(YEndPoints);

   UNPROTECT(numprotected);

   return Ans;
}

