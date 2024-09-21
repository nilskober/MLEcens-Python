/* reduc.c */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <R.h>
#include <Rinternals.h>

#include "MLEcens.h"

/* Algorithm to find maximal intersections: the areas of possible mass support.
 * The code follows the pseudo code given in Maathuis (2005, JCGS),  
 * except for the fact that all arrays start from zero here. 
 */
SEXP ReductionStepForR(SEXP RR, SEXP BB, SEXP hm, SEXP cm)
{
   int n = nrows(RR);
   int *h  = R_Calloc(2*n, int);   /* value of height map, initialized at zero */
   int *e  = R_Calloc(2*n, int);   /* value of last entered rectangle */
   int *rx = R_Calloc(2*n, int);   /* indices of rectangles corresponding to 
                                  *   sorted X endpoints */
   int *ry = R_Calloc(2*n, int);   /* indices of rectangles corresponding to 
                                  *   sorted Y endpoints */
   int *lb = R_Calloc(2*n, int);   /* left/right boundary (1/0) */

   /* canonical observation rectangles */
   SCanonRect *CanonObsRectangles    = R_Calloc(n, SCanonRect);  

   /* canonical maximal intersections */
   int MaxNrMIs = n;  /* starting value, will be increased if necessary */
   SCanonRect *CanonMaxIntersections = R_Calloc(MaxNrMIs, SCanonRect); 

   /* variables for input */
   double *pRR = REAL(RR);
   int    *pBB = INTEGER(BB);
   int    ihm = *INTEGER(hm);
   int    icm = *INTEGER(cm);
   int    BBexplicit = (isMatrix(BB) && nrows(BB)==n);
   int    LengthBB = length(BB);

   /* variables for output/answer: 
    * MI=Maximal intersections, HM=heightmap, CM=clique matrix. */
   SEXP AnsRects, AnsBounds=NULL, AnsHM=NULL, AnsCM=NULL, ListNames, Ans;
   double *pAnsRects;
   int *pAnsBounds=NULL, *pAnsHM=NULL, *pAnsCM=NULL;
   char *names[4] = {"rects","bounds","hm","cm"};
   SIntPoint point;

   int b, m;
   int i, j, k, numprotected=0;

   VerifyInputRectangles(RR, BB);

   if (ihm!=0 && ihm!=1)
      error("invalid third argument\n");
   if (icm!=0 && icm!=1)
      error("invalid fourth argument\n");

   if (ihm)
   {
      PROTECT(AnsHM = allocMatrix(INTSXP, 2*n+1, 2*n+1));
      numprotected++;
      pAnsHM = INTEGER(AnsHM);
   }
      
   /* 1. Transform to canonical rectangles */
   RealToCanonical(n, pRR, pBB, CanonObsRectangles, rx, ry, lb, LengthBB);

   /* 3. m := 0 */
   m = 0;

   /* 4. and 5. set e0..2n-1 to -1 (array h was set to zero at allocation) */
   for (i=0; i<2*n; i++)
      e[i] = -1;

   /* make column 2*n+1 of pAnsHM contain zeroes */
   if (ihm){
      pAnsHM[2*n*(2*n+1)] = 0;
      MemCopy(&pAnsHM[2*n*(2*n+1)+1],h,2*n);
   }

   /* 6. for j = 0 to 2n-1 do */
   for (j=0; j<2*n; j++)
   {
      if (ihm){
         pAnsHM[j*(2*n+1)] = 0;  /* make column 0 of pAnsHM contain zeroes */
         MemCopy(&pAnsHM[j*(2*n+1)+1],h,2*n);
      }

      /* 7. if rx_j is a left boundary */
      if (lb[j])
      {
         /* 8. for k = y_{1,rx_j} to (y_{2,rx_j} - 1) do */
         for (k=CanonObsRectangles[rx[j]].y1; k<CanonObsRectangles[rx[j]].y2; k++)
         {
            /* 9. h_k := h_{k} + 1; e_k := rx_j */
            h[k] += 1;
            e[k] = rx[j];  /* index of last entered rectangle */
         }
      }
      /* 10. else, i.e., if rx_j is a right boundary */
      else 
      {
         /* 11. b := y_{1,rx_j}: bottom coordinate of maximal intersection, 
          * -1 blocks output 
          */
         b = CanonObsRectangles[rx[j]].y1;

         /* look for local maxima in rows y_{1,rx_j} to (y_{2,rx_j} - 2)
          * 12. for k = ( y_{1,rx_j} ) to ( y_{2,rx_j} - 2 ) do 
          */
         for (k=CanonObsRectangles[rx[j]].y1; k<CanonObsRectangles[rx[j]].y2-1; k++)
         {
            /* 13. if ( (h_{k+1}<h_k) and (b>-1) ) then */
            if ( (h[k+1]<h[k]) && (b>-1) )
            {
               /* 14. if ( e_i>-1 for i=b,..,k); */
               for (i=b; (i<=k && e[i]>-1); i++)
			   {}  /* empty loop to advance i */
               if (i>k)
               {
                  if (m == MaxNrMIs) 
                  {
                     MaxNrMIs = 2*MaxNrMIs;
                     CanonMaxIntersections = R_Realloc(CanonMaxIntersections, MaxNrMIs, SCanonRect);   
                  }

                  /* 15. A_m := (x_{1,e_k}, j, b, k+1); m := m+1 */
                  CanonMaxIntersections[m].x1 = CanonObsRectangles[e[k]].x1;
                  CanonMaxIntersections[m].x2 = j;
                  CanonMaxIntersections[m].y1 = b;
                  CanonMaxIntersections[m].y2 = k+1;
                  m += 1;
                     
                  /* 16. e_{b} := -1 */
                  e[b] = -1;
               }  
               /* 17. b := -1 */
               b = -1;
            }
            /* 18. if (h_{k+1}>h_k) then */
            if (h[k+1]>h[k])
            {
               /* 19. b := k+1 */
               b = k+1;
            }
         }

         /* look for local maximum in row y_{2,rx_j} -1
          * 20. k := y_{2,rx_j} - 1 
          */
         k = CanonObsRectangles[rx[j]].y2 - 1;

         /* 21. if (b>-1) */ 
         if (b>-1)
         {
            /* 22. if (e_i>-1 for i=b,...,k) */
            for (i=b; (i<=k && e[i]>-1); i++)
			{} /* empty loop to advance i */
            if (i>k)
            {
               if (m == MaxNrMIs) 
               {
                  MaxNrMIs = 2*MaxNrMIs;
                  CanonMaxIntersections = R_Realloc(CanonMaxIntersections, MaxNrMIs, SCanonRect); 
               }

               /* 23. m := m + 1;  Am := (x_{1,e_k}, j, b, k+1) */
               CanonMaxIntersections[m].x1 = CanonObsRectangles[e[k]].x1;
               CanonMaxIntersections[m].x2 = j;
               CanonMaxIntersections[m].y1 = b;
               CanonMaxIntersections[m].y2 = k+1;         
               m += 1;
                  
               /* 24. e_{b} := -1 */
               e[b] = -1;
            }   
         }

         /* 25. for k = (y_{1,rx_j}) to (y_{2,rx_j} - 1) do */
         for (k=CanonObsRectangles[rx[j]].y1; k<CanonObsRectangles[rx[j]].y2; k++)
         {
            /* 26. hk := hk - 1 */
            h[k]--;
         }
      }
   }

   if (icm)
   {
      PROTECT(AnsCM = allocMatrix(INTSXP, m, n));
      numprotected++;
      pAnsCM = INTEGER(AnsCM);
      for (i=0; i<n; i++)
      {
         for (j=0; j<m; j++)
         {
            point.x = CanonMaxIntersections[j].x2;
            point.y = CanonMaxIntersections[j].y2;
            pAnsCM[j + i*m] = IsInRectangle(&point, &CanonObsRectangles[i]);
         }
      }
   }

   /* allocate space for output */
   PROTECT(AnsRects = allocMatrix(REALSXP, m, 4));   
   numprotected++;
   pAnsRects = REAL(AnsRects);

   if (BBexplicit)
   {
      PROTECT(AnsBounds = allocMatrix(INTSXP, m, 4));  
      numprotected++;
      pAnsBounds = INTEGER(AnsBounds);
   }

   /* 27. transform the canonical maximal intersections A_1..A_m 
          back to the original coordinates */
   CanonicalToReal(CanonMaxIntersections, m, m, NULL, n, pRR, pBB, BBexplicit, rx, ry, 
                   pAnsRects, pAnsBounds);
 
   /* Output list containing the transformed rectangles and their boundaries */
   PROTECT(ListNames = allocVector(STRSXP, 4));  
   numprotected++;
   for(i=0; i<4; i++)
      SET_STRING_ELT(ListNames,i,mkChar(names[i]));
 
   PROTECT(Ans = allocVector(VECSXP, 4));              
   numprotected++;
   SET_VECTOR_ELT(Ans,0,AnsRects);            

   if(BBexplicit)
      SET_VECTOR_ELT(Ans,1,AnsBounds);    
   else
      SET_VECTOR_ELT(Ans,1,BB);

   if(ihm)
      SET_VECTOR_ELT(Ans,2,AnsHM);
   
   if (icm)
      SET_VECTOR_ELT(Ans,3,AnsCM);

   setAttrib(Ans,R_NamesSymbol,ListNames);   /* attach vector names */

   UNPROTECT(numprotected);

   /* release memory of all intermediate data structures */
   R_Free(h);
   R_Free(e);
   R_Free(rx);
   R_Free(ry);
   R_Free(lb);
   R_Free(CanonMaxIntersections);
   R_Free(CanonObsRectangles);

   return Ans;
}

/* vim:set et ts=3 sw=3: */

