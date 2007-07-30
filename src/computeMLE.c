/* computeMLE.c */

/* All arrays start at zero */

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>

#include <R.h>
#include <Rinternals.h>

#include "MLEcens.h"

int FenchelViol(int m, int *ind, int mm, double nabla[], double tol,
                int *indexmin, double *minimum, double *maxnorm);
void	ComputeProbabilities(int ndata, SCanonRect *R, int m, SIntPoint *t, 
                           double alpha[], double P[]);
void	ComputeNabla(int ndata, SCanonRect *R, double P[], int m, SIntPoint *t, 
                   double tol, double nabla[]);
double Phi(int ndata, double P[], int m, double alpha[], double tol);
int ArmijoViol1(int m, double eps, double phi_old, double phi_new, 
                double nabla_old[], double alpha_old[], double alpha_new[]);
int ArmijoViol2(int m, double eps, double phi_old, double phi_new, double nabla_old[], 
                double alpha_old[], double alpha_new[]);
void 	Sort_t_alpha(int n, SIntPoint *t, double alpha[], int ind[], int *sorted,
                   int *i_dummy, double *d_dummy);
void	IterationStepIQM(int ndata, SCanonRect *R, double P[], double *w, 
                       int ndata2, SIntPoint *tt, int *m1, SIntPoint *t, 
                       double alpha[], double alpha_new[], int ind[], int *i_dummy_mm,
                       int *i_dummy_2mm, double *d_dummy_mm, int *indexmin1, 
                       int *iteration_inner, double *minimum1);
void	CreateValidMasses(int ndata, SCanonRect *R, double *w, int *m1, SIntPoint *t, 
                        double alpha[], double alpha_new[], int ind[], int* i_dummy_mm,
                        double *d_dummy_mm);
void ComputeWeightsIQM(int ndata, double P[], double tol, double *w);
double ComputeMinimumIQM(int ndata, SCanonRect *R, double P[], double *w, int ndata2, 
                         SIntPoint *tt, int m, int *index); 
void ComputeW(int ndata, SCanonRect *R, double *w, int m, SIntPoint *t, 
              double *W);
void ComputeAlphasIQM(int ndata, SCanonRect *R, double *w, int m, SIntPoint *t, 
                      double alpha[], int* i_dummy_mm, double *d_dummy_mm);
void HeightMapAlgorithmCanonical(int n, SCanonRect *CanonObsRectangles, int *rx, 
                                 int *lb, SCanonRect **CanonMaxIntersections, int *nrMaximalIntersections);
void MLE_IQM(int ndata, SCanonRect *R, int mm, SIntPoint *tt, 
             int *m1, SIntPoint *t, int *ind, double alpha_m[], 
             double tol, int MaxNumIterationsOuter, int MaxNumIterationsInner,
             double *criterion_function, double *sum_alpha, int *converged); 

             /* Function computes whether Fenchel conditions for outer loop are satisfied
             * (up to a tolerance 'tol'), see Maathuis (2003, page 49, eq (5.4). 
             * If these conditions are satisfied, we have found the MLE
             * 
             * Input: m = number of maximal intersections with positive mass
             * ind = indices of maximal intersections with positive mass
             * t[i] = tt[ind[i]], i=0,..m-1 are all maximal intersections with positive mass.
             * mm = number of maximal intersections
             * nabla = vector of derivatives of the object function 
             * (see Maathuis (2003, page 49, eq (5.5))
             * tol = specified tolerance with which Fenchel conditions must be satisfied.
             * 
             * Output: indexmin = index where nabla[i] is most negative (in 0,..,mm-1).
             * minimum = min(nabla[i]), i=0,..,mm-1.
             * maxnorm = max(abs(nabla[ind[i]]), i=0,..m.
*/
int FenchelViol(int m, int *ind, int mm, double nabla[], double tol,
                int *indexmin, double *minimum, double *maxnorm)
{
   double	min, max;
   int	 i, imin;
   
   min = 0; /* minimum of nabla[i], i=1,..,mm */
   imin = 0; /* argmin of nabla[i], i=1,.,,mm */
   max = 0; /* max of absolute value of nabla[i] for i st t[i] has positive mass */
   
   for (i=0; i<mm; i++)
   {
      if (nabla[i]<min)
      {
         min = nabla[i];
         imin = i;
      }
   }
   
   if (m>0)
   {
      for (i=0; i<m; i++)
      {
         if (fabs(nabla[ind[i]])>max)
            max = fabs(nabla[ind[i]]);
      }
   }
   
   *minimum = min;
   *indexmin = imin;
   *maxnorm = max;
   
   return (min<-tol || max>tol);
}


/* Computes P_F(R_i), probability masses in the observation rectangles
*
* Input: ndata = number of observation rectangles
*	 R	 = observation rectangles (in canonical form)
* m = number of maximal intersections 
* t = upper right corners of maximal intersections with positive mass 
* (in canonical form)
* alpha = alpha[i] is probability mass corresponding to t[i]
* 
* Output: P = P[i] is probability mass in observation rectangle R[i]
*/
void	ComputeProbabilities(int ndata, SCanonRect *R, int m, SIntPoint *t, 
                           double alpha[], double P[])
{
   int i, k;
   
   for (i=0; i<ndata; i++)
   {
      P[i] = 0;
      for (k=0; k<m; k++)
      {
         if (IsInRectangle(&t[k],&R[i]))
            P[i] += alpha[k]; 
      }
   }
}


/* Computes nabla: vector of partial derivatives of object function w.r.t. alpha[j],
* see Maathuis (2003, page 49, eq (5.5))
*
* Input: ndata = number of observation rectangles
* R = observation rectangles (in canonical form)
* P = P[i] is probability mass in observation rectangle R[i]
* m = number of maximal intersections with positive mass
* t = upper right corners of maximal intersections with positive mass
* (in canonical form)
* 
* Output: nabla = vector of partial derivatives of object function w.r.t. 
alpha[j], j=0..m-1, see Maathuis (2003, page 49, eq (5.5)).
*/
void	ComputeNabla(int ndata, SCanonRect *R, double P[], int m, SIntPoint *t, 
                   double tol, double nabla[])
{
   int j, k;
   double sum;
   
   for (k=0; k<m; k++)
   {
      sum = 0;
      for (j=0; j<ndata; j++)
      {
         if (IsInRectangle(&t[k],&R[j]))
         {
            if (P[j]>tol)
               sum += 1/P[j];
            else
               sum += 1/tol;
         }
      }
      nabla[k] = 1.0 - sum/ndata;
   }
}


/* Computes phi, the function that must be minimized: 
* - (1/n) sum_{i=1}^n log P_F(R_i) + sum alpha_j - 1,
* see Maathuis (2003, page 49, eq (5.2))
*
* Input: ndata = number of observation rectangles 
* P	 = P[i] is probability mass in observation rectangle R[i]
* m = number of maximal intersections with positive mass
* alpha = alpha[i] is probability mass corresponding to maximal intersection t[i]
* 
* Output: the value of phi under this iterate.
*/
double Phi(int ndata, double P[], int m, double alpha[], double tol)
{
   int		i;
   double	sum1=0, sum2=0;
   
   for (i=0; i<ndata; i++)
   {
      if (P[i]>tol)
         sum1 -= log(P[i]);
      else
         sum1 -= log(tol);
   }
   
   for (i=0; i<m; i++)
      sum2 += alpha[i];
			
   return sum1/ndata + sum2 - 1;
}


int ArmijoViol1(int m, double eps, double phi_old, double phi_new, 
                double nabla_old[], double alpha_old[], double alpha_new[])
{
   int	 i;
   double inprod=0;
   
   for (i=0; i<m; i++)
      inprod += nabla_old[i]*(alpha_new[i] - alpha_old[i]);
   
   return (phi_new < phi_old + (1-eps)*inprod);
}


int ArmijoViol2(int m, double eps, double phi_old, double phi_new, double nabla_old[], 
                double alpha_old[], double alpha_new[])
{
   int	 i;
   double inprod=0;
   
   for (i=0; i<m; i++)
      inprod += nabla_old[i]*(alpha_new[i]-alpha_old[i]);
   
   return (phi_new > phi_old + eps*inprod);
}


/* Helper function for Sort_t_alpha() */
static int SortInt(const void *elem1, const void *elem2)
{
   int *a = (int*) elem1;
   int *b = (int*) elem2;
   
   if (a[0] < b[0])
      return -1;
   
   return (a[0]>b[0]);
}


/* Sort array ind[] so that it is increasing,
* and order arrays t[] and alpha[] accordingly */
void Sort_t_alpha(int m, SIntPoint *t, double alpha[], int ind[], int *sorted, 
                  int *i_dummy, double *d_dummy)
{
   int j;
   
   /* copy ind[0]..ind[m-1], the values to be sorted, in sorted[2*j], j=0..m-1
   * copy indices 0..m-1 in sorted[2*j+1], j=0..m-1 */
   for (j=0; j<m; j++)
   {
      sorted[2*j] = ind[j];
      sorted[2*j+1] = j;
   }
   
   /* sort array with qsort algorithm */
   qsort(sorted, m, 2*sizeof(int), SortInt);
   
   /* use indices in sorted[2*j+1],j=0..m-1 to sort arrays t[][0],t[][1] and alpha[] */
   for (j=0; j<m; j++) i_dummy[j] = t[j].x;
   for (j=0; j<m; j++) t[j].x = i_dummy[sorted[2*j+1]];
   for (j=0; j<m; j++) i_dummy[j] = t[j].y;
   for (j=0; j<m; j++) t[j].y = i_dummy[sorted[2*j+1]];
   for (j=0; j<m; j++) d_dummy[j] = alpha[j];
   for (j=0; j<m; j++) alpha[j] = d_dummy[sorted[2*j+1]];
   
   /* copy sorted values in ind[j],j=0..m-1 */
   for (j=0; j<m; j++) 
      ind[j] = sorted[2*j];
}


/* Iteration step for the LS problem 
* 
* Input: ndata = number of observation rectangles
* R = observation rectangles (in canonical form)
* w = w[i]=1/max(tol, P[i]), where P[i] is probability mass in 
* observation rectangle R[i].
* mm = number of maximal intersections. 
* tt = upper right corners of maximal intersections (in canonical form)
* m = number of maximal intersections with positive mass.
*		 t = upper right corners of maximal intersections with positive mass 
* (in canonical form).
*		 alpha = alpha[i] is probability mass corresponding to t[i]
* alpha_new = dummy vector for candidate vector alpha
* (this vector is allocated in MLE_IQM() to limit number of memory allocations).
* ind = idices of maximal intersections with positive mass
* (t[i]=tt[ind[i]], i=0,..m-1, are all maximal intersections with positive mass.
* i_dummy_mm, i_dummy_2mm, d_dummy_mm = dummy vectors, allocated in MLE_IQM()
* to limit number of memory allocations.
* indexmin1 = index of minimum of derivatives of criterion function for LS part
* (see Maathuis (2003, page 51, eq (5.10)).
* iteration_inner = number of performed iterations in the current inner loop. 
*
* Output: P	 = updated mass in observation rectangles.
* m = updated number of maximal intersections with positive mass. 
* t = updated maximal intersections with positive mass.
*		 ind = updated indices of maximal intersections with positive mass.
*		 alpha = alpha[i] is updated probability mass corresponding to t[i].
* indexmin1 = index of updated minimum of derivatives of criterion function for LS part
* (see Maathuis (2003, page 51, eq (5.10)).
* iteration_inner is increased by one. 
* minimum1 = updated minimum of derivatives of criterion function for LS part 
* (see Maathuis (2003, page 51, eq (5.10)).
*/
void	IterationStepIQM (int ndata, SCanonRect *R, double P[], double *w, 
                        int mm, SIntPoint *tt, int *m1, SIntPoint *t, 
                        double alpha[], double alpha_new[], int ind[], 
                        int *i_dummy_mm, int *i_dummy_2mm, double *d_dummy_mm, 
                        int *indexmin1, int *iteration_inner, double *minimum1)
{
   int m, indexmin;
   
   m = *m1;
   indexmin = *indexmin1;
   
   /* add masspoint with index indexmin */
   ind[m] = indexmin;
   t[m].x = tt[indexmin].x;
   t[m].y = tt[indexmin].y;
   alpha[m] = 0;
   m++;
   
   /* Sort array ind[] so that it is increasing,
   * and order arrays t[] and alpha[] accordingly */
   if (m>1)
      Sort_t_alpha(m,t,alpha,ind,i_dummy_2mm,i_dummy_mm,d_dummy_mm);
   
   /* compute alpha_new, based on vectors R, w and t */
   ComputeAlphasIQM(ndata,R,w,m,t,alpha_new,i_dummy_mm,d_dummy_mm);
   
   if (m==1)
      alpha[0] = alpha_new[0];
   else	
      CreateValidMasses(ndata,R,w,&m,t,alpha,alpha_new,ind,i_dummy_mm,d_dummy_mm);
   
   ComputeProbabilities(ndata,R,m,t,alpha,P);
   *minimum1 = ComputeMinimumIQM(ndata,R,P,w,mm,tt,m,&indexmin);
   
   *m1 = m;
   *indexmin1 = indexmin;
   (*iteration_inner)++;
}

/* Create valid mass vector from old vector alpha and candidate vector alpha_new
* We check whether all alpha_new's are nonnegative. If this is not the case 
* we remove mass the point for which alpha[i]/(alpha[i]-alpha_new[i]) 
* is smallest. We continue until there is no such violation. 
* In the end, the new values are copied into alpha, and the new length is m1. 
* The vectors t and ind are adjusted as well.
*
* Input: ndata = number of observation rectangles
* R = observation rectangles (in canonical form)
* w = w[i]=1/max(tol,P[i]), where P[i] is probability mass in observation rectangle R[i]
* m1 = number of maximal intersections with positive mass. 
* t = maximal intersections with positive mass. 
* alpha = alpha[i] is mass corresponding to t[i], i=0,..m1-1.
* alpha_new = candidate value for mass corresponding to t[i], i=0,..,m1-1.
* ind = indices of maximal intersections with positive mass. 
* t[i]=tt[ind[i]], i=0,..,m1-1 are all maximal intersections with positive mass. 
* i_dummy_mm, d_dummy_mm = dummy vectors needed in function ComputeAlphasIQM
* (they are allocated in main() to limit number of memory allocations).
*
* Output: alpha = updated mass vector.
* m1 = updated number of maximal intersections with positive mass.
* t = updated upper right corners of maximal intersections with positive mass.
* ind = updated indices of maximal intersections with positive mass.
*/
void	CreateValidMasses (int ndata, SCanonRect *R, double *w, int *m1, SIntPoint *t, 
                         double alpha[], double alpha_new[], int ind[], 
                         int *i_dummy_mm, double *d_dummy_mm)
{
   double min, min1, crit, crit0;
   int i, k;
   int m, imin;
   
   m = *m1;
   
   crit0 = 0; 
   min = 0; /* minimum of alpha_new */
   min1 = 1.0; /* minimum of alpha[i]/(alpha[i]-alpha_new[i]) */
   imin = 0;
   
   for (i=0; i<m; i++)
   {
      if (alpha_new[i]<0)
      {
         crit = alpha[i]/(alpha[i]-alpha_new[i]);
         if (crit<min1)
         {
            min1 = crit;
            imin = i;
         }
         
         if (alpha_new[i]<min)
            min = alpha_new[i];
      }
   }
   
   while (min<0)
   {
   /* create new alpha1 from linear combination of alpha and alpha_new
      * note that, by definition, the new alpha_new[imin]=0 */
      for (k=0; k<m; k++)
         alpha_new[k] = alpha[k] + min1*(alpha_new[k]-alpha[k]);
      
         /* we now remove index imin
         * if imin==m-1, we can just do m--, and we are done
         * if imin<m-1, then we first need to move the entries that are larger 
      * than imin each one place to the left */
      if (imin < m-1)
      {
         for (k=imin; k<m-1; k++)
         {
            t[k].x = t[k+1].x;
            t[k].y = t[k+1].y;
            alpha_new[k] = alpha_new[k+1];
            ind[k] = ind[k+1];
         }
      }
      m--;
      
      MemCopy(alpha,alpha_new,m);
      ComputeAlphasIQM(ndata,R,w,m,t,alpha_new,i_dummy_mm,d_dummy_mm);
      
      min1 = 1.0;
      min = 0;
      
      for (i=0; i<m; i++)
      {
         if (alpha_new[i] < 0)
         {
            crit = alpha[i]/(alpha[i]-alpha_new[i]);
            if (crit < min1)
            {
               min1 = crit;
               imin = i;
            }
            if (alpha_new[i] < min)
               min = alpha_new[i];
         }
      }
   }
   MemCopy(alpha,alpha_new,m);
   *m1=m;
}

/* Compute weights w[i] = 1/max(tol, P[i]), where P[i] is probability mass in 
* observation rectangle R[i].
* Input: ndata = number of observation rectangles. 
* P = P[i] is probability mass in observation rectangle R[i].
* tol = tolerance to prevent dividing by zero.
* 
* Output: w = w[i] = 1/max(tol,P[i]).
*/
void ComputeWeightsIQM (int ndata, double P[], double tol, double *w)
{
   int	i;
   
   for (i=0; i<ndata; i++)
   {
      if (P[i]>tol)
         w[i] = 1.0/P[i];
      else
         w[i] = 1.0/tol;
   }
}


/* The following procedure computes the minimum of the derivatives of the 
* criterion function for the LS part (see Maathuis (2003, page 51, eq (5.10)).
* At the end of the iteration loop, this minimum should be 
* nonnegative or bigger than some allowed negative tolerance. 
*
* Input: ndata = number of observation rectangles
* R = observation rectangles (in canonical form)
* 	 P = probability mass in each observation rectangle (not used if m==0)
*		 w = w[i] = 1/max(tol,P[i])
*		 mm = number of maximal intersections
*		 tt = upper right corners of maximal intersections (in canonical form)
*		 m = number of maximal intersections that get positive mass
*
* Output: index = index (in 0,..mm-1) at which minimum occurs.
* return value = minimum of the derivatives of the criterion function 
*						 for the least squares part.
*/
double ComputeMinimumIQM (int ndata, SCanonRect *R, double P[], double *w, int mm, 
                          SIntPoint *tt, int m, int *indexmin)
{
   int j, k;
   double sum, min=0;
   
   if (m==0)
   {
      for (k=0; k<mm; k++)
      {
         sum = 0;
         for (j=0; j<ndata; j++)
         {
            if (IsInRectangle(&tt[k],&R[j]))
               sum -= w[j];
         }
         sum *= 2;
         
         if (sum<min)
         {
            min = sum;
            *indexmin = k;
         }
      }
      return min;
   }
   
   /* else, if m>0: */
   for (k=0; k<mm; k++)
   {
      sum = 0;
      for (j=0; j<ndata; j++)
      {
         if (IsInRectangle(&tt[k],&R[j]))
         {
            sum += P[j]*SQR(w[j]) - 2*w[j];
         }
      }
      sum = sum/ndata + 1;
      
      if (sum<min)
      {
         min = sum;
         *indexmin = k;
      }
   }
   return min;
}


/* Computation of the matrix W that is used in the function ComputeAlphasIQM().
*
* Input: ndata = number of observation rectangles.
* R = observation rectangles (in canonical form).
* w = w[i]=1/max(P[i],tol), where P[i] is the mass in 
* observation rectangle R[i].
*		 m = number of maximal intersections with positive mass.
*		 t = upper right corners of maximal intersections with positive mass
*		 (in canonical form).
* Output: W = matrix W used in the function ComputeAlphasIQM().
*				 W is a symmetric matrix. We only store its upper triangle 
*				 (including the diagonal), column by column. To translate between 
*				 the full matrix a[i][j], i=0..m-1, j=0..m-1 and the compactly stored 
*				 matrix W[i], i=0..m*(m+1)/2, use a[i][j] = W[i+j*(j+1)/2], for i<=j.
*/
void ComputeW(int ndata, SCanonRect *R, double *w, int m, SIntPoint *t, double *W)
{
   int i, j, k, lengthW;
   
   lengthW = m*(m+1)/2;
   
   for (i=0; i<lengthW; i++)
      W[i] = 0;
   
   for (i=0; i<m; i++)
   {
      for (j=0; j<ndata; j++)
      {
         if (IsInRectangle(&t[i],&R[j]))
         {
            for (k=i; k<m; k++)
            {
               if (IsInRectangle(&t[k],&R[j]))
                  W[i+k*(k+1)/2] += SQR(w[j]); /* add to element a[i][k] */
            }
         }
      }
   }
   
   for (i=0; i<lengthW; i++)
      W[i] = W[i]/ndata;
}


/* Computation of the coefficients alpha that satisfy the equality part of 
* Maathuis (2003, page 50, eq (5.9)) for all alpha's with positive mass. 
* The alpha's are a solution of a linear system of equations with 
* a symmetric (but not positive definite) matrix of coefficients. It is solved using 
* the function SolveSymmetricLinearSystem(),which uses a LAPACK routine. 
*
* Input: ndata = number of observation rectangles 
* R = observation rectangles (in canonical form)
* w = w[i]=1/max(tol,P[i]), where P[i] is mass in observation rectangle R[i]
* m = number of maximal intersections with positive mass
* t = upper right corners of maximal intersections with positive mass
* i_dummy_mm = dummy array needed in function SolveSymmetricLinearSystem
* (this array is allocated in MLE_IQM() to limit number of memory allocations)
* d_dummy_mm = right hand side of linear system.
* (this array is allocated in MLE_IQM() to limit number of memory allocations)
*
* Output: alpha = alpha[i] is new probability mass corresponding to t[i].
*/
void ComputeAlphasIQM (int ndata, SCanonRect *R, double *w, int m, SIntPoint *t, 
                       double alpha[], int *i_dummy_mm, double *d_dummy_mm)
{
   int i, j;
   double *W;
   
   W = Calloc(m*(m+1)/2, double);
   
   /* compute matrix a, based on vectors R, w and t */
   ComputeW(ndata,R,w,m,t,W);
   
   for (i=0; i<m; i++)
   {
      d_dummy_mm[i]=0;
      for (j=0; j<ndata; j++)
      {
         if (IsInRectangle(&t[i],&R[j]))
            d_dummy_mm[i] += w[j];
      }
   }
   
   for (i=0; i<m; i++)
      d_dummy_mm[i] = 2*d_dummy_mm[i]/ndata - 1;
   
   SolveSymmetricLinearSystem(W,m,d_dummy_mm,1,i_dummy_mm);
   
   MemCopy(alpha,d_dummy_mm,m);
   Free(W);
}

/* Algorithm to find maximal intersections: the areas of possible mass support.
* The code follows the pseudo code given in the JCGS paper, 
* except for the fact that all arrays start from zero here. 
* Input:
* n = number of observation rectangles
* CanonObsRects = canonical observation rectangles, stored as SCanonRect
* rx: vector of length 2*n containing the indices of the rectangles 
* corresponding to the sorted x-coordinates
* lb: vector of length 2*n indicating the left/right boundaries of the 
* sorted x-coordinates (1=left, 0=right)
* Output: 
* CanonMaxIntersections = canonical maximal intersections, stored as *ScanonRect
* m = number of maximal intersections
*/
void HeightMapAlgorithmCanonical(int n, SCanonRect *CanonObsRectangles, int *rx, 
                                 int *lb, SCanonRect **CanonMaxInt, int *nrMaximalIntersections)
{
   /* allocate space for canonical maximal intersections */
   int MaxNrMIs = n; /* starting size for nr of CanonMaxIntersections */
   
   int *h = Calloc(2*n, int); /* value of height map, initialized at zero */
   int *e = Calloc(2*n, int); /* value of last entered rectangle */
   
   int b, m;
   int i, j, k;
   
   SCanonRect *CanonMaxIntersections = *CanonMaxInt;
   
   CanonMaxIntersections = Calloc(MaxNrMIs, SCanonRect); /* starting size, will be increased if necessary */
   
   /* 3. m := 0 */
   m = 0;
   
   /* 4. and 5. set e0..2n-1 to -1 (array h was set to zero at allocation) */
   for (i=0; i<2*n; i++)
      e[i] = -1;
   
   /* 6. for j = 0 to 2n-1 do */
   for (j=0; j<2*n; j++)
   {
      /* 7. if rx_j is a left boundary */
      if (lb[j])
      {
         /* 8. for k = y_{1,rx_j} to (y_{2,rx_j} - 1) do */
         for (k=CanonObsRectangles[rx[j]].y1; k<CanonObsRectangles[rx[j]].y2; k++)
         {
            /* 9. h_k := h_{k} + 1; e_k := rx_j */
            h[k] += 1;
            e[k] = rx[j]; /* index of last entered rectangle */
         }
      }
      /* 10. else, i.e., if rx_j is a right boundary */
      else 
      {
      /* 11. b := y_{1,rx_j}: bottom coordinate of maximal intersection, 
      * -1 blocks output. 
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
                  ;		/* empty loop */

               if (i>k)
               {
                  if (m == MaxNrMIs)
                  {
                     MaxNrMIs = 2*MaxNrMIs;
                     CanonMaxIntersections = Realloc(CanonMaxIntersections, MaxNrMIs, SCanonRect);	
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
            for (i=b; (i<=k && e[i]>-1); i++);
            if (i>k)
            {
               if (m == MaxNrMIs) 
               {
                  MaxNrMIs = 2*MaxNrMIs;
                  CanonMaxIntersections = Realloc(CanonMaxIntersections, MaxNrMIs, SCanonRect);	
               }
               /* 23. m := m + 1; Am := (x_{1,e_k}, j, b, k+1) */
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
   
   /* 28. return A_1..A_m */
   *nrMaximalIntersections = m;
   *CanonMaxInt = CanonMaxIntersections; 
   
   /* release memory of all intermediate data structures */
   Free(h);
   Free(e);
}


void MLE_IQM(int ndata, SCanonRect *R, int mm, SIntPoint *tt, 
             int *m1, SIntPoint *t, int *ind, double alpha_m[], 
             double tol, int MaxNumIterationsOuter, int MaxNumIterationsInner,
             double *criterion_function, double *sum_alpha, int *converged) 
{
   int i, m, iteration_outer, iteration_inner, indexmin, armijo_viol1, armijo_viol2;
   int *i_dummy_mm, *i_dummy_2mm;
   
   double norm, phi_old, phi, partialsum, eps, l, s, minimum, sum;
   double *P, *P_old, *P_temp, *w, *nabla, *d_dummy_mm;
   double *alpha_m_temp, *alpha_old, *alpha, *alpha_temp;
   
   P = Calloc(ndata, double);
   P_old = Calloc(ndata, double);
   P_temp = Calloc(ndata, double);
   w = Calloc(ndata, double);
   
   nabla = Calloc(mm, double);
   alpha_m_temp = Calloc(mm, double);
   alpha_old = Calloc(mm, double);
   alpha = Calloc(mm, double);
   alpha_temp = Calloc(mm, double);
   
   /* dummy variables, allocated here so that they aren't allocated and freed inside a loop */
   i_dummy_mm = Calloc(mm, int); /* dummy vector of integers of length mm */
   i_dummy_2mm = Calloc(2*mm, int); /* dummy vector of integers of length 2*mm */
   d_dummy_mm = Calloc(mm, double); /* dummy vector of doubles of length mm */
   iteration_inner = 0;
   iteration_outer = 0;
   eps = 0.1;
   
   /* start with uniform mass in all mass points */
   for (i=0; i<mm; i++)
      alpha_old[i] = 1.0/mm; 
   ComputeProbabilities(ndata,R,mm,tt,alpha_old,P_old);
   ComputeWeightsIQM(ndata,P_old,tol,w);
   
   m = 0;
   
   /* Outer iteration loop. Solution is obtained when fenchel conditions are 
   * satisfied. */
   while (iteration_outer==0 || (iteration_outer<MaxNumIterationsOuter && 
      FenchelViol(m,ind,mm,nabla,tol,&indexmin,&partialsum,&norm)))
   {
      iteration_inner = 0;
      
      phi_old = Phi(ndata,P_old,mm,alpha_old,tol);
      
      /* print interation and phi_old for outside loop */
      /* printf("\nOuter loop iteration %d, phi=%lf, m=%d\n",iteration_outer,phi_old,m); */
      
      minimum = ComputeMinimumIQM(ndata,R,P,w,mm,tt,m,&indexmin);
      
      /* Inner iteration loop: solve LS problem obtained by quadratic approximation
      * around alpha_old. As long as minimum<-tol, solution has not been obtained. */
      while (minimum<-tol && iteration_inner<MaxNumIterationsInner)
      {
         IterationStepIQM(ndata,R,P,w,mm,tt,&m,t,alpha_m,alpha_m_temp,
            ind,i_dummy_mm,i_dummy_2mm,d_dummy_mm,
            &indexmin,&iteration_inner,&minimum);
         /* P, m, t, alpha_m, ind, indexmin, iteration and minimum are updated */
         /* printf("Inner loop %d, min=%lf\n",iteration_inner,minimum); */
      }
      
      phi = Phi(ndata,P,m,alpha_m,tol);
      
      /* put vector alpha_m of lenght m into vector alpha of length mm */
      for (i=0; i<mm; i++)
         alpha[i] = 0; 
      for (i=0; i<m; i++)
         alpha[ind[i]] = alpha_m[i]; 
      
      /* determine correct stepsize using Armijo's rule (see, e.g., Jongbloed 1998) */
      armijo_viol1 = ArmijoViol1(mm,eps,phi_old,phi,nabla,alpha_old,alpha);
      armijo_viol2 = ArmijoViol2(mm,eps,phi_old,phi,nabla,alpha_old,alpha);
      
      if (!(armijo_viol2)) 
      {
      /* upate alpha_old by taking a step into the direction of alpha, 
         * of size 0.1*(alpha-alpha_old) */
         Convex_Comb(mm,0.1,alpha_old,alpha,alpha_old); 
         
         /* update P and P_old accordingly */
         Convex_Comb(ndata,0.1,P_old,P,P_old); 
      }
      else															 
      {
         MemCopy(alpha_temp,alpha,mm); 
         MemCopy(P_temp,P,ndata);
         
         l = 1.0;
         s = 0.5;
         
         while (s>tol && (armijo_viol1 || armijo_viol2)) 
         {
            if (armijo_viol1) 
               l = l+s;
            else
               l = l-s;
            
            Convex_Comb(mm,(1-l),alpha_old,alpha_temp,alpha); 
            Convex_Comb(ndata,(1-l),P_old,P_temp,P); 
            
            phi = Phi(ndata,P,mm,alpha,tol); 
            armijo_viol1 = ArmijoViol1(mm,eps,phi_old,phi,nabla,alpha_old,alpha);
            armijo_viol2 = ArmijoViol2(mm,eps,phi_old,phi,nabla,alpha_old,alpha);
            
            s = s/2;
         }
         /* update alpha_old by taking a step into the direction of alpha, 
         * of size 0.1*(alpha-alpha_old) */
         Convex_Comb(mm,0.1,alpha_old,alpha,alpha_old); 
         
         /* update P and P_old accordingly */
         Convex_Comb(ndata,0.1,P_old,P,P_old); 
      }
      ComputeWeightsIQM(ndata,P_old,tol,w);
      ComputeAlphasIQM(ndata,R,w,m,t,alpha_m_temp,i_dummy_mm,d_dummy_mm); 
      CreateValidMasses(ndata,R,w,&m,t,alpha_m,alpha_m_temp,ind,i_dummy_mm,d_dummy_mm); 
      ComputeProbabilities(ndata,R,m,t,alpha_m,P); 
      ComputeNabla(ndata,R,P,mm,tt,tol,nabla);
      
      iteration_outer++;
   }
   
   sum=0;
   for (i=0; i<m; i++)
      sum += alpha_m[i];
   
   phi = Phi(ndata,P,m,alpha_m,tol);
   *criterion_function = phi;
   *m1 = m;
   *sum_alpha = sum;
   *converged = !FenchelViol(m,ind,mm,nabla,tol,&indexmin,&partialsum,&norm);
   
   Free(P);
   Free(P_old);
   Free(P_temp);
   Free(w);
   Free(nabla);
   Free(alpha_m_temp);
   Free(alpha_old);
   Free(alpha);
   Free(alpha_temp);
   Free(i_dummy_mm);
   Free(i_dummy_2mm);
   Free(d_dummy_mm);
}

SEXP ComputeMLEForR(SEXP RR, SEXP BB, SEXP MaxIterInner, SEXP MaxIterOuter, SEXP tol)
{
   /* Variables for input */
   double *pRR = REAL(RR);
   int n = nrows(RR);
   int *pBB = INTEGER(BB);
   int LengthBB = length(BB);
   int BBexplicit = (LengthBB==4*n);
   int iMaxIterInner = *INTEGER(MaxIterInner);
   int iMaxIterOuter = *INTEGER(MaxIterOuter);
   double dtol = *REAL(tol);
   
   /* Variables for output */
   SEXP Ans, AnsAlpha, AnsRects, AnsBounds=NULL, AnsConverge, AnsLogl, ListNames;
   double *pAnsAlpha, *pAnsRects, *pAnsLogl, phi, sum_alpha;
   int *pAnsBounds=NULL, *pAnsConverge, iconv;
   char *names[5] = {"p","rects","bounds","conv","llh"};
   
   /* Internal variables */
   SIntPoint *tt, *t;
   int *ind, i, mm, m, numprotected=0;
   double *alpha;
   SCanonRect *CanonMaxIntersections; 
   /* CanonMaxIntersections will be allocated in HeightMapAlgorithmCanonical */
   
   /* Allocate space */
   SCanonRect *CanonObsRectangles = Calloc(n, SCanonRect);
   int *rx = Calloc(2*n, int);
   int *ry = Calloc(2*n, int);
   int *lb = Calloc(2*n, int);
   
   VerifyInputRectangles(RR, BB);
   
   RealToCanonical(n,pRR,pBB,CanonObsRectangles,rx,ry,lb,LengthBB);
   
   /* Do HeightMap algorithm */
   HeightMapAlgorithmCanonical(n,CanonObsRectangles,rx,lb,&CanonMaxIntersections,&mm);
   
   /* Allocate space for input and output of optimization step */
   tt = Calloc(mm, SIntPoint);
   t = Calloc(mm, SIntPoint);
   ind = Calloc(mm, int); 
   alpha = Calloc(mm, double);
   
   /* Copy upper right corners of maximal intersections into the array tt. 
   * (It is sufficient to only consider the upper right corners for the 
   * optimization step, see Maathuis (2003, page 31, Lemma 4.1)).
   */
   for (i=0; i<mm; i++)
   {
      tt[i].x = CanonMaxIntersections[i].x2;
      tt[i].y = CanonMaxIntersections[i].y2;
   }
   
   /* Do optimization */
   MLE_IQM(n,CanonObsRectangles,mm,tt,&m,t,ind,alpha,dtol,iMaxIterOuter,
      iMaxIterInner,&phi,&sum_alpha,&iconv);
   
   if (!iconv)
      warning("the algorithm did not converge");
   if (sum_alpha>1+dtol || sum_alpha<1-dtol)
      warning("total probability mass is not equal to one");
   
   PROTECT(AnsLogl = allocVector(REALSXP,1));
   numprotected++;
   pAnsLogl = REAL(AnsLogl);
   
   PROTECT(AnsConverge = allocVector(INTSXP,1));
   numprotected++;
   pAnsConverge = INTEGER(AnsConverge);
   
   pAnsLogl[0] = -n*(phi-sum_alpha+1);
   pAnsConverge[0] = iconv;
   
   PROTECT(AnsAlpha = allocVector(REALSXP, m));
   numprotected++;
   pAnsAlpha = REAL(AnsAlpha);
   MemCopy(pAnsAlpha, alpha, m);
   
   PROTECT(AnsRects = allocMatrix(REALSXP, m, 4));
   numprotected++;
   pAnsRects = REAL(AnsRects);
   
   if (BBexplicit)
   {
      PROTECT(AnsBounds = allocMatrix(INTSXP, m, 4));
      numprotected++;
      pAnsBounds = INTEGER(AnsBounds);
   }
   
   /* Transform canonical rectangles back to original coordinates */
   CanonicalToReal(CanonMaxIntersections,mm,m,ind,n,pRR,pBB,BBexplicit,rx,ry,
      pAnsRects,pAnsBounds);
   
   /* Output list */
   PROTECT(ListNames = allocVector(STRSXP, 5)); 
   numprotected++;
   for(i=0; i<5; i++)
      SET_STRING_ELT(ListNames,i,mkChar(names[i]));
   
   PROTECT(Ans = allocVector(VECSXP, 5)); 
   numprotected++;
   
   SET_VECTOR_ELT(Ans,0,AnsAlpha);
   SET_VECTOR_ELT(Ans,1,AnsRects);
   
   if (BBexplicit)
      SET_VECTOR_ELT(Ans,2,AnsBounds);
   else
      SET_VECTOR_ELT(Ans,2,BB);
   
   SET_VECTOR_ELT(Ans,3,AnsConverge);
   SET_VECTOR_ELT(Ans,4,AnsLogl);
   
   setAttrib(Ans,R_NamesSymbol,ListNames); /* attach vector names */
   
   Free(CanonObsRectangles);
   Free(CanonMaxIntersections);
   Free(rx);
   Free(ry);
   Free(lb);
   Free(tt);
   Free(t);
   Free(ind);
   Free(alpha);
   
   UNPROTECT(numprotected);
   
   return Ans;
}

/* vim:set et ts=3 sw=3: */
