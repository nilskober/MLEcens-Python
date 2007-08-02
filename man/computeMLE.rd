\name{computeMLE}
\alias{computeMLE}
\title{Compute the MLE for bivariate censored data}
\description{
      This function computes the MLE for bivariate censored data. To be 
      more precise, we compute the MLE for the bivariate 
      distribution of (X,Y) in the following situation: realizations of (X,Y) 
      cannot be observed directly; instead, we observe a set of rectangles 
      (that we call 'observation rectangles') that are known to contain the 
      unobservable realizations of (X,Y).
}
\usage{computeMLE(R, B=c(0,1), max.inner=10, max.outer=1000, tol=1e-10)}
\arguments{
   \item{R}{
       A nx4 matrix of observation rectangles. Each row corresponds to
       a rectangle, represented as (x1,x2,y1,y2). The point (x1,y1) is the 
       lower left corner of the rectangle and the point (x2,y2) is the upper 
       right corner of the rectangle. 
   }
   \item{B}{
       This describes the boundaries of the rectangles (0=open or 1=closed). 
       It can be specified in three ways:
       \item A nx4 matrix containing 0's and 1's. Each row corresponds to a
       rectangle and is denoted as (cx1, cx2, cy1, cy2), where cx1 denotes the 
       boundary type of x1, cx2 denotes the boundary type of x2, etc. 
       \item A vector (cx1, cx2, cy1, cy2) containing 0's and 1's. This 
       representation can be used if all rectangles have the same type of 
       boundaries.
       \item A vector (c1, c2) containing 0's and 1's. This representation can 
       be used if all x and y intervals have the same type of boundaries. 
       c1 denotes the boundary type of x1 and y1, and c2 denotes the boundary 
       type of x2 and y2.

      The default value is c(0,1).
   }
   \item{max.inner}{Maximum number of iterations for the 
      inner iteration loop (see section 'Details').}
   \item{max.outer}{Maximum number of iterations for the 
      outer loop (see section 'Details').}
   \item{tol}{Tolerance up to which the necessary and sufficient conditions
      of the MLE must be satisfied (see section 'Details'). }
}
\details{
   Let PF(Ri) be the probability mass in observation rectangle Ri 
   under distribution F. Then the MLE maximizes 
   log PF(R1)+...+log(PF(Rn)) over the space of all bivariate distribution 
   functions F. The MLE can only assign mass to a finite number of disjoint
   sets, called maximal intersections, and the MLE is indifferent to the 
   distribution of mass within these sets. Hence, the computation of the MLE 
   can be split into two steps:
   a reduction step and an optimization step. In the reduction step, the 
   maximal intersections are computed. Next, in the optimization step, it is 
   determined how much probability mass should be assigned to each of these 
   areas. 

   The function \code{\link{computeMLE}} uses the height map algorithm
   of Maathuis (2005) for the reduction step (see \code{\link{reduc}})
   For the optimization step, it uses a 
   combination of sequential quadratic programming (outer iteration loop) and 
   the support reduction algorithm of Groeneboom, Jongbloed and Wellner (2007) 
   (inner iteration loop). It terminates when the necessary and sufficient 
   conditions for the MLE (see Maathuis (2003, page 49, eq (5.4)))  
   are satisfied up to a tolerance \code{tol}, or when 
   the maximum number of iterations is reached, whichever comes first. 
   We have found that it works well in practice to set \code{max.inner} low;
   hence it's default value is 10. 

   The MLE is typically non-unique, in two ways: (1) 
  The MLE is indifferent to the distribution of mass within the maximal
   intersections (called representational non-uniqueness by Gentleman and 
   Vandal (2002); (2) 
The amounts of mass assigned to the maximal intersections may be 
   non-unique (called mixture non-uniqueness by Gentleman and Vandal (2002)). 

   Hence, the algorithm \code{\link{computeMLE}} returns \emph{a} MLE. One 
   can deal with representational 
   non-uniqueness by creating an upper bound (assign all mass to the lower
   left corners of the maximal intersections) and a lower bound (assign
   all mass to the upper right corners of the maximal intersections) of the MLE
   (see also \code{\link{plotCDF1}}, \code{\link{plotCDF2}}, 
   \code{\link{plotDens1}}, \code{\link{plotDens2}}). The  
   algorithm does not (yet) allow to account for mixture non-uniqueness. 
}
\value{A list containing the following elements:
   \item{p}{Vector of length m with positive probability masses.}
   \item{rects}{A mx4 matrix containing the maximal intersections that
       correspond to the positive probability masses \code{p}.}
   \item{bounds}{Boundaries of \code{rects}, 
specified in the same way as \code{B}.}
   \item{conv}{Boolean indicating if the algorithm has converged.} 
   \item{llh}{Value of the log likelihood log(P(R1))+...+log(P(Rn)).}
}
\references{
   Groeneboom, Jongbloed and Wellner (2007). The support reduction 
   algorithm for computing nonparametric function estimates in mixture
   models. Submitted.

   Gentleman and Vandal (2002). Nonparametric estimation of the bivariate CDF 
   for arbitrarily censored data. \emph{Canadian Journal of Statistics}
   \bold{30} 557-571.

   M.H. Maathuis (2003). Nonparametric maximum likelihood estimation 
   for bivariate censored data. Master's thesis, Delft University of 
   Technology, The Netherlands. Available at 
   http://www.stat.washington.edu/marloes/papers.

   M.H. Maathuis (2005). Reduction algorithm for the NPMLE for
   the distribution function of bivariate interval censored data.
   \emph{Journal of Computational and Graphical Statistics} \bold{14}
   252--262.
}
\author{Marloes Maathuis: \email{marloes@u.washington.edu}. Part of the 
 code for the optimization step is adapted from code that was written by Piet 
 Groeneboom.}
\seealso{\code{\link{reduc}}, \code{\link{plotCDF1}}, \code{\link{plotCDF2}},
         \code{\link{plotDens1}}, \code{\link{plotDens2}}}
\examples{
# Load example data:
data(ex)
mle <- computeMLE(ex)
par(mfrow=c(2,2))

# Bivariate density plots of the MLE:

#   The colors represent the density=p/(area of maximal intersection)
plotDens2(mle, xlim=range(ex[,1:2]), ylim=range(ex[,3:4]), 
 main="Bivariate density plot of the MLE")
plotRects(ex, add=TRUE)

#   Alternative: numbers represent the mass p in the maximal intersections
plotDens2(mle, xlim=range(ex[,1:2]), ylim=range(ex[,3:4]), 
 col="lightgray", main="Bivariate density plot of the MLE", 
 key=FALSE, numbers=TRUE)
plotRects(ex, add=TRUE)

# Univariate density plots of the MLE:

#   Plot univariate density for X
plotDens1(mle, margin=1, xlim=range(ex[,1:2]), 
 main="Marginal density plot, 
 x-margin", xlab="x", ylab=expression(f[X](x)))

#   Plot univariate density for Y
plotDens1(mle, margin=2, xlim=range(ex[,3:4]), 
 main="Marginal density plot, 
 y-margin", xlab="y", ylab=expression(f[Y](y))) 

# Bivariate CDF plots of the MLE:

#   Plot lower bound for representational non-uniqueness
plotCDF2(mle, xlim=c(min(ex[,1])-1,max(ex[,2])+1), 
 ylim=c(min(ex[,3])-1, max(ex[,4])+1), bound="l", n.key=4,
 main="Bivariate CDF plot of the MLE,
 lower bound")

#   Add observation rectangles and shaded maximal intersections
plotRects(ex, add=TRUE) 
plotRects(mle$rects, density=20, border=NA, add=TRUE) 

#   Plot upper bound for representational non-uniqueness
plotCDF2(mle, xlim=c(min(ex[,1])-1,max(ex[,2])+1), 
 ylim=c(min(ex[,3])-1, max(ex[,4])+1), bound="u", n.key=4,
 main="Bivariate CDF plot of the MLE,
 upper bound")

#   Add observation rectangles and shaded maximal intersections
plotRects(ex, add=TRUE)
plotRects(mle$rects, density=20, border=NA, add=TRUE)

# Marginal CDF plots of the MLE:

#   Plot marginal CDF for X
plotCDF1(mle, margin=1, xlim=c(min(ex[,1])-1,max(ex[,2])+1), 
 bound="b", xlab="x", ylab="P(X<=x)", main="MLE for P(X<=x)")

#   Plot marginal CDF for Y
plotCDF1(mle, margin=2, xlim=c(min(ex[,3])-1,max(ex[,4])+1), 
 bound="b", xlab="y", ylab="P(Y<=y)", main="MLE for P(Y<=y)")
}

\keyword{survival}
\keyword{nonparametric}
\keyword{optimize}
\keyword{iteration}
\concept{nonparametric maximum likelihood estimator}
\concept{censored data}

