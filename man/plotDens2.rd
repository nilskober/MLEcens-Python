\name{plotDens2}
\alias{plotDens2}
\title{Create a bivariate density plot of the MLE}
\description{This function creates a bivariate density plot 
of the MLE. It plots the maximal intersections that
get positive mass under the MLE, filled with colors that represent the 
estimated density of F under the assumption that the mass is 
distributed uniformly over the maximal intersections. 
In other words, if the MLE assigns mass p
to a maximal intersection, then the color represents 
p/(area of maximal intersection). 
}
\usage{
plotDens2(mle, col=gray(seq(.9,.3,len=30)), border=NA,  
          key=TRUE, n.key=10, round.key=2, 
          numbers=FALSE, round.numbers=2, cex.numbers=0.6, 
          xlim=NULL, ylim=NULL, zlim=NULL, breaks=NULL, 
          xlab="", ylab="", main="", sub="")
}
\arguments{
   \item{mle}{List with elements 'p' and 'rects', as outputted by 
       \code{\link{computeMLE}}.}
   \item{col}{Color vector used to represent the values of the density, and
to fill the maximal intersections. The default value is 
\code{col=gray(seq(.9,.3,len=30))}.}
   \item{border}{Color for rectangle borders of the maximal intersections. 
The default value NA means that the borders are omitted.}
   \item{key}{Logical, indicating if a color key should be drawn. The 
              default value is TRUE.}
   \item{n.key}{Approximate number of tickmarks for the color key.
                The default value is 10.}
   \item{round.key}{Number of decimals used for the labels of the color key.
                    The default value is 2.}
   \item{numbers}{Logical, indicating if the total amounts of mass \code{mle$p}
                 should be printed
                 in the maximal intersections. The default value is FALSE.}
   \item{round.numbers}{Number of decimals used for \code{numbers}. The 
                        default value is 2.}
   \item{cex.numbers}{Numerical value giving the amount by which text 
  size of \code{numbers} should be scaled relative to the default. The default 
value is 0.6.}
   \item{xlim, ylim}{Ranges for the plotted x and y values, defaulting to 
     the ranges of the x- and y-coordinates of the relevant corners
     of the maximal intersections.}
   \item{zlim}{The minimum and maximum values of the density for which colors 
should be plotted, defaulting to the range of the finite values of the 
density. Each of the given colors will be used to color an equispaced 
interval of this range. The midpoints of the intervals cover the range, 
so that values just outside the range will be plotted (see the 
documentation of \code{\link{image}}). This parameter is not used if 
\code{breaks} is specified.}
   \item{breaks}{Numeric vector with break points for the colors, satisfying
      length(breaks)=length(col)+1. This parameter overrides \code{zlim}.
   }
   \item{xlab, ylab}{Labels for the x- and y-axis. The default values are
  empty.}
   \item{main}{Title of the plot. The default value is empty.}
   \item{sub}{Sub title of the plot. The default value is empty.}
}
\value{No value is returned.}
\details{In many cases we assign specific values to represent +/- infinity and
(see, e.g., \code{\link{actg181}}). Note that these values 
determine the size of maximal intersections that extend to +/- infinity, 
and hence they also determine the value of the density at such maximal 
intersections. The value of the density at such maximal intersections is
therefore meaningless.}
\author{Marloes Maathuis: \email{marloes@u.washington.edu}}
\seealso{\code{\link{computeMLE}}}
\examples{
# Load example data:
data(ex)

# Compute the MLE:
mle <- computeMLE(ex)

# Bivariate density plots of the MLE: 
#   The colors represent the density=p/(area of maximal intersection)
par(mfrow=c(1,1))
plotDens2(mle, xlim=range(ex[,1:2]), ylim=range(ex[,3:4]), 
 main="Plot of the MLE. Colors represent the density.")
plotRects(ex, add=TRUE)

#   Alternative: numbers show the amount of mass in each maximal intersection
plotDens2(mle, col="lightgray", xlim=range(ex[,1:2]), 
 ylim=range(ex[,3:4]), numbers=TRUE, key=FALSE, 
 main="Plot of the MLE")
plotRects(ex, add=TRUE)
}
\keyword{hplot}
\keyword{dplot}
\concept{nonparametric maximum likelihood estimator}
\concept{censored data}

