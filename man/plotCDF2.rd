\name{plotCDF2}
\alias{plotCDF2}
\title{Create a bivariate CDF (or survival function) plot of the MLE}
\description{This function plots the MLE for the bivariate CDF of (X,Y) (or
the bivariate survival function). 
The value of the estimate at the point (x,y) is computed 
by summing all probability mass of the MLE that falls in the region 
(-infinity,x] x (-infinity,y]. The plot uses colors/shades to represent
the value of the MLE, and is generated using the function 
\code{\link{image}}. 
}
\usage{
plotCDF2(mle, bound, col=gray(seq(0.9,0.3,len=30)), surv=FALSE, 
         key=TRUE, n.key=10, round.key=2, cex.key=0.6, xlim=NULL, 
         ylim=NULL, zlim=NULL, breaks=NULL, xlab="", ylab="", 
         main="", sub="")
}
\arguments{
   \item{mle}{List with elements 'p' and 'rects', as outputted by 
       \code{\link{computeMLE}}.}
   \item{bound}{Parameter taking the values "u" and "l". It 
indicates how  representational non-uniqueness of
the MLE should be handled. 
Option "u" (upper) indicates an upper bound, obtained by assigning all mass 
to the lower left corners of the maximal intersections. Option "l" (lower) 
indicates a lower bound, obtained by assigning all mass to the upper right 
corners of the maximal intersections.}
   \item{col}{Color vector used to represent the values of the 
MLE. The default value is \code{gray(seq(0.9,0.3,len=30))}.}
   \item{surv}{Logical. If FALSE, the bivariate CDF P(X<=x, Y<=y) is plotted.
 If TRUE, the bivariate survival function 
 1-P(X<=x,Y<=y) is plotted.  The default value is FALSE.
}
   \item{key}{Logical, indicating if a color key should be drawn. The 
              default value is TRUE.}
   \item{n.key}{Approximate number of tickmarks for the color key.
                The default value is 10.}
   \item{round.key}{Number of decimals used for the labels of the color key.
                    The default value is 2.}
   \item{cex.key}{Numerical value giving the amount by which text 
  in the key should be scaled relative to the default. The default value is 
0.6.}
   \item{xlim, ylim}{Ranges for the plotted x and y values, defaulting to 
     the ranges of the x- and y-coordinates of the relevant corners
     of the maximal intersections.}
   \item{zlim}{The minimum and maximum values of the estimate for which colors 
should be 
plotted, defaulting to the range of the finite values of the estimate. 
Each of the given colors will be used to color an equispaced interval of 
this range. The midpoints of the intervals cover the range, so that values 
just outside the range will be plotted (see the documentation of 
\code{\link{image}}). This parameter is not used if \code{breaks} is 
specified.}
   \item{breaks}{Numeric vector with break points for the colors, satisfying
      length(breaks)=length(col)+1. This parameter overrides \code{zlim}.
   }
   \item{xlab, ylab}{Labels for the x- and y-axis. The default values are 
empty.}
   \item{main}{Title of the plot. The default value is empty.}
   \item{sub}{Sub title of the plot. The default value is empty.}
}
\value{No value is returned.}
\author{Marloes Maathuis: \email{maathuis@stat.math.ethz.ch}}
\seealso{\code{\link{computeMLE}}}
\examples{
# Load example data:
data(ex)

# Compute the MLE:
mle <- computeMLE(ex)

### Bivariate CDF plot of the MLE

# Plot lower bound for representational non-uniqueness
par(mfrow=c(1,1))
plotCDF2(mle, xlim=c(min(ex[,1])-1,max(ex[,2])+1), 
 ylim=c(min(ex[,3])-1, max(ex[,4])+1), bound="l", n.key=4,
 main="Bivariate CDF plot of the MLE,
 lower bound")

# Add observation rectangles and shaded maximal intersections
plotRects(ex, add=TRUE) 
plotRects(mle$rects, density=20, border=NA, add=TRUE) 

# Plot upper bound for representational non-uniqueness
plotCDF2(mle, xlim=c(min(ex[,1])-1,max(ex[,2])+1), 
 ylim=c(min(ex[,3])-1, max(ex[,4])+1), bound="u", n.key=4,
 main="Bivariate CDF plot of the MLE,
 upper bound")

# Add observation rectangles and shaded maximal intersections
plotRects(ex, add=TRUE)
plotRects(mle$rects, density=20, border=NA, add=TRUE)
}
\keyword{hplot}
\keyword{dplot}
\concept{nonparametric maximum likelihood estimator}
\concept{censored data}

