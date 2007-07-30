\name{plotDens1}
\alias{plotDens1}
\title{Create a univariate density plot of the MLE}
\description{This function creates a univariate density plot 
of the MLE. To be precise, it can plot the density (d/dx) P(X<=x, a<Y<=b)
as a function of x, and (d/dy) P(a<X<=b, Y<=y) as a function of y, 
under the assumption that the mass is distributed uniformly over the 
maximal intersections. 
}
\usage{
plotDens1(mle, margin, int=NULL, col=1, lty=1, add=FALSE,
          xlim=NULL, ylim=NULL, xlab="", ylab="", main="", sub="")
}
\arguments{
   \item{mle}{List with elements 'p' and 'rects', as outputted by 
       \code{\link{computeMLE}}.}
   \item{margin}{Indicates which margin should be plotted: 1 = x-margin, 
     2 = y-margin. So if margin=1, the MLE for (d/dx) P(X<=x, a<Y<=b) is 
plotted, and if margin=2, then the MLE for (d/dy) P(a<X<=b, Y<=y) is plotted. 
   }
   \item{int}{This indicates the range of interest of the variable that was 
       \emph{not} chosen in \code{margin}.
       If \code{int} is specified, it should be of the form c(a,b), with a<b. 
       If margin=1, the MLE for (d/dx) P(X<=x, a<Y<=b) is plotted as a 
       function of x. If 
       margin=2, the MLE for (d/dy) P(a<X<=b, Y<=y) is plotted as a function 
       of y. This parameter defaults to (-infinity,infinity), yielding plots 
       of the marginal density of X and Y.  
   }
   \item{col}{Line color. The default value is 1="black".}
   \item{lty}{Line type. The default value is 1="solid".}
    \item{add}{Logical, indicating if the lines should be added to an existing
      plot. The default value is FALSE.
   }
   \item{xlim}{Range for the horizontal axis, defaulting 
    to the range of x-coordinates (if margin=1) or y-coordinates (if margin=2)
    of the relevant corners of maximal interesctions. 
   }
   \item{ylim}{Range for the vertical axis, defaulting to the range of 
    values of the estimate.}
   \item{xlab,ylab}{Labels of the x- and y-axis. The default values are empty.}
   \item{main}{Title of the plot.}
   \item{sub}{Sub title of the plot.}

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

# Bivariate density plot of the MLE:
#   Numbers represent the mass p in the maximal intersections
par(mfrow=c(2,2))
plotDens2(mle, xlim=range(ex[,1:2]), ylim=range(ex[,3:4]), 
 col="lightgray", main="Bivariate density plot of the MLE", 
 key=FALSE, numbers=TRUE)
plotRects(ex, add=TRUE)

# Univariate density plots of the MLE:

#   Plot of the marginal density of Y
plotDens1(mle, margin=2, xlim=range(ex[,3:4]), 
 main="Marginal density plot, 
 y-margin", xlab="y", ylab=expression(f[Y](y))) 

#   Plot of the marginal density of X 
plotDens1(mle, margin=1, xlim=range(ex[,1:2]), 
 main="Marginal density plot, 
 x-margin", xlab="x", ylab=expression(f[X](x)))
}
\keyword{hplot}
\keyword{dplot}
\keyword{aplot}
\concept{nonparametric maximum likelihood estimator}
\concept{censored data}

