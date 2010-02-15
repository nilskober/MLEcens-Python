\name{plotCDF1}
\alias{plotCDF1}
\title{Create a marginal CDF (or survival function) plot of the MLE}
\description{This function plots the MLE for a marginal (sub)-CDF (or 
survival function) for 
one of the two variables of interest. To be precise, it can plot the 
MLE for P(X<=x, a<Y<=b) (or 1-P(X<=x, a<Y<=b)) as a function of x, 
and the MLE for P(a<X<=b, Y<=y) (or 1-P(a<X<=b, Y<=y))
as a function of y, where a and b may take the values -infinity and 
infinity. The values of these estimates are computed by summing all 
probability mass of the MLE that falls in the regions (-infinity,x] x (a,b] 
and (a,b] x (-infinity,y], respectively.
}
\usage{
plotCDF1(mle, margin, bound="b", int=NULL, surv=FALSE,
         add=FALSE, col=1, lty=1, xlim=NULL, 
         ylim=NULL, xlab="", ylab="", main="", sub="")
}
\arguments{
   \item{mle}{List with elements 'p' and 'rects', as outputted by 
       \code{\link{computeMLE}}.}
   \item{margin}{Indicates which margin should be plotted: 1 = x-margin, 
     2 = y-margin. So if margin=1, the MLE for P(X<=x, a<Y<=b) is plotted, 
     and if margin=2, then the MLE for P(a<X<=b, Y<=y) is plotted. 
   }
   \item{bound}{Parameter taking the values "u", "l" and "b". It 
indicates how representational non-uniqueness of
the MLE should be handled. Option "u" (upper) indicates an upper bound, 
obtained by assigning all mass to the 
lower left corners of the maximal intersections. Option "l" (lower) indicates a 
lower bound, obtained by assigning all mass to the upper right corners of 
the maximal 
intersections. Option "b" (both) indicates that both the upper and the lower
bound should be plotted. The default value is "b".
}
   \item{int}{This indicates the range of interest of the variable that was 
       \emph{not} chosen in \code{margin}.
       If \code{int} is specified, it should be of the form c(a,b), with a<b. 
       If margin=1, the MLE for P(X<=x, a<Y<=b) is plotted as a 
       function of x. If 
       margin=2, the MLE for P(a<X<=b, Y<=y) is plotted as a function 
       of y. This parameter defaults to (-infinity,infinity), yielding plots 
       of the estimates for P(X<=x) and P(Y<=y). 
   }
   \item{surv}{Logical. The default value is FALSE. 
    If TRUE, the function 
    1-P(X<=x, a<Y<=b) is plotted instead of P(X<=x, a<Y<=b), and 
    the function
    1-P(a<X<=b, Y<=y) is plotted instead of P(a<X<=b, Y<y).
   }
   \item{add}{Logical, indicating if the lines should be added to an existing
      plot. The default value is FALSE.
   }
   \item{col}{Line color. The default value is 1="black".}
   \item{lty}{Line type. The default value is 1="solid".}
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
\author{Marloes Maathuis: \email{maathuis@stat.math.ethz.ch}}
\seealso{\code{\link{computeMLE}}}
\examples{
# Load example data:
data(ex)

# Compute the MLE:
mle <- computeMLE(ex)

# Plot marginal CDF for X
par(mfrow=c(2,2))
plotCDF1(mle, margin=1, xlim=c(min(ex[,1])-1,max(ex[,2])+1), 
 bound="b", xlab="x", ylab="P(X<=x)", main="MLE for P(X<=x)")

# Plot marginal survival function for X
plotCDF1(mle, margin=1, surv=TRUE, xlim=c(min(ex[,1])-1,max(ex[,2])+1), 
 bound="b", xlab="x", ylab="P(X>x)", main="MLE for P(X>x)")

# Plot marginal CDF for Y
plotCDF1(mle, margin=2, xlim=c(min(ex[,3])-1,max(ex[,4])+1), 
 bound="b", xlab="y", ylab="P(Y<=y)", main="MLE for P(Y<=y)")

# Plot marginal survival function for Y
plotCDF1(mle, margin=2, surv=TRUE, xlim=c(min(ex[,3])-1,max(ex[,4])+1), 
 bound="b", xlab="y", ylab="P(Y>y)", main="MLE for P(Y>y)")
}
\keyword{hplot}
\keyword{dplot}
\keyword{aplot}
\concept{nonparametric maximum likelihood estimator}
\concept{censored data}

