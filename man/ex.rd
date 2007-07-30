\name{ex}
\docType{data}
\alias{ex}
\title{Example data set (artificial)}
\description{
  Example data set with artificial bivariate interval censored data. These
  data are used to illustrate various functions of this R-package.  
}
\usage{data(ex)}
\format{A matrix containing 6 rows and 4 columns. Each row    
(x1,x2,y1,y2) represents a rectangle that is known to contain 
the unobservable realization of the variables of interest (X,Y).
The point (x1,y1) is the lower left corner of the rectangle and
(x2,y2) is the upper right corner of the rectangle.
}
\examples{
# Load the data
data(ex)

# Plot the rectangles
par(mfrow=c(1,1))
plotRects(ex)
}
\keyword{datasets}
\concept{censored data}

