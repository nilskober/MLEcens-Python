\name{plotCM}
\alias{plotCM}
\title{Plot a clique matrix}
\description{This function can be used to make an image of a clique matrix 
(or any other 0-1 matrix). It is basically just the function 
\code{\link{image}}, with some pre-defined settings.
}
\usage{
plotCM(cm, col=c("white","black"), at.x=NULL, at.y=NULL, 
       xlab="Observation rectangles", ylab="Maximal intersections", 
       main="", sub="")
}
\arguments{
   \item{cm}{A mxn clique matrix, where m is the number of maximal 
       intersections, and n is the number of observation rectangles. 
       The (i,j)th element is 1 if the ith maximal intersection 
       is contained in the jth observation rectangle, and it is 0 otherwise. 
   }
   \item{col}{Colors to be used. The default value is c("white","black").}
   \item{at.x}{The points at which tick-marks are to be drawn along the 
       x-axis. The default value is NULL, meaning that tickmark locations 
       are computed automatically. See also the documentation of 
       \code{\link{axis}}.}
   \item{at.y}{The points at which tick-marks are to be drawn along the 
       y-axis. The default value is NULL, meaning that tickmark locations 
       are computed automatically. See also the documentation of 
       \code{\link{axis}}.}
   \item{xlab}{Label for the x-axis. The default value is "Observation 
     rectangles".}
   \item{ylab}{Label for the y-axis. The default value is "Maximal 
     intersections".}
   \item{main}{Title of the plot. The default value is empty.}
   \item{sub}{Sub-title of the plot. The default value is empty.}
}
\value{No value is returned.}
\author{Marloes Maathuis: \email{marloes@u.washington.edu}}
\seealso{\code{\link{reduc}}}
\examples{
# Load example data and plot observation rectangles
data(ex)
par(mfrow=c(2,1))
plotRects(ex,main="Rectangles and maximal intersections")

# Perform reduction step and plot maximal intersections (shaded)
res<-reduc(ex, cm=TRUE)
plotRects(res$rects, density=15, border=NA, add=TRUE)

# Plot clique matrix
plotCM(res$cm, main="Clique matrix")
}
\keyword{hplot}
\concept{nonparametric maximum likelihood estimator}
\concept{censored data}
\concept{parameter reduction}
\concept{clique matrix}
