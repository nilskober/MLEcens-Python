\name{plotHM}
\alias{plotHM}
\title{Plot a height map}
\description{This function can be used to plot a 'height map' of a set
of rectangles in a new plot, or to add it to an existing plot. The value of 
the heightmap at a point equals the number of rectangles that overlap at 
this point.}
\usage{
plotHM(hm, R, grid=TRUE, grid.lty=3, grid.col="lightgray", key=TRUE, 
       n.key=10, cex.key=0.6, numbers=FALSE, col=terrain.colors(max(hm)+1), 
       xlim=NULL, ylim=NULL, xlab="", ylab="", main="", sub="")
}
\arguments{
   \item{hm}{A (2n+1)x(2n+1) matrix with the values of the heightmap, 
    outputted by the function \code{\link{reduc}}.}
   \item{R}{A nx4 matrix with the real or canonical observation rectangles. 
         Each row corresponds to
       a rectangle, represented as (x1,x2,y1,y2). The point (x1,y1) is the 
       lower left corner of the rectangle and the point (x2,y2) is the upper 
       right corner of the rectangle. }
   \item{grid}{Logical, indicating if a grid should be drawn. The default
    value is TRUE.}
   \item{grid.lty}{Line type of the grid lines. The default value is 3=dotted.}
   \item{grid.col}{Color of the grid lines. The default value is light gray.}
   \item{key}{Logical, indicating if a color key should be drawn in 
    the right margin of the plot. The default value is TRUE.}
   \item{n.key}{Approximate number of tickmarks for the color key. 
      The default value is 10.}
   \item{cex.key}{Numerical value giving the amount by which text 
  in the key should be scaled relative to the default. The default value is 
  0.6.}
   \item{numbers}{Logical, indicating if numbers should be plotted in 
   the grid cells that indicate the values of the height map.}
   \item{col}{Color vector used to represent the values of the height map. 
   The length of \code{col} must equal \code{max(hm)+1}. }
   \item{xlim}{Range for the plotted x values, defaulting to 
   \code{c(min(x-coordinates)-1, max(x-coordinates)+1).} }
   \item{ylim}{Range for the plotted y values, defaulting to
   \code{c(min(x-coordinates)-1, max(x-coordinates)+1)}.}
   \item{xlab, ylab}{Labels of the x and y axis. The default values are empty.}
   \item{main}{Title of the plot. Default value is empty.}
   \item{sub}{Sub title of the plot. Default value is empty.}
}
\details{We chose to create a thin color key that fits
  in the margin of the plot. In this way, the 
  plotting margins and plotting region do not have to be adjusted, so that
  other elements can be easily added to the plot later on (like observation 
  rectangles or maximal intersections).
}
\value{No value is returned.}
\author{Marloes Maathuis: \email{marloes@u.washington.edu}}
\seealso{\code{\link{reduc}}, \code{\link{plotRects}}, 
  \code{\link{palette}}
}
\examples{
# Load example data
data(ex)

# Perform reduction step
res <- reduc(ex, hm=TRUE)

# Plot the height map:
par(mfrow=c(1,1))
plotHM(res$hm, ex, main="Height map")

# Add observation rectangles in black:
plotRects(ex, add=TRUE, border="black")

# Add shaded maximal intersections:
plotRects(res$rects, add=TRUE, border=NA, density=15)

# Compute heightmap of the canonical rectangles:
canon <- real2canon(ex)
res2 <- reduc(canon, hm=TRUE)
# Note that res$hm and res2$hm are identical. So we only need to change
# the x- and y-coordinates of the height map.

# Plot height map of the canonical rectangles
plotHM(res$hm, canon, key=FALSE, numbers=TRUE, main="Canonical height map")

# Add canonical rectangles in black:
plotRects(canon, add=TRUE, border="black")

# Add canonical maximal intersections (local maxima of height map) in red:
plotRects(res2$rects, add=TRUE, border="red")
}
\keyword{hplot}
\concept{nonparametric maximum likelihood estimator}
\concept{censored data}
\concept{parameter reduction}
