\name{plotRects}
\alias{plotRects}
\title{Plot a set of rectangles}
\description{This function can be used to plot a set of rectangles in a 
   new plot, or to add them to an existing plot. It is basically the 
   function \code{\link{rect}}, with some pre-defined settings. 
}
\usage{
plotRects(R, grid=FALSE, grid.lty=3, grid.col="lightgray", density=NULL, 
          angle=45, col=NA, border=rainbow(nrow(R)), lty=1, lwd=1, 
          add=FALSE, xlim=NULL, ylim=NULL, xlab="", ylab="", main="", sub="")
}
\arguments{
   \item{R}{A nx4 matrix of rectangles. Each row corresponds to an 
       rectangle, represented as (x1,x2,y1,y2). The point (x1,y1) is 
       the lower left corner of the rectangle and the point (x2,y2) is the 
       upper right corner of the rectangle.
   }
   \item{grid}{Logical, indicating if a grid should be drawn. 
The default value is FALSE.}
   \item{grid.lty}{Line type of the grid lines. The default value is 3=dotted.}
   \item{grid.col}{Line color of the grid lines. The default value is light 
gray.}
   \item{density}{Density
     of shading lines for the rectangles, in lines per inch. 
The default value is NULL, meaning that 
     no shading lines are drawn. A zero value of density means no shading 
     lines whereas negative values (and NA) suppress shading (and so allow 
     color filling).
   }
   \item{angle}{Angle (in 
     degrees) of the shading lines.
   }
   \item{col}{Color(s) to 
     fill or shade the rectangles with. The default NA (or also NULL) means 
     do not fill, i.e., draw transparent rectangles, unless density is 
     specified.
   }
   \item{border}{Color for 
     rectangle borders. Use border = NA to omit borders. If there are shading
     lines, border = TRUE means use the same color for the border as for the 
     shading lines. The default value is rainbow(n): a vector of n contiguous 
     colors from the palette \code{\link{rainbow}}. 
   }
   \item{lty}{Line type for 
      borders and shading of the rectangles. The default value is 1="solid".} 
   \item{lwd}{Line width for 
      borders and shading of the rectangles. The default value is 1.}
   \item{add}{Logical, indicating if the rectangles should be added 
to an existing plot. The default value is FALSE. 
   } 
   \item{xlim}{Range of the x-axis. The default value is the range of 
x-coordinates of the rectangles.}
   \item{ylim}{Range of the y-axis. The default value is the range of 
y-coordinates of the rectagnles.}
   \item{xlab, ylab}{Labels of the x and y axis. The default values are empty.}
   \item{main}{Title of the plot. The default value is empty.}
   \item{sub}{Sub title of the plot. The default value is empty.}
}
\value{No value is returned.}
\author{Marloes Maathuis: \email{marloes@u.washington.edu}}
\seealso{\code{\link{rect}}, \code{\link{palette}}}
\examples{
n <- 10
x <- c(0:(n-1))
R <- cbind(x,x+3,x,x+3)  # first rectangle is (0,3)x(0,3), 
                         # second rectangle is (1,4)x(1,4), etc...
par(mfrow=c(1,1))
plotRects(R,main="Example")
}
\keyword{hplot}
\keyword{aplot}
\concept{nonparametric maximum likelihood estimator}
\concept{censored data}

