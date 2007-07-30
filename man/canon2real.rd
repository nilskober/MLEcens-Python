\name{canon2real}
\alias{canon2real}
\title{Transform (intersections of) canonical rectangles back to their original coordinates}
\description{
      This function transforms a set of (intersections of) canonical 
rectangles (see \code{\link{real2canon}} for a definition) back to 
their original coordinates. It performs the inverse 
operation of the function \code{\link{real2canon}}.
}
\usage{canon2real(Rcanon, R, B = c(0,1))}
\arguments{
   \item{Rcanon}{
       A mx4 matrix of (intersections of) canonical rectangles that are 
to be transformed back to their original coordinates.  
Each row corresponds to a rectangle, represented as (x1,x2,y1,y2). 
The point (x1,y1) is the lower left corner of the rectangle and 
(x2,y2) is the upper right corner of the rectangle. 
   }
   \item{R}{
       A nx4 matrix with the coordinates of the original rectangles. 
Each row corresponds to a rectangle, represented as (x1,x2,y1,y2). 
   }
   \item{B}{
       This describes the boundaries of the original rectangles 
(0=open or 1=closed). It can be specified in three ways:
       \item A nx4 matrix containing 0's and 1's. Each row corresponds to a
rectangle and is denoted as (cx1, cx2, cy1, cy2), where cx1 denotes the 
boundary type of x1, cx2 denotes the boundary type of x2, etc. 
       \item A vector (cx1, cx2, cy1, cy2) containing 0's and 1's. This 
representation can be used if all rectangles have the same type of boundaries. 
       \item A vector (c1, c2) containing 0's and 1's. This representation 
can be used if all x and y intervals have the same type of boundaries. c1 
denotes the boundary type of x1 and y1, and c2 denotes the boundary type of
x2 and y2.
   }
}
\details{
   The functions \code{\link{real2canon}} and \code{\link{canon2real}} are carried out 
   automatically in C-code as part of the functions \code{\link{reduc}} and \code{\link{computeMLE}}. 
   We chose to make the functions available separately as well, in order to illustrate our
   algorithm for computing the MLE. 

   As a first step in the computation of the MLE, we 
   transform rectangles into canonical rectangles, using \code{\link{real2canon}}). 
   This is useful for two reasons. Firstly, it forces us in the 
   very beginning to deal with possible ties and with the fact whether endpoints are 
   open or closed. As a consequence, we do not have to account for ties and open or 
   closed endpoints in the actual computation of the MLE. Secondly, it is convenient to work
   with the integer coordinates of the canonical rectangles in the computation of the MLE. 
   After all computations are done, we transform the canonical rectangles back to their original 
   coordinates, using \code{\link{canon2real}}. For more details, see Maathuis (2005, Section 2.1).
}
\value{
  A list with the following elements:
  \item{rects}{A mx4 matrix giving the original coordinates of the 
input rectangles. Each row (x1,x2,y1,y2) represents a rectangle. 
  }
  \item{bounds}{This describes the boundaries \code{rects}. 
  It is given in the same format as \code{B}.
  }
}
\references{
   M.H. Maathuis (2005). Reduction algorithm for the NPMLE for
   the distribution function of bivariate interval censored data.
   \emph{Journal of Computational and Graphical Statistics} \bold{14}
   252--262.
}
\author{Marloes Maathuis: \email{marloes@u.washington.edu}}
\seealso{\code{\link{real2canon}}}
\examples{
# An example with 3 arbitrarily chosen observation rectangles
R <- rbind(c(3.5, 4.2, 3.3, 9.1),    # first rectangle
           c(4.2, 4.9, 3, 4.5),      # second rectangle
           c(3.8, 5.1, 8.1, 9.5))    # third rectangle

# Plot the rectangles
par(mfrow=c(2,2))
plotRects(R, lwd=2, main="Original rectangles")

# Transform rectangles into canonical rectangles
res1 <- real2canon(R, c(0,1))   
plotRects(res1, grid=TRUE, lwd=2, main="Canonical rectangles")

# Transform canonical rectangles back to original coordinates   
res2 <- canon2real(res1, R, c(0,1))
plotRects(res2$rects, lwd=2, main="Original rectangles")
   
# Only transform rectangle (2,3)x(4,5), which is the  
#   the intersection of the canonical rectangles R1 and R3. 
#   The result is the intersection of the original rectangles R1 and R3. 
R.1.3 <- matrix(c(2,3,4,5),nrow=1)
res3 <- canon2real(R.1.3, R, c(0,1))
res3$rects
 
# Note that the algorithm keeps track of the boundaries of the rectangles:
B <- rbind(c(1,0,1,0),
           c(1,1,1,1),
           c(0,1,0,1))
res4 <- canon2real(R.1.3, R, B)
res4$bounds    
}
\keyword{manip}
\concept{nonparametric maximum likelihood estimator}
\concept{censored data}

