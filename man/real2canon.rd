\name{real2canon}
\alias{real2canon}
\title{Transform a set of rectangles into canonical rectangles}
\description{
      This function transforms a set of rectangles R1,..,Rn into canonical 
      rectangles R1',..,Rn' with the following properties:
      \item R1',..Rn' have the same intersection structure as the original rectangles, i.e., 
      Ri' and Rj' intersect if and only if Ri and Rj intersect.
      \item All x-coordinates of R1',..,Rn' are distinct, and take values in {1,..,2n}.
      \item All y-coordinates of R1',..,Rn' are distinct, and take values in {1,..,2n}.

      The function \code{\link{real2canon}} performs the inverse operation of the function 
      \code{\link{canon2real}}.
}
\usage{real2canon(R,B=c(0,1))}
\arguments{
   \item{R}{
       A nx4 real matrix of rectangles. Each row corresponds to
       a rectangle, represented as (x1,x2,y1,y2). The point (x1,y1) is the 
       lower left corner of the rectangle and the point (x2,y2) is the upper right
       corner of the rectangle. 
   }
   \item{B}{
       This describes the boundaries of the rectangles (0=open or 1=closed). 
       It can be specified in three ways:
       \item A nx4 matrix containing 0's and 1's. Each row corresponds to a
       rectangle, and is denoted as (cx1, cx2, cy1, cy2). Here cx1 denotes the boundary 
       type of x1, cx2 denotes the boundary type of x2, etc.
       \item A vector (cx1, cx2, cy1, cy2) containing 0's and 1's. This representation 
       can be used if all rectangles have the same type of boundaries.
       \item A vector (c1, c2) containing 0's and 1's. This representation can be 
       used if all x and y intervals have the same type of boundaries. 
       c1 denotes the boundary type of x1 and y1, and c2 denotes the boundary type of 
       x2 and y2.
       The default value is c(0,1).
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
\value{A nx4 matrix of canonical observation rectangles. Each row 
      (x1,x2,y1,y2) represents a rectangle.
  }
}
\references{
   M.H. Maathuis (2005). Reduction algorithm for the NPMLE for
   the distribution function of bivariate interval censored data.
   \emph{Journal of Computational and Graphical Statistics} \bold{14}
   252--262.
}
\author{Marloes Maathuis: \email{marloes@u.washington.edu}}
\seealso{\code{\link{canon2real}}}
\examples{
# An example with 3 arbitrarily chosen observation rectangles
R <- rbind(c(3.5, 4.2, 3.3, 9.1),   # first rectangle
           c(4.2, 4.9, 3, 4.5),     # second rectangle
           c(3.8, 5.1, 8.1, 9.5))   # third rectangle

# Plot the rectangles
par(mfrow=c(2,2))
plotRects(R, lwd=2, main="Original rectangles")

# Transform rectangles to canonical rectangles. Since the 
#   boundaries of R1 and R2 coincide, it matters which boundaries
#   we define to be open or closed. 

# With boundary structure c(0,1), R1 and R2 do *not* overlap:
res1 <- real2canon(R, c(0,1))    
plotRects(res1, grid=TRUE, lwd=2, main="Canonical rectangles
boundary c(0,1)")

# But with boundary structure c(1,1), R1 and R2 *do* overlap:
res2<-real2canon(R, c(1,1))
plotRects(res2, grid=TRUE, lwd=2, main="Canonical rectangles
boundary c(1,1)")
}
\keyword{manip}
\concept{nonparametric maximum likelihood estimator}
\concept{censored data}

