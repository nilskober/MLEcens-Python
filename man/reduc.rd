\name{reduc}
\alias{reduc}
\title{Determine areas of possible mass support of the MLE}
\description{
      The MLE for censored data can only assign mass to a finite set distinct 
      regions, called maximal intersections. The function \code{\link{reduc}} 
      computes 
      these areas, using the height map algorithm described in Maathuis (2005).
      In addition to the maximal intersections, the function can also output 
      the height map and the clique matrix.  
}
\usage{reduc(R,B=c(0,1),hm=FALSE,cm=FALSE)}
\arguments{
   \item{R}{
       A nx4 matrix containing the observation rectangles. Each row 
       corresponds to
       a rectangle, represented as (x1,x2,y1,y2). The point (x1,y1) is the 
       lower left corner of the rectangle and the point (x2,y2) is the upper 
       right corner of the rectangle. 
   }
   \item{B}{
       This describes the boundaries of the rectangles (0=open or 1=closed). 
       It can be specified in three ways:
       (1) A nx4 matrix containing 0's and 1's. Each row corresponds to a
       rectangle, and is denoted as (cx1, cx2, cy1, cy2). Here cx1 denotes the 
       boundary type of x1, cx2 denotes the boundary type of x2, etc. 
       (2) A vector (cx1, cx2, cy1, cy2) containing 0's and 1's. This 
       representation can be used if all rectangles have the same type of 
       boundaries.
       (3) A vector (c1, c2) containing 0's and 1's. This representation can 
       be used if all x and y intervals have the same type of boundaries. 
       c1 denotes the boundary type of x1 and y1, and c2 denotes the boundary 
       type of x2 and y2.

      The default value is c(0,1).
   }
   \item{hm}{
      Logical, indicating if the heightmap must be outputted. 
      The default value is FALSE. The height map 
      is a (2n+1)x(2n+1) matrix, where n is the number of observation rectangles. 
      Its values give the number of rectangles that 
      overlap at any given point.  
   }
   \item{cm}{
      Logical, indicating if the clique matrix must be outputted. 
      The default value is FALSE. The clique matrix
      is a mxn matrix, where m is the number of maximal intersections and n
      is the number of observation rectangles.  
      The (i,j)th element of the matrix is 1 if the ith maximal intersection 
      is contained in the jth observation rectangle, and it is 0 otherwise.  
   }
}
\details{
   The computation of the MLE for censored data can be split into two steps:
   a reduction step and an optimization step. In the reduction step, the areas
   of possible mass support are computed. Next, in the optimization step, it is 
   determined how much probability mass should be assigned to each of these 
   areas. The function \code{\link{reduc}} can be used for the reduction step. 
   It is carried out automatically as part of the function 
   \code{\link{computeMLE}}.
 
   The time and space complexity of the function \code{\link{reduc}} depend 
   on the parameters \code{hm} and \code{cm}.
   If \code{hm}=FALSE and \code{cm}=FALSE, then the algorithm is O(n^2)
   in time and O(n) in memory space. If \code{hm}=TRUE and \code{cm}=FALSE, 
   then the algorithm 
   is O(n^2) in both time and space. If \code{cm}=TRUE, then the algorithm is 
   O(n^3) in time and space. 

   The function \code{\link{reduc}} 
   uses the height map algorithm
   of Maathuis (2005). It first converts the observation rectangles to 
   canonical rectangles. Next, it computes the local maxima of the height map, 
   and then it convert these back to the original coordinates. 
   This process can be mimicked by hand as 
   follows: (1) use \code{\link{real2canon}} to convert the 
   rectangles to canonical rectangles; (2) use \code{\link{reduc}}
   to find the canonical maximal intersections 
   (local maxima of the height map of the canonical rectangles); (3)
   use \code{\link{canon2real}}
   to convert the canonical maximal intersections back to the original 
   coordinates. 
}
\value{A list containing the following elements:
   \item{rects}{A mx4 matrix of maximal intersections. Each row 
      (x1,x2,y1,y2) represents a maximal intersection, i.e., an area where
      the MLE can possibly assign mass.}
   \item{bounds}{This describes the boundaries of \code{rects}. It is 
   given in the same format as \code{B}.}
   \item{hm}{(Optional) A (2n)x(2n) matrix containing the height map. 
   Its values represent the number of rectangles that overlap at any 
   given point.
   }
   \item{cm}{(Optional) A mxn matrix containing the clique matrix. 
The (i,j)th element of the matrix is 
1 if the ith maximal intersection is contained in the jth observation rectangle, and it is 0 otherwise. 
   }
}
\references{
   M.H. Maathuis (2005). Reduction algorithm for the NPMLE for
   the distribution function of bivariate interval censored data.
   \emph{Journal of Computational and Graphical Statistics} \bold{14}
   252--262.
}
\author{Marloes Maathuis: \email{maathuis@stat.math.ethz.ch}}
\seealso{\code{\link{real2canon}}, \code{\link{canon2real}}, 
\code{\link{computeMLE}},
 \code{\link{plotCM}}, \code{\link{plotHM}}}
\examples{
# Load example data:
data(ex)
par(mfrow=c(1,1))

# Plot the observation rectangles
plotRects(ex,main="Example")

# Perform the reduction step
res<-reduc(ex, hm=TRUE, cm=TRUE)

# Shade the maximal intersections
plotRects(res$rects, density=15, add=TRUE, border=NA)

# Plot the height map, together with the observation 
# rectangles (in black) and the maximal intersections (shaded)
plotHM(res$hm, ex)
plotRects(ex, add=TRUE, border="black")
plotRects(res$rects, add=TRUE, border=NA, density=15)

# Print the clique matrix 
res$cm

# Make a plot of the clique matrix (useful for large data sets)
plotCM(res$cm)
}

\keyword{survival}
\keyword{nonparametric}
\keyword{manip}
\keyword{optimize}
\concept{nonparametric maximum likelihood estimator}
\concept{censored data}
\concept{parameter reduction}
