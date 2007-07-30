\name{actg181Mod}
\docType{data}
\alias{actg181Mod}
\title{Modified data from the Aids Clinical Trials Group protocol ACTG 181}
\description{
  An example data set with bivariate interval censored data. 
  The data come from the AIDS Clinical Trials Group protocol ACTG 181, 
  and contain information on 
  the time (in months) to shedding of cytomegalovirus (CMV) 
  in the urine and blood and the time (in months) 
  to colonization of mycobacterium 
  avium complex (MAC) in the sputum and stool (Betensky and Finkelstein, 1999). 

  The format of the data has been modified 
  to allow for easy plotting with the functions \code{\link{plotHM}}, 
  \code{\link{plotDens1}} and \code{\link{plotDens2}} (see section 'Format').   
}
\usage{data(actg181Mod)}
\format{A matrix containing 204 rows and 4 columns. Each row    
  (x1,x2,y1,y2) corresponds to a subject in the study, and 
  represents the rectangle that is known to contain the unobservable 
  times of CMV shedding (x) and MAC colonization (y) of this person:
  x1<=x<=x2 and y1<=y<=y2. The times are 
  given in months. We use the values +/- 100 to represent +/- infinity.
  
  In order to allow easy plotting with \code{\link{plotHM}}, 
  \code{\link{plotDens1}} and \code{\link{plotDens2}}, the x- and y- intervals
  were modified as follows:
  [x1,x2] was changed into [x1-0.5, x2+0.5]
  and [y1,y2] was changed into [y1-0.5,y2+0.5]. 
}
\details{Extracted from Betensky and Finkelstein (1999): 
The data describe 204 of the 232 subjects in the study
who were tested for CMV shedding and MAC colonization at least once
during the trial, and did not have a prior CMV or MAC diagnosis. 
Tests were performed during clinic visits, scheduled at regular monthly
intervals. For patients who did not miss any clinic visits, 
the time of event was
recorded as the month that the first positive test occurred, resulting 
in discrete failure time data. For patients who missed some visits, and 
who were detected to be positive directly following one or more missed
visits, the event time was recorded as having occurred in a time interval, 
resulting in discrete interval censored failure time data. All visit times
were rounded to the closest quarter.

One should use closed boundaries (B=c(1,1,1,1)) in order to reproduce the 
results of Betensky and Finkelstein (1999). In that case the probability
masses of the MLE that we find are exactly equal to those given in Table 
IV of Betensky and Finkelstein, but there are some discrepancies in 
the maximal intersections (compare rows 1, 2, 7, 8 and 10 of their table IV).
}
\source{Betensky and Finkelstein (1999). A non-parametric maximum 
   likelihood estimator for bivariate interval censored data. 
   \emph{Statistics in Medicine} \bold{18} 3089-3100.}
\seealso{\code{\link{actg181}}}
\examples{
# Load the data
data(actg181Mod)

# Compute the MLE
mle <- computeMLE(R=actg181Mod, B=c(1,1,1,1))

# Create density plots
par(mfrow=c(2,2))

#   Bivariate density plot
plotDens2(mle, main="Bivariate density",
 xlab="time to CMV shedding (months)", 
 ylab="time to MAC colonization (months)")

#   Marginal density plot for time to MAC colonization
plotDens1(mle, margin=2, main="Density for time
 to MAC colonization", xlab="t (months)", 
 ylab="density")

#   Marginal density plot for time to CMV shedding
plotDens1(mle, margin=1, main="Density for time 
 to CMV shedding", xlab="t (months)", 
 ylab="density")

# Note that many maximal intersections extend to 
# infinity, and hence the value of the density is
# not very meaningful. 
}
\keyword{datasets}
\concept{censored data}

