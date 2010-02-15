\name{actg181}
\docType{data}
\alias{actg181}
\title{Data from the Aids Clinical Trials Group protocol ACTG 181}
\description{
  An example data set with bivariate interval censored data. 
  The data come from the AIDS Clinical Trials Group protocol ACTG 181, 
  and contain information on 
  the time (in months) to shedding of cytomegalovirus (CMV) 
  in the urine and blood and the time (in months) 
  to colonization of mycobacterium 
  avium complex (MAC) in the sputum and stool (Betensky and Finkelstein, 1999). 
}
\usage{data(actg181)}
\format{A matrix containing 204 rows and 4 columns. Each row    
(x1,x2,y1,y2) corresponds to a subject in the study, and 
represents the rectangle that is known to contain the unobservable
times of CMV shedding (x) and MAC colonization (y) of this person:
x1<=x<=x2 and y1<=y<=y2. The times are 
given in months. We use the values +/- 100 to represent +/- infinity.
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
\seealso{\code{\link{actg181Mod}}}
\examples{
# Load the data
data(actg181)

# Compute the MLE
mle <- computeMLE(R=actg181, B=c(1,1,1,1))

# Create CDF plots of the MLE: 
# (Maximal intersections are denoted in red)
par(mfrow=c(2,2))

# Lower bound for bivariate CDF
plotCDF2(mle, bound="l", xlim=c(-1,101), ylim=c(-1,101),
 n.key=5, main="Bivariate CDF (lower bound)", 
 xlab="time to CMV shedding (months)", 
 ylab="time to MAC colonization (months)")
plotRects(mle$rects, border="red", add=TRUE)

# Upper bound for bivariate CDF
plotCDF2(mle, bound="u", xlim=c(-1,101), ylim=c(-1,101),
 n.key=5, main="Bivariate CDF (upper bound)", 
 xlab="time to CMV shedding (months)", 
 ylab="time to MAC colonization (months)")
plotRects(mle$rects, border="red", add=TRUE)

# Marginal CDF for X
plotCDF1(mle, margin=1, xlim=c(0,90), 
 main="CDF for time to CMV shedding",   
 xlab="t (months)", ylab="P(time to CMV shedding <= t)")

# Marginal CDF for Y
plotCDF1(mle, margin=2, xlim=c(0,90), 
 main="CDF for time to MAC colonization", 
 xlab="t (months)",  ylab="P(time to MAC colonization <= t)")

# Note that the difference between the upper and lower bound 
# of the MLE (because of representational non-uniqueness)
# is large, especially for the time to MAC colonization. 
}
\keyword{datasets}
\concept{censored data}

