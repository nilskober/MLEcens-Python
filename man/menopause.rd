\name{menopause}
\docType{data}
\alias{menopause}
\title{Menopause data}
\description{
  An example data set with interval censored data and competing risks. 
  The data come from Cycle I of the 
  Health Examination Survey of the National Center for Health Statistics, 
  and contain information on the menopausal status of 2423 women 
  (MacMahon and Worcestor, 1966). 
}
\usage{data(menopause)}
\format{A matrix containing 2423 rows and 4 columns. Each row    
(x1,x2,y1,y2) corresponds to a subject in the study. The interval (x1,x2] 
contains the unobservable age of menopause X. The interval [y1,y2]
contains the type of menopause Y, where Y=1 represents operative
menopause and Y=2 represents natural menopause. We use the value 100
to represent infinity.  
}
\details{The Health Examination Survey used a nationwide 
  probability sample of people 
  between age 18 and 79 from the United States civilian, noninstitutional
  population. The participants were asked to complete a self-administered   
  questionnaire. The sample contained 4211 females, of whom 3581 completed
  the questionnaire. We restrict attention to the age range 25-59 years. 
  Furthermore, seven women who were less than 35 years of age and reported 
  having had a natural menopause were excluded as being an error
  or abnormal. The remaining data set contains information on 2423 women.
 
  MacMahon and Worcestor (1966) found that there was 
  marked terminal digit clustering in the response of this question, 
  especially for women who had a natural menopause. 
  Therefore, Krailo and Pike (1983) decided to only consider the menopausal
  status of women at the time of the questionnaire, yielding current 
  status data on the age of menopause with two competing 
  risks: operative menopause and natural menopause. 
}
\source{MacMahon and Worcestor (1966). Age at menopause, United States 
1960 - 1962. \emph{National Center for Health Statistics. Vital and Health 
Statistics}, volume 11, number 19.}
\references{
Krailo and Pike (1983). Estimation of the distribution of age at 
natural menopause from prevalence data. \emph{American Journal of 
Epidemiology} \bold{117} 356-361.
}
\seealso{\code{\link{menopauseMod}}}
\examples{
# Load the data
data(menopause)

# Compute the MLE
mle <- computeMLE(R=menopause, B=c(0,1,1,1))

# Plot first sub-distribution function P(X<=x, 0.5<Y<=1.5) = P(X<=x, Y=1)
par(mfrow=c(1,1))
plotCDF1(mle, margin=1, bound="b", int=c(0.5,1.5), col="red", ylim=c(0,1),
 xlab="x", main="P(X<=x, Y=k), k=1,2")

# Plot second sub-distribution function P(X<=x, 1.5<Y<=2.5) = P(X<=x, Y=2)
plotCDF1(mle, margin=1, bound="b", int=c(1.5,2.5), col="black", add=TRUE)

# Add legend
legend(0,1,c("k=1: operative","k=2: natural"), col=c("red","black"), lty=1)

# Plot marginal distribution of the failure cause Y
plotCDF1(mle, margin=2, bound="u", col="black", xlim=c(0,3), 
 xlab="y", main="P(Y<=y)")
}
\keyword{datasets}
\concept{censored data}
\concept{competing risk}

