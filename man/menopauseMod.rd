\name{menopauseMod}
\docType{data}
\alias{menopauseMod}
\title{Modified menopause data}
\description{
  Example data set with interval censored data and competing risks. 
  The data come from Cycle I of the 
  Health Examination Survey of the National Center for Health Statistics, 
  and contain information on the menopausal status of 2423 women 
  (MacMahon and Worcestor, 1966). 

  The format of the data has been modified 
  to allow for easy plotting with the functions \code{\link{plotHM}} and 
  \code{\link{plotDens2}} (see section 'Format').  
}
\usage{data(menopauseMod)}
\format{A matrix containing 2423 rows and 4 columns. Each row    
(x1,x2,y1,y2) corresponds to a subject in the study. The interval (x1,x2] 
contains the unobservable age of menopause X. The interval [y1,y2]
contains the type of menopause Y, where Y=1 represents operative
menopause and Y=2 represents natural menopause. We use the value 100 
to represent infinity. 

In order to allow easy plotting with \code{\link{plotHM}} and 
\code{\link{plotDens2}}, 
the y-intervals were modified as follows: [y1,y2] was changed into
[y1-0.25, y2+0.25]. 
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
  status data on the time of menopause with two competing 
  risks: operative menopause and natural menopause. 
}
\source{MacMahon and Worcestor (1966). Age at menopause, United States 
1960 - 1962. \emph{National Center for Health Statistics. Vital and Health 
Statistics}, volume 11, number 19. 
}
\references{
Krailo and Pike (1983). Estimation of the distribution of age at 
natural menopause from prevalence data. \emph{American Journal of 
Epidemiology} \bold{117} 356-361.
}
\seealso{\code{\link{menopause}}}
\examples{
# Load the data
data(menopauseMod)

# Compute the MLE
mle <- computeMLE(menopauseMod)

# Create density plot
par(mfrow=c(1,1))
plotDens2(mle, xlim=c(0,100), border="black", xlab="age in years", 
 ylab="cause of menopause (1=operative, 2=natural)", 
 main="Density plot of the MLE for the menopause data")
}
\keyword{datasets}
\concept{censored data}
\concept{competing risk}
