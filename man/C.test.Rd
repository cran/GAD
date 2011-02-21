\name{C.test}
\alias{C.test}
\title{Cochran's C test of homogeneity of variances}
\description{
Performs a Cochran's test of the null hypothesis that the largest variance in several sampled variances are the same.
}
\usage{
C.test(object)
}
\arguments{
  \item{object}{An object of class lm, containing the specified design.}
}
\details{
The test statistic is a ratio that relates the largest variance to the sum of the sampled variances.
}
\value{
A list of class htest containing the following components:
\item{statistic}{Cochran's C test statistic}
\item{p-value}{The p-value of the test}
\item{alternative}{A character string describing the alternative hypothesis}
\item{method}{The character string Cochran test of homogeneity of variances}
\item{data.name}{A character string giving the name of the lm object}
\item{estimate}{Sample estimates of variances}
}
\author{
Leonardo Sandrini-Neto (\email{leonardosandrini@gmail.com})
}
\seealso{\code{\link{gad}}}
\examples{
  library(GAD)
  data(rohlf95)
  cg <- as.fixed(rohlf95$cages)
  mq <- as.random(rohlf95$mosquito) 
  model <- lm(wing ~ cg + mq\%in\%cg, data = rohlf95)
  C.test(model)
}
\keyword{htest}

