\name{rohlf95}
\alias{rohlf95}
\docType{data}
\title{
Dataset - Mosquitos' wing data colleted by Rohlf and cited in Sokal & Rohlf (1995)
}
\description{
Three different types of cage are tested on the growth of \emph{Aedes intrudens}, a kind of mosquito pupae. In each one, four mosquitos are added and its wings measured twice. There are 24 observations (3 cages X 4 jars X 2 measures).
}
\usage{data(rohlf95)}
\format{
  A data frame with 24 observations on the following 4 variables.
  \describe{
    \item{\code{cages}}{a fixed factor with levels \code{cage1} \code{cage2} \code{cage3}}
    \item{\code{mosquito}}{a random factor with levels \code{m1} \code{m2} \code{m3} \code{m4} nested in cages}
    \item{\code{measure}}{sample size}
    \item{\code{wing}}{response variable}
  }
}
\references{
Sokal, R.R., Rohlf, F.J. 1995. \emph{Biometry: the principles and practice of statistics in biological research.}  3rd edition. W. H. Freeman and Co. New York. 887 pp. 
}
\examples{
data(rohlf95)
}
\keyword{datasets}
