\name{snails}
\alias{snails}
\docType{data}
\title{
Dataset - Growth rates of snails on large boulders on different rock shores.
}
\description{
This design was extracted from Underwood (1997) but data are artifitial. Snails were transplanted from origin to different shores. Several boulders were used on each shore. Cages with snail of each origin on each boulder were replicated. All factors (origin, shore, boulder and cage) are random.
}
\usage{data(snails)}
\format{
  A data frame with 240 observations on the following 6 variables.
  \describe{
    \item{\code{origin}}{a random factor with levels \code{O1} \code{O2}}
    \item{\code{shore}}{a random factor with levels \code{S1} \code{S2} \code{S3} \code{S4} orthogonal to origin}
    \item{\code{boulder}}{a random factor with levels \code{B1} \code{B2} \code{B3} nested in shore}
    \item{\code{cage}}{a random factor with levels \code{C1} \code{C2} nested in the combination of boulder and origin}
    \item{\code{replicate}}{Sample size}
    \item{\code{growth}}{response variable}
  }
}
\references{
Underwood, A.J. 1997. \emph{Experiments in Ecology: Their Logical Design and Interpretation Using Analysis of Variance.} Cambridge University Press, Cambridge.
}
\examples{
data(snails)
}
\keyword{datasets}
