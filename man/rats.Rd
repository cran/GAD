\name{rats}
\alias{rats}
\docType{data}
\title{
Dataset - Glycogen content of rat livers.
}
\description{
Duplicate readings were made on each of three preparations of rat livers from each of two rats for three different treatments (Sokal & Rohlf, 1995).
}
\usage{data(rats)}
\format{
  A data frame with 36 observations on the following 4 variables.
  \describe{
    \item{\code{treat}}{a fixed factor}
    \item{\code{rat}}{a random factor nested in treat}
    \item{\code{liver}}{sample size}
    \item{\code{glycog}}{response variable}
  }
}
\references{
Sokal, R. R. and F. J. Rohlf. 1995. \emph{Biometry: the principles and practice of statistics in biological research.}  3rd edition. W. H. Freeman and Co. New York. 887 pp. 
}
\examples{
data(rats)
}
\keyword{datasets}
