\name{as.random}
\alias{as.random}
\title{Encodes a vector as a "random factor"}
\description{
Assigns a class "random" to a vector
}
\usage{
as.random(x)
}
\arguments{
  \item{x}{a vector of data}
}
\details{
The function works the same way as \emph{as.factor}, but assigns an additional class informing that it is a random factor.
}
\value{
as.factor returns an object of class "factor" and "random". 
}
\author{
Leonardo Sandrini-Neto (\email{leonardosandrini@gmail.com})
}
\seealso{
\code{\link{as.fixed}}
}
\examples{
library(GAD)
data(rohlf95)
CG <- as.fixed(rohlf95$cages)
MQ <-as.random(rohlf95$mosquito) 
}
\keyword{design}

