\name{as.fixed}
\alias{as.fixed}
\title{Encodes a vector as a "fixed factor"}
\description{
Assigns a class "fixed" to a vector
}
\usage{
as.fixed(x)
}
\arguments{
  \item{x}{a vector of data}
}
\details{
The function works the same way as \emph{as.factor}, but assigns an additional class informing that it is a fixed factor.
}
\value{
as.factor returns an object of class "factor" and "fixed". 
}
\author{
Leonardo Sandrini-Neto (\email{leonardosandrini@gmail.com})
}
\seealso{
\code{\link{as.random}}
}
\examples{
library(GAD)
data(rohlf95)
CG <- as.fixed(rohlf95$cages)
MQ <-as.random(rohlf95$mosquito) 
}
\keyword{design}

