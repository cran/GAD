\name{is.fixed}
\alias{is.fixed}
\title{Tests if a factor is fixed}
\description{
This function works the same way of \emph{is.factor}.
}
\usage{is.fixed(x)}
\arguments{
  \item{x}{a vector of data}
}
\value{
is.fixed returns "TRUE" or "FALSE" depending on whether its argument is a fixed factor or not}
\author{
Leonardo Sandrini-Neto (\email{leonardosandrini@gmail.com})
}
\seealso{
\code{\link{is.random}}
}
\examples{
library(GAD)
data(rohlf95)
CG <- as.fixed(rohlf95$cages)
MQ <- as.random(rohlf95$mosquito) 
is.fixed(CG)
is.random(MQ)
}
\keyword{design}

