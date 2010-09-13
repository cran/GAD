\name{estimates}
\alias{estimates}
\title{Estimates of an ANOVA design}
\description{
This function is used to construct the mean squares estimates of an ANOVA design, considering the complications imposed by nested/orthogonal and fixed/random factors.
}
\usage{
estimates(object)
}
\arguments{
  \item{object}{an object of class lm, containing the specified design with random and/or fixed factors}
}
\details{
Determines what each mean square estimates in an ANOVA design by a set of procedures originally described by Cornfield and Tukey (1956). This version is a modification proposed by Underwood (1997), which does not allow for the use of fixed nested factors. The steps involve the construction of a table of multipliers with a row for each source of variation and a column for each term in the model that is not an interaction. The mean square estimates for each source of variation is obtained by determining which components belong to each mean square and what is their magnitude. This enables the recognition of appropriate F-ratios.
}
\value{
A list of length 3, containing the table of multipliers ($tm), the mean squares estimates ($mse) and the F-ratio versus ($f.versus) for the model.
}
\references{
Cornfield, J. and Tukey, J.W. (1956). Average values of mean squares in factorials. \emph{Annals of Mathematical Statistics}, \bold{27}, 907-949.\cr
\cr
Underwood, A.J., 1997. \emph{Experiments in Ecology: Their Logical Design and Interpretation Using Analysis of Variance.} Cambridge University Press, Cambridge.
}
\author{
  Leonardo Sandrini-Neto (\email{leonardosandrini@gmail.com})
}
\seealso{
\code{\link{gad}}, \code{\link{estimates}}
}
\examples{
#Example 1
library(GAD)
data(rohlf95)
CG <- as.fixed(rohlf95$cages)
MQ <-as.random(rohlf95$mosquito) 
model_R <- lm(wing ~ CG + CG\%in\%MQ, data = rohlf95)
estimates(model_R)
##
##
#Example 2
data(chapman86)
O <- as.random(chapman86$origin)
S <- as.random(chapman86$shore)
B <- as.random(chapman86$boulder)
C <- as.random(chapman86$cage)
model_C <- lm(growth ~ O + S + O*S + B\%in\%S + O*(B\%in\%S) + C\%in\%(O*(B\%in\%S)),
              data = chapman86)
estimates(model_C)
}
\keyword{design}

