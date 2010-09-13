\name{GAD-package}
\alias{GAD-package}
\docType{package}
\title{Analysis of variance from general principles}
\description{
This package contains functions for the analysis of any complex ANOVA models with any combination of orthogonal/nested and fixed/random factors, based on general principles described by Underwood (1997).}
\details{
\tabular{ll}{
Package: \tab GAD\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2010-09-14\cr
License: \tab GPL (>= 2.0)\cr
}
}
\author{
Leonardo Sandrini-Neto (\email{leonardosandrini@gmail.com})\cr
Mauricio G. Camargo (\email{camargo.ufpr@gmail.com})\cr
\cr
Maintainer: leonardosandrini@gmail.com
}
\references{
Cornfield, J. and Tukey, J.W. (1956). Average values of mean squares in factorials. \emph{Annals of Mathematical Statistics}, \bold{27}, 907-949.\cr
\cr
Sokal, R. R. and F. J. Rohlf. 1995. \emph{Biometry: the principles and practice of statistics in biological research.}  3rd edition. W. H. Freeman and Co. New York. 887 pp.\cr 
\cr
Underwood, A.J., 1997. \emph{Experiments in Ecology: Their Logical Design and Interpretation Using Analysis of Variance.} Cambridge University Press, Cambridge.\cr
}
\keyword{ package }
\seealso{
\code{\link{gad}}, \code{\link{estimates}}, \code{\link{C.test}}, \code{\link{snk.test}}
}
\examples{
#Example 1
library(GAD)
data(rohlf95)
CG <- as.fixed(rohlf95$cages)
MQ <-as.random(rohlf95$mosquito) 
model <- lm(wing ~ CG + CG\%in\%MQ, data = rohlf95)
gad(model)
##
##
#Example 2
data(rats)
names(rats)
TR <- as.fixed(rats$treat)
RA <- as.random(rats$rat)
LI <- as.random(rats$liver)
model <- lm(glycog ~ TR + RA\%in\%TR + LI\%in\%RA\%in\%TR, data=rats)
gad(model)
##
##
#Example 3
data(chapman86)
O <- as.random(chapman86$origin)
S <- as.random(chapman86$shore)
B <- as.random(chapman86$boulder)
C <- as.random(chapman86$cage)
model <- lm(growth ~ O + S + O*S + B\%in\%S + O*(B\%in\%S) + C\%in\%(O*(B\%in\%S)), data = chapman86)
gad(model)
}
\keyword{package}
\keyword{design}

