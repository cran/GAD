\name{gad}
\alias{gad}
\title{General ANOVA Design}
\description{Fits a general ANOVA design with any combination of orthogonal/nested and fixed/random factors through function \code{\link{estimates}}}
\usage{gad(object)}
\arguments{
  \item{object}{an object of class lm, containing the specified design with random and/or fixed factors}
}
\details{
Function \emph{gad} returns an analysis of variance table using the \code{\link{estimates}} function to identify the appropriate F-ratios and consequently p-values for any complex model of orthogonal or nested, fixed or random factors as described by Underwood(1997).}
\value{An object of class "anova" inheriting from class "data.frame"}
\references{
Underwood, A.J. 1997. \emph{Experiments in Ecology: Their Logical Design and Interpretation Using Analysis of Variance.} Cambridge University Press, Cambridge.
}
\author{
Leonardo Sandrini-Neto (\email{leonardosandrini@gmail.com})
}
\seealso{\code{\link{estimates}}}
\examples{
#Example 1
library(GAD)
data(rohlf95)
CG <- as.fixed(rohlf95$cages)
MQ <- as.random(rohlf95$mosquito) 
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
data(snails)
O <- as.random(snails$origin)
S <- as.random(snails$shore)
B <- as.random(snails$boulder)
C <- as.random(snails$cage)
model <- lm(growth ~ O + S + O*S + B\%in\%S + O*(B\%in\%S) + C\%in\%(O*(B\%in\%S)), 
            data = snails)
gad(model)
}
\keyword{design}

