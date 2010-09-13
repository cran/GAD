\name{snk.test}
\alias{snk.test}
\title{Student-Newman-Keuls (SNK) procedure}
\description{
This function perforns a SNK post-hoc test of means on the factors of a chosen term of the model, comparing among levels of one factor within each level of other factor or combination of factors.
}
\usage{snk.test(object, term, among = NULL, within = NULL)}
\arguments{
  \item{object}{An object of class lm, containing the specified design.}
  \item{term}{Term of the model to be analysed. Use \code{\link{estimates}} to see the right form to inform it.}
  \item{among}{Specifies the factor which levels will be compared among. Need to be specified if the term to be analysed envolves more than one factor.}  
  \item{within}{Specifies the factor or combination of factors that will be compared within level among.}
}
\details{SNK is a stepwise procedure for hypothesis testing. First the sample means are sorted, then the pairwise studentized range (q) is calculated by dividing the  differences between means by the standard error, which is based upon the average variance of the two sample.}
\value{
  A list containing the standard error, the degree of freedom and pairwise comparisons among levels of one factor within each level of other(s) factor(s).
}
\references{
Underwood, A.J. 1997. \emph{Experiments in Ecology: Their Logical Design and Interpretation Using Analysis of Variance.} Cambridge University Press, Cambridge.
}
\author{
  Mauricio G. Camargo (\email{camargo.ufpr@gmail.com})\cr
  Leonardo Sandrini-Neto (\email{leonardosandrini@gmail.com})
}
\seealso{
\code{\link{gad}}, \code{\link{estimates}}
}
\examples{
library(GAD)
data(rohlf95)
CG <- as.fixed(rohlf95$cages)
MQ <-as.random(rohlf95$mosquito) 
model <- lm(wing ~ CG + CG\%in\%MQ, data = rohlf95)
gad(model)
##Check estimates to see model structure
estimates(model)
snk.test(model,term='CG:MQ',among='MQ',within='CG')
##
##
##Example using chapman86
data(chapman86)
O <- as.random(chapman86$origin)
S <- as.random(chapman86$shore)
B <- as.random(chapman86$boulder)
C <- as.random(chapman86$cage)
model <- lm(growth ~ O + S + O*S + B\%in\%S + O*(B\%in\%S) + C\%in\%(O*(B\%in\%S)), 
            data = chapman86)
gad(model)
##Check estimates to see model structure
estimates(model)
snk.test(model,term='O')
snk.test(model,term='O:S',among='S',within='O')
#if term O:S:B were significant, we could try 
snk.test(model,term='O:S:B',among='B',within='O:S')
}
\keyword{htest}

