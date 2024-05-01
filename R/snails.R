#' Growth rates of snails on large boulders on different rock shores
#'
#' This design was extracted from Underwood (1997), but data are artificial.
#' Snails were transplanted from origin to different shores. Several boulders
#' were used on each shore. Cages with snails of each origin on each boulder
#' were replicated. All factors (origin, shore, boulder and cage) are random.
#'
#' \itemize{
#'   \item origin: a random factor with two levels (\code{O1, O2})
#'   \item shore: a random factor with four levels (\code{S1, S2, S3, S4}) orthogonal to origin
#'   \item boulder: a random factor with three levels (\code{B1, B2, B3}) nested in shore
#'   \item cage: a random factor with two levels (\code{C1, C2}) nested in the combination of boulder and origin
#'   \item replicate: sample size
#'   \item growth: response variable
#' }
#'
#' @docType data
#' @keywords datasets
#' @name snails
#' @usage data(snails)
#' @format A data frame with 240 rows and 6 variables
#' @references Underwood, A.J. 1997. \emph{Experiments in Ecology}: Their Logical Design and Interpretation Using Analysis of Variance. Cambridge University Press, Cambridge.
#' @examples
#' library(GAD)
#' data(snails)
#' snails
NULL
