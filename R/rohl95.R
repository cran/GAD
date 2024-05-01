#' Mosquitos' wing data colleted by Rohlf and cited in Sokal & Rohlf (1995)
#'
#' Three different types of cage are tested on the growth of \emph{Aedes intrudens}, a kind of
#' mosquito pupae. In each one, four mosquitos are added and its wings measured twice. There are
#' 24 observations (3 cages x 4 jars x 2 measures).
#'
#' \itemize{
#'   \item cages: a fixed factor with three levels (\code{cage 1, cage2, cage3})
#'   \item mosquito: a random factor with four levels (\code{m1, m2, m3, m4}) nested in cages
#'   \item measure: sample size
#'   \item wing: response variable
#' }
#'
#' @docType data
#' @keywords datasets
#' @name rohlf95
#' @usage data(rohlf95)
#' @format A data frame with 24 rows and 4 variables
#' @references Sokal, R.R., Rohlf, F.J. 1995. \emph{Biometry}: The Principles and Practice of Statistics in Biological Research. 3rd Edition. W.H. Freeman and Co., New York.
#' @examples
#' library(GAD)
#' data(rohlf95)
#' rohlf95
NULL
