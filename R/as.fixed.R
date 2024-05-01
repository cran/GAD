#' Encodes a vector as a "fixed factor"
#'
#' Assigns a class \code{"fixed"} to a vector
#'
#' @param x a vector of data, usually a nominal variable.
#' @details The function works the same way as \code{\link{as.factor}}, but assigns an additional class informing that it is a fixed factor.
#' @return Function \code{as.fixed} returns an object of class \code{"factor"} and \code{"fixed"}.
#' @author Leonardo Sandrini-Neto (\email{leonardosandrini@@ufpr.br})
#' @seealso \code{\link{as.random}}
#' @examples
#' library(GAD)
#' data(rohlf95)
#' CG <- as.fixed(rohlf95$cages)
#' MQ <- as.random(rohlf95$mosquito)
#' class(CG)
#' class(MQ)
#' @export
as.fixed <-
function(x){
  f <- factor(x)
  class(f) <- c("factor", "fixed")
  f
}
