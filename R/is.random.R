#' Check if a factor is random
#'
#' This function works the same way as \code{\link{is.factor}}.
#'
#' @param x a vector of data, usually a nominal variable encoded as a \code{"random factor"} using \code{\link{as.random}}.
#' @return Function \code{is.random} returns \code{TRUE} or \code{FALSE} depending on whether its argument is a random factor or not.
#' @author Leonardo Sandrini-Neto (\email{leonardosandrini@@ufpr.br})
#' @seealso \code{\link{is.fixed}}, \code{\link{as.random}}
#' @examples
#' library(GAD)
#' data(rohlf95)
#' CG <- as.fixed(rohlf95$cages)
#' MQ <- as.random(rohlf95$mosquito)
#' is.fixed(CG)
#' is.random(MQ)
#' @export
is.random <-
function(x){
  inherits(x, "random")
}
