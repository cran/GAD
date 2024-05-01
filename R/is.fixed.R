#' Check if a factor is fixed
#'
#' This function works the same way as \code{\link{is.factor}}.
#'
#' @param x a vector of data, usually a nominal variable encoded as a \code{"fixed factor"} using \code{\link{as.fixed}}.
#' @return Function \code{is.fixed} returns \code{TRUE} or \code{FALSE} depending on whether its argument is a fixed factor or not.
#' @author Leonardo Sandrini-Neto (\email{leonardosandrini@@ufpr.br})
#' @seealso \code{\link{is.random}}, \code{\link{as.fixed}}
#' @examples
#' library(GAD)
#' data(rohlf95)
#' CG <- as.fixed(rohlf95$cages)
#' MQ <- as.random(rohlf95$mosquito)
#' is.fixed(CG)
#' is.random(MQ)
#' @export
is.fixed <-
function(x){
  inherits(x, "fixed")
}
