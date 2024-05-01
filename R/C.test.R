#' Cochran's C test of homogeneity of variances
#'
#' Performs a Cochran's test of the null hypothesis that the largest variance in several sampled variances are the same.
#'
#' @param object an object of class "\code{\link{lm}}", containing the specified design.
#' @details The test statistic is a ratio that relates the largest variance to the sum of the sampled variances.
#' @return A list of class htest containing the following components:
#' @return \item{statistic}{Cochran's C test statistic}
#' @return \item{p-value}{The p-value of the test}
#' @return \item{alternative}{A character string describing the alternative hypothesis}
#' @return \item{method}{The character string Cochran test of homogeneity of variances}
#' @return \item{data.name}{A character string giving the name of the lm object}
#' @return \item{estimate}{Sample estimates of variances}
#' @author Leonardo Sandrini-Neto (\email{leonardosandrini@@ufpr.br})
#' @seealso \code{\link{gad}}
#' @examples
#' library(GAD)
#' data(rohlf95)
#' CG <- as.fixed(rohlf95$cages)
#' MQ <- as.random(rohlf95$mosquito)
#' model <- lm(wing ~ CG + MQ%in%CG, data = rohlf95)
#' C.test(model)
#' @importFrom stats pf
#' @export
C.test <-
function(object){
  model <- deparse(substitute(object))
  by.factor <- as.factor(1:object$rank)
  n <- length(object$model[,1])/object$rank
  k <- object$rank
  var <- tapply(object$model[,1], rep(1:k, each = n), var)
  int <- interaction(object$model[,-1], lex.order = TRUE)
  f.int <- factor(int, levels = unique(int))
  names(var) <- levels(f.int)
  mean <- tapply(object$model[,1], rep(1:k, each = n), mean)
  C <- max(var)/sum(var)
  group <- names(var)[which(var == max(var))]
  method <- "Cochran test of homogeneity of variances"
  alt <- paste("Group", group, "has outlying variance")
  f <- (1/C - 1)/(k - 1)
  p <- 1 - pf(f, (n - 1) * (k - 1), (n - 1)) * k
  pval <- 1 - p
  pval[pval < 0] <- 0
  pval[pval > 1] <- 1
  result <- list(statistic = c(C = C), parameter = c(n = n, k = k), alternative = alt, p.value = pval, method = method, estimate = round(var, 4), mean = mean, var = var, data.names = model)
  class(result) <- "htest"
  return(result)
}
