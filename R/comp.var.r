#' Components of variation
#'
#' This function calculates components of variation for fixed and random factors.
#'
#' @param object an object of class "\code{\link{lm}}", containing the specified design with random and/or fixed factors.
#' @param anova.tab an object containing the results returned by \code{gad} or \code{pooling} functions. See examples below.
#' @details {
#' This function calculates components of variation for any combination of orthogonal/nested and fixed/random factors. For pooled terms,
#' \code{comp.var} function seeks for the denominator of the removed term and keep its type (fixed or random). \cr
#'
#' Note that there are differences on the interpretation of the results between fixed and random factors. For fixed factors, the components of
#' variation calculated are just the sum of squares divided by the degrees of freedom and the hypothesis only concern those levels that were
#' included in the model. On the other hand, for random factors the components of variation calculated are truly variance components and are
#' interpretable as measures of variability of a population of levels, which are randomly sampled. However, for most of studies that aims to
#' compare the amount of variation attributed to each term, the estimates of components of variation for both fixed and random factors are important
#' and directly comparable (Anderson et al., 2008). \cr
#'
#' Eventually, the estimates of components of variation for some terms in the model can be negative. The term in question generally has a large p-value
#' and its estimative is usually set to zero, since a negative estimate is illogical. An alternative to this problem is remove the term using
#' \code{\link{pooling}} function and re-analyse the model (Fletcher and Underwood, 2002).
#' }
#' @return {
#' A list of length 2, containing the mean squares estimates (\code{$mse}) and the table of components of variation (\code{$comp.var}) for the model. \cr
#'
#' The \code{$comp.var} table contains four columns: “Type” shows whether terms are fixed or random; “Estimate” shows the estimate of component of variation
#' (in squared units); “Square root” shows the square root of the “Estimate” column (original unit); and “Percentage” shows the percentage of variability
#' attributable to each term (both fixed and random factors).
#' }
#' @author Eliandro Gilbert (\email{eliandrogilbert@@gmail.com})
#' @references Anderson, M.J., Gorley, R.N., Clarke, K.R. 2008. \emph{PERMANOVA+ for PRIMER}: Guide to Software and Statistical Methods. PRIMER-E: Plymouth, UK.
#' @references Fletcher, D.J., Underwood, A.J. 2002. How to cope with negative estimates of components of variance in ecological field studies. Journal of Experimental Marine Biology and Ecology 273, 89-95.
#' @seealso \code{\link{gad}}, \code{\link{pooling}}, \code{\link{estimates}}
#' @examples
#' library(GAD)
#' data(crabs)
#' head(crabs)
#'
#' Re <- as.random(crabs$Region)   # random factor
#' Lo <- as.random(crabs$Location) # random factor nested in Region
#' Si <- as.random(crabs$Site)     # random factor nested in Location
#'
#' model <- lm(Density ~ Re + Lo%in%Re + Si%in%Lo%in%Re, data = crabs)
#'
#' C.test(model)                     # Checking homogeneity of variances
#' estimates(model)                  # Checking model structure
#' model.tab <- gad(model)           # Storing the result of ANOVA on a new object
#' model.tab                         # Checking ANOVA results
#' comp.var(model, anova.tab = model.tab) # Calculating the components of variations
#' @export
comp.var <-
function(object, anova.tab = NULL){
  tab <- anova.tab
  if (is.null(tab)){
    stop("Data cannot be null")
  }
  if(is.null(tab$anova)){
    stop("Something is wrong here. The argument anova.tab argument must store the result of gad or pooling functions")
  }
  f.versus <- tab$f.versus
  anova.res <- tab$anova
  quasi.f <- 0
  if(is.null(f.versus$Num.Df)){
    quasi.f <- FALSE
  }
  if(!is.null(f.versus$Num.Df)){
    quasi.f <- TRUE
  }
  mse <- tab$mse
  if(is.null(tab$mse)){
    a <- estimates(object, quasi.f = quasi.f)
    mse <- a$mse
    f.versus <- a$f.versus
  }
  table <- anova.res[, 1:4]
  cvar <- data.frame(table)
  cvar[,] <- NA
  names(cvar) <- c("Type", "Estimate", "Square root", "Percentage")
  nomes <- matrix(ncol=2,nrow=nrow(f.versus))
  nomes[, 2] <- rownames(f.versus)
  mod <- gsub(pattern = "\\*", replacement = " + ",  x = mse)
  mod <- strsplit(mod,  " \\+ ")
  mse2 <- list(1)
  mse21 <- list(1)
  for(i in 1:nrow(mse)){
    mse2[[i]] <- mod[[i]][seq(from = 2, to = length(mod[[i]]),by = 2)]
    mse21[[i]] <- mod[[i]][seq(from = 1, to = length(mod[[i]]),by = 2)]
    }
  for(i in 1:nrow(nomes)){
    for(j in 1:length(mse2[[i]])){
      if(mse2[[i]][j] == nomes[i,2]){
        nomes[i,1] <- mse21[[i]][j]
        }
      }
    }
  mse.final <- matrix(ncol = 1, nrow = length(mse2))
  colnames(mse.final) <- "Mean square estimates"
  rownames(mse.final) <- rownames(mse)
  for(i in 1:nrow(mse.final)) {
    a <- mse2[[i]]
    for(j in 1:length(mse2[[i]])){
      a[j] <- paste(mse21[[i]][j], "*", mse2[[i]][j], sep = "")
      }
    mse.final[i, 1] <- paste(a, collapse = " + ")
    }
  cvar[nrow(cvar), 2] <- as.numeric(table$Mean[nrow(table)])
  nome <- rownames(cvar)[nrow(cvar)]
  for(i in nrow(f.versus):1){
    teste <- strsplit(f.versus$Denominator[i], " \\+ ")
    if(length(teste[[1]]) > 1){
      if (quasi.f == TRUE){
        num <- strsplit(f.versus$Numerator[i], " \\+ ")
        den <- strsplit(f.versus$Denominator[i], " \\+ ")
        ms.num <- ((table$Mean[which(rownames(table) == num[[1]][1])])  + (table$Mean[which(rownames(table) == num[[1]][2])]))
        ms.den <- ((table$Mean[which(rownames(table) == den[[1]][1])])  + (table$Mean[which(rownames(table) == den[[1]][2])]))
        cvar[i,2] <- as.numeric(((ms.num)- (ms.den))/ as.numeric(nomes[i,1]))
        } else {
        cvar[i,2] <- NA
        }
      }
    if(f.versus$Denominator[i] == "No test"){
      cvar[i,2] <- NA
      }
    if(f.versus$Denominator[i] == "Residuals" | f.versus$Denominator[i] == nome) {
      cvar[i,2] <- as.numeric(((table$Mean[i]) - (table$Mean[nrow(table)]))/ as.numeric(nomes[i,1]))
      }
    if((f.versus$Denominator[i] != "Residuals") & (length(teste[[1]]) == 1) & (f.versus$Denominator[i] != "No test") & (f.versus$Denominator[i] != nome)){
      cvar[i,2] <- as.numeric(((table$Mean[i]) - (table$Mean[which(rownames(table) == f.versus$Denominator[i])]))/ as.numeric(nomes[i,1]))
      }
    }
  raiz <- vector(length = nrow(cvar))
  for(i in 1:nrow(cvar)){
    if(is.na(cvar[i,2])){
      cvar[i,3] <- NA
      raiz[i] <- 0
      next
      }
    if(cvar[i,2] < 0){
      cvar[i,3] <- NA
      raiz[i] <- 0
      next
      }
    if(!is.na(cvar[i,2]) & cvar[i,2] > 0){
      cvar[i,3] <- cvar[i,2]^0.5
      raiz[i] <- cvar[i,2]^0.5
      }
    }
  soma <- sum(raiz)
  for(i in 1:nrow(cvar)){
    cvar[i,4] <- ((raiz[i]*100)/soma)
  }
  cvar[,2:4]<-round(cvar[,2:4],digits = 3)
  tab.est<-estimates(object)
  for(i in 1:nrow(cvar)){
    if (rownames(cvar)[i] == "Residuals"){
      cvar[i,1] <- "Random"
      next
    }
    termo <- rownames(cvar)[i]
    teste <- strsplit(termo, "\\:")
    a <- numeric(length = 1)
    tipo <- character(length = length(teste[[1]]))
    if(termo == "Pooled1" | termo == "Pooled2" | termo == "Pooled3" |
       termo == "Pooled4" | termo == "Pooled5" | termo == "Pooled6" |
       termo == "Pooled7" | termo == "Pooled8" | termo == "Pooled9" |
       termo == "Pooled10"){
      ptab <- tab$pool.table
      termo <- ptab[(rownames(cvar)[i]),2]
      while (all(rownames(tab.est$tm) != termo) & termo != "Residuals"){
        termo <- ptab[termo,2]
      }
      if(termo == "Residuals"){
        cvar[i,1] <- "Random"
        next
      }
      teste <- strsplit(termo, "\\:")
    }
    for(j in 1:length(teste[[1]])){
      a <- tab.est$tm[which(rownames(tab.est$tm) == termo) , which(colnames(tab.est$tm) == teste[[1]][j])]
      if(a > 0){
        tipo[j] <- "Random"
      } else {
        tipo[j]<-"Fixed"
      }
    }
    if(any(tipo == "Random")){
      cvar[i,1] <- "Random"
    }
    if(all(tipo == "Fixed")){
      cvar[i,1] <- "Fixed"
    }
  }
  for(i in 1:nrow(cvar)){
    for(j in 1:ncol(cvar)){
      if(is.na(cvar[i,j])){
        cvar[i,j] <- ""
      }
    }
  }
  aviso <- NA
  if(any(cvar[,2] < 0)){
    aviso <- "There is one or more negative estimates of components of variation in your model. Perhaps you should consider pooling these terms"
  }
  if(is.na(aviso)){
    comp.var <- list(mse = mse.final, comp.var = cvar)
  } else {
    comp.var <- list(mse = mse.final, comp.var = cvar, note = aviso)
  }
  return(comp.var)
}
