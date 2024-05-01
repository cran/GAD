#' Estimates of an analysis of variance design
#'
#' This function is used to construct the mean squares estimates of an ANOVA design, considering the complications imposed by nested/orthogonal and fixed/random factors.
#'
#' @param object an object of class "\code{\link{lm}}", containing the specified design with random and/or fixed factors.
#' @param quasi.f logical, indicating whether to use quasi F-ratio when there is no single error term appropriate in the analysis. Default to \code{FALSE}.
#' @details Determines what each mean square estimates in an ANOVA design by a set of procedures originally described by Cornfield and Tukey (1956). This version is a modification proposed by Underwood (1997), which does not allow for the use of fixed nested factors. The steps involve the construction of a table of multipliers with a row for each source of variation and a column for each term in the model that is not an interaction. The mean square estimates for each source of variation is obtained by determining which components belong to each mean square and what is their magnitude. This enables the recognition of appropriate F-ratios.
#' @return A list containing the table of multipliers (\code{$tm}), the mean squares estimates (\code{$mse}) and the F-ratio versus (\code{$f.versus}) for the model.
#' @author Leonardo Sandrini-Neto (\email{leonardosandrini@@ufpr.br})
#' @author Eliandro Gilbert (\email{eliandrogilbert@@gmail.com})
#' @references Cornfield, J., Tukey, J.W. 1956. Average values of mean squares in factorials. Annals of Mathematical Statistics 27, 907-949.
#' @references Underwood, A.J. 1997. \emph{Experiments in Ecology}: Their Logical Design and Interpretation Using Analysis of Variance. Cambridge University Press, Cambridge.
#' @seealso \code{\link{gad}}
#' @examples
#' # Example 1
#' library(GAD)
#' data(rohlf95)
#' CG <- as.fixed(rohlf95$cages)
#' MQ <- as.random(rohlf95$mosquito)
#' model <- lm(wing ~ CG + MQ%in%CG, data = rohlf95)
#' estimates(model)
#'
#' # Example 2
#' data(rats)
#' names(rats)
#' TR <- as.fixed(rats$treat)
#' RA <- as.random(rats$rat)
#' LI <- as.random(rats$liver)
#' model2 <- lm(glycog ~ TR + RA%in%TR + LI%in%RA%in%TR, data = rats)
#' estimates(model2)
#'
#' # Example 3
#' data(snails)
#' O <- as.random(snails$origin)
#' S <- as.random(snails$shore)
#' B <- as.random(snails$boulder)
#' C <- as.random(snails$cage)
#' model3 <- lm(growth ~ O + S + O*S + B%in%S + O*(B%in%S) + C%in%(O*(B%in%S)), data = snails)
#' estimates(model3) # 'no test' for shore
#' estimates(model3, quasi.f = TRUE) # suitable test for shore
#' @import matrixStats
#' @importFrom stats var
#' @export
estimates <-
  function (object, quasi.f = FALSE){
    balance <- var(as.vector(table(object$model[, 2:(length(object$x) + 1)])))
    if (balance > 0) {
      stop("Design unbalanced! This function can only handle balanced designs.\n")
    }
    tm <- attr(object$terms, "factors")
    tm <- tm[-1, , drop = FALSE]
    tm <- t(tm)
    tm2 <- tm
    for (i in 1:length(rownames(tm))){
      rownames(tm2)[i]<-as.character(x=i)
    }
    nomes <- matrix(ncol = 3,nrow = nrow(tm))
    nomes[,1] <- rownames(tm2)
    nomes[,2] <- rownames(tm)
    tm.class <- tm
    for (i in 2:length(object$model)) {
      colnames(tm.class)[i - 1] <- class(object$model[, i])[2]
    }
    nest.fixed <- subset(tm.class, rowMaxs(tm.class) > 1)
    if (length(nest.fixed) != 0) {
      nest.logic <- matrix(nrow = nrow(nest.fixed), ncol = 1)
      for (i in 1:nrow(nest.fixed)) {
        nest.logic[i] <- length(which(nest.fixed[i, ] == 1)) == 1
      }
      nest.fixed <- nest.fixed[nest.logic, , drop = FALSE]
      for (i in 1:nrow(nest.fixed)) {
        if (names(subset(nest.fixed[i, ], nest.fixed[i, ] == 1)) == "fixed") {
          stop("This function does not allow for the use of fixed nested factors")
        }
      }
    }
    factors <- object$model[, -1]
    f.levels <- numeric(ncol(tm))
    if (length(f.levels) == 1) {
      f.levels <- nlevels(factors)
    }
    else {
      for (j in 1:ncol(tm)) {
        f.levels[j] = nlevels(factors[, j])
      }
    }
    tm.final <- matrix(nrow = nrow(tm), ncol = ncol(tm))
    rownames(tm.final) <- rownames(tm)
    colnames(tm.final) <- colnames(tm)
    for (j in 1:ncol(tm.class)) {
      for (i in 1:nrow(tm.class)) {
        if (tm.class[i, j] == 1) {
          if (colnames(tm.class)[j] == "fixed")
            tm.final[i, j] = 0
          if (colnames(tm.class)[j] == "random")
            tm.final[i, j] = 1
        }
        if (tm.class[i, j] > 1)
          tm.final[i, j] = 1
        if (tm.class[i, j] == 0)
          tm.final[i, j] = f.levels[j]
      }
    }
    subs <- tm != 0
    tr.list <- vector("list", length = nrow(tm))
    for (i in 1:nrow(tm)) {
      tr.list[[i]] <- tm[, subs[i, ], drop = FALSE]
    }
    tr.logic <- vector("list", length = nrow(tm))
    for (i in 1:nrow(tm)) {
      tr.logic[[i]] <- rowMins(tr.list[[i]]) > 0
    }
    tr.subs <- vector("list", length = nrow(tm))
    for (i in 1:nrow(tm)) {
      tr.subs[[i]] <- tr.list[[i]][tr.logic[[i]], , drop = FALSE]
    }
    tr.table <- vector("list", length = nrow(tm))
    for (i in 1:nrow(tm)) {
      tr.table[[i]] <- tm.final[rownames(tr.subs[[i]]), -match(colnames(tr.list[[i]]),colnames(tm)), drop = FALSE]
    }
    tr.comps <- vector("list", length = nrow(tm))
    for (i in 1:nrow(tm)) {
      if (ncol(tr.table[[i]]) == 0) {
        tr.comps[[i]] <- tr.table[[i]]
      }
      else {
        tr.comps[[i]] <- tr.table[[i]][rowMins(tr.table[[i]]) > 0, , drop = FALSE]
      }
    }
    mse.list <- vector("list", length = nrow(tm))
    for (i in 1:nrow(tm)) {
      mse.list[[i]] <- rev(rownames(tr.comps[[i]]))
    }
    mse.H0.list <- vector("list", length = nrow(tm))
    for (i in 1:nrow(tm)) {
      if (length(mse.list[[i]]) == 1) {
        mse.H0.list[[i]] <- "Residuals"
      }
      else {
        mse.H0.list[[i]] <- mse.list[[i]][which(rownames(tm)[[i]] != mse.list[[i]])]
      }
    }
    if (quasi.f == TRUE | quasi.f == T){
      list <- mse.list
      names(list)<-rownames(tm)
      for (i in 1:length(mse.list)){
        for (j in 1:length(mse.list[[i]])){
          list[[i]][j] <- (nomes[,1][which(nomes[,2] == mse.list[[i]][j])])
        }
      }
    }
    mse.tab <- matrix(ncol = 1, nrow = nrow(tm))
    rownames(mse.tab) <- rownames(tm)
    for (i in 1:nrow(tm)) {
      mse.tab[i, ] <- paste(mse.list[[i]], collapse = " + ")
    }
    mse.H0.tab <- matrix(ncol = 1, nrow = nrow(tm))
    rownames(mse.H0.tab) <- rownames(tm)
    for (i in 1:nrow(tm)) {
      mse.H0.tab[i, ] <- paste(mse.H0.list[[i]], collapse = " + ")
    }
    f.versus <- matrix(ncol = 1, nrow = nrow(tm))
    rownames(f.versus) <- rownames(tm)
    colnames(f.versus) <- "F-ratio versus"
    pseudof <- matrix(ncol=2, nrow=nrow(tm))
    pseudof <- as.data.frame(pseudof)
    names(pseudof) <- c("Numerator", "Denominator")
    rownames(pseudof) <- rownames(tm)
    for (i in 1:nrow(tm)) {
      if (length(rownames(mse.tab)[which(mse.H0.tab[i] == mse.tab)]) == 0) {
        if (quasi.f == TRUE | quasi.f == T){
          a <- (list[i])
          for (j in nrow(tm):i){
            b <- (list[j])
            if (names(a) == names(b)) next
            for (k in nrow(tm):i) {
              d <- (list[k])
              for (l in nrow(tm):i) {
                e <- (list[l])
                if (names(d) == names(e)) next
                if(length(a[[1]]) != 1){
                  if (((d[[1]][length(d[[1]])]) != (a[[1]][length(a[[1]])-1])) & ((e[[1]][length(e[[1]])]) != (a[[1]][length(a[[1]])-1]))) next
                }
                ms1 <- c(b[[1]], a[[1]][-length(a[[1]])])
                ms2 <- c(e[[1]], d[[1]])
                sort1 <- sort(ms1)
                sort2 <- sort(ms2)
                if(!identical(sort1, sort2)) next
                if(identical(sort1, sort2)){
                  pseudof$Numerator[i] <- paste(names(b), names(a), sep = " + ")
                  pseudof$Denominator[i] <- paste(names(e), names(d), sep = " + ")
                  break
                }
              }
            }
          }
        }
        f.versus[i] <- "No test"
      } else {
        f.versus[i] <- rownames(mse.tab)[which(mse.H0.tab[i] == mse.tab)]
        if (quasi.f == TRUE | quasi.f == T){
          pseudof$Numerator[i] <- rownames(tm)[i]
          pseudof$Denominator[i] <- rownames(mse.tab)[which(mse.H0.tab[i] == mse.tab)]
        }
      }
      if (mse.H0.tab[i] == "Residuals") {
        f.versus[i] <- "Residuals"
        if (quasi.f == TRUE | quasi.f == T){
          pseudof$Numerator[i] <- rownames(tm)[i]
          pseudof$Denominator[i] <- "Residuals"
        }
      }
    }
    mse <- matrix(ncol = 1, nrow = nrow(tm) + 1)
    colnames(mse) <- "Mean square estimates"
    rownames(mse) <- c(rownames(tm), "Residuals")
    for (i in 1:nrow(tm)) {
      mse[i, ] <- paste(c("Residuals", rev(mse.tab[[i]])), collapse = " + ")
    }
    mse[nrow(mse), ] <- "Residuals"
    n <- nrow(object$model)/object$rank
    Res <- 1
    tm.res <- cbind(tm.final, n)
    tm.res <- rbind(tm.res, Res)
    nomes[,2] <- rownames(f.versus)
    for(i in 1:nrow(f.versus)){
      a <- strsplit(rownames(f.versus)[i],  "\\:")
      b <- object$model[match(a[[1]],colnames(object$model))]
      c <- interaction(b, lex.order = TRUE)
      res <- length(c[which(c == levels(c)[1])])
      nomes[i,3] <- res
    }
    mse.final <- mse
    mod <- strsplit(mse.final,  " \\+ ")
    for(i in 1:length(mod)){
      mod[[i]][1] <- paste("1*", mod[[i]][1], sep = "")
      g <- mod[[i]]
      for(j in 2:length(mod[[i]])){
        if (length(g) == 1) next
        g[j] <- paste(c((nomes[,3][which(nomes[,2] == g[j])]), "*", g[j]), collapse = "")
      }
      h <- paste(g, collapse = " + ")
      mod[[i]] <- paste(h, collapse = " + ")
      mse.final[i] <- mod[[i]]
    }
    if (quasi.f == FALSE | quasi.f == F){
      pseudof$Numerator <- rownames(f.versus)
      pseudof$Denominator <- f.versus[,1]
    }
    estimates <- list(tm = tm.res, mse = mse.final, f.versus = pseudof)
    return(estimates)
  }
