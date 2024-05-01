#' Student-Newman-Keuls (SNK) procedure
#'
#' This function performs a SNK post-hoc test of means on the factors of a chosen term of the model, comparing among levels of one factor within each level of other factor or combination of factors.
#'
#' @param object an object of class "\code{\link{lm}}" containing the specified design.
#' @param term	term of the model to be analysed. Argument \code{term} can be a main effect, an interaction or a nested factor, which must be specified using "quotes". Use \code{\link{estimates}} to see the right form to inform it. See examples below.
#' @param among specifies the factor which levels will be compared among. Need to be specified if the term to be analysed envolves more than one factor.
#' @param within specifies the factor or combination of factors that will be compared within level among.
#' @param anova.tab an object containing the results returned by \code{\link{gad}} or \code{\link{pooling}} functions. See examples below.
#' @details SNK is a stepwise procedure for hypothesis testing. First the sample means are sorted, then the pairwise studentized range (q) is calculated by dividing the differences between means by the standard error, which is based upon the average variance of the two samples.
#' @return A list containing the standard error, the degrees of freedom and pairwise comparisons among levels of one factor within each level of other(s) factor(s).
#' @author Maurício Camargo (\email{mauricio.camargo@@furg.br})
#' @author Eliandro Gilbert (\email{eliandrogilbert@@gmail.com})
#' @author Leonardo Sandrini-Neto (\email{leonardosandrini@@ufpr.br})
#' @seealso \code{\link{gad}}, \code{\link{estimates}}
#' @examples
#' library(GAD)
#'
#' # Example 1
#' data(rohlf95)
#' CG <- as.fixed(rohlf95$cages)     # a fixed factor
#' MQ <- as.random(rohlf95$mosquito) # a random factor nested in cages
#'
#' model <- lm(wing ~ CG + CG%in%MQ, data = rohlf95)
#' model.tab <- gad(model) # storing ANOVA table in an object
#' model.tab               # checking ANOVA results
#' estimates(model)        # checking model structure
#'
#' # Comparison among levels of mosquito ("MQ") within each level of cage ("CG")
#' snk.test(model, term = "CG:MQ", among = "CG", within = "MQ", anova.tab = model.tab)
#'
#' # Example 2
#' data(snails)
#' O <- as.random(snails$origin)   # a random factor
#' S <- as.random(snails$shore)    # a random factor orthogonal to origin
#' B <- as.random(snails$boulder)  # a random factor nested in shore
#' C <- as.random(snails$cage)     # a random factor nested in the combination of boulder and origin

#' model2 <- lm(growth ~ O + S + O*S + B%in%S + O*(B%in%S) + C%in%(O*(B%in%S)), data = snails)
#' model2.tab <- gad(model2, quasi.f = FALSE) # storing ANOVA table in an object
#' model2.tab                                 # checking ANOVA results
#' estimates(model2, quasi.f = FALSE)         # checking model structure
#'
#' # Comparison among levels of "origin"
#' snk.test(model2, term = "O", anova.tab = model2.tab)
#' # Comparison among levels of "shore" within each level of "origin"
#' snk.test(model2, term = "O:S", among = "S", within = "O", anova.tab = model2.tab)
#' # If term "O:S:B" were significant, we could try
#' snk.test(model2, term = "O:S:B", among = "B", within = "O:S", anova.tab = model2.tab)
#' @importFrom stats ptukey
#' @export
snk.test <-
  function(object, term = NULL, among = NULL, within = NULL, anova.tab = NULL){
    tab <- anova.tab
    if(is.null(tab)){
      stop("The argument anova.tab cannot be empty")
    }
    if(is.null(tab$anova)){
      stop("Something is wrong here! The argument anova.tab must store the result of gad or pooling functions")
    }

    if(is.null(term)){
      stop("You must specify a term of the linear model")
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
      if(is.null(f.versus)){
        f.versus <- a$f.versus
        }
      }
    nomes <- matrix(ncol = 2, nrow = nrow(mse))
    nomes[,2] <- rownames(mse)
    mod <- gsub(pattern = "\\*", replacement = " + ",  x = mse)
    mod <- strsplit(mod,  " \\+ ")
    mse2 <- list(1)
    for(i in 1:nrow(mse)){
      mse2[[i]] <- mod[[i]][seq(from = 1, to = length(mod[[i]]),by = 2)]
    }
    names(mse2) <- rownames(mse)
    for(i in 1:nrow(nomes)){
      nomes[i,1] <- mse2[[i]][length(mse2[[i]])]
    }
    teste <- strsplit(f.versus$Denominator[match(term, rownames(f.versus))], " \\+ ")
    if(quasi.f == TRUE){
      if(length(teste[[1]]) > 1) {
        se.num <- ((anova.res$Mean[which(rownames(anova.res) == teste[[1]][1])]) + (anova.res$Mean[which(rownames(anova.res) == teste[[1]][2])]))
        n.den <- strsplit(f.versus$Numerator[match(term, rownames(f.versus))], " \\+ ")
        se.den <- (as.numeric(nomes[,1][which(n.den[[1]][1] == nomes[,2])]) + as.numeric(nomes[,1][which(n.den[[1]][2] == nomes[,2])]))
        se <- sqrt(se.num/se.den)
        Df <- f.versus$Den.Df[match(term, rownames(f.versus))]
        }
      }
    if(quasi.f == FALSE | length(teste[[1]]) == 1){
      se.num <- anova.res$Mean[which(rownames(anova.res) == f.versus$Denominator[match(term, rownames(f.versus))])]
      se.den <- as.numeric(nomes[,1][which(term == nomes[,2])])
      se <- sqrt(se.num/se.den)
      if (length(se) == 0) stop ("F-ratio versus = No test \nStandard error for this set of means cannot be calculated!")
      Df <- anova.res$Df[which(rownames(anova.res) == f.versus$Denominator[match(term, rownames(f.versus))])]
      }
    if(length(among) == 0){
      set.means <- tapply(object$model[,1], subset(object$model, select = term), mean)
      a <- length(set.means)
      rank.means <- sort(set.means)
      m.means <- matrix(NA, ncol = a, nrow = a-1)
      m.final <- matrix("", ncol = a, nrow = a-1)
      i <- 0; g <- a:2
      for(i in 1:nrow(m.means)){
        for(j in seq(i)){
          m.means[i,j] <- rank.means[a-i+j]-rank.means[j]
        }
      }
      Q <- m.means/se
      pvalue <- Q
      for(i in 1:nrow(m.means)){
        for(j in 1:ncol(m.means)){
          pvalue[i,j] <- ptukey(Q[i,j], nmeans = g[i], df = Df, lower.tail = FALSE)
        }
      }
      p <- pvalue
      for(i in 1:nrow(pvalue)){
        for(j in seq(i)){
          if(pvalue[i,j] <= 0.001) p[i,j] <- "***"
          else if (pvalue[i,j] <= 0.01) p[i,j] <- "**"
          else if (pvalue[i,j] <= 0.05) p[i,j] <- "*"
          else if (pvalue[i,j] > 0.05) p[i,j] <- "ns"
        }
      }
      ns <- which(p == "ns", arr.ind = TRUE)
      if(nrow(ns) == 0){
        pos <- ns
      } else {
        pos <- subset(ns, ns[,1] != max(ns[,1]))
      }
      if(nrow(pos) != 0){
        sublist <- vector("list", length = nrow(pos))
        for(k in 1:nrow(pos)){
          sublist[[k]] <- p[pos[k,][1]:nrow(p), pos[k,][2]:(pos[k,][2] + length(pos[k,][1]:nrow(p)) - 1)]
          for(i in 2:nrow(sublist[[k]])){
            for(j in seq(i)){
              sublist[[k]][i,j] <- "x"
              p[pos[k,][1]:nrow(p), pos[k,][2]:(pos[k,][2] + length(pos[k,][1]:nrow(p)) - 1)] <- sublist[[k]]
            }
          }
        }
      }
      for(i in 1:nrow(p)){
        for(j in seq(i)){
          if(p[i,j] == "***") m.final[i,j] <- paste(paste(c(a-i+j, j), collapse = "-"), "***")
          else if(p[i,j] == "**") m.final[i,j] <- paste(paste(c(a-i+j, j), collapse = "-"), "**")
          else if(p[i,j] == "*") m.final[i,j] <- paste(paste(c(a-i+j, j), collapse = "-"), "*")
          else if(p[i,j] == "ns") m.final[i,j] <- paste(paste(c(a-i+j, j), collapse = "-"), "ns")
          else if(p[i,j] == "x") m.final[i,j] <- paste(paste(c(a-i+j, j), collapse = "-"), "x")
        }
      }
      result <- data.frame(m.final, stringsAsFactors = FALSE)
      result <- rbind(as.character(""), result)
      result <- rbind(round(rank.means, 4), result)
      result <- rbind(as.character(1:ncol(m.final)), result)
      colnames(result) <- names(rank.means)
      rownames(result) <- c("Rank order:", "Ranked means:", "Comparisons:", 1:nrow(m.means))
      p.line <- p[nrow(p),]
      p.line <- p.line[-length(p.line)]
      p.res <- p.line
      for(i in 1:length(p.line )){
        if(p.line[i] == "ns" | p.line[i] == "x") { p.res[i] <- "="
        } else {
          p.res[i] <- "<"
        }
      }
      res <- character(length = (length(names(rank.means))*2)-1)
      res[seq(from = 1, to = length(res), by = 2)] <- names(rank.means)
      res[seq(from = 2, to = length(res), by = 2)] <- p.res
    } else {
      within.split <- strsplit(within, ":")
      within.int <- interaction(object$model[, within.split[[1]]], lex.order = TRUE)
      new.model <- cbind(object$model, within.int)
      set <- split(new.model, f = factor(within.int, levels = unique(within.int)))
      set.means <- set
      for(i in 1:length(set)){
        set.means[[i]] <- tapply(set[[i]][,1], subset(set[[i]], select = among), mean)
      }
      a <- length(set.means[[1]])
      rank.means <- set.means
      for(i in 1:length(set.means)){
        rank.means[[i]] <- sort(set.means[[i]])
      }
      m.means <- matrix(NA, ncol = a, nrow = a-1)
      m.final <- matrix("", ncol = a, nrow = a-1)
      i <- 0; g <- a:2
      l.means <- rank.means
      for(i in 1:length(l.means)){
        l.means[[i]] <- m.means
      }
      for(i in 2:a-1){
        for(j in 1:i){
          for(k in 1:length(l.means)){
            l.means[[k]][i,j] <- as.vector(rank.means[[k]])[a-i+j]-as.vector(rank.means[[k]])[j]
          }
        }
      }
      Q <- l.means
      for(i in 1:nrow(m.means)){
        for(j in 1:ncol(m.means)){
          for(k in 1:length(Q)){
            Q[[k]][i,j] <- l.means[[k]][i,j]/se
          }
        }
      }
      pvalue <- Q
      for(i in 1:nrow(m.means)){
        for(j in 1:ncol(m.means)){
          for(k in 1:length(Q)){
            pvalue[[k]][i,j] <- ptukey(Q[[k]][i,j], nmeans = g[i], df = Df, lower.tail = FALSE)
          }
        }
      }
      p <- pvalue
      for(i in 1:length(p)){
        p[[i]] <- m.final
      }
      for(i in 2:a-1){
        for(j in 1:i){
          for(k in 1:length(p)){
            if(pvalue[[k]][i,j] <= 0.001) p[[k]][i,j] <- "***"
            else if(pvalue[[k]][i,j] <= 0.01) p[[k]][i,j] <- "**"
            else if(pvalue[[k]][i,j] < 0.05) p[[k]][i,j] <- "*"
            else if(pvalue[[k]][i,j] > 0.05) p[[k]][i,j] <- "ns"
          }
        }
      }
      ns <- p
      for(i in 1:length(p)){
        ns[[i]] <- which(p[[i]] == "ns", arr.ind = TRUE)
      }
      pos <- ns
      for(i in 1:length(ns)){
        if(nrow(ns[[i]]) == 0){
          pos[[i]] <- ns[[i]]
        } else {
          pos[[i]] <- subset(ns[[i]], ns[[i]][,1] != max(ns[[i]][,1]))
        }
      }
      l.final <- pvalue
      for(i in 1:length(l.final)){
        l.final[[i]] <- m.final
      }
      for(l in 1:length(pos)){
        if(nrow(pos[[l]]) != 0){
          sublist <- vector("list", length = nrow(pos[[l]]))
          for(k in 1:nrow(pos[[l]])){
            sublist[[k]] <- p[[l]][pos[[l]][k,][1]:nrow(p[[l]]), pos[[l]][k,][2]:(pos[[l]][k,][2] + length(pos[[l]][k,][1]:nrow(p[[l]])) - 1)]
            for(i in 2:nrow(sublist[[k]])){
              for(j in seq(i)){
                sublist[[k]][i,j] <- "x"
                p[[l]][pos[[l]][k,][1]:nrow(p[[l]]), pos[[l]][k,][2]:(pos[[l]][k,][2] + length(pos[[l]][k,][1]:nrow(p[[l]])) - 1)] <- sublist[[k]]
              }
            }
          }
        }
      }
      for(i in 2:a-1){
        for(j in 1:i){
          for(k in 1:length(l.final)) {
            if(p[[k]][i,j] == "***") l.final[[k]][i,j] <- paste(paste(c(a-i+j, j), collapse = "-"), "***")
            else if(p[[k]][i,j] == "**") l.final[[k]][i,j] <- paste(paste(c(a-i+j, j), collapse = "-"), "**")
            else if(p[[k]][i,j] == "*") l.final[[k]][i,j] <- paste(paste(c(a-i+j, j), collapse = "-"), "*")
            else if(p[[k]][i,j] == "ns") l.final[[k]][i,j] <- paste(paste(c(a-i+j, j), collapse = "-"), "ns")
            else if(p[[k]][i,j] == "x") l.final[[k]][i,j] <- paste(paste(c(a-i+j, j), collapse = "-"), "x")
          }
        }
      }
      result <- l.final
      for(i in 1:length(l.final)){
        result[[i]] <- data.frame(l.final[[i]], stringsAsFactors = FALSE)
        result[[i]] <- rbind(as.character(""), result[[i]])
        result[[i]] <- rbind(round(rank.means[[i]], 4), result[[i]])
        result[[i]] <- rbind(1:ncol(m.means), result[[i]])
        colnames(result[[i]]) <- names(rank.means[[i]])
        rownames(result[[i]]) <- c("Rank order:", "Ranked means:", "Comparisons:", 1:nrow(m.means))
      }
      p.line <- vector("list", length = length(p))
      for(l in 1:length(p)){
        p.line[[l]] <- p[[l]][nrow(p[[l]]),]
        p.line[[l]] <- p.line[[l]][-length(p.line[[l]])]
      }
      p.res <- p.line
      for(l in 1:length(p)){
        for(i in 1:length(p.line[[l]])){
          if(p.line[[l]][i] == "ns" | p.line[[l]][i] == "x") { p.res[[l]][i] <- "="
          } else {
            p.res[[l]][i] <- "<"
          }
        }
      }
      res <- vector("list", length = length(p))
      for(l in 1:length(p)){
        res[[l]] <- character(length = (length(names(rank.means[[l]]))*2)-1)
        res[[l]][seq(from = 1, to = length(res[[l]]), by = 2)] <- names(rank.means[[l]])
        res[[l]][seq(from = 2, to = length(res[[l]]), by = 2)] <- p.res[[l]]
      }
    }
    cat("Student-Newman-Keuls test for:", term, "\n")
    cat("\nStandard error =", round(se, 4))
    cat("\nDf =", Df, "\n")
    if(length(among) == 0){
      print(result, right = FALSE)
      cat("Summary:", res, "\n")
    } else {
      cat("\nPairwise comparisons among levels of:", among,
          "\nwithin each level of:", within, "\n")
      for(i in 1:length(result)){
        cat("\nLevel:", names(result)[[i]], "\n")
        print(result[[i]], right = FALSE)
        cat("Summary:", res[[i]], "\n")
      }
    }
    cat("---")
    cat("\nSignif. codes: ***p < 0.001; **p < 0.01; *p < 0.05; ns, p > 0.05\n")
  }
