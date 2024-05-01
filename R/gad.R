#' General analysis of variance design
#'
#' Fits a general ANOVA design with any combination of orthogonal/nested and fixed/random factors through function \code{\link{estimates}}.
#'
#' @param object an object of class lm, containing the specified design with random and/or fixed factors.
#' @param quasi.f logical, indicating whether to use quasi F-ratio when there is no single error term appropriate in the analysis. Default to \code{FALSE}.
#' @details Function \code{gad} returns an analysis of variance table using \code{\link{estimates}} to identify the appropriate F-ratios and consequently p-values for any complex model of orthogonal or nested, fixed or random factors as described by Underwood(1997).
#' @return A "\code{list}" containing an object of class \code{"anova"} inheriting from class \code{"data.frame"}.
#' @author Leonardo Sandrini-Neto (\email{leonardosandrini@@ufpr.br})
#' @references Underwood, A.J. 1997. \emph{Experiments in Ecology}: Their Logical Design and Interpretation Using Analysis of Variance. Cambridge University Press, Cambridge.
#' @seealso \code{\link{estimates}}
#' @examples
#' # Example 1
#' library(GAD)
#' data(rohlf95)
#' CG <- as.fixed(rohlf95$cages)
#' MQ <- as.random(rohlf95$mosquito)
#' model <- lm(wing ~ CG + MQ%in%CG, data = rohlf95)
#' model.tab <- gad(model)
#' model.tab
#'
#' # Example 2
#' data(rats)
#' TR <- as.fixed(rats$treat)
#' RA <- as.random(rats$rat)
#' LI <- as.random(rats$liver)
#' model2 <- lm(glycog ~ TR + RA%in%TR + LI%in%RA%in%TR, data = rats)
#' model2.tab <- gad(model2)
#' model2.tab
#'
#' # Example 3
#' data(snails)
#' O <- as.random(snails$origin)
#' S <- as.random(snails$shore)
#' B <- as.random(snails$boulder)
#' C <- as.random(snails$cage)
#' model3 <- lm(growth ~ O + S + O*S + B%in%S + O*(B%in%S) + C%in%(O*(B%in%S)), data = snails)
#' model3.tab <- gad(model3)
#' model3.tab # 'no test' for shore
#' model3.tab2 <- gad(model3, quasi.f = TRUE)
#' model3.tab2 # suitable test for shore
#' @importFrom stats anova
#' @importFrom stats formula
#' @importFrom stats pf
#' @export
gad <-
function (object, quasi.f=FALSE)
{
  anova.res <- anova(object) # anova normal do modelo, com F errado
  table <- anova.res[, 1:3] # aproveita somente df, SS e MS
  est<-estimates(object, quasi.f=quasi.f) # chama o estimates pro modelo, com ou sem quasi.f
  f.versus <- est$f.versus #extrai o data.frame f.versus
  f.versus$Num.Df<-numeric(nrow(f.versus)) # cria as colunas de graus de liberdade
  f.versus$Den.Df<-numeric(nrow(f.versus)) # para o numerador e denominador

  f <- numeric(nrow(table))
  for (i in 1:nrow(f.versus)) {
    if (f.versus$Denominator[i] == "No test") {
      f[i]<-NA
      } # quando for 'no test' nao tem F value e quasi.f sera FALSE
    teste<-strsplit(f.versus$Denominator[i], " \\+ ") # tenta separar o denominador em dois (pra ver se foi usado quasi.f)
    if (length(teste[[1]]) > 1){ # se dividir eh pq foi usado quasi.f
      num<-strsplit(f.versus$Numerator[i], " \\+ ") # divide numerador e denominador
      den<-strsplit(f.versus$Denominator[i], " \\+ ") # e faz uma caralhada de calculos

      ms.num<-((table$Mean[which(rownames(table) == num[[1]][1])])  + (table$Mean[which(rownames(table) == num[[1]][2])]))
      num.1<- (((table$Mean[which(rownames(table) == num[[1]][1])])  + (table$Mean[which(rownames(table) == num[[1]][2])]))^2)
      num.2<- (((table$Mean[which(rownames(table) == num[[1]][1])])^2) / table$Df[which(rownames(table) == num[[1]][1])])
      num.3<- (((table$Mean[which(rownames(table) == num[[1]][2])])^2) / table$Df[which(rownames(table) == num[[1]][2])])
      df.num<- ( num.1 / (num.2 + num.3) )

      ms.den<-((table$Mean[which(rownames(table) == den[[1]][1])])  + (table$Mean[which(rownames(table) == den[[1]][2])]))
      den.1<- (((table$Mean[which(rownames(table) == den[[1]][1])])  + (table$Mean[which(rownames(table) == den[[1]][2])]))^2)
      den.2<- (((table$Mean[which(rownames(table) == den[[1]][1])])^2) / table$Df[which(rownames(table) == den[[1]][1])])
      den.3<- (((table$Mean[which(rownames(table) == den[[1]][2])])^2) / table$Df[which(rownames(table) == den[[1]][2])])
      df.den<- ( den.1 / (den.2 + den.3) )

      f.versus$Num.Df[i]<-round(df.num, digits = 2) # armazena os graus de liberdade aproximados do numerador
      f.versus$Den.Df[i]<-round(df.den, digits = 2) # armazena os graus de liberdade aproximados do denominador

      f[i]<- (ms.num/ms.den) # calcula e armazena o F aproximado
      }
    if (f.versus$Denominator[i] == "Residuals") {
      f[i] <- table$Mean[i]/table$Mean[nrow(table)]
      } # caso o denominador seja o residuo
    if ((f.versus$Denominator[i] != "Residuals") & (length(teste[[1]]) == 1) & (f.versus$Denominator[i] != "No test")){
      f[i] <- table$Mean[i]/table$Mean[which(rownames(table)==f.versus$Denominator[i])]
      } # se o denominador nao for o residuo, nao for um quasi.f test, e nao for 'no test'
    } # calcula o F usando o denominador apropriado da tabela de f.versus
  P <- numeric(nrow(table)) # depois parte para o p value
  for (i in 1:nrow(f.versus)) {
    if (f.versus$Denominator[i] == "No test") {
      P[i]<-NA # se o denominador do termo for 'no test' nao tem calculo do p value
    }
    teste<-strsplit(f.versus$Denominator[i], " \\+ ")
    if (length(teste[[1]]) > 1){
      P[i] <- pf(as.numeric(f[i]), f.versus$Num.Df[i], f.versus$Den.Df[i], lower.tail = FALSE)
    } # se o denominador for uma combinacao de MS utilizando quasi.f usa os graus de liberdade aproximados
    if (f.versus$Denominator[i] == "Residuals") {
      P[i] <- pf(as.numeric(f[i]), table$Df[i], object$df.residual,lower.tail = FALSE)
      f.versus$Num.Df[i]<-signif(table$Df[i],3)
      f.versus$Den.Df[i]<-signif(object$df.residual,3)
    }# se o denominador for o residuo calcula normalmente e armazena os graus de liberdade nas colunas do data.frame f.versus
    if ((f.versus$Denominator[i] != "Residuals") & (length(teste[[1]]) == 1) & (f.versus$Denominator[i] != "No test")){
      P[i] <- pf(as.numeric(f[i]), table$Df[i], table$Df[which(rownames(table)==f.versus$Denominator[i])], lower.tail = FALSE)
      f.versus$Num.Df[i]<-signif(table$Df[i],3)
      f.versus$Den.Df[i]<-signif(table$Df[which(rownames(table)==f.versus$Denominator[i])],3)
    }# se o denominador nao for 'no test', nem quasi.f test e nem o residuo calcula normalmente e armazena os graus de liberdade nas colunas do data.frame f.versus
  }
  anova.table <- data.frame(table, f, P) # junta tudo em uma nova tabela de anova
  anova.table[length(P), 4:5] <- NA
  colnames(anova.table) <- colnames(anova.res)
  rownames(anova.table) <- c(rownames(f.versus), "Residuals")
  if(quasi.f==TRUE | quasi.f == T){
    aviso<-character(length=nrow(f.versus)) # caso tenha sido usado quasi.f=TRUE sera gerado um aviso

    for(i in 1:nrow(f.versus)){ # que ira indicar quais termos tem valores de F aproximados
      teste<-strsplit(f.versus$Denominator[i], " \\+ ")
      if (length(teste[[1]]) > 1){
        aviso[i]<-rownames(f.versus)[i]
        }
      }
    aviso<-aviso[which(aviso != "")]

    if (length(aviso) == 1){ # se somente 1 'no test' foi resolvido com quasi.f
      aviso.final<- paste("The F-ratio of", aviso, "is approximate. See 'help(estimates)' for further details\n")
      }

    if (length(aviso) > 1){# se mais que 1 'no test' foi resolvido com quasi.f
      aviso1<-paste(aviso[1:(length(aviso)-1)], collapse=", ")
      aviso.final<- paste("The F-ratios of", aviso1, "and", aviso[length(aviso)], "are approximations. See 'help(estimates)' for further details\n")
      }

    if (length(aviso) == 0){ # se todos os termos tinham um F versus apropriado e mesmo assim quasi.f tenha sido utilizado como TRUE
      aviso.final<- "All of your terms has an appropriated F-ratio, Quasi F-test is not necessary\n"
      }

    tab.anova<-structure(anova.table, heading = c("Analysis of Variance Table: Quasi F-ratios\n",
                                                  aviso.final,
                                                  paste("Response:", deparse(formula(object)[[2L]]))),
                         class = c("anova", "data.frame"))
    result <- list(f.versus = f.versus, anova = tab.anova) # A tabela de F versus sera mostrada apenas com quasi.f=TRUE
    } # nao faz sentido mostrar quando quasi.f=FALSE
  else{ # caso quasi.f=FALSE
    tab.anova<-structure(anova.table, heading = c("Analysis of Variance Table\n",
                                                  paste("Response:", deparse(formula(object)[[2L]]))), class = c("anova", "data.frame"))
    result <- list(anova = tab.anova)
    }
  return(result)
  }
