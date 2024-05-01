#' Post-hoc pooling
#'
#' Performs a \emph{post-hoc} pooling by combining or completely excluding terms from linear models
#'
#' @param object an object of class "\code{\link{lm}}", containing the specified design with random and/or fixed factors.
#' @param term the term which will be removed from model.
#' @param method {method for removing a term from the model. Could be \code{method = "eliminate"} for completely exclude the term from model,
#'                or \code{method = "pool"} for pool the selected term with its appropriated F-ratio. Default to \code{method = "pool"}.}
#' @param anova.tab an object containing the results returned by \code{gad} or \code{pooling} functions. See examples below.
#'
#' @details {
#' \emph{Post-hoc} pooling is a procedure to remove terms from a model. It might be done by several reasons:
#' (i) lack of evidence against the null hypothesis of that term; (ii) a negative estimate of that term's component of variation (Fletcher and Underwood, 2002);
#' (iii) the hypothesis of interest can not be tested until some terms are excluded from the model (Anderson et al., 2008).
#' According to literature the term's p-value should exceed 0.25 before removing it (Underwood, 1997).\cr
#'
#' There are two different methods to remove a term from the model, determinated by \code{method} argument. When \code{method = "eliminate"}
#' the chosen term is completely excluded from the model and its sum of squares and degrees of freedom are pooled with the residual sum of squares and
#' degrees of freedom, as if the selected term had never been part of the model. When \code{method = "pool"} the chosen term's sum of squares and
#' degrees of freedom are pooled with its denominator's sum of squares and degrees of freedom. The removal of terms using \code{method = "pool"} will be
#' appropriated for most of situations (Anderson et al., 2008).\cr
#'
#' Note that removing a term has consequences for the construction of F-ratios (or quasi F-ratios), p-values and the estimation of components for the
#' remaining terms, so should be done wisely. When there is more than one term which might be removed from the model (which p-value exceed 0.25), it is
#' recommended to begin with the one having the smallest mean square (Anderson et. al, 2008).\cr
#'
#' Function \code{pooling} removes one term at once. After the removal of the term of interest, one should re-assess whether or not more terms should be
#' removed. If it is the case, the output of \code{pooling} function should be stored in a new object and the function should be run again, using this new
#' object in the \code{data} argument. This can be done successively. The way of \code{pooling} function does the analysis, step-by-step and storing the result
#' in a new object at each step, gives the user total control of what happens and makes it easier return to the previous results.
#' }
#' @return {A list of length 4, containing the table of pooled terms (\code{$pool.table}), the mean squares estimates (\code{$mse}),
#'          the F-ratio versus (\code{$f.versus}) and the result of the analysis of variance (\code{$anova}).}
#' @author Eliandro Gilbert (\email{eliandrogilbert@@gmail.com})
#' @references Anderson, M.J., Gorley, R.N., Clarke, K.R. 2008. \emph{PERMANOVA+ for PRIMER}: Guide to Software and Statistical Methods. PRIMER-E: Plymouth, UK.
#' @references Fletcher, D.J., Underwood, A.J. 2002. How to cope with negative estimates of components of variance in ecological field studies. Journal of Experimental Marine Biology and Ecology 273, 89-95.
#' @references Underwood, A.J. 1997. \emph{Experiments in Ecology}: Their Logical Design and Interpretation Using Analysis of Variance. Cambridge University Press, Cambridge.
#' @seealso \code{\link{estimates}}, \code{\link{gad}}, \code{\link{comp.var}}
#' @examples
#' library(GAD)
#' data(snails)
#' O <- as.random(snails$origin)   # a random factor
#' S <- as.random(snails$shore)    # a random factor orthogonal to origin
#' B <- as.random(snails$boulder)  # a random factor nested in shore
#' C <- as.random(snails$cage)     # a random factor nested in the combination of boulder and origin
#'
#' model <- lm(growth ~ O + S + O*S + B%in%S + O*(B%in%S) + C%in%(O*(B%in%S)),data = snails)
#' estimates(model, quasi.f = FALSE)  # 'no test' for shore
#' gad(model, quasi.f = FALSE)        # no results for shore term
#' estimates(model, quasi.f = TRUE)   # suitable test for shore
#' gad(model, quasi.f = TRUE)         # test result for shore
#'
#' # An alternative of using linear combinations of mean squares is the pooling function.
#' model.tab <- gad(model, quasi.f = FALSE) # stores the result of ANOVA on a new object
#' pooling(model, term = "S:B", method = "pool", anova.tab = model.tab)      # pooling terms
#' pooling(model, term = "S:B", method = "eliminate", anova.tab = model.tab) # or eliminating terms
#' @importFrom stats pf
#' @export
pooling <-
function(object, term = NULL, method = "pool", anova.tab = NULL){
  tab <- anova.tab
  if (is.null(tab)){
    stop ("Data cannot be null")
    }
  if (is.null(tab$anova)){
    stop ("Something is wrong here. The argument anova.tab must store the result of gad or pooling functions")
    }

  if (is.null(term)){
    stop ("You must specify which term should be removed")
  }

  f.versus<-tab$f.versus
  anova.res<-tab$anova
  quasi.f<-0
  if (is.null(f.versus$Num.Df)){
    quasi.f<-FALSE
  }
  if (!is.null(f.versus$Num.Df)){
    quasi.f<-TRUE
  }

  mse<-tab$mse
  if (is.null(tab$mse)){
    a<-estimates(object)
    mse<-a$mse
  }

  nomes<-matrix(ncol=2,nrow=nrow(mse))
  nomes[,2]<-rownames(mse)
  for (i in 1:nrow(nomes)){
    nomes[i,1]<-as.character(x=i)
  }

  mod<-gsub(pattern = "\\*", replacement = " + ",  x = mse)
  mod<-strsplit(mod,  " \\+ ")
  mse2<-list(1)
  mse21<-list(1)
  for (i in 1:nrow(mse)){
    mse2[[i]]<-mod[[i]][seq(from = 2, to=length(mod[[i]]),by = 2)]
    mse21[[i]]<-mod[[i]][seq(from = 1, to=length(mod[[i]]),by = 2)]
  }
  names(mse2)<-rownames(mse)
  names(mse21)<-rownames(mse)

  names(mse2)<-nomes[,1]
  names(mse21)<-nomes[,1]
  termo<-term

  if (all(rownames(mse) != termo)){
    texto1<-"does not exist in your model"
    aviso<-paste(termo, texto1, sep=" ")
    stop (aviso)
  }
  termo2<-nomes[,1][which(nomes[,2]==termo)]

  for (i in 1:length(mse2)){
    for (j in 1:length(mse2[[i]])){
      mse2[[i]][j]<-(nomes[,1][which(nomes[,2] == mse2[[i]][j])])
    }
  }

  for(i in 1:length(mse2)){
    b<-gsub(pattern = termo2, replacement = "",  x = mse2[[i]])
    for (j in 1:length(mse2[[i]])){
      if(b[j] == ""){
        mse21[[i]][j]<-"_"
      }
    }
    c<-gsub(pattern = "_", replacement = "",  x = mse21[[i]])
    mse2[[i]]<-b[which(b != "")]
    mse21[[i]]<-c[which(c != "")]
  }
  termo3 <- mse2[which(names(mse2) == termo2)]

  if (method == "pool"){
    pooled <-  list(NA)
    for (j in 1:length(mse2)){
      a<-paste(mse2[[j]],collapse=" + ")
      b<-paste(termo3[[1]],collapse=" + ")
      if(names(mse2)[[j]] == names (termo3)[[1]])next
      if (a == b){
        pooled[[1]]<- mse2[[j]]
        names(pooled)<-names(mse2)[[j]]
      }
    }
    if (is.na(pooled)){
      pooled<-mse2[[length(mse2)]]
      names(pooled)<-names(mse2)[[length(mse2)]]
    } # se nao encontrar na mse, o termo sera o residuo
    names(pooled)<-nomes[,2][which(nomes[,1]==names(pooled))]
  }

  mse3<-mse2
  names(mse3)<-nomes[,2][which(nomes[,1] == names(mse3))]
  names(mse21)<- names(mse3)
  mse3<-mse3[which(names(mse3) != termo)] # exclui a linha do termo selecionado da mse
  mse21<-mse21[which(names(mse21) != termo)]

  ################# TESTEEEEE
  if (method == "pool"){
    if (!is.null(tab$mse)){#se tab$mse existir, data ja armazena resultado de pooling anterior
      y<-rownames(tab$pool.table)[which(grepl(pattern="Pooled", x = rownames(tab$pool.table)))]
      if(length(y) != 0){#caso ja exista algum Pooledx
        a<-(gsub(pattern = "Pooled", replacement = "",  x = y))
        nome<-paste("Pooled",c(as.numeric(max(a))+1),sep="")
        names(mse3)[which(names(mse3) == names(pooled))] <- nome
        names(mse21)[which(names(mse21) == names(pooled))] <- nome
        }#cria um Pooledx+1
      if (length(y) == 0){#caso nao exista um Pooledx em pool.table
        nome<-paste("Pooled",1,sep ="")
        names(mse3)[which(names(mse3) == names(pooled))] <- nome
        names(mse21)[which(names(mse21) == names(pooled))] <- nome
        }#cria o Pooled1
      }
    if (is.null(tab$mse)){ #se tab$mse nao existir, data armazena resultado do gad
      nome<-paste("Pooled",1,sep ="")
      names(mse3)[which(names(mse3) == names(pooled))] <- nome
      names(mse21)[which(names(mse21) == names(pooled))] <- nome
      } #primeiro pooling = Pooled1
    }

  #################

#  if (method == "pool"){
#    for (i in 1:length(mse)){
#      y<-rownames(mse)[which(rownames(mse) == "Pooled1" |
#        rownames(mse) == "Pooled2" | rownames(mse) == "Pooled3" |
#        rownames(mse) == "Pooled4" | rownames(mse) == "Pooled5" |
#        rownames(mse) == "Pooled6" | rownames(mse) == "Pooled7" |
#        rownames(mse) == "Pooled8" | rownames(mse) == "Pooled9" |
#        rownames(mse) == "Pooled10")]
#    }
#    if(length(y)!=0){
#      a<-as.numeric(gsub(pattern = "Pooled", replacement = "",  x = y))
#      nome<-paste("Pooled",c(as.numeric(max(a))+1),sep="")
#      names(mse3)[which(names(mse3) == names(pooled))] <- nome
#      names(mse21)[which(names(mse21) == names(pooled))] <- nome
#    } # caso jah exista algum termo 'Pooledx' no modelo (fruto de um pooling anterior) sera adicionado um valor maior ao nome do novo 'Pooledx'
#    x<-format(textcnt(rownames(mse), n=6L, useBytes=T))
#    y<-x[which(rownames(x)=="Pooled"),1]
#    if(length(y)!=0){
#      a<-nomes[which(grepl(pattern="Pooled", x = rownames(mse))),2]
#      a<-as.numeric(gsub(pattern = "Pooled", replacement = "",  x = a))
#      nome<-paste("Pooled",c(max(a)+1),sep="")
#      names(mse3)[which(names(mse3) == names(pooled))] <- nome
#      names(mse21)[which(names(mse21) == names(pooled))] <- nome
#    } # caso jah exista algum termo 'Pooledx' no modelo (fruto de um pooling anterior) sera adicionado um valor maior ao nome do novo 'Pooledx'
#    if (length(y)==0){
#      nome<-paste("Pooled",1,sep ="")
#      names(mse3)[which(names(mse3) == names(pooled))] <- nome
#      names(mse21)[which(names(mse21) == names(pooled))] <- nome
#      }
#    }# se nao existir nenhum termo 'Pooledx' no modelo sera gerado o 'Pooled1' no lugar do termo que recebeu os 'df' e 'SS'

  f.versus<-matrix(ncol=4,nrow=length(mse3))
  f.versus<-as.data.frame(f.versus)
  names(f.versus)<-c("Numerator","Denominator","Num.Df","Den.Df")
  rownames(f.versus)<-names(mse3)
  f.versus$Numerator<-names(mse3) # recria a tabela f.versus

  for (i in 1:length(mse3)){ # constroi a coluna de denominadores
    a<-(mse3[i])
    if (length(a[[1]])==2){
      f.versus$Denominator[i]<-rownames(f.versus)[nrow(f.versus)]
      next
    }
    for (j in length(mse3):1){
      b<-(mse3[j])
      if (length(b[[1]]) >= length(a[[1]]))next
      ms1<-c(a[[1]])
      ms2<-c(b[[1]])

      c<-paste(ms1[-length(ms1)],collapse=" + ")
      d<-paste(ms2, collapse=" + ")

      if (c==d){
        f.versus$Denominator[i]<-names(b)
      }
    }
  }

  for (z in 1:length(mse3)){
    if (is.na(f.versus$Denominator)[z]){
      f.versus$Denominator[z]<-"No test"
    } # se nao encontrar um denominador adequado eh 'no test' pro termo

    if(quasi.f==TRUE){ # se quasi.f for TRUE faz o mesmo procedimento da funcao estimates para encontrar uma combinacao de MS adequada
      if (f.versus$Denominator[z] == "No test"){
        num1<-mse3[z]
        for (j in nrow(f.versus):z){
          num2<-mse3[j]
          for (k in nrow(f.versus):z){
            den1<-(mse3[k])
            for (l in nrow(f.versus):z){
              den2<-(mse3[l])
              if (names(den1) == names(den2))next
              if(length(num1[[1]]) != 1){
                if (((den1[[1]][length(den1[[1]])]) != (num1[[1]][length(num1[[1]])-1])) & ((den2[[1]][length(den2[[1]])]) != (num1[[1]][length(num1[[1]])-1])))next
              }
              ms1<-c(num2[[1]], num1[[1]][-length(num1[[1]])])# junta todos os elementos do numerador, menos o proprio termo
              ms2<-c(den2[[1]], den1[[1]])# junta todos os elementos do denominador

              sort1<-sort(ms1) # ordena os termos
              sort2<-sort(ms2) # ordena os termos
              if(!identical(sort1,sort2))next # se forem diferentes pula pra proxima interacao
              if(identical(sort1,sort2)){ # SEEE os dois modelos forem iguais
                f.versus$Numerator[z]<-paste(names(num2),names(num1),sep =" + ")
                f.versus$Denominator[z]<-paste(names(den2),names(den1),sep =" + ")
                break # para tudo
              }
            }
          }
        }
      }
    }
  } # o negocio eh encardido...
  mse4<-mse3
  for (i in 1:length(mse3)){
    for (j in 1:length(mse3[[i]])){
      mse4[[i]][j]<-(nomes[,2][which(nomes[,1] == mse3[[i]][j])])
    }
  }# substitui os numeros pelos termos originais

  if (method == "pool"){ # aqui vai substituir o nome do termo que recebeu os 'df' e 'ss' por 'Pooledx'
    for (i in 1:length(mse3)){
      for (j in 1:length(mse3[[i]])){
        if (mse4[[i]][j] == names(pooled)){
          mse4[[i]][j]<-nome
#          y<-rownames(mse)[which(grepl(pattern="Pooled", x = rownames(mse)))]
#          if(length(y) != 0){
#            a<-(gsub(pattern = "Pooled", replacement = "",  x = y))
#            nome<-paste("Pooled",c(as.numeric(max(a))+1),sep="")
#            mse4[[i]][j]<-nome
#          }
#          if (length(y) == 0){
#            nome<-paste("Pooled",1,sep ="")
#            mse4[[i]][j]<-nome
#          }
        }
      }
    }
    } # ta quase na metade....

  mse5 <- matrix(ncol = 1, nrow = length(mse4))
  colnames(mse5) <- "Mean square estimates"
  rownames(mse5) <- names(mse4)
  for (i in 1:nrow(mse5)) {
    a<-mse4[[i]]
    for (j in 1:length(mse4[[i]])){
      a[j]<-paste(mse21[[i]][j],"*",mse4[[i]][j],sep="")
    }
    mse5[i,1] <- paste(a, collapse = " + ")
  } # junta novamente os termos com seus respectivos numeros de niveis

  if (method == "pool"){ # cria a tabela 'pool.tab' que vai armazenar as informacoes
    pool.tab<-matrix(nrow=1, ncol=2) # de quem foi removido, foi juntado com qual termo
    pool.tab[1,1]<-termo # quem foi eliminado, etc.
    pool.tab[1,2]<-names(pooled)

    if (!is.null(tab$pool.tab)){ # se a tabela 'pool.tab' jah existir vai ser adicionada uma nova linha
      pool.tab<-rbind(tab$pool.tab,pool.tab)
      rownames(pool.tab)[nrow(pool.tab)]<-nome
    }
    if (is.null(tab$pool.tab)){
      rownames(pool.tab)<-"Pooled1"
      colnames(pool.tab)<-c("Term", "Pooled with")
      }
    }

  if (method == "eliminate"){
    pool.tab<-matrix(nrow=1, ncol=2)
    pool.tab[1,1]<-termo
    pool.tab[1,2]<-"Residuals"

    if (!is.null(tab$pool.tab)){
      pool.tab<-rbind(tab$pool.tab,pool.tab)
      rownames(pool.tab)[nrow(pool.tab)]<-"Eliminated"
    }
    if (is.null(tab$pool.tab)){
      rownames(pool.tab)<-"Eliminated"
      colnames(pool.tab)<-c("Term", "Pooled with")
      }
    } # aqui finaliza a tabela 'pool.tab'

  if (is.na(anova.res$Pr[which(rownames(anova.res) == termo)])){
    stop ("You can not apply post hoc pooling in a term which does not has an appropriate F-ratio")
  } # identifica se o pooling estaria sendo realizado num termo sem F versus apropriado

  nota<-NA
  if (anova.res$Pr[which(rownames(anova.res) == termo)] < 0.25){
    texto<-"has p value < 0.25. You should reconsider removing this term of the analysis"
    nota<-paste(termo, texto, sep=" ")
  } # verifica se o p valor do termo selecionado eh inferior a 0.25, se for vai dar um aviso

  table <- anova.res[, 1:3]
  table2<-table
  table2<-table2[which(rownames(table2) != termo),]

  if (method == "pool"){ # aqui vai somar 'df' e 'SS' do termo selecionado ao F versus apropriado (que passa a se chamar 'Pooledx')
    rownames(table2)[which(rownames(table2) == pool.tab[nrow(pool.tab),2])] <- rownames(pool.tab)[nrow(pool.tab)]

    ss.termo<-table$Sum[which(rownames(table)==pool.tab[nrow(pool.tab),1])]
    ss.pooled<-table$Sum[which(rownames(table)==pool.tab[nrow(pool.tab),2])]
    df.termo<-table$Df[which(rownames(table)==pool.tab[nrow(pool.tab),1])]
    df.pooled<-table$Df[which(rownames(table)==pool.tab[nrow(pool.tab),2])]
    SSpooled<-(ss.termo + ss.pooled)
    Dfpooled<-(df.termo + df.pooled)
    MSpooled<-( SSpooled/Dfpooled )

    table2[rownames(pool.tab)[nrow(pool.tab)],"Df"]<-Dfpooled
    table2[rownames(pool.tab)[nrow(pool.tab)],"Sum Sq"]<-SSpooled
    table2[rownames(pool.tab)[nrow(pool.tab)],"Mean Sq"]<-MSpooled
    }
  if (method == "eliminate"){ # se o metodo for de eliminacao o 'df' e 'SS' vai para o residuo
    ss.termo<-table$Sum[which(rownames(table) == termo)]
    ss.res<-table2$Sum[nrow(table2)]
    df.termo<-table$Df[which(rownames(table) == termo)]
    df.res<-table2$Df[nrow(table2)]
    SSres<-(ss.termo + ss.res)
    Dfres<-(df.termo + df.res)
    MSres<-( SSres/Dfres )

    table2[nrow(table2),"Df"]<-Dfres
    table2[nrow(table2),"Sum Sq"]<-SSres
    table2[nrow(table2),"Mean Sq"]<-MSres
    }

  f.versus<-f.versus[-nrow(f.versus),]

  f <- numeric(nrow(table2)) # F value
  for (i in 1:nrow(f.versus)){
    if(quasi.f==TRUE){ # mesmo procedimento do gad para qdo quasi.f=TRUE
      teste<-strsplit(f.versus$Denominator[i], " \\+ ")
      if (length(teste[[1]]) > 1){
        num<-strsplit(f.versus$Numerator[i], " \\+ ")
        den<-strsplit(f.versus$Denominator[i], " \\+ ")

        ms.num<-((table2$Mean[which(rownames(table2) == num[[1]][1])])  + (table2$Mean[which(rownames(table2) == num[[1]][2])]))
        num.1<- (((table2$Mean[which(rownames(table2) == num[[1]][1])])  + (table2$Mean[which(rownames(table2) == num[[1]][2])]))^2)
        num.2<- (((table2$Mean[which(rownames(table2) == num[[1]][1])])^2) / table2$Df[which(rownames(table2) == num[[1]][1])])
        num.3<- (((table2$Mean[which(rownames(table2) == num[[1]][2])])^2) / table2$Df[which(rownames(table2) == num[[1]][2])])
        df.num<- ( num.1 / (num.2 + num.3) )

        ms.den<-((table2$Mean[which(rownames(table2) == den[[1]][1])])  + (table2$Mean[which(rownames(table2) == den[[1]][2])]))
        den.1<- (((table2$Mean[which(rownames(table2) == den[[1]][1])])  + (table2$Mean[which(rownames(table2) == den[[1]][2])]))^2)
        den.2<- (((table2$Mean[which(rownames(table2) == den[[1]][1])])^2) / table2$Df[which(rownames(table2) == den[[1]][1])])
        den.3<- (((table2$Mean[which(rownames(table2) == den[[1]][2])])^2) / table2$Df[which(rownames(table2) == den[[1]][2])])
        df.den<- ( den.1 / (den.2 + den.3) )

        f.versus$Num.Df[i]<-round(df.num, digits = 2)
        f.versus$Den.Df[i]<-round(df.den, digits = 2)

        f[i]<- (ms.num/ms.den)
      }
      if (f.versus$Denominator[i] == names(mse3[length(mse3)])){
        f[i] <- table2$Mean[i]/table2$Mean[nrow(table2)]
      }# se o denominador for o residuo

      if ((length(teste[[1]]) == 1) & (f.versus$Denominator[i] != names(mse3[length(mse3)]))){
        f[i] <- table2$Mean[i]/table2$Mean[which(rownames(table2) == f.versus$Denominator[i])]
      }# quando o denominador tiver o F versus apropriado
    }

    if(quasi.f == FALSE){ # quando quasi.f=FALSE
      if (f.versus$Denominator[i] == "No test"){
        f[i] <- NA
      }
      if (f.versus$Denominator[i] == names(mse3[length(mse3)])){
        f[i] <- table2$Mean[i]/table2$Mean[nrow(table2)]
      }
      if ((f.versus$Denominator[i] != names(mse3[length(mse3)])) & (f.versus$Denominator[i] != "No test")){
        f[i] <- table2$Mean[i]/table2$Mean[which(rownames(table2) == f.versus$Denominator[i])]
      }
    }
  }

  P <- numeric(nrow(table2)) # P value
  for (i in 1:nrow(f.versus)){
    if(quasi.f==TRUE){ # mesmo procedimento do gad
      teste<-strsplit(f.versus$Denominator[i], " \\+ ")
      if (length(teste[[1]]) > 1){
        P[i] <- pf(as.numeric(f[i]), f.versus$Num.Df[i], f.versus$Den.Df[i], lower.tail = FALSE)
      }
      if (f.versus$Denominator[i] == names(mse3[length(mse3)])){
        P[i] <- pf(as.numeric(f[i]), table2$Df[i], table2$Df[which(rownames(table2) == f.versus$Denominator[i])],lower.tail = FALSE)
        f.versus$Num.Df[i]<-signif(table2$Df[i],3)
        f.versus$Den.Df[i]<-signif(table2$Df[which(rownames(table2) == f.versus$Denominator[i])],3)
      }
      if ((length(teste[[1]]) == 1) & (f.versus$Denominator[i] != names(mse3[length(mse3)]))){
        P[i] <- pf(as.numeric(f[i]), table2$Df[i], table2$Df[which(rownames(table2) == f.versus$Denominator[i])], lower.tail = FALSE)
        f.versus$Num.Df[i]<-signif(table2$Df[i],3)
        f.versus$Den.Df[i]<-signif(table2$Df[which(rownames(table2) == f.versus$Denominator[i])],3)
      }
    }

    if(quasi.f == FALSE){
      if (f.versus$Denominator[i] == "No test"){
        P[i] <- NA
      }
      if (f.versus$Denominator[i] == names(mse3[length(mse3)])){
        P[i] <- pf(as.numeric(f[i]), table2$Df[i], object$df.residual,lower.tail = FALSE)
      }
      if ((f.versus$Denominator[i] != names(mse3[length(mse3)])) & (f.versus$Denominator[i] != "No test")){
        P[i] <- pf(as.numeric(f[i]), table2$Df[i], table2$Df[which(rownames(table2) == f.versus$Denominator[i])], lower.tail = FALSE)
      }
    }
  }

  anova.table <- data.frame(table2, f, P) # junta tudo em uma nova tabela de anova
  anova.table[length(P), 4:5] <- NA
  colnames(anova.table) <- colnames(anova.res)
  rownames(anova.table) <- rownames(table2)
  if (quasi.f == FALSE){
    f.versus<-f.versus[,1:2] # se quasi.f=FALSE nao faz sentido mostrar os graus de liberdade da tabela de F versus
    tab.anova<-structure(anova.table, heading = c("Analysis of Variance Table\n",
                                                  paste("Response:", deparse(formula(object)[[2L]]))),
                         class = c("anova", "data.frame"))
  }
  if (quasi.f == TRUE){ # quando quasi.f=TRUE o usuario sera notificado sobre quais termos tem sua razao F aproximada
    aviso<-character(length=nrow(f.versus))
    for(i in 1:nrow(f.versus)){
      teste<-strsplit(f.versus$Denominator[i], " \\+ ")
      if (length(teste[[1]]) > 1){
        aviso[i]<-rownames(f.versus)[i]
      }
    }
    aviso<-aviso[which(aviso != "")]
    if (length(aviso) == 1){
      aviso.final<- paste("The F-ratio of", aviso, "is approximate. See 'help(gad)' for further details\n")
    }
    if (length(aviso) > 1){
      aviso1<-paste(aviso[1:(length(aviso)-1)], collapse=", ")
      aviso.final<- paste("The F-ratios of", aviso1, "and", aviso[length(aviso)], "are approximations. See 'help(gad)' for further details\n")
    }
    if (length(aviso) == 0){
      aviso.final<- "All of your terms has an appropriated F-ratio, Quasi F-test is no longer necessary\n"
    }
    tab.anova<-structure(anova.table, heading = c("Analysis of Variance Table: Quasi F-ratios\n",
                                                  aviso.final,
                                                  paste("Response:", deparse(formula(object)[[2L]]))),
                         class = c("anova", "data.frame"))
  }
  if (is.na(nota)){ # se o p valor do termo removido for menor que 0.25
    result <- list(pool.table = pool.tab, mse = mse5, f.versus = f.versus, anova = tab.anova)
  }
  else{ # se for maior nao tem nota de aviso
    result <- list(pool.table = pool.tab, mse = mse5, f.versus = f.versus, anova = tab.anova, note = nota)
  }
  return(result)
}
