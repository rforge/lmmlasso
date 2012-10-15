summary.glmmlasso <-
function(object,...)
 {
   ## 1. part
    table1 <- round(c(object$aic,object$bic,object$logLik,object$deviance,
                      object$objective),1)
    names(table1) <- c("AIC","BIC","logLik","deviance","objective")
    cat("Generalized linear mixed model fit by the Laplace approximation","\n")
    cat("family =",object$family,"; ntot =",object$ntot,"; N =",object$N,
        "; p =",object$p,"; lambda =",round(object$lambda,3),"\n")
    print(table1)

    ## 2a. part
    cat("","\n")
    cat("Random effects:\n") 
    if(object$covStruct == "Sym"){
      m <- fill.mat(object$theta, object$stot)
      ##
      cov  <- tcrossprod(m)
      sd   <- sqrt(diag(cov)) ## standard deviations
      corr <- cov2cor(cov)    ## correlation matrix
      tableCov  <- round(cov,5)
      tableSd   <- round(sd,5)
      matrixCorr <- round(corr,5)
      ##colnames(matrixVar) <- c("Variance","Std.Dev.")
      cat("\nCovariance Matrix:\n")
      print(tableCov)
      cat("\nCorrelation Matrix:\n")
      print(matrixCorr)
      cat("\nStandard Deviations:\n")
      print(tableSd)
      out <- list(cov = tableCov, corr = matrixCorr, sd = tableSd)
    }else{
      tableVar <- round(c(object$theta^2),5)
      tableSd <- round(c(object$theta),5)
      matrixVar <- cbind(tableVar,tableSd)
      colnames(matrixVar) <- c("Variance","Std.Dev.")
      ##if (length(object$ranInd)>1) rownames(matrixVar) <-
      ##c("(Intercept)",paste("X",object$ranInd[-1],sep=""))  
      ##else rownames(matrixVar) <- c("(Intercept)")
      print(matrixVar)
      out <- list(var = matrixVar)
    }
    
    ## 3. part
    cat("","\n")
    cat("Fixed effects:","\n")
    penalty <- rep("",length(object$coefficients))
    penalty[object$unpenalized] <- "(n)"
    table3 <- data.frame(round(object$coefficients,5),penalty)

    xnames <- rownames(object$data$X)

    if (is.null(xnames)){
      rownames(table3) <- c("(Intercept)",paste("X", 2:length(object$fixef),
                                                    sep=""))
    }else{
      rownames(table3) <- c("(Intercept)",rownames(object$data$X)[-1]) 	
    }
    
    table3 <- table3[object$coefficients!=0,]
    colnames(table3) <- c("Estimate"," ")
    if (dim(table3)[1]==1) rownames(table3) <- c("(Intercept)")
    cat("|active set|=",sum(object$fixef!=0),"\n")
    print(table3)

    ## 4th part
    cat("","\n")
    if (object$family=="binomial") cat("Misclassification Error:", object$gof,
          "\n")
    if (object$family=="poisson") cat("Lack-of-fit:",object$gof,"\n")
    cat("Number of iterations:",object$nIter,"\n")

    invisible(out)
 }

print.glmmlasso <-
function(x,...)
  {

    ## 1. part
    table1 <- round(c(x$aic,x$bic,x$logLik,x$deviance,x$objective),1)
    names(table1) <- c("AIC","BIC","logLik","deviance","objective")
    cat("Generalized linear mixed model fit by the Laplace approximation","\n")
    cat("family =",x$family,"; ntot =",x$ntot,"; N =",x$N,"; p =",x$p,
        "; lambda =",round(x$lambda,3),"\n")
    print(table1)

    ## 2a. part
    cat("","\n")
    cat("Random effects:\n")
    if(x$covStruct == "Sym"){
      m <- fill.mat(x$theta, x$stot)
      ##
      cov  <- tcrossprod(m)
      sd   <- sqrt(diag(cov)) ## standard deviations
      corr <- cov2cor(cov)    ## correlation matrix
      tableCov  <- round(cov,5)
      tableSd   <- round(sd,5)
      matrixCorr <- round(corr,5)
      ##colnames(matrixVar) <- c("Variance","Std.Dev.")
      cat("\nCovariance Matrix:\n")
      print(tableCov)
      cat("\nCorrelation Matrix:\n")
      print(matrixCorr)
      cat("\nStandard Deviations:\n")
      print(tableSd)
    }else{
      tableVar <- round(c(x$theta^2),5)
      tableSd <- round(c(x$theta),5)
      matrixVar <- cbind(tableVar,tableSd)
      colnames(matrixVar) <- c("Variance","Std.Dev.")
      if (length(x$ranInd)>1) rownames(matrixVar) <-
        c("(Intercept)",paste("X",x$ranInd[-1],sep=""))
      else rownames(matrixVar) <- c("(Intercept)")
      print(matrixVar)
    }
    
    ## 3. part
    cat("","\n")
    cat("Fixed effects:","\n")
    cat("|active set|=",sum(x$fixef!=0),"\n")
  }

plot.glmmlasso <-
function(x,...)
 {
   par(mfrow=c(2,2))

   ## Tukey-Anscombe plots
   plot(x$workResid~x$eta,col=x$data$group,main="Tukey-Anscombe Plot",
        xlab="eta",ylab="working residuals")
   abline(h=0,col="grey")

   if (x$family=="binomial")
    {
     plot(x$respResid~x$mu,col=x$data$group,main="Tukey-Anscombe Plot",
          xlab="mu",ylab="response residuals")
     abline(h=0,col="grey")
    } else
   {
     plot(x$data$y~x$fitted,col=x$data$group,main="Tukey-Anscombe Plot",
          xlab="Predicted values",ylab="Observed values")
   }

   ## QQ-plot of the random effects
   qqnorm(x$ranef,main="QQ-Plot of the random effects")
   qqline(x$ranef)

   ## histogram of the fixed effects
   fixEf <- as.vector(x$fixef)
   hist(fixEf,main=paste("|active set| = ",sum(fixEf!=0)),xlab="fixed effects")
   rug(fixEf,lwd=3)
 }


