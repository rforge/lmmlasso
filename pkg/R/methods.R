summary.lmmlasso <-
function(object,...)
 {
    if (object$stopped)
      {
        cat("Since the algorithm stopped due to","\n")
        cat("|active.set|>=min(p,ntot)","\n")
        cat("no summary is available.","\n")
        cat("Set stopSaturated=FALSE or increase lambda!","\n")
      } else {
    # 1. part
    t.1 <- round(c(object$aic,object$bic,object$logLik,object$deviance),1)
    names(t.1) <- c("AIC","BIC","logLik","deviance")
    cat("Model fitted by",object$method,"for","lambda =",round(object$lambda,3),":\n")
    print(t.1)

    # 2a. part
    cat("","\n")
    cat("Random effects:",object$pdMat,"\n")
    t.var <- c(diag(object$Psi),object$sigma^2)
    t.sd <- c(sqrt(diag(object$Psi)),object$sigma)
    m.var <- cbind(t.var,t.sd)
    colnames(m.var) <- c("Variance","Std.Dev.")
    if (dim(object$data$z)[[2]]>1) rownames(m.var) <- c("(Intercept)",paste("X", 2:dim(object$data$z)[[2]],sep=""),"Residual")
    else rownames(m.var) <-c("(Intercept)","Residual")
    print(m.var)

    if (object$pdMat=="pdSym")
      {
        cat("","\n")
        cat("Random effects: Correlations","\n")
        m.cor <- object$corPsi
        if (dim(object$data$z)[[2]]>1) colnames(m.cor) <- rownames(m.cor) <- c("(Intercept)",paste("X", 2:dim(object$data$z)[[2]],sep=""))
        print(m.cor)
      }
    
    # 3. part
    cat("","\n")
    cat("Fixed effects:","\n")
    penalty <- rep("",length(object$coefficients))
    penalty[object$nonpen] <- "(n)"
    t.3 <- data.frame(object$coefficients,penalty)

    xnames <- colnames(object$data$x)

    if (is.null(xnames))
     {
       rownames(t.3) <- c("(Intercept)",paste("X", 2:length(object$coefficients),sep=""))
     } else
     {
       rownames(t.3) <- c("(Intercept)",colnames(object$data$x)[-1]) 	
     }
    t.3 <- t.3[object$coefficients!=0,]
    colnames(t.3) <- c("Estimate"," ")
    cat("|active set|=",sum(object$coefficients!=0),"\n")
    print(t.3)

    # 4th part
    cat("","\n")
    cat("Number of iterations:",object$counter,"\n")
  }
 }

print.lmmlasso <-
function(x,...) {
    # 0. part
    print(x$call)
  
    # 2a. part
    cat("","\n")
    cat("Random effects:",x$pdMat,"\n")
    t.var <- c(diag(x$Psi),x$sigma^2)
    t.sd <- c(sqrt(diag(x$Psi)),x$sigma)
    m.var <- cbind(t.var,t.sd)
    colnames(m.var) <- c("Variance","Std.Dev.")
    if (dim(x$data$z)[[2]]>1) rownames(m.var) <- c("(Intercept)",paste("X", 2:dim(x$data$z)[[2]],sep=""),"Residual")
    else rownames(m.var) <-c("(Intercept)","Residual")
    print(m.var)
    
    # 3. part
    cat("","\n")
    cat("Fixed effects:","\n")
    cat("|active set|=",sum(x$coefficients!=0),"\n")

}

plot.lmmlasso <-
function(x,...)
{

 par(mfrow=c(3,2))

 # Tukey-Anscombe plot
 plot(x$residuals~x$fitted.values,col=x$data$grp,main="Tukey-Anscombe Plot",xlab="fitted values",ylab="raw residuals")
 abline(h=0,col="grey")

 # QQ-Plot of the residuals
 qqnorm(x$residuals,col=x$data$grp,main="QQ-Plot of the residuals")
 qqline(x$residuals)

 # QQ-plot of the random effects
 t.ranef <- unlist(x$random,recursive=FALSE)
 qqnorm(t.ranef,main="QQ-Plot of the random effects",col=rep(1:dim(x$data$z)[[2]],length(levels(x$data$grp))))
 qqline(t.ranef)

 # boxplots of the random effects
 boxplot(t.ranef~as.factor(rep(1:dim(x$data$z)[[2]],length(levels(x$data$grp)))),border=1:dim(x$data$z)[[2]],
        main="Random effects by effects",xlab="random effects index")
 abline(h=0,col="grey",lty=2)

 # fitted vs. measured y
 plot(x$fitted,x$data$y,xlab="fitted values",ylab="response variable",main="response ~ fitted")
 abline(a=0,b=1)

 # histogram of the fixed effects
 hist(x$coefficients,main=paste("|active set| = ",sum(x$coefficients!=0)),xlab="fixed effects")
 rug(x$coefficients,lwd=3)

}

