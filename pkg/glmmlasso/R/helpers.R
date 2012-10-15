ObjFctPirls <-
function(pirls,beta,weights,penalized,lambda){
  2*pirls$Qu + sum(2*determinant(pirls$Chol)$modulus) +
    lambda*sum(abs(beta[penalized,drop=FALSE])/weights[penalized])
}

armijo <-
function(data,beta,theta,pirls,k,beta_k,lambda_k,direction,deriv,lambda,weights,
         fct0,penalized,Wsqrt,convergence=convergence,tZL,pen,
         control,minArmijo,Xb0,family,Chol)
{
   
  # Initialisierungen vor den Schritten
  
  eta0 <- pirls$eta
  u.cp <- pirls$cputilde
  beta.new <- beta
  lsabeta0 <- lambda*sum(abs(beta[penalized,drop=FALSE])/weights[penalized])
  increment <- direction*deriv$gradient + control$gamm*direction^2*
    deriv$hessian + pen*lambda_k*(abs(beta_k+direction)-abs(beta_k))
  tXk <- data$X[k,]
    
  # Verändern der Schrittweite
  for (l in minArmijo:control$maxArmijo)
    {
      alpha <- control$a_init*control$delta^l
      beta.newk <- beta_k + alpha*direction
 
      delta <- tXk*alpha*direction
      eta.new <- eta0 + delta
      lsabeta <- lsabeta0 - pen*lambda_k*(abs(beta_k) - abs(beta.newk))
      fct1 <- armijoFct(y=data$y,tZL=tZL,u.cp=u.cp,Wsqrt=Wsqrt,eta=eta.new,
                        lsabeta=lsabeta,family=family,Chol=Chol)

      add <- alpha*control$rho*increment
      if (is.na(fct1)) cat("l=",l,"fct1",fct1,"fct0",fct0,"add",add,"\n")
      
      if (fct1 <= fct0 + add)
       {
       	 beta.new[k] <- beta.newk
         beta <- beta.new
	 Xb <- Xb0 + delta
         fct <- fct1
         break
       }
       if (l==control$maxArmijo)
        {
          convergence <- convergence+2
          #cat("Armijo not successful","\n")
	  Xb <- Xb0
          fct <- fct0
        }
    }
   #cat("l=",l,"\n")

  return(list(beta=beta,fct=sum(fct),l=l,convergence=convergence,Xb=Xb,
              u=pirls$utilde))
  
}

armijoExact <-
function(data,tX,beta,theta,u,k,beta_k,direction,lambda_k,deriv,lambda,weights,
         fct0,penalized,Wsqrt,convergence=convergence,tZL,pen,
         control,minArmijo,Xb0,family,Chol)
{

if (minArmijo!=control$maxArmijo)
 {
   #cat("1...\n")
  # Initialisierungen vor den Schritten
  beta.new <- beta
  lsabeta <- lambda*sum(abs(beta[penalized,drop=FALSE])/weights[penalized])
  increment <- direction*deriv$gradient + control$gamm*direction^2*
    deriv$hessian + pen*lambda_k*(abs(beta_k+direction)-abs(beta_k))
  #cat("increment:",increment,"\n")
  tXk <- data$X[k,]
  u.init <- u

  # Verändern der Schrittweite
  for (l in minArmijo:control$maxArmijo)
    {
      #cat("l=",l,"\n")
      # new beta_k
      alpha <- control$a_init*control$delta^l
      beta.newk <- beta_k + alpha*direction
      # new Xbeta
      delta <- tXk*alpha*direction
      Xb.new <- Xb0 + delta
      beta.new[k] <- beta.newk
      # new u
      pirls.new <- pirls(u=u.init,data=data,beta=beta.new,Wsqrt=Wsqrt,
                         Xb=Xb.new,tZL=tZL,family=family,Chol=Chol,
                         tol=control$tol4)
      u.new <- u.init <- pirls.new$utilde
      u.new.cp <- pirls.new$cputilde
      # new eta
      eta.new <- pirls.new$eta
      lsabeta <- lsabeta - pen*lambda_k*(abs(beta_k) - abs(beta.newk)) 
      fct1 <- armijoFctExact(pirls=pirls.new,lsabeta=lsabeta)

      add <- alpha*control$rho*increment      
      if (fct1 <= fct0 + add)
       {
         beta <- beta.new
	 Xb <- Xb.new
         u <- u.new
         fct <- fct1
         break
       }
       if (l==control$maxArmijo)
        {
          convergence <- convergence+2
          #cat("Armijo not successful","\n")
	  Xb <- Xb0
          fct <- fct0
        }
    }
  #cat("l=",l,"\n")
 } else { convergence <- convergence+2 ; Xb <- Xb0 ; fct <- fct0 ;
          l <- minArmijo}
  
  return(list(beta=beta,fct=sum(fct),l=l,convergence=convergence,Xb=Xb,u=u))
  
}

armijoFct <-
function(y,tZL,u.cp,Wsqrt,eta,lsabeta,family,Chol)
 {

  if (family=="binomial")
   {
     exp.eta <- exp(eta)
     Wsqrt@x <- sqrt(exp.eta)/(1+exp.eta)
   }
  if (family=="poisson")
   {
     exp.eta <- exp(eta)
     Wsqrt@x <- sqrt(exp.eta)
   }
  
  C1 <- Cholesky(tcrossprod(tcrossprod(tZL,Wsqrt)),Imult=1,perm=FALSE)

  #C1 <- update(Chol,tcrossprod(tcrossprod(tZL,Wsqrt)),mult=1,perm=FALSE)
  
  log.det.L2 <- 2*determinant(C1)$modulus

  fct <- 2*glmPart(y=y,eta=eta,family=family) + log.det.L2 + u.cp + lsabeta
  
  return(sum(fct))
}

armijoFctExact <-
function(pirls,lsabeta)
 {
   
   2*pirls$Qu + sum(2*determinant(pirls$Chol)$modulus) + lsabeta 
     
 }

derivBeta <-
function(Xk,Xk2,y,mu,tZL,lower,upper,W1k,unit,Chol,family,CInv=NULL,diagIndex)
 {

  if (family=="binomial")
    {
      mu1 <- mu*(1-mu)
      W1k@x <- Xk*mu1*(1-2*mu)
    }
  if (family=="poisson")
    {
      mu1 <- mu
      W1k@x <- Xk*mu
    }

  #M1 <- solve(Chol,unit)
  M2 <- tcrossprod(tcrossprod(tZL,W1k),tZL)
  M3 <- CInv%*%M2

  gradient <- -2*(crossprod(Xk,(y-mu))) +  sum(M3@x[diagIndex])

  hessian <- min(max(2*crossprod(Xk2,mu1),lower),upper)
  
  return(list(gradient=gradient,hessian=hessian))
  
 }

derivBetaExact <-
function(Xk,Xk2,y,mu,tZL,lower,upper,W1k,W2k,Wu,unit,Chol,family,CInv=NULL,
         diagIndex,u,diag2Index)
 {
  if (family=="binomial")
    {
      mu1 <- mu*(1-mu)
      Wu@x <- sqrt(mu1)
      xkmu <- Xk*mu1
      mu2 <- mu1*(1-2*mu)
    }
  if (family=="poisson")
    {
      mu1 <- mu
      Wu@x <- sqrt(mu1)
      xkmu <- Xk*mu1
      mu2 <- mu
    }
  
  gradFu <- tcrossprod(tcrossprod(tZL,Wu)) + unit
  
  gradFb <- tZL%*%xkmu
  gradub <- -1*solve(gradFu,unit)%*%gradFb 
  
  tZLgradub <- crossprod(tZL,gradub)@x
  XktZLgradub <- Xk + tZLgradub
  W1k@x <- XktZLgradub*mu2
 

  M2 <- tcrossprod(tcrossprod(tZL,W1k),tZL)
  M3 <- CInv%*%M2

  gradient <-  -2*sum(crossprod(XktZLgradub,(y-mu))) +
    sum(M3@x[diagIndex]) + 2*sum(crossprod(u,gradub))
  
  hessian <- min(max(2*crossprod(Xk2,mu1),lower),upper)
  
  return(list(gradient=gradient,hessian=hessian))
  
 }

desDir <-
function(deriv,lambda_k,beta_k,penalized_k) {
  gradient <- deriv$gradient
  hessian <- deriv$hessian
  if (penalized_k) desdir <- median(c((lambda_k-gradient)/hessian,-beta_k,
                                      (-lambda_k-gradient)/hessian))
  else desdir <- -gradient/hessian

  return(desdir)
}

devianceFct <-
function(y,mu,u,Wsqrt,tZL,fct,family,lambda,beta,weights,penalized)
{
  if (family=="binomial")
   {
     deviance <- fct - lambda*sum(abs(beta[penalized,drop=FALSE])/
                                  weights[penalized])
   }
  
  if (family=="poisson")
   {
     Wsqrt@x <- sqrt(mu)

     C <- Cholesky(tcrossprod(tcrossprod(tZL,Wsqrt)),Imult=1,perm=FALSE)
     ldL2 <-  2*determinant(C)$modulus
     
     usqr <-  crossprod(u)
  
     index <- y!=0
  
     disc <- 2*sum(y[index]*log(y[index]/mu[index]))-2*sum(y-mu)

     deviance <- sum(disc + ldL2 + usqr)
   }
  
  return(deviance)
}

fittedValues <-
function(eta,family)
 {
   if (family=="binomial")
    {
      mu <- exp(eta)/(1+exp(eta))
      fitted <- as.numeric(eta > 0)
      
    }
   if (family=="poisson")
    {
      fitted <- mu <- exp(eta)
    }
  return(list(mu=mu,fitted=fitted))
 }

glmPart <-
function(y,eta,family)
 {
  if (family=="binomial") return(-sum(y*eta-log(1+exp(eta))))
  
  if (family=="poisson")  return(sum(lfactorial(y))-sum(y*eta-exp(eta)))
   
 }

pirls <-
function(u,data,beta,tol,Wsqrt,theta=NULL,Xb=NULL,tZL=NULL,family,inverse=FALSE,
         unit=NULL,CInv=NULL,Chol)
 {
   if (missing(Xb)) Xb <- t(data$X)%*%beta 
   if (missing(tZL)) tZL <- theta*data$Z
   y <- data$y
   W <- G <- Wsqrt
   uStart <- u
   
   counter <- 0
   #cat("...","\n")

   cptZLu <- crossprod(tZL,u)@x
   eta.new <- Xb + cptZLu
   
   repeat 
    {
      counter <- counter+1
      eta <- eta.old <- eta.new
      if (family=="binomial")
       { 
         exp.eta <- exp(eta)
         exp.eta.1 <- 1+exp.eta
         mu <- exp.eta/exp.eta.1
         W@x <- mu/exp.eta.1
         Wsqrt@x <- sqrt(exp.eta)/exp.eta.1
         G@x <- (exp.eta.1*exp.eta.1)/exp.eta
       }
      if (family=="poisson")
       {
         mu <- exp(eta)
         W@x <- mu
         Wsqrt@x <- sqrt(mu)
         G@x <- 1/mu
       }

      z <- cptZLu + G%*%(y-mu) # das ist wahrscheinlich extrem teuer
      
      C <- Cholesky(tcrossprod(tcrossprod(tZL,Wsqrt)),Imult=1,perm=FALSE)
      #C <- update(Chol,tcrossprod(tcrossprod(tZL,Wsqrt)),mult=1,perm=FALSE)
      if (!inverse) u <- solve(C,tcrossprod(tZL,W)%*%z)
      else {CInv <- solve(C,unit) ; u <- CInv%*%tcrossprod(tZL,W)%*%z }

      cptZLu <- crossprod(tZL,u)@x

      eta.new <- Xb + cptZLu
 
      if (sum(sqrt(abs(crossprod(eta.old-eta.new)/crossprod(eta.old))))<
          tol|(counter>50)) break
    }

   if (counter>50) u <- uStart

   cpu <- sum(crossprod(u))
   
   Qu <- glmPart(y=data$y,eta=eta.new,family=family) + 0.5*cpu
   ##Qu <- -sum(y*eta.new-log(1+exp(eta.new))) + 0.5*cpu

   if (family=="binomial")
    {
      exp.eta.new <- exp(eta.new)
      mu <- exp.eta.new/(1+exp.eta.new)
    }
   if (family=="poisson")
    {
      mu <- exp(eta.new)
    }

   #cat("counter:",counter,"\n")

   return(list(utilde=u,Xb=Xb,tZL=tZL,Qu=Qu,Chol=C,mu=mu,eta=eta.new,
               Wsqrt=Wsqrt,cputilde=cpu,CInv=CInv))
   
 }

## NEW ##############################################################
thetaFctExactpdSym <- function(thetal,data,Xb,beta,tZ,u,add,Wsqrt,L,index,N,
                               family,tol4,Chol=NULL)
  ## L:      Wurzel der Kovarianzmatrix in Vektorform
  ## index:  zu aenderndes Element der Matrix/Vektors
  ## thetal: Wert des Elements (1-dim)
 {
  ## Wurzel der Kovarianzmatrix ändern

  L@x[index] <- thetal 

  tZL <- crossprod(L,tZ) ## Multiplikation mit L, check Transponierung...??????
  
  pirlsOutput <- pirls(u=u,data=data,beta=beta,Wsqrt=Wsqrt,Xb=Xb,tZL=tZL,
                       family=family,Chol=Chol,tol=tol4)
  
  log.det.L2 <- 2*determinant(pirlsOutput$Chol)$modulus

  fct <- 2*pirlsOutput$Qu + log.det.L2 + add
}
###################################################################

thetaFctExactpdDiag <-
function(thetal,data,Xb,beta,tZ,u,add,Wsqrt,L,index,N,family,tol4,Chol=NULL)
 {
  L@x[index] <- rep(thetal,N)
  tZL <- crossprod(L,tZ)

  pirlsOutput <- pirls(u=u,data=data,beta=beta,Wsqrt=Wsqrt,Xb=Xb,tZL=tZL,
                       family=family,Chol=Chol,tol=tol4)
  
  log.det.L2 <- 2*determinant(pirlsOutput$Chol)$modulus # 0.001

  fct <- 2*pirlsOutput$Qu + log.det.L2 + add
}

thetaFctExactpdIdent <-
function(theta,data,Xb,beta,u,add,Wsqrt,family,tol4,Chol=NULL)
 { 
  tZL <- theta*data$Z
  pirlsOutput <- pirls(u=u,data=data,beta=beta,Wsqrt=Wsqrt,Xb=Xb,tZL=tZL,
                       family=family,Chol=Chol,tol=tol4)

  log.det.L2 <- 2*determinant(pirlsOutput$Chol)$modulus

  fct <- 2*pirlsOutput$Qu + log.det.L2 + add
}


thetaOptExactpdSym <-
function(theta,data,beta,u,lambda,weights,penalized,Wsqrt,Xb,family,L,stot,
         control,ind,N,Chol=NULL)
{

  add <- lambda*sum(abs(beta[penalized,drop=FALSE])/weights[penalized])
  for (s in 1:length(theta)){
      if (control$CovOpt=="nlminb")
        {
          optRes <- nlminb(theta[s],thetaFctExactpdSym,data=data,Xb=Xb,
                           beta=beta,tZ=data$Z,u=u,add=add,
                           Wsqrt=Wsqrt,L=L,index=ind[[s]],N=N,
                           family=family,Chol=Chol,tol4=control$tol4)
          if (abs(optRes$par)<control$thres)
            {
              theta[s] <- 0
              L@x[ind[[s]]] <- 0
            } else {
              theta[s] <- optRes$par
              L@x[ind[[s]]] <- optRes$par
            }
        }
    }
  tZL <- crossprod(L,data$Z)
  
  return(list(theta=theta,fct=optRes$objective,tZL=tZL,L=L))
}


thetaOptExactpdDiag <-
function(theta,data,beta,u,lambda,weights,penalized,Wsqrt,Xb,family,L,stot,
         control,ind,N,Chol=NULL)
{

  add <- lambda*sum(abs(beta[penalized,drop=FALSE])/weights[penalized])
  
  for (s in 1:stot)
   {
     index <- which(ind==s)
     if (control$CovOpt=="nlminb")
       {
         optRes <- nlminb(theta[s],thetaFctExactpdDiag,data=data,Xb=Xb,
                          beta=beta,tZ=data$Z,u=u,add=add,
                           Wsqrt=Wsqrt,L=L,index=index,N=N,family=family,
                          Chol=Chol,tol4=control$tol4,lower=0)
         if (abs(optRes$par)<control$thres)
         {
           theta[s] <- 0
           L@x[index] <-  rep(theta[s],N)
         } else {theta[s] <- optRes$par ; L@x[index] <-  rep(theta[s],N)}
       }

   }
  
  tZL <- crossprod(L,data$Z)
  
  return(list(theta=theta,fct=optRes$objective,tZL=tZL,L=L))
}

thetaOptExactpdIdent <-
function(theta,data,beta,u,lambda,weights,fctstart,penalized,Wsqrt,Xb,
         family,Chol=NULL,control)
{

  add <- lambda*sum(abs(beta[penalized,drop=FALSE])/weights[penalized])

  if (control$CovOpt=="nlminb")
   {
     optRes <- nlminb(theta,thetaFctExactpdIdent,data=data,Xb=Xb,beta=beta,
                      u=u,add=add,Wsqrt=Wsqrt,family=family,Chol=Chol,
                      tol4=control$tol4,lower=0)
     theta <- optRes$par
   }

  #if (optRes$objective>fctstart) cat("Non-decreasing objective function!","\n")
  
  tZL <- theta*data$Z
  
  return(list(theta=theta,fct=sum(optRes$objective),tZL=tZL))
  
}

fill.mat <- function(vec, dim)
{
  ## Purpose:
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Lukas Meier, Date:  9 Jul 2012, 10:05

  m <- matrix(0, nrow = dim, ncol = dim)
  diag(m) <- vec[1:dim]
  count <- dim + 1
  ## fill row-wise, needs to be done more elegant...
  for(i in 1:(dim-1)){
    for(j in (i+1):dim){
      m[i,j] <- vec[count]
      count <- count + 1
    }
  }
  m
}



