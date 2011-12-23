glmmlasso <- function(y,...)
  UseMethod("glmmlasso")

glmmlasso.default <- function(y,group,X,Z=NULL,family=c("binomial","poisson"),covStruct=c("Identity","Diagonal"),lambda,weights=NULL,
                              coefInit=NULL,exactStep=FALSE,exactDeriv=FALSE,random=NULL,unpenalized=1:stot,ranInd=1:stot,
                              control=glmmlassoControl(family=family),...)
{
  
# generate the Z matrix via formula term
if (missing(random)&missing(Z))
 {
   Z <- sparse.model.matrix(~-1+group)
   cat("A random-intercept model is fitted.","\n")
 }

if (!missing(random))
 {
  if (!is.character(random)) stop("random must be a character")
  Z <- t(glmer(as.formula(paste("y~1+",random)),data=data.frame(X,y=y,group=group),family=family,nAGQ=1)@Zt)
 }

# transpose and save the data in an appropriate way
data <- list() ; data$y <- y ; data$group <- group ; data$X <- t(X) ; data$Z <- t(Z)

# --- Introductory checks and reformulations ---
# ----------------------------------------------

family <- match.arg(family)
covStruct <- match.arg(covStruct)

# do some checks
if (missing(data)) stop("Missing values are not allowed.")
if (!is.matrix(data$X)) stop("x has to be a matrix")
if (!is(data$Z,"Matrix")) stop("Z has to be a sparse Matrix")
if (!is.numeric(y)) stop("y has to be of type 'numeric'")
if (ncol(data$X)!=length(y)) stop("Length of X and y differ.")
if (any(data$X[1,]!=rep(1,dim(data$X)[[2]]))) stop("first column is not the intercept")
if (length(levels(group))==1) stop("Only one group. No covariance parameters!")
if (length(lambda)!=1) stop("lambda can only be one number")
if (lambda<0) stop("regularization parameter must be non-negative")
if ((family=="binomial")&(any(y<0)|any(1<y)|(length(levels(factor(y)))>2))) stop("y values must be 0 or 1")
if ((family=="poisson")&(any(abs(y-round(y)) > .Machine$double.eps^0.5))) stop("y values must be integers")
if ((family=="poisson")&any(y<0)) stop("y values must be positive")

ptm <- proc.time()
p <- dim(data$X)[[1]]

# weights (NA: not penalized, 0: drop from the analysis)
if (missing(weights))
 {
   weights <- rep(1,p)
 } else
    {
      if (!length(weights)==p) stop("Weights vector has not length p")

      unpenalized <- which(is.na(weights))
      if (any(weights[-unpenalized]<0)) stop("Weights must be positive")
      remove <- weights==0
      remove[unpenalized] <- FALSE

      xnames <- rownames(data$X)
      if (is.null(xnames)) xnames <- rownames(data$X) <- c(paste("X", 1:p,sep=""))

      data$X <- data$X[!remove,,drop=FALSE]
      rownames(data$X) <- xnames[!remove]
      
      weights <- weights[!remove]
      unpenalized <- which(is.na(weights))
   }

# crucial allocations and definitions
group <- factor(group)
p <- dim(data$X)[[1]]
N <- length(levels(group))
q <- data$Z@Dim[1] 
ntot <- length(y)
stot <- q/N

tX2 <- data$X^2

penalized <- c(1:p)[-unpenalized]
if (!(1%in%unpenalized)) cat("Intercept is penalized?!?","\n")

if (control$verbose>=2) control$fctSave <- TRUE

Wsqrt <- .symDiagonal(n=ntot,x=seq(1,ntot))
unit <- .symDiagonal(n=q)

LStart <- .symDiagonal(n=q,x=rep.int(1,q))
ind <- rep(1:stot,each=N)

# --- Determine or calculate the starting values ---
# --------------------------------------------------

if (missing(coefInit))
 {
  # ... for the fixed effects (using cv-optimal tuning parameter \lambda)
  set.seed(control$seed)
  initCv <- cv.glmnet(x=t(data$X)[,-1,drop=FALSE],y=y,family=family,alpha=1) # weights berucksichtigen !!!
  initGlmnet <- glmnet(x=t(data$X)[,-1,drop=FALSE],y=y,family=family,alpha=1)
  betaStart <- as.vector(predict(initGlmnet,s=initCv$lambda.1se,type="coef"))

  # ... for Xb
  XbStart <- as.vector(crossprod(data$X,betaStart))

  # ... for the naive random effects
  uStart <- rep.int(0,q)

  # ... for Chol
  tZLInit <- crossprod(LStart,data$Z)
  CholStart <- Cholesky(tcrossprod(tcrossprod(tZLInit,Wsqrt)),Imult=1,perm=FALSE)
  
  # ... for the covariance parameters
  if (covStruct=="Identity")
   {
    newtonStart <- thetaOptExactpdIdent(theta=1,data=data,beta=betaStart,u=uStart,
                        lambda=lambda,weights=weights,fctstart=.Machine$integer.max,penalized=penalized,Wsqrt=Wsqrt,
				Xb=XbStart,family=family,Chol=CholStart,control=control)
    thetaStart <- newtonStart$theta
   } else if (covStruct=="Diagonal")
  {
    newtonStart <- thetaOptExactpdDiag(theta=rep(1,stot),data=data,beta=betaStart,u=uStart,
                        lambda=lambda,weights=weights,penalized=penalized,Wsqrt=Wsqrt,
			Xb=XbStart,L=LStart,stot=stot,control=control,ind=ind,N=N,family=family,Chol=CholStart)
    thetaStart <- newtonStart$theta
    LStart <- newtonStart$L
  }
  
  # ... for the advised random effects
  pirlsStart <- pirls(u=uStart,data=data,beta=betaStart,Wsqrt=Wsqrt,Xb=XbStart,theta=thetaStart,tZL=newtonStart$tZL,family=family,
                      inverse=TRUE,unit=unit,Chol=CholStart,tol=control$tol4)
  uStart <- pirlsStart$utilde

  # ... for tZL
  tZLStart <- pirlsStart$tZL
  
 } else
 {
  if ((p!=length(coefInit[[1]]))|(q!=length(coefInit[[3]])))
    stop("The length of the starting values do not coincide with the model you want to fit.")
   
  betaStart <- coefInit[[1]]
  XbStart <- as.vector(crossprod(data$X,betaStart))

  uStart <- coefInit[[3]]

  if (covStruct=="Identity")
   {
     thetaStart <- mean(coefInit[[2]])
     tZLInit <- thetaStart*data$Z
   } else if (covStruct=="Diagonal")
   {
     thetaStart <- coefInit[[2]]
     for (s in 1:stot)
       {
         index <- which(ind==s)
         LStart@x[index] <- rep(thetaStart[s],N)
         tZLInit <- crossprod(LStart,data$Z)
       }
   }
  
  CholStart <- Cholesky(tcrossprod(tcrossprod(tZLInit,Wsqrt)),Imult=1,perm=FALSE)
  pirlsStart <- pirls(u=uStart,data=data,beta=betaStart,Wsqrt=Wsqrt,Xb=XbStart,theta=thetaStart,tZL=tZLInit,family=family,
                      inverse=TRUE,unit=unit,Chol=CholStart,tol=control$tol4)
  uStart <- pirlsStart$utilde
  tZLStart <- pirlsStart$tZL # entspricht tZLInit
 }

# diagIndex
tZLW1kZL <- pirlsStart$CInv%*%tcrossprod(tcrossprod(tZLStart,Wsqrt),tZLStart)
diagIndex <- match(diag(tZLW1kZL),tZLW1kZL@x)

# diagIndex2
if (exactDeriv)
 {
   tZLDiagZL <- tcrossprod(tcrossprod(tZLStart,Wsqrt),tZLStart)+.symDiagonal(n=q,x=seq(1,q)/q)
   diag2Index <-  match(diag(tZLDiagZL),tZLDiagZL@x)
 }


# ... for the function value
fctStart <- ObjFctPirls(pirls=pirlsStart,beta=betaStart,weights=weights,penalized=penalized,lambda=lambda)

if (control$verbose>=2) print(fctStart)
if (control$fctSave) fctEval <- fctStart else fctEval <- nfctEval <- NULL

# --- GCD Iterations ---
# ----------------------

# some necessary allocations
convergence <- 0
doAll <- FALSE
pirlsEval <- FALSE
penalizedIter <- logical(p) ; penalizedIter[penalized] <- TRUE
lIter <- gradIter <- hessIter <- numeric(p)

nIter <- 0
nIterIn <- 0
nIterPirls <- 0

uIter <- uStart
betaIter <- betaStart
thetaIter <- thetaStart
fctIter <- fctStart
XbIter <- XbStart
tZLIter <- tZLStart
pirlsIter <- pirlsStart
if (covStruct=="Diagonal") LIter <- LStart
gradientIter <- numeric(p)
unit <- .symDiagonal(n=q)

convPar <- sum(crossprod(betaIter))
convFct <- convLem <- sum(fctStart)
convCov <- sum(crossprod(thetaIter))

while ( (control$maxIter > nIter) &  (convFct > control$tol1 | convPar > control$tol2 | convCov > control$tol3 | convLem > 10^(-3) | !doAll ) )
{

 nIter <- nIter + 1
 if (control$verbose>=1) print(nIter)
 betaIterOld <- betaIter
 fctIterOld <- fctIter
 thetaIterOld <- thetaIter

 if (nIterIn==0 | nIterIn>=control$number)
  {
   doAll <- TRUE
   activeSet <- 1:p
   nIterIn <- 1    
  } else
  {
   doAll <- FALSE
   activeSet <- which(betaIter!=0)
   nIterIn <- nIterIn+1
  }
 
 # 1) --- Calculating the blocks wrt beta ---

 for (k in activeSet) 
  { # START: cycle through the active set

   # 1) a) PIRLS + Laplace step
   if ((pirlsEval))#&(k==1))
    {
       nIterPirls <- nIterPirls+1
       pirlsIter <- pirls(u=uIter,data=data,beta=betaIter,Wsqrt=Wsqrt,Xb=XbIter,tZL=tZLIter,family=family,inverse=TRUE,
                          unit=unit,Chol=CholStart,tol=control$tol4)
       uIter <- pirlsIter$utilde
       fctIter <- ObjFctPirls(pirls=pirlsIter,beta=betaIter,weights=weights,penalized=penalized,lambda=lambda)
     if (control$verbose>=2) cat(fctIter,"\n")
     if (control$fctSave) fctEval <- c(fctEval,fctIter)
    }

   # 1) b) quadratic approximation + Armijo step
   if (exactDeriv)
     {
       derivIter <- derivBetaExact(Xk=data$X[k,],Xk2=tX2[k,],y=y,mu=pirlsIter$mu,tZL=tZLIter,lower=control$lower,upper=control$upper,
                          W1k=Wsqrt,W2k=Wsqrt,Wu=Wsqrt,unit=unit,Chol=pirlsIter$Chol,family=family,CInv=pirlsIter$CInv,diagIndex=diagIndex,u=uIter,
                                   diag2Index=diag2Index)       
     } else
     {
       derivIter <- derivBeta(Xk=data$X[k,],Xk2=tX2[k,],y=y,mu=pirlsIter$mu,tZL=tZLIter,lower=control$lower,upper=control$upper,
                          W1k=Wsqrt,unit=unit,Chol=pirlsIter$Chol,family=family,CInv=pirlsIter$CInv,diagIndex=diagIndex)
     }


   betaIter_k <- betaIter[k]
   lambda_k <- ifelse(penalizedIter[k],lambda/weights[k],0)
   gradientIter[k] <- derivIter$gradient+(penalizedIter[k])*lambda_k*sign(betaIter_k)
   
   directionIter <- desDir(deriv=derivIter,lambda_k=lambda_k,beta_k=betaIter_k,penalized_k=penalizedIter[k])

   gradIter[k] <- directionIter ; hessIter[k] <- derivIter$hessian
   
   if (directionIter!=0)
    {
     if (exactStep|k%in%unpenalized)
       {
         
         armijoIter <- armijoExact(data=data,beta=betaIter,beta_k=betaIter_k,lambda_k=lambda_k,theta=thetaIter,u=uIter,convergence=convergence,
                           k=k,direction=directionIter,deriv=derivIter,lambda=lambda,weights=weights,fct0=fctIter,penalized=penalized,
                           pen=penalizedIter[k],Wsqrt=Wsqrt,tZL=tZLIter,Xb0=XbIter,control=control,minArmijo=lIter[k],family=family,Chol=CholStart)
       } else
       {
         armijoIter <- armijo(data=data,beta=betaIter,pirls=pirlsIter,beta_k=betaIter_k,lambda_k=lambda_k,theta=thetaIter,convergence=convergence,
                           k=k,direction=directionIter,deriv=derivIter,lambda=lambda,weights=weights,fct0=fctIter,penalized=penalized,
                              pen=penalizedIter[k],Wsqrt=Wsqrt,tZL=tZLIter,Xb0=XbIter,control=control,minArmijo=lIter[k],family=family,Chol=pirlsIter$Chol)
       }
   
     lIter[k] <- ifelse(control$min.armijo,armijoIter$l,0)
     betaIter <- armijoIter$beta
     fctIter <- armijoIter$fct
     convergence <- armijoIter$convergence
     XbIter <- armijoIter$Xb
    
     if (control$verbose>=2) cat(fctIter,"\n")
     if (control$fctSave) fctEval <- c(fctEval,fctIter)
    
     pirlsEval <- TRUE
    } else pirlsEval <- FALSE

  } # END: cycle through the active set

 # 2) --- Calculating the blocks wrt theta ---
   
 # 2) a) PIRLS + Laplace step (kann evt. weggelassen werden)
 if (pirlsEval)
  {
   nIterPirls <- nIterPirls+1
   pirlsIter <- pirls(u=uIter,data=data,beta=betaIter,Wsqrt=Wsqrt,Xb=XbIter,tZL=tZLIter,family=family,Chol=CholStart,tol=control$tol4)
   uIter <- pirlsIter$utilde
   fctIter <- ObjFctPirls(pirls=pirlsIter,beta=betaIter,weights=weights,penalized=penalized,lambda=lambda)
   if (control$verbose>=2) {cat("--------","\n") ; print(fctIter)}
   if (control$fctSave)  {fctEval <- c(fctEval,fctIter)}
  } else
  {
   if (control$verbose>=2) cat("--------","\n")
  }                                     
 
 # 2) b) Newton step
 
 if (covStruct=="Identity")
  {
    newtonIter <- thetaOptExactpdIdent(theta=thetaIter,data=data,beta=betaIter,u=uIter,
                  lambda=lambda,weights=weights,fctstart=fctIter,penalized=penalized,Wsqrt=Wsqrt,Xb=XbIter,
                               family=family,Chol=CholStart,control=control)
  } else if (covStruct=="Diagonal")
  {
    newtonIter <- thetaOptExactpdDiag(theta=thetaIter,data=data,beta=betaIter,u=uIter,
                     lambda=lambda,weights=weights,penalized=penalized,Wsqrt=Wsqrt,Xb=XbIter,
                     L=LIter,stot=stot,control=control,ind=ind,N=N,family=family,Chol=CholStart)
    thetaIter <- newtonIter$theta
    LIter <- newtonIter$L
  }
 pirlsIter <- pirls(u=uIter,data=data,beta=betaIter,Wsqrt=Wsqrt,Xb=XbIter,tZL=newtonIter$tZL,family=family,inverse=TRUE,
                    unit=unit,Chol=CholStart,tol=control$tol4)
 uIter <- pirlsIter$utilde
 pirlsEval <- FALSE
     
 thetaIter <- newtonIter$theta
 fctIter <- newtonIter$fct
 tZLIter <- newtonIter$tZL
 if (control$verbose>=2) {print(fctIter) ; cat("--------","\n")}
 if (control$fctSave) {fctEval <- c(fctEval,fctIter)}

 
 # 3) --- Convergence analysis ---
 convPar <- sum(sqrt(crossprod(betaIter-betaIterOld))/(1+sqrt(crossprod(betaIter))))
 convFct <- abs(sum((fctIterOld-fctIter)/(1+abs(fctIter))))
 convCov <- sum(sqrt(crossprod(thetaIter-thetaIterOld))/(1+sqrt(crossprod(thetaIter))))
 convLem <- 0

 if (control$verbose>=3) cat("convPar:",convPar,"convFct:",convFct,"convCov:",convCov,"convLem:",convLem,"\n")

 if ((convFct <= control$tol1) & (convPar <= control$tol2) & (convCov <= control$tol3) & (convLem<=10^(-3))) nIterIn <- 0
 
} # END while

if (control$maxIter==nIter) {cat("maxIter reached","\n") ; convergence <- convergence+1}

relComp <- gradientIter[penalized][c(betaIter!=0)[penalized]]
maxGrad <- ifelse((exactStep|(length(relComp)<1)),Inf,max(abs(gradientIter[penalized][c(betaIter!=0)[penalized]])))

if ((control$verbose>=2)&(maxGrad>control$gradTol)) cat("Gradient tolerance exceeded.","\n")

#  --- final calculations ---
# ---------------------------

# random effects, sorted per effect
ranef <- as.vector(uIter)

# linear predictor eta
eta <- as.vector(crossprod(data$X,betaIter) + crossprod(tZLIter,uIter)@x)

# mu and fitted value
fitValues <- fittedValues(eta,family)
mu <- fitValues$mu
fitted <- fitValues$fitted

if (family=="binomial") { tabl <- ftable(y~fitted) ; gof <- (ntot-sum(diag(tabl)))/ntot }
if (family=="poisson") { gof <- sum((y-mu)^2/mu) }

# conditional variances
if (family=="binomial") Var <- mu*(1-mu) else Var <- mu

# residuals
respResid <- resid <- y - mu
pearResid <- resid/sqrt(Var)
workResid <- respResid/Var

if (family=="binomial")
 {
   di <- -2*(y*eta+log(1-mu))
 } else if (family=="poisson")
 {
   di <- numeric(ntot) ;  index <- y!=0 ; di[index] <- 2*(y[index]*log(y[index]/mu[index])-(y[index]-mu[index]))
 }

devResid <- sign(y-mu)*sqrt(di)

# varia
cpu <- proc.time() - ptm
activeSet <- which(betaIter!=0)
data$X <- t(data$X) ; data$Z <- t(data$Z)

# --- Summary Information ---
# ---------------------------

npar <- sum(betaIter!=0) + length(thetaIter)
logLik <- -1/2*(fctIter - lambda*sum(abs(betaIter[penalized,drop=FALSE])/weights[penalized]))
deviance <- devianceFct(y=y,mu=mu,u=uIter,Wsqrt=Wsqrt,tZL=tZLIter,fct=fctIter,family=family,lambda=lambda,
                        beta=betaIter,weights=weights,penalized=penalized)

aic <- deviance + 2*npar
bic <- deviance + log(ntot)*npar

if (covStruct=="Identity") {thetaIter <- rep(thetaIter,stot) ; thetaStart <- rep(thetaStart,stot)}

list.out <- list(fixef=betaIter,coefficients=as.vector(betaIter),theta=thetaIter,ranef=ranef,u=as.vector(uIter),objective=fctIter,
                 logLik=logLik,deviance=deviance,aic=aic,bic=bic,activeSet=activeSet,eta=eta,mu=mu,fitted=fitted,lambda=lambda,
                 weights=weights,data=data,family=family,ntot=ntot,p=p,N=N,unpenalized=unpenalized,ranInd=ranInd,exactStep=exactStep,exactDeriv=exactDeriv,
                 coefInit=list(beta=betaStart,theta=thetaStart,u=uStart),coefOut=list(beta=betaIter,theta=thetaIter,u=uIter),    
                 convergence=convergence,nIter=nIter,nIterPirls=nIterPirls,nfctEval=length(fctEval),fctEval=fctEval,
                 gradient=gradientIter,maxGrad=maxGrad,maxArmijo=lIter,control=control,         
                 resid=resid,pearResid=pearResid,respResid=respResid,workResid=workResid,devResid=devResid,Var=Var,di=di,
                 gof=gof,cpu=cpu)

list.out
structure(list.out,class="glmmlasso")
  
}

