glmmlassoControl <-
function(family,verbose=0,maxIter=200,number=0,CovOpt=c("nlminb"),fctSave=TRUE,
         a_init=1,delta=0.5,rho=0.1,gamm=0,lower=10^(-6),
         upper=ifelse(family=="binomial",10^5,10^3),
         seed=418,maxArmijo=20,min.armijo=TRUE,thres=10^(-4),
         tol1=10^(-6),tol2=10^(-6),tol3=10^(-3),tol4=10^(-8),gradTol=10^(-3))
  {
  
    list(verbose=verbose,maxIter=maxIter,number=number,a_init=a_init,
         delta=delta,rho=rho,gamm=gamm,lower=lower,upper=upper,seed=seed,
         maxArmijo=maxArmijo,min.armijo=min.armijo,thres=thres,tol1=tol1,
         tol2=tol2,tol3=tol3,tol4=tol4,gradTol=gradTol,
         CovOpt=CovOpt,fctSave=fctSave)
  }

