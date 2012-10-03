\name{glmmlasso}
\alias{glmmlasso}
\alias{glmmlasso.default}
\title{Function to fit high-dimensional generalized linear mixed models.}
\description{Fits the solution for a high-dimensional generalized linear
  mixed model.}

\usage{
glmmlasso(y,...)

\method{glmmlasso}{default}(y, group, X, Z = NULL, family = c("binomial", "poisson"),
 covStruct = c("Identity", "Diagonal", "Sym"), lambda, weights = NULL,
 coefInit = NULL, exactStep = FALSE, exactDeriv = FALSE, random = NULL,
 unpenalized = 1:stot,ranInd = 1:stot,
 control = glmmlassoControl(family = family), ...)
}

\arguments{
  \item{y}{response variable of length n.}
  \item{group}{grouping variable of length n.}
  \item{X}{fixed-effects matrix as an n x p matrix. An intercept has to be included in the first column as (1,...,1).}
  \item{Z}{random-effects matrix as an n x q matrix, must be in sparse Matrix format (see package Matrix).}
  \item{family}{GLM family. Currently only "binomial" and "poisson" are implemented.}
  \item{covStruct}{Covariance Stucture to be used. "Identity" fits
  \eqn{\bm{\Sigma}_{\theta}=\theta^2\bm{1}_q}, i.e. one single
  covariance parameter for all random effects. "Diagonal" fits a
  diagonal matrix for \eqn{\bm{\Sigma}_{\theta}}, i.e.  for each random
  effect a different covariance parameter}
  \item{lambda}{non-negative regularization parameter}
  \item{weights}{weights for the fixed-effects covariates: NA means no
    penalization, 0 means drop this covariate ; if given, the argument
    unpenalized is ignored. By default each covariate has weight 1}
  \item{coefInit}{starting values. They must be of the same length as the model to be fitted, i.e. the number of variables p, the number of random-effects q and the number of covariance parameters stot must coincide. Otherwise a warning is issued.}
  \item{exactStep}{logical. Should the Armijo step include the update of
  the random effects u to ensure that the objective function strictly decreases?}
  \item{exactDeriv}{logical. Should the exact derivate be calculated or
  the derivative for fixed random effects u.}
  \item{random}{expression for the random-effects structure for
  non-correlated random effects of the form "(1|group)+(0+X2|group)+(0+X3|group)". It is used only for generating the corresponding Z matrix,
  and it dominates the Z matrix, i.e. if random and a Z matrix is given, the Z matrix corresponding to random is used.}
  \item{unpenalized}{indices as subset of \eqn{\{1,...,p\}} indicating
  which fixed-effects covariates are not subject to the
  \eqn{\ell_1}-penalty. Ignored if weights is given.}
  \item{ranInd}{indices of the random effects as subset of
    \eqn{\{1,...,p\}}. Only used for
    \code{summary.glmmlasso}. They need to be specified if the random-effects covariates do not corespond to the first columns in X.}
  \item{control}{control parameters for the algorithm, see \code{lmmlassoControl} for the details}
  \item{\dots}{not used.}
}
\details{
  All the details of the algorithm can be found in the article.

  Concerning logLik, deviance, aic and bic, we have to be very careful when
  comparing them with other generalized linear mixed model
  functions. If we study a low-dimensional data example and set
  \eqn{\lambda=0}, the log-likelihood is calculated as given in the
  paper. Deviance, aic and bic are computed such that they coincide with
  the results from \code{glmer} from the \code{lme4} package. The latter
  does not employ the standard definitions of deviance and
  log-likelihood function value. Nevertheless, the differences only
  depends on constants, and not the parameters.
}
\value{A \code{glmmlasso} object is returned, for which
 \code{coef},\code{resid}, \code{fitted}, \code{logLik}
  \code{print}, \code{summary}, \code{plot} methods exist.
  \item{fixef}{fixed-effects parameter beta}
  \item{coefficients}{fixed-effects parameter beta}
  \item{theta}{covariance parameter estimates}
  \item{ranef}{random effects in sparse vector format}
  \item{u}{random effects in dense vector format}
  \item{objective}{Value of the objective function corresponding to the estimates}
  \item{logLik}{value of the log-likelihood function. See details.}
  \item{deviance}{value of the deviance function. See details.}
  \item{aic}{AIC. See details.}
  \item{bic}{BIC. See details.}
  \item{activeSet}{Indices of the non-zero fixed-effects coefficients.}
  \item{eta}{The linear predictor at the current values.}
  \item{mu}{The conditional expectation of the response at the current values.}
  \item{fitted}{The fitted values at the current values.}
  \item{lambda}{non-negative regularization parameter}
  \item{weights}{weights (possible adapted to the argument weights)}
  \item{data}{data. List with y, group, X and Z}
  \item{family}{GLM family used.}
  \item{ntot}{total number of observations}
  \item{p}{number of fixed-effects parameters}
  \item{N}{number of groups/clusters}
  \item{unpenalized}{indices of the non-penalized fixed-effects covariates}
  \item{ranInd}{indices of the random effects as subset of \eqn{\{1,...,p\}}.}
  \item{exactStep}{logical. If the Armijo step includes the update of the random effects u or not.}
  \item{exactDeriv}{logical. If the exact derivate has been calculated
  or if the derivative for fixed random effects u has been calculated.}
  \item{coefInit}{starting values used in the algorithm for beta, theta and u}
  \item{coefOut}{list with estimates in the form required for the argument coefInit}
  \item{convergence}{integer giving convergence information. Each time
  maxArmijo was reached, convergence is increased by 2. If maxIter was
  reached, convergence is increased by 1.}
  \item{nIter}{number of outer iterations.}
  \item{nIterPirls}{number of pirls evaluation within the outer
  iteration. Pirls-Evaluation within the Armijo steps are not counted.}
  \item{nfctEval}{number of function evaluation. See value fctEval.}
  \item{fctEval}{vector of all function values calculated during the
  algorithm. It may be interesting if studying the convergence behaviour
  of the algorithm. Only if argument fctSave=TRUE}
  \item{gradient}{gradient of the objective function with respect to the
  fixed-effects coefficients}
  \item{maxGrad}{maximal value of the gradient which has to be close to
  zero for convergence.}
  \item{maxArmijo}{the maximal value of l used in each fixed-effects component.}
  \item{control}{see \code{lmmlassoControl}}
  \item{resid}{response residuals. See McCullagh and Nelder (1989).}
  \item{pearResid}{pearson residuals. See McCullagh and Nelder (1989).}
  \item{respResid}{response residuals. See McCullagh and Nelder (1989).}
  \item{workResid}{working residuals. See McCullagh and Nelder (1989).}
  \item{devResid}{deviance residuals. See McCullagh and Nelder (1989).}
  \item{Var}{conditional variance of the response at the current
  values. See McCullagh and Nelder (1989).}
  \item{di}{contribution of each observation to the deviance. See McCullagh and Nelder (1989).}
  \item{gof}{goodness-of-fit criterion. For family="binomial", it is the
  in-sample missclassification rate. For family="poisson", it is the
  Pearson X2 statistic. See McCullagh and Nelder (1989).}
  \item{cpu}{cpu time needed for the algorithm.}
}

\references{
P. McCullagh and J. A. Nelder (1989), Generalized Linear Models,
Chapman-Hall, London.
}

\author{Juerg Schelldorfer, \email{schell@stat.math.ethz.ch}}

\examples{
# (1) Use glmmlasso on the xerop data set
data(xerop)
fit1 <- glmmlasso(y=xerop$y,group=xerop$group,X=xerop$X,Z=xerop$Z,
         family="binomial",covStruct="Identity",lambda=30)
summary(fit1)
plot(fit1)

# (2) Use the glmmlasso on a small simulated data set
set.seed(142)

N <- 40           # number of groups
p <- 6            # number of covariates (including intercept)
q <- 2            # number of random effect covariates
ni <- rep(10,N)    # observations per group
ntot <- sum(ni)   # total number of observations

group <- factor(rep(1:N,times=ni)) # grouping variable

beta <- c(0,1,-1,1,-1,0,0) # fixed-effects coefficients
X <- cbind(1,matrix(rnorm(ntot*p),nrow=ntot)) # fixed-effects design matrix
Z <- t(glmer(rbinom(ntot,1,runif(ntot,0.3,1))~1+(1|group)+(0+X2|group),
      data=data.frame(X),family="binomial")@Zt) # random-effects design matrix

bi <- c(rnorm(N,0,2),rnorm(N,0,1))

eta <- X\%*\%beta + Z\%*\%bi
mu <- exp(eta)/(1+exp(eta))
y <- rbinom(ntot,1,mu@x)

# correct random effects structure
fit2 <- glmmlasso(y=y,group=group,X=X,Z=Z,family="binomial",lambda=50,
         covStruct="Diagonal")
summary(fit2)
plot(fit2)

# wrong random effects structure
fit3 <- glmmlasso(y=y,group=group,X=X,family="binomial",lambda=50,
         covStruct="Diagonal",random="(1|group) + (0+X2|group) + (0+X3|group)")
summary(fit3)
plot(fit3)
}
\keyword{models}
\keyword{regression}
