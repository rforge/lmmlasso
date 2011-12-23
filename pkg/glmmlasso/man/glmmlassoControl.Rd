\name{glmmlassoControl}
\alias{glmmlassoControl}
\title{Options for the glmmlasso algorithm}
\description{Definition of various kinds of options in the algorithm.}
\usage{
glmmlassoControl(family, verbose = 0, maxIter = 200, number = 0,
CovOpt=c("nlminb"), fctSave = TRUE, a_init = 1, delta = 0.5,
rho = 0.1,gamm = 0, lower = 10^(-6),
upper = ifelse(family == "binomial", 10^5,10^3), seed = 418,
maxArmijo = 20, min.armijo = TRUE, thres = 10^(-4),
tol1 = 10^(-6), tol2 = 10^(-6), tol3 = 10^(-3), tol4 = 10^(-8),
gradTol = 10^(-3))
}
\arguments{
  \item{family}{a GLM family. Currently implemented are "binomial" (default) and "poisson".}
  \item{verbose}{integer. 0 prints no output, 1 prints the outer
    iteration step, 2 prints the current function value, 3 prints the values of the convergence criteria}
  \item{maxIter}{maximum number of (outer) iterations}
  \item{number}{integer. Determines the active set algorithm. The zero
    fixed-effects coefficients are only updated each number
    iteration. Use \eqn{0 \le number \le 10}.}
  \item{CovOpt}{character string indicating which covariance parameter optimizer to use. Currently, only "nlminb" is implemented}
  \item{fctSave}{Should all evaluation of the objective function be stored? It may help to identify the convergence pattern of the algorithm.}
  \item{a_init}{\eqn{\alpha_{init}} in the Armijo step.}
  \item{delta}{\eqn{\delta} in the Armijo step.}
  \item{rho}{\eqn{\rho} in the Armijo step.}
  \item{gamm}{\eqn{\gamma} in the Armijo step.}
  \item{lower}{lower bound for the Hessian}
  \item{upper}{upper bound for the Hessian}
  \item{seed}{set.seed in order to choose the same starting value in the cross-validation for the fixed effects}
  \item{maxArmijo}{maximum number of steps to be chosen in the Armijo
    step. If the maximum is reached, the algorithm continues with optimizing the next coordinate.}
  \item{min.armijo}{logical. If TRUE, the smallest l in the Armijo step
    is increased, as suggested in Tseng and Yun (2009). Otherwise l always starts with 0.}
  \item{thres}{if a variance or covariance parameter has smaller absolute value than thres, the parameter is set to exactly zero,}
  \item{tol1}{convergence tolerance for the relative change in the function value}
  \item{tol2}{convergence tolerance for the relative change in the fixed-effects parameters}
  \item{tol3}{convergence tolerance for the relative change in the covariance parameters}
  \item{tol4}{convergence tolerance in the PIRLS algorithm}
  \item{gradTol}{the tolerance for the gradient accepted without giving a warning}
}
\details{For the Armijo step parameters, see Bertsekas (2003).}

\references{Dimitri P. Bertsekas (2003) \emph{Nonlinear Programming},
  Athena Scientific.}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}
