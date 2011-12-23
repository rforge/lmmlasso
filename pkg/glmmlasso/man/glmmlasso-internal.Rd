\name{glmmlasso-internal}
\alias{ObjFctPirls}
\alias{armijo}
\alias{armijoExact}
\alias{armijoFct}
\alias{armijoFctExact}
\alias{derivBeta}
\alias{derivBetaExact}
\alias{desDir}
\alias{devianceFct}
\alias{glmPart}
\alias{fittedValues}
\alias{pirls}
\alias{thetaFctExactpdDiag}
\alias{thetaFctExactpdIdent}
\alias{thetaOptExactpdDiag}
\alias{thetaOptExactpdIdent}
\title{Internal glmmlasso objects}
\description{Internal glmmlasso objects.}
\usage{
ObjFctPirls(pirls, beta, weights, penalized, lambda)
armijo(data, beta, theta, pirls, k, beta_k, lambda_k, direction, deriv,
lambda, weights, fct0, penalized, Wsqrt, convergence = convergence, tZL,
pen, control, minArmijo, Xb0, family, Chol)
armijoExact(data, tX, beta, theta, u, k, beta_k, direction, lambda_k, deriv, lambda, weights, fct0, penalized, Wsqrt, convergence = convergence, tZL, pen, control, minArmijo, Xb0, family, Chol)
armijoFct(y, tZL, u.cp, Wsqrt, eta, lsabeta, family, Chol)
armijoFctExact(pirls, lsabeta)
derivBeta(Xk, Xk2, y, mu, tZL, lower, upper, W1k, unit, Chol, family,CInv = NULL, diagIndex)
derivBetaExact(Xk, Xk2, y, mu, tZL, lower, upper, W1k, W2k, Wu, unit, Chol, family, CInv = NULL, diagIndex, u, diag2Index)
desDir(deriv, lambda_k, beta_k, penalized_k)
devianceFct(y, mu, u, Wsqrt, tZL, fct, family, lambda, beta, weights,penalized)
glmPart(y, eta, family)
fittedValues(eta, family)
pirls(u, data, beta, tol, Wsqrt, theta = NULL, Xb = NULL, tZL = NULL,
family, inverse = FALSE, unit = NULL, CInv = NULL, Chol)
thetaFctExactpdDiag(thetal, data, Xb, beta, tZ, u, add, Wsqrt, L, index,
N, family, tol4, Chol = NULL)
thetaFctExactpdIdent(theta, data, Xb, beta, u, add, Wsqrt, family, tol4,
Chol = NULL)
thetaOptExactpdDiag(theta, data, beta, u, lambda, weights, penalized, Wsqrt, Xb, family, L, stot, control, ind, N, Chol = NULL)
thetaOptExactpdIdent(theta, data, beta, u, lambda, weights, fctstart,
penalized, Wsqrt, Xb, family, Chol = NULL, control)
}
\details{These are not to be called by the user.}
\keyword{internal}