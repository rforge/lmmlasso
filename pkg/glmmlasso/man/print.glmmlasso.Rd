\name{print.glmmlasso}
\alias{print.glmmlasso}
\title{Print a short summary of a lmmlasso object.}
\description{Prints a short summary of a \code{glmmlasso} object
  comprising information about the variance components parameters and
  the number of nonzero fixed-effects coefficients.}
\usage{
\method{print}{glmmlasso}(x, ...)
}

\arguments{
  \item{x}{a \code{glmmlasso} object}
  \item{\dots}{not used}
}
\seealso{
\code{summary.glmmlasso}
}

\examples{
data(xerop)
fit <- glmmlasso(y=xerop$y,group=xerop$group,X=xerop$X,Z=xerop$Z,
        family="binomial",covStruct="Identity",lambda=30)
print(fit)
}

\keyword{methods}
