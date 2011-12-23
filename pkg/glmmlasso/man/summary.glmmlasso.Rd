\name{summary.glmmlasso}
\alias{summary.glmmlasso}
\title{Summarize a glmmlasso object}
\description{Providing an elaborate summary of a \code{glmmlasso} object.}
\usage{
\method{summary}{glmmlasso}(object, ...)
}

\arguments{
  \item{object}{a \code{glmmlasso} object}
  \item{\dots}{not used.}
}
\details{This functions shows a detailed summary of a \code{glmmlasso}
  object. In the fixed-effects part, (n) right from a fixed-effects
  coefficient means that this coefficient was not subject to penalization.}

\examples{
data(xerop)
fit <- glmmlasso(y=xerop$y,group=xerop$group,X=xerop$X,Z=xerop$Z,
        family="binomial",covStruct="Identity",lambda=30)
summary(fit)
}

\keyword{methods}

