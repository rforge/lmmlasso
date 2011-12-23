\name{plot.glmmlasso}
\alias{plot.glmmlasso}
\title{Diagnostic Plots for a lmmlasso object}
\description{Plots four diagnostic plots for checking the model
  assumptions and supporting model selection for a \code{glmmlasso} object
}
\usage{
\method{plot}{glmmlasso}(x, ...)
}

\arguments{
  \item{x}{a \code{lmmlasso} object}
  \item{\dots}{not used.}
}
\details{
\code{plot.glmmlasso} shows four diagnostic plots which support
  checking the model assumption, model fit and may give hints for
  another model. 1) The first plot depicts the Tukey-Anscombe plot on
  the predictor scale. Points with the same
  color belong to the same group. 2) Plot depending on the GLM
  family. For family="poisson", the fitted values against the observed
  values is shown. For family="binomial", the Tukey-Anscombe plot on the response
  scale is depicted. 3) QQ-Plot of the
  predicted random effects. Be careful with the interpretation since the
  random effects have not been standardized. The color shows which
  points belong to the same random-effects covariate. 4) A histogram of
  the fixed-effects coefficients. For the interpretation of the
  Tukey-Anscombe plot in GLMs, see Faraway (2006).
}

\references{ Julian J. Faraway (2006) \emph{Extending the Linear Model with R},
  Chapman and Hall/CRC.
}

\examples{
data(xerop)
fit <- glmmlasso(y=xerop$y,group=xerop$group,X=xerop$X,Z=xerop$Z,
        family="binomial",covStruct="Identity",lambda=30)
plot(fit)
}

\keyword{hplot}
