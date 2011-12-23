\name{xerop}
\alias{xerop}
\docType{data}
\title{Dataset of Xerophthalmia}
\description{This is a subset of the Xerophthalmia data described in
  Diggle et al. (2002) and Zeger and Karim (1991).
}
\usage{data(xerop)}
\format{
  \describe{
   A list with the following four components.
   \item{y}{Binary response variable. If the child suffers from xerophthalmia.}
   \item{grp}{Grouping variable comprising the child id.}
   \item{X}{Fixed-effect design matrix. The first column is the intercept,
   then age, xero, cos, sin, sex, height and stunt. The
   covariates are all standardized with mean 0 and variance 1.}
   \item{Z}{Random-effects design matrix for a random-intercept model.}
  }
}
\details{
A detailed description of the covariates can be found in Diggle
et. al. and Zeger and Karim (1991).
}
\source{
\url{http://faculty.washington.edu/heagerty/Books/AnalysisLongitudinal/xerop.data}
}
\references{
Peter J. Diggle, Patrick Heagerty, Kung-Yee Liang and Scott L. Zeger
(2002), \emph{Analysis of Longitudinal Data}, Oxford University Press.

Scott L. Zeger and M. Rezaul Karim (1991), \emph{Generalized Linear
  Models With Random Effects; A Gibbs Sampling Approach}, Journal of the American Statistical Association, Vol. 86, No. 413 (Mar., 1991), pp. 79-86.
}

\examples{
data(xerop)
}
\keyword{datasets}
