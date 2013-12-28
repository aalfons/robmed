\name{covHuber}
\alias{covHuber}
\alias{print.covHuber}
\title{Huber M-estimator of location and scatter}
\usage{
  covHuber(x, control = covHuber.control(...), ...)
}
\arguments{
  \item{x}{a numeric matrix or data frame.}

  \item{control}{a list of tuning parameters as generated
  by \code{\link{covHuber.control}}.}

  \item{\dots}{additional arguments can be used to specify
  tuning parameters directly instead of via
  \code{control}.}
}
\value{
  An object of class \code{"covHuber"} with the following
  components:

  \item{center}{a numeric vector containing the location
  vector estimate.}

  \item{cov}{a numeric matrix containing the scatter matrix
  estimate.}

  \item{prob}{numeric; probability for the quantile of the
  \eqn{\chi^{2}}{chi-squared} distribution used as cutoff
  point in the Huber weight function.}

  \item{weights}{a numeric vector containing the relative
  robustness weights for the observations.}

  \item{tau}{numeric; correction for Fisher consistency
  under multivariate normal distributions.}

  \item{converged}{a logical indicating whether the
  iterative reweighting algorithm converged.}

  \item{iterations}{an integer giving the number of
  iterations required to obtain the solution.}
}
\description{
  Compute a Huber M-estimator of location and scatter,
  which is reasonably robust for a small number of
  variables.
}
\author{
  Andreas Alfons
}
\references{
  Huber, P. J. (1981) \emph{Robust statistics}. John Wiley
  & Sons.

  Zu, J. and Yuan, K.-H. (2010) Local influence and robust
  procedures for mediation analysis. \emph{Multivariate
  Behavioral Research}, \bold{45}(1), 1--44.
}
\seealso{
  \code{\link{covHuber.control}},
  \code{\link{testMediation}}, \code{\link{fitMediation}}
}
\keyword{multivariate}

