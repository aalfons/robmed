\name{covHuber.control}
\alias{covHuber.control}
\title{Tuning parameters for (skipped) Huber M-estimation of location and scatter}
\usage{
  covHuber.control(prob, skipped = FALSE, tol = 1e-06,
    maxIterations = 200)
}
\arguments{
  \item{prob}{numeric; probability for the quantile of the
  \eqn{\chi^{2}}{chi-squared} distribution to be used as
  cutoff point in the (skipped) Huber weight function
  (defaults to 0.95 for the Huber weight function and 0.975
  for the skipped Huber weight function).}

  \item{skipped}{a logical indicating whether the Huber
  weight function or the skipped Huber weight function
  should be used.}

  \item{tol}{a small positive numeric value to be used to
  determine convergence of the iteratively reweighted
  algorithm.}

  \item{maxIterations}{an integer giving the maximum number
  of iterations in the iteratively reweighted algorithm.}
}
\value{
  A list with components corresponding to the arguments.
}
\description{
  Obtain a list with tuning paramters for
  \code{\link{covHuber}}.
}
\author{
  Andreas Alfons
}
\references{
  Huber, P. J. (1981) \emph{Robust statistics}. John Wiley
  & Sons.
}
\seealso{
  \code{\link{covHuber}}
}
\keyword{multivariate}

