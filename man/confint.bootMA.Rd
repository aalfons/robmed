\name{confint.bootMA}
\alias{confint.bootMA}
\alias{confint.sobelMA}
\title{Confidence intervals for (robust) mediation analysis}
\usage{
  \method{confint}{bootMA} (object, parm = NULL,
    level = NULL, ...)

  \method{confint}{sobelMA} (object, parm = NULL,
    level = 0.95, ...)
}
\arguments{
  \item{object}{an object of class \code{"bootMA"} or
  \code{"sobelMA"} containing results from (robust)
  mediation analysis, as returned by
  \code{\link{mediate}}.}

  \item{parm}{an integer, character or logical vector
  specifying the coefficients for which to extract or
  compute confidence intervals, or \code{NULL} to extract
  or compute confidence intervals for all coefficients.}

  \item{level}{for the \code{"bootMA"} method, this is
  ignored and the confidence level of the bootstrap
  confidence interval for the indirect effect is used.  For
  the \code{"sobelMA"} method, the confidence level of the
  confidence intervals to be computed.  The default is to
  compute 95\% confidence intervals.}

  \item{\dots}{additional arguments are currently ignored.}
}
\value{
  A numeric matrix containing the requested confidence
  intervals.
}
\description{
  Extract or compute confidence intervals for coefficients
  from (robust) mediation analysis.
}
\author{
  Andreas Alfons
}
\seealso{
  \code{\link{mediate}}, \code{\link[=coef.testMA]{coef}}
}
\keyword{utilities}

