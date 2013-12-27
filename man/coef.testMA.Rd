\name{coef.testMA}
\alias{coef.testMA}
\title{Coefficients in (robust) mediation analysis}
\usage{
  \method{coef}{testMA} (object, parm = NULL, ...)
}
\arguments{
  \item{object}{an object of class \code{"bootMA"} or
  \code{"sobelMA"} containing results from (robust)
  mediation analysis, as returned by
  \code{\link{mediate}}.}

  \item{parm}{an integer, character or logical vector
  specifying the coefficients to be extracted, or
  \code{NULL} to extract all coefficients.}

  \item{\dots}{additional arguments are currently ignored.}
}
\value{
  A numeric vetor containing the requested coefficients.
}
\description{
  Extract coefficients from regression models computed in
  (robust) mediation analysis.
}
\author{
  Andreas Alfons
}
\seealso{
  \code{\link{mediate}},
  \code{\link[=confint.bootMA]{confint}}
}
\keyword{utilities}

