\name{mediate}
\alias{mediate}
\alias{summary.bootMA}
\alias{summary.sobelMA}
\title{(Robust) mediation analysis}
\usage{
  mediate(x, y, m, method = c("boot", "sobel"), R = 5000,
    alpha = 0.05, robust = TRUE, transform = FALSE,
    control, ...)
}
\arguments{
  \item{x}{the explanatory variable.}

  \item{y}{the response variable.}

  \item{m}{the mediator variable.}

  \item{method}{the type of test to be performed.  Possible
  values are \code{"boot"} for the bootstrap or
  \code{"sobel"} for Sobel's test.}

  \item{R}{the number of bootstrap replicates.}

  \item{alpha}{the significance level of the test.}

  \item{\dots}{additional arguments to be passed to
  \code{\link[boot]{boot}}.}
}
\value{
  An object containing the necessary information for
  mediation analysis.  A \code{summary} method to produce
  summaries is available.
}
\description{
  Perform (robust) mediation analysis via a (fast and
  robust) bootstrap test or Sobel's test.
}
\author{
  Andreas Alfons
}
\seealso{
  \code{\link[robustbase]{lmrob}}, \code{\link[stats]{lm}},
  \code{\link[boot]{boot}}
}
\keyword{multivariate}

