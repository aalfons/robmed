\name{fortify.bootMA}
\alias{fortify.bootMA}
\alias{fortify.sobelMA}
\title{Convert (robust) mediation analysis results into a data frame for plotting}
\usage{
  \method{fortify}{bootMA} (model, data,
    method = c("dot", "density"), parm = c("c", "ab"), ...)

  \method{fortify}{sobelMA} (model, data,
    method = c("dot", "density"), parm = c("c", "ab"),
    level = 0.95, ...)
}
\arguments{
  \item{model}{an object of class \code{"bootMA"} or
  \code{"sobelMA"} containing results from (robust)
  mediation analysis, as returned by
  \code{\link{mediate}}.}

  \item{data}{for the \code{"bootMA"} method, this is
  currently ignored.  For the \code{"sobelMA"} method, this
  is an optional numeric vector containing the
  \eqn{x}-values at which to evaluate the assumed normal
  density from Sobel's test (only used in case of a density
  plot).  The default is to take 100 equally spaced points
  between the estimated indirect effect \eqn{\pm}{+/-}
  three times the standard error according to Sobel's
  formula.}

  \item{method}{a character string specifying for which
  plot to construct the data frame.  Possible values are
  \code{"dot"} for a dot plot of selected coefficients, or
  \code{"density"} for a density plot of the indirect
  effect.}

  \item{parm}{a character string specifying the
  coefficients to be included in a dot plot.  The default
  is to include the direct and the indirect effect.}

  \item{level}{numeric; the confidence level of the
  confidence intervals from Sobel's test to be included in
  a dot plot.  The default is to include 95\% confidence
  intervals.}

  \item{\dots}{additional arguments are currently ignored.}
}
\value{
  A data frame containing the necessary data for the
  selected plot, as well as additional information stored
  in attributes.
}
\description{
  Supplement the estimated coefficients with other useful
  information for informative visualization of the (robust)
  mediation analysis results.  It is thereby possible to
  construct data frames for dot plots of selected
  coefficients, as well as density plots of the indirect
  effect.
}
\author{
  Andreas Alfons
}
\seealso{
  \code{\link{mediate}}, \code{\link[=mediatePlot]{plot}}
}
\keyword{utilities}

