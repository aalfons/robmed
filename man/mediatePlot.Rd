\name{mediatePlot}
\alias{autoplot.bootMA}
\alias{autoplot.sobelMA}
\alias{mediatePlot}
\alias{mediatePlot.bootMA}
\alias{mediatePlot.default}
\alias{mediatePlot.list}
\alias{mediatePlot.sobelMA}
\alias{plot.bootMA}
\alias{plot.sobelMA}
\title{Plot (robust) mediation analysis results}
\usage{
  mediatePlot(object, ...)

  \method{mediatePlot}{bootMA} (object,
    method = c("dot", "density"), parm = c("c", "ab"), ...)

  \method{mediatePlot}{sobelMA} (object, data,
    method = c("dot", "density"), parm = c("c", "ab"),
    level = 0.95, ...)

  \method{mediatePlot}{list} (object, data, ...)

  \method{mediatePlot}{default} (object,
    mapping = attr(object, "mapping"),
    facets = attr(object, "facets"), ...)

  \method{autoplot}{bootMA} (object, ...)

  \method{autoplot}{sobelMA} (object, ...)

  \method{plot}{bootMA} (x, ...)

  \method{plot}{sobelMA} (x, ...)
}
\arguments{
  \item{object,x}{an object of class \code{"bootMA"} or
  \code{"sobelMA"} containing results from (robust)
  mediation analysis, as returned by \code{\link{mediate}}.
  For \code{mediatePlot}, a list of such objects may be
  supplied as well.}

  \item{data}{an optional numeric vector containing the
  \eqn{x}-values at which to evaluate the assumed normal
  density from Sobel's test (only used in case of a density
  plot).  The default is to take 100 equally spaced points
  between the estimated indirect effect \eqn{\pm}{+/-}
  three times the standard error according to Sobel's
  formula.}

  \item{method}{a character string specifying which plot to
  produce. Possible values are \code{"dot"} for a dot plot
  of selected coefficients, or \code{"density"} for a
  density plot of the indirect effect.}

  \item{parm}{a character string specifying the
  coefficients to be included in a dot plot.  The default
  is to include the direct and the indirect effect.}

  \item{level}{numeric; the confidence level of the
  confidence intervals from Sobel's test to be included in
  a dot plot.  The default is to include 95\% confidence
  intervals.}

  \item{mapping}{an aesthetic mapping to override the
  default behavior (see \code{\link[ggplot2]{aes}} or
  \code{\link[ggplot2]{aes_string}}).}

  \item{facets}{a faceting formula to override the default
  behavior (only used in case of a dot plot).  If supplied,
  \code{\link[ggplot2]{facet_wrap}} or
  \code{\link[ggplot2]{facet_grid}} is called depending on
  whether the formula is one-sided or two-sided.}

  \item{\dots}{additional arguments to be passed to and
  from methods.}
}
\value{
  An object of class \code{"ggplot"} (see
  \code{\link[ggplot2]{ggplot}}).
}
\description{
  Produce dot plots of selected coefficients from
  regression models computed in (robust) mediation
  analysis, or density plots of the indirect effect.
}
\author{
  Andreas Alfons
}
\seealso{
  \code{\link{mediate}},
  \code{\link[=fortify.bootMA]{fortify}}
}
\keyword{hplot}

