% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.R
\name{summary.test_mediation}
\alias{summary.test_mediation}
\alias{summary.boot_test_mediation}
\alias{summary.sobel_test_mediation}
\title{Summary of results from (robust) mediation analysis}
\usage{
\method{summary}{boot_test_mediation}(object, type = c("boot", "data"), plot = TRUE, ...)

\method{summary}{sobel_test_mediation}(object, ...)
}
\arguments{
\item{object}{an object inheriting from class
\code{"\link{test_mediation}"} containing results from (robust) mediation
analysis.}

\item{type}{a character string specifying how to summarize the effects
other than the indirect effect(s).  Possible values are \code{"boot"} (the
default) to compute significance tests using the normal approximation of the
bootstrap distribution (i.e., to assume a normal distribution of the
corresponding effect with the standard deviation computed from the bootstrap
replicates), or \code{"data"} to compute significance tests via
statistical theory based on the original data (e.g., t-tests if the
coefficients are estimated via regression).  Note that this is only relevant
for mediation analysis via a bootstrap test, where significance of the
indirect effect is always assessed via a percentile-based confidence
interval due to the asymmetry of its distribution.}

\item{plot}{a logical indicating whether to include a diagnostic plot of
robust regression weights (see \code{\link{weight_plot}()}).  This is only
used for mediation analysis objects fitted with the robust MM-estimator (see
\code{\link{test_mediation}()}).  Note that the diagnostic plot is only
shown when the returned object is printed in order to maintain a clear
separation between computing results and printing/plotting them.}

\item{\dots}{additional arguments are currently ignored.}
}
\value{
An object of class \code{"summary_test_mediation"} with the
following components:
\item{object}{the \code{object} passed to the \code{summary} method, which
contains the results from testing the indirect effect(s).}
\item{summary}{an object containing all necessary information to summarize
the effects other than the indirect effect(s).}
\item{plot}{if applicable, an object inheriting from class
\code{"\link[ggplot2]{ggplot}"} containing the diagnostic plot.}
}
\description{
Summarize results from (robust) mediation analysis for proper interpretation.
}
\examples{
data("BSG2014")

## seed to be used for the random number generator
seed <- 20211117

## simple mediation
# set seed of the random number generator
set.seed(seed)
# The results in Alfons et al. (2022a) were obtained with an
# older version of the random number generator.  To reproduce
# those results, uncomment the two lines below.
# RNGversion("3.5.3")
# set.seed(20150601)
# perform mediation analysis
boot_simple <- test_mediation(TeamCommitment ~
                                m(TaskConflict) +
                                  ValueDiversity,
                              data = BSG2014)
summary(boot_simple)
# the diagnostic plot is not shown when the summary is
# computed, only when the resulting object is printed
summary_simple <- summary(boot_simple)  # does not show plot
summary_simple                          # shows output and plot

\donttest{
## serial multiple mediators
# set seed of the random number generator
set.seed(seed)
# perform mediation analysis
boot_serial <- test_mediation(TeamScore ~
                                serial_m(TaskConflict,
                                         TeamCommitment) +
                                ValueDiversity,
                              data = BSG2014)
summary(boot_serial)

## parallel multiple mediators and control variables
# set seed of the random number generator
set.seed(seed)
# perform mediation analysis
boot_parallel <- test_mediation(TeamPerformance ~
                                  parallel_m(ProceduralJustice,
                                             InteractionalJustice) +
                                  SharedLeadership +
                                  covariates(AgeDiversity,
                                             GenderDiversity),
                                data = BSG2014)
summary(boot_parallel)
}

}
\references{
Alfons, A., Ates, N.Y. and Groenen, P.J.F. (2022a) A Robust Bootstrap Test
for Mediation Analysis.  \emph{Organizational Research Methods},
\bold{25}(3), 591--617.  doi:10.1177/1094428121999096.

Alfons, A., Ates, N.Y. and Groenen, P.J.F. (2022b) Robust Mediation Analysis:
The \R Package \pkg{robmed}.  \emph{Journal of Statistical Software},
\bold{103}(13), 1--45.  doi:10.18637/jss.v103.i13.
}
\seealso{
\code{\link{test_mediation}()}, \code{\link{weight_plot}()}
}
\author{
Andreas Alfons
}
\keyword{utilities}
