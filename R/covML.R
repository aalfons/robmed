# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

#' Maximum likelihood estimator of mean vector and covariance matrix
#'
#' Compute the maximum likelihood estimator of the mean vector and the
#' covariance matrix.
#'
#' @aliases print.covML
#'
#' @param x  a numeric matrix or data frame.
#' @param \dots  additional arguments are currently ignored.
#'
#' @return An object of class \code{"covML"} with the following components:
#' \item{center}{a numeric vector containing the mean vector estimate.}
#' \item{cov}{a numeric matrix containing the covariance matrix estimate.}
#' \item{n}{an integer giving the number of observations.}
#'
#' @author Andreas Alfons
#'
#' @references
#' Zu, J. and Yuan, K.-H. (2010) Local influence and robust procedures for
#' mediation analysis. \emph{Multivariate Behavioral Research}, \bold{45}(1),
#' 1--44.
#'
#' @seealso \code{\link{testMediation}}, \code{\link{fitMediation}}
#'
#' @keywords multivariate
#'
#' @export

covML <- function(x, ...) {
  # initializations
  x <- as.matrix(x)
  n <- nrow(x)
  # compute MLEs of mean vector and covariance matrix
  mu <- colMeans(x)
  Sigma <- crossprod(sweep(x, 2, mu, check.margin=FALSE)) / n
  # return results
  result <- list(center=mu, cov=Sigma, n=n)
  class(result) <- "covML"
  result
}
