# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

#' Maximum likelihood estimator of mean vector and covariance matrix
#'
#' Compute the maximum likelihood estimator of the mean vector and the
#' covariance matrix.
#'
#' @aliases print.cov_ML
#'
#' @param x  a numeric matrix or data frame.
#' @param \dots  additional arguments are currently ignored.
#'
#' @return An object of class \code{"cov_ML"} with the following components:
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
#' @seealso \code{\link{test_mediation}()}, \code{\link{fit_mediation}()}
#'
#' @examples
#' data("BSG2014")
#'
#' # define variables
#' x <- "ValueDiversity"
#' y <- "TeamCommitment"
#' m <- "TaskConflict"
#'
#' # compute Huber M-estimator
#' cov_ML(BSG2014[, c(x, y, m)])
#'
#' @keywords multivariate
#'
#' @export

cov_ML <- function(x, ...) {
  # initializations
  x <- as.matrix(x)
  n <- nrow(x)
  # compute MLEs of mean vector and covariance matrix
  mu <- colMeans(x)
  Sigma <- crossprod(sweep(x, 2, mu, check.margin=FALSE)) / n
  # return results
  result <- list(center=mu, cov=Sigma, n=n)
  class(result) <- "cov_ML"
  result
}
