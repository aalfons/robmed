# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

## extract number of observations from Huber M-estimator of location and scatter
#' @export
#' @import stats
nobs.cov_Huber <- function(object, ...) {
  weights <- weights(object, type="relative")
  length(weights)
}

## extract number of observations from MLE of mean vector and covariance matrix
#' @export
#' @import stats
nobs.cov_ML <- function(object, ...) object$n
