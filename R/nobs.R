# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

## extract number of observations from Huber M-estimator of location and scatter
#' @export
#' @import stats
nobs.covHuber <- function(object, ...) {
  weights <- weights(object, type="relative")
  length(weights)
}

## extract number of observations from MLE of mean vector and covariance matrix
#' @export
#' @import stats
nobs.covML <- function(object, ...) object$n
