# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

## extract number of observations from Huber-type M-estimator of 
## location and scatter
#' @S3method nobs covHuber
#' @import stats
nobs.covHuber <- function(object, ...) {
  weights <- weights(object, type="relative")
  length(weights)
}

## extract number of observations from MLE of location and scatter
nobs.covMLE <- function(object, ...) object$n
