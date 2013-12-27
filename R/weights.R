# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

## extract robustness weights from Huber-type M-estimator
#' @S3method weights covHuber
weights.covHuber <- function(object, type = c("consistent", "relative"), ...) {
  # initializations
  type <- match.arg(type)
  # extract weights
  weights <- object$weights
  if(type == "consistent") weights <- weights / sqrt(object$tau)
  weights
}
