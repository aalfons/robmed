# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

covMLE <- function(x, ...) {
  # initializations
  x <- as.matrix(x)
  n <- nrow(x)
  # compute MLEs of mean vector and covariance matrix
  mu <- colMeans(x)
  Sigma <- crossprod(sweep(x, 2, mu, check.margin=FALSE)) / n
  # return results
  result <- list(center=mu, cov=Sigma, n=n)
  class(result) <- "covMLE"
  result
}
