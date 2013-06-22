# ------------------------------------
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ------------------------------------

#' @export
covHuber <- function(x, control = covHuber.control(...), ...) {
  ## initialization
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  # check control arguments
  defaults <- formals(covHuber.control)
  prob <- control$prob
  if(!is.finite(prob)) prob <- defaults$prob
  else if(prob < 0) prob <- 0
  else if(prob > 1) prob <- 1
  tol <- control$tol
  if(!is.finite(tol) || tol < 0) tol <- 0
  maxIterations <- control$maxIterations
  if(!is.finite(maxIterations) || maxIterations < 0) maxIterations <- 0
  # define tuning parameters 
  # tau is chosen such that E[chi_p^2 u_2(chi_p^2)] = p (i.e., such that Sigma 
  # is asymptotically unbiased)
  kappa <- 1 - prob
  r <- sqrt(qchisq(prob, df=p))  # cutoff point
  tau <- (-(r^p * exp(-r^2/2))/(2^(p/2-1)*gamma(p/2))+p*(prob)+r^2*kappa)/p
  # compute starting values for location vector and scatter matrix
  mu <- colMeans(x)
  Sigma <- cov(x)
  ## perform iterative reweighting
  i <- 0
  continue <- TRUE
  while(continue && i < maxIterations) {
    oldMu <- mu
    oldSigma <- Sigma
    # compute weights based on mahalanobis distances
    d <- sqrt(mahalanobis(x, center=mu, cov=Sigma))  # mahalanobis distances
    u <- ifelse(d > r, r/d, 1)
    # update location vector and scatter matrix
    mu <- apply(x, 2, weighted.mean, w=u)
    Sigma <- crossprod(u * sweep(x, 2, mu, check.margin=FALSE)) / (n*tau)
    # check convergence
    i <- i + 1
    continue <- max(abs(mu-oldMu)) > tol || max(abs(Sigma-oldSigma)) > tol
  }
  ## return results
  if(i == 0) u <- rep.int(1, n)
  if(i == maxIterations && continue) {
    warning(sprintf("no convergence in %d iterations", maxIterations))
  }
  list(mu=mu, Sigma=Sigma, weights=u/sqrt(tau), nIterations=i)
}

## utility function to define control object for covHuber()
#' @export
covHuber.control <- function(prob = 0.95, tol = 1e-06, maxIterations = 200) {
  prob <- rep(as.numeric(prob), length.out=1)
  tol <- rep(as.numeric(tol), length.out=1)
  maxIterations <- rep(as.integer(maxIterations), length.out=1)
  list(prob=prob, tol=tol, maxIterations=maxIterations)
}
