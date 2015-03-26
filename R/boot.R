# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

# simple bootstrap with prespecified bootstrap samples 
# (implemented for use in simulation studies)
localBoot <- function(data, statistic, R, indices = NULL, ...) {
  # initializations
  matchedCall <- match.call()
  matchedCall[[1]] <- as.name("boot")
  n <- NROW(data)
  if (is.null(n) || n == 0) stop("no data in bootstrap")
  # prepare bootstrap samples
  if(is.null(indices)) {
    R <- round(rep(R, length.out=1))
    if(!isTRUE(R > 0)) stop("invalid number of replications in bootstrap")
    cat("generating bootstrap samples...\n")
    indices <- replicate(R, sample(n, replace=TRUE))
  } else {
    cat("using prespecified bootstrap samples...\n")
    indices <- as.matrix(indices)
    R <- ncol(indices)
  }
  # compute statistic on original data
  t0 <- statistic(data, seq_len(n), ...)
  # perform the bootstrap
  t <- lapply(seq_len(R), function(r) statistic(data, indices[, r], ...))
  t <- do.call(rbind, t)
  # return results
  out <- list(t0=t0, t=t, R=R, data=data, seed=NULL, statistic=statistic, 
              sim="ordinary", call=matchedCall, stype="i", 
              strata=rep.int(1, n), weights=rep.int(1/n, n))
  out$indices <- t(indices)
  class(out) <- "boot"
  out
}

boot.array <- function(boot.out, indices=FALSE) {
  out <- boot.out$indices
  if(!indices) out <- boot:::freq.array(out)
  out
}
