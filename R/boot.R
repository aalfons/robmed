# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

# # simple bootstrap with prespecified bootstrap samples 
# # (implemented for use in simulation studies)
# localBoot <- function(data, statistic, indices, ...) {
#   # initializations
#   matchedCall <- match.call()
#   matchedCall[[1]] <- as.name("boot")
#   n <- NROW(data)
#   if (is.null(n) || n == 0) stop("no data in bootstrap")
#   # prepare bootstrap samples
#   indices <- as.matrix(indices)
#   R <- ncol(indices)
#   # compute statistic on original data
#   t0 <- statistic(data, seq_len(n), ...)
#   # perform the bootstrap
#   t <- lapply(seq_len(R), function(r, ...) statistic(data, indices[, r], ...))
#   t <- do.call(rbind, t)
#   # return results
#   out <- list(t0=t0, t=t, R=R, data=data, seed=NULL, statistic=statistic, 
#               sim="ordinary", call=matchedCall, stype="i", 
#               strata=rep.int(1, n), weights=rep.int(1/n, n))
#   class(out) <- "boot"
#   out
# }
