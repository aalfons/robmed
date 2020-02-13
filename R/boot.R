# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------


#' Draw bootstrap samples
#'
#' Draw bootstrap samples to be used for (fast and robust) bootstrap tests
#' for mediation analysis.  Note that this function is intended for use in
#' simulation studies, it is not expected to be called by the user.
#'
#' @param n  an integer giving the number of observations in the original data
#' set.
#' @param R  an integer giving the number of bootstrap samples to be generated.
#'
#' @return An object of class \code{"boot_samples"} with the following
#' components:
#' \item{indices}{an integer matrix in which each column contains the indices
#' of the corresponding bootstrap sample.}
#' \item{seed}{the state of the random number generator before the bootstrap
#' samples were drawn}
#'
#' @author Andreas Alfons
#'
#' @seealso \code{\link{test_mediation}()}
#'
#' @examples
#' \donttest{
#' # control parameters
#' n <- 100
#' a <- b <- c <- 0.2
#'
#' # generate data
#' set.seed(20150326)
#' x <- rnorm(n)
#' m <- a * x + rnorm(n)
#' y <- b * m + c * x + rnorm(n)
#' simulated_data <- data.frame(x, y, m)
#'
#' # perform boostrap tests
#' indices <- boot_samples(n, R = 5000)
#' standard_boot <- test_mediation(simulated_data,
#'                                 x = "x", y = "y", m = "m",
#'                                 robust = FALSE,
#'                                 indices = indices)
#' summary(standard_boot)
#' robust_boot <- test_mediation(simulated_data,
#'                               x = "x", y = "y", m = "m",
#'                               robust = TRUE,
#'                               indices = indices)
#' summary(robust_boot)
#' }
#'
#' @keywords utilities
#'
#' @export

boot_samples <- function(n, R) {
  # retrieve current seed of the random number generator
  seed <- .Random.seed
  # draw bootstrap samples the same way as package 'boot'
  indices <- sample.int(n, n*R, replace = TRUE)
  dim(indices) <- c(R, n)
  # transpose to have bootstrap samples in columns and return indices and seed
  # as object of class "boot_samples"
  out <- list(indices = t(indices), seed = seed)
  class(out) <- "boot_samples"
  out
}


## Function boot() from package 'boot' doesn't allow to supply prespecified
## bootstrap samples, but that functionality is useful for simulation studies.
## This is a wrapper function that implements a simple nonparametric bootstrap
## if indices of bootstrap samples are supplied, otherwise function boot() is
## called.
##
## Argument 'indices' should be an object of class "boot_samples" as generated
## by function boot_samples().

local_boot <- function(data, statistic, R, indices = NULL, ...,
                       # the following arguments of function boot() from
                       # package 'boot' should be ignored
                       sim, stype, L, m, ran.gen, mle) {
  # check if indices of bootstrap samples are supplied
  if (is.null(indices)) {
    # call function boot() from package 'boot'
    boot(data = data, statistic = statistic, R = R, ...)
  } else {
    # initializations
    matched_call <- match.call()
    matched_call[[1]] <- as.name("boot")
    n <- NROW(data)
    if (is.null(n) || n == 0) stop("no data to bootstrap")
    # extract information from object for bootstrap samples
    seed <- indices$seed
    indices <- as.matrix(indices$indices)
    R <- ncol(indices)
    # compute statistic on original data
    t0 <- statistic(data, seq_len(n), ...)
    # perform the bootstrap
    t <- lapply(seq_len(R), function(r) statistic(data, indices[, r], ...))
    t <- do.call(rbind, t)
    # return results
    out <- list(t0 = t0, t = t, R = R, data = data, seed = seed,
                statistic = statistic, sim = "ordinary", call = matched_call,
                stype = "i", strata = rep.int(1, n), weights = rep.int(1/n, n))
    # out$indices <- t(indices)  # for testing
    class(out) <- "boot"
    out
  }
}


### simple version of function boot.array() from package 'boot'
#
# boot.array <- function(boot.out, indices = FALSE) {
#   out <- boot.out$indices
#   if(!indices) out <- freq.array(out)
#   out
# }
