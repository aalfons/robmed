# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------


#' Draw bootstrap samples
#'
#' Draw bootstrap samples to be used for (fast-and-robust) bootstrap tests
#' for mediation analysis.  Note that this function is intended for use in
#' simulation studies by experienced users.
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
#' a <- b <- c <- 0.4
#'
#' # generate data
#' set.seed(20200309)
#' x <- rnorm(n)
#' m <- a * x + rnorm(n)
#' y <- b * m + c * x + rnorm(n)
#' simulated_data <- data.frame(x, y, m)
#'
#' # perform boostrap tests
#' indices <- boot_samples(n, R = 5000)
#' robust_boot <- test_mediation(simulated_data,
#'                               x = "x", y = "y", m = "m",
#'                               robust = TRUE,
#'                               indices = indices)
#' summary(robust_boot)
#' ols_boot <- test_mediation(simulated_data,
#'                            x = "x", y = "y", m = "m",
#'                            robust = FALSE,
#'                            indices = indices)
#' summary(ols_boot)
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


# ## wrapper function for boot() that ignores unused arguments, but allows
# ## arguments for parallel computing to be passed down
# local_boot <- function(..., sim, stype, L, m, ran.gen, mle) boot(...)

# ## simple version of function boot.array() from package 'boot'
# boot.array <- function(boot.out, indices = FALSE) {
#   out <- boot.out$indices
#   if(!indices) out <- freq.array(out)
#   out
# }


## internal function to compute confidence intervals from bootstrap results
# For regression fits, only the bootstrap replicates of the regression
# coefficients are stored, and the bootstrap replicates of the effects in the
# mediation model are subsequently extracted.  Hence the former object of class
# "boot" needs to be used as a donor for computing confidence intervals.
# t0 ....... point estimates on the original sample
# t ........ a matrix of bootstrap replicates
# object ... an object of class "boot" to be used as donor
# index .... integer vector of columns for which to compute confidence interval
boot_ci <- function(t0, t, object, parm = NULL, ...) {
  # copy original point estimates and bootstrap replicates to donor object
  object$t0 <- t0
  object$t <- t
  # this shouldn't be necessary, but just in case: removing the function for
  # the bootstrap replicates ensures that the empirical influence function
  # for BCa confidence intervals is not computed with another resampling
  # approach (jackknife) that uses the wrong function
  object$statistic <- NULL
  # compute confidence intervals
  if (is.null(parm)) parm <- seq_along(t0)
  if (length(parm) == 1L) ci <- extract_ci(parm, object = object, ...)
  else {
    # extract all confidence intervals
    ci <- lapply(parm, extract_ci, object = object, ...)
    ci <- do.call(rbind, ci)
    # copy row names from point estimates on original data
    rownames(ci) <- names(t0)
  }
  # return confidence intervals
  ci
}

## internal function compute a confidence interval from bootstrap results
# parm ..... integer index of a single column of bootstrap replicates
# object ... object of class "boot" containing bootstrap replicates
extract_ci <- function(parm = 1L, object, alternative = "twosided",
                       level = 0.95, type = "bca", ...) {
  # initializations
  which <- if(type == "perc") "percent" else type
  # extract confidence interval
  if (level == 0) {
    ci <- rep.int(mean(object$t[, parm], na.rm = TRUE), 2L)
  } else if (level == 1) {
    ci <- c(-Inf, Inf)
  } else if (alternative == "twosided") {
    ci <- boot.ci(object, conf = level, type = type,
                  index = parm)[[which]][4:5]
  } else {
    alpha <- 1 - level
    ci <- boot.ci(object, conf = 1-2*alpha, type = type,
                  index = parm)[[which]][4:5]
    if (alternative == "less") ci[1] <- -Inf
    else ci[2] <- Inf
  }
  # add names for confidence bounds and return confidence interval
  names(ci) <- c("Lower", "Upper")
  ci
}


## internal function compute p-values based on bootstrap confidence intervals
# For regression fits, only the bootstrap replicates of the regression
# coefficients are stored, and the bootstrap replicates of the effects in the
# mediation model are subsequently extracted.  Hence the former object of class
# "boot" needs to be used as a donor for computing confidence intervals.
# t0 ....... point estimates on the original sample
# t ........ a matrix of bootstrap replicates
# object ... an object of class "boot" to be used as donor
# index .... integer vector of columns for which to compute p-values
boot_p_value <- function(t0, t, object, parm = NULL, ...) {
  # copy original point estimates and bootstrap replicates to donor object
  object$t0 <- t0
  object$t <- t
  # this shouldn't be necessary, but just in case: removing the function for
  # the bootstrap replicates ensures that the empirical influence function
  # for BCa confidence intervals is not computed with another resampling
  # approach (jackknife) that uses the wrong function
  object$statistic <- NULL
  # compute p-values
  if (is.null(parm)) parm <- seq_along(t0)
  if (length(parm) == 1L) {
    p_values <- extract_p_value(parm, object = object, ...)
  } else {
    # extract all p-values
    p_values <- sapply(parm, extract_p_value, object = object, ...,
                       USE.NAMES = FALSE)
    # copy names from point estimates on original data
    names(p_values) <- names(t0)
  }
  # return p-values
  p_values
}

## internal function compute a p-value based on a bootstrap confidence interval
# parm ..... integer index of a single column of bootstrap replicates
# object ... object of class "boot" containing bootstrap replicates
extract_p_value <- function(parm = 1L, object, digits = 4L,
                            alternative = "twosided",
                            type = "bca", ...) {
  # set lower bound of significance level to 0
  lower <- 0
  # loop over the number of digits and determine the corresponding digit after
  # the comma of the p-value
  for (digit in seq_len(digits)) {
    # set step size
    step <- 1 / 10^digit
    # reset the significance level to the lower bound as we continue from there
    # with a smaller stepsize
    alpha <- lower
    # there is no rejection at the lower bound, so increase significance level
    # until there is rejection
    reject <- FALSE
    while (!reject) {
      # update lower bound and significance level
      lower <- alpha
      alpha <- alpha + step
      # retest at current significance level and extract confidence interval
      ci <- extract_ci(parm, object = object, alternative = alternative,
                       level = 1 - alpha, type = type)
      # reject if 0 is not in the confidence interval
      reject <- prod(ci) > 0
    }
  }
  # return smallest significance level where 0 is not in the confidence interval
  alpha
}
