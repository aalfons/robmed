# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------


#' Confidence intervals from (robust) mediation analysis
#'
#' Extract or compute confidence intervals for effects in (robust)
#' mediation analysis.
#'
#' @name confint.test_mediation
#'
#' @param object  an object inheriting from class \code{"\link{test_mediation}"}
#' containing results from (robust) mediation analysis.
#' @param parm  an integer, character or logical vector specifying the paths
#' for which to extract or compute confidence intervals, or \code{NULL}
#' to extract or compute confidence intervals for all coefficients.  In
#' case of a character vector, possible values are \code{"a"}, \code{"b"},
#' \code{"d"} (only serial multiple mediator models), \code{"total"},
#' \code{"direct"}, and \code{"indirect"}.
#' @param level  for the \code{"boot_test_mediation"} method, this is ignored
#' and the confidence level of the bootstrap confidence interval for the
#' indirect effect is used.  For the other methods, the confidence level of
#' the confidence intervals to be computed.  The default is to compute 95\%
#' confidence intervals.
#' @param type  a character string specifying how to compute the confidence
#' interval of the effects other than the indirect effect(s).  Possible values
#' are \code{"boot"} (the default) to compute bootstrap confidence intervals
#' using the normal approximation (i.e., to assume a normal distribution of the
#' corresponding effect with the standard deviation computed from the bootstrap
#' replicates), or \code{"data"} to compute confidence intervals via
#' statistical theory based on the original data (e.g., based on a
#' t-distribution if the coefficients are estimated via regression).  Note that
#' this is only relevant for mediation analysis via a bootstrap test, where
#' the confidence interval of the indirect effect is always computed via a
#' percentile-based method due to the asymmetry of its distribution.
#' @param \dots  additional arguments are currently ignored.
#'
#' @return A numeric matrix containing the requested confidence intervals.
#'
#' @author Andreas Alfons
#'
#' @seealso
#' \code{\link{test_mediation}()}, \code{\link[=coef.test_mediation]{coef}()},
#' \code{\link{p_value}()}, \code{\link[boot]{boot.ci}()}
#'
#' @examples
#' data("BSG2014")
#'
#' # run fast-and-robust bootstrap test
#' robust_boot <- test_mediation(BSG2014,
#'                               x = "ValueDiversity",
#'                               y = "TeamCommitment",
#'                               m = "TaskConflict",
#'                               robust = TRUE)
#' confint(robust_boot, type = "boot")
#'
#' # run OLS bootstrap test
#' ols_boot <- test_mediation(BSG2014,
#'                            x = "ValueDiversity",
#'                            y = "TeamCommitment",
#'                            m = "TaskConflict",
#'                            robust = FALSE)
#' confint(ols_boot, type = "data")
#'
#' @keywords utilities

NULL


#' @rdname confint.test_mediation
#' @method confint boot_test_mediation
#' @export

## argument 'level' is ignored
confint.boot_test_mediation <- function(object, parm = NULL, level = NULL,
                                        type = c("boot", "data"), ...) {
  # confidence intervals of other effects
  type <- match.arg(type)
  if (type == "boot") {
    ci_list <- get_ci_list(object$fit, level = object$level, boot = object$reps)
  } else ci_list <- get_ci_list(object$fit, level = object$level)
  # add confidence intervals for indirect effects
  ci_list$indirect <- object$ci
  # if requested, take subset of confidence intervals
  if (!is.null(parm)) ci_list <- ci_list[check_parm(parm)]
  # stack confidence intervals on top of each other
  ci <- do.call(rbind, ci_list)
  # define row names
  keep <- names(ci_list)
  rn <- get_effect_names(effects = object[keep])
  # define column names for lower and upper bound
  cn <- get_ci_names(object$level, alternative = object$alternative)
  # add row and column names and return confidence intervals
  dimnames(ci) <- list(rn, cn)
  ci
}


#' @rdname confint.test_mediation
#' @method confint sobel_test_mediation
#' @export

confint.sobel_test_mediation <- function(object, parm = NULL, level = 0.95,
                                         ...) {
  # initializations
  level <- rep(as.numeric(level), length.out = 1L)
  if(is.na(level) || level < 0 || level > 1) level <- formals()$level
  # confidence intervals of other effects
  ci_list <- get_ci_list(object$fit, level = level)
  # add confidence interval for indirect effect
  ci_list$indirect <- confint_z(object$fit$indirect, object$se, level = level,
                                alternative = object$alternative)
  # if requested, take subset of confidence intervals
  if (!is.null(parm)) ci_list <- ci_list[check_parm(parm)]
  # stack confidence intervals on top of each other
  ci <- do.call(rbind, ci_list)
  # define row names
  keep <- names(ci_list)
  rn <- get_effect_names(effects = object$fit[keep])
  # define column names for lower and upper bound
  cn <- get_ci_names(level, alternative = object$alternative)
  # add row and column names and return confidence intervals
  dimnames(ci) <- list(rn, cn)
  ci
}


# there is no confint() method for median regression results
#' @export
confint.rq <- function(object, parm = NULL, level = 0.95, ...) {
  # compute the usual summary and extract coefficient matrix
  summary <- summary(object, se = "iid")
  coef_mat <- coef(summary)
  coef_names <- rownames(coef_mat)
  # extract point estimates and standard errors
  coef <- coef_mat[, 1L]
  se <- coef_mat[, 2L]
  # check parameters to extract
  if (missing(parm) || is.null(parm)) parm <- coef_names
  else if (is.numeric(parm)) parm <- coef_names[parm]
  # significance level and quantile of the t distribution
  a <- (1 - level) / 2
  a <- c(a, 1 - a)
  q <- qt(a, df = summary$rdf)
  # column names for output should contain percentages
  cn <- get_ci_names(level)
  # construct confidence interval
  ci <- array(NA_real_, dim = c(length(parm), 2L), dimnames = list(parm, cn))
  ci[] <- coef[parm] + se[parm] %o% q
  ci
}


## internal function to compute confidence intervals for estimated effects
## other than the indirect effect

get_ci_list <- function(object, level = 0.95, ...) UseMethod("get_ci_list")

get_ci_list.reg_fit_mediation <- function(object, level = 0.95,
                                          boot = NULL, ...) {
  # initializations
  have_serial <- object$model == "serial"
  # either extract confidence intervals from regression models or compute
  # bootstrap confidence intervals using the normal approximation
  if (is.null(boot)) {
    # further initializations
    alpha <- 1 - level
    x <- object$x
    p_x <- length(x)
    m <- object$m
    p_m <- length(m)
    have_yx <- !is.null(object$fit_yx)
    # compute confidence intervals of effects for a path
    if (p_m == 1L) confint_a <- confint(object$fit_mx, parm = x, level = level)
    else {
      confint_a <- lapply(object$fit_mx, confint, parm = x, level = level)
      confint_a <- do.call(rbind, confint_a)
    }
    # compute confidence intervals of effects for b path
    confint_b <- confint(object$fit_ymx, parm = m, level = level)
    # compute confidence intervals of total effects
    if (have_yx) {
      # extract confidence intervals from regression model
      confint_total <- confint(object$fit_yx, parm = x, level = level)
    } else {
      # confidence intervals not available
      confint_total <- matrix(NA_real_, nrow = p_x, ncol = 2L)
    }
    # compute confidence intervals of direct effects
    confint_direct <- confint(object$fit_ymx, parm = x, level = level)
    # construct list of confidence intervals
    if (have_serial) {
      # compute confidence intervals of effects for d path
      j_list <- lapply(seq_len(p_m-1L), function(j) 1L + seq_len(j))
      confint_d <- mapply(confint, object$fit_mx[-1L], parm = j_list,
                          SIMPLIFY = FALSE, USE.NAMES = FALSE)
      confint_d <- do.call(rbind, confint_d)
      # construct list
      ci_list <- list(a = confint_a, b = confint_b, d = confint_d,
                      total = confint_total, direct = confint_direct)
    } else {
      ci_list <- list(a = confint_a, b = confint_b, total = confint_total,
                      direct = confint_direct)
    }
  } else {
    # extract bootstrap replicates for effects other than the indirect effect
    keep <- c("a", "b", if (have_serial) "d", "total", "direct")
    boot_list <- extract_boot(object, boot = boot)[keep]
    # construct list of confidence intervals for the different effects
    ci_list <- lapply(boot_list, function(boot) {
      # compute means and standard errors of bootstrap replicates
      estimates <- colMeans(boot, na.rm = TRUE)
      se <- apply(boot, 2L, sd, na.rm = TRUE)
      # compute confidence intervals
      tmp <- mapply(confint_z, mean = estimates, sd = se,
                    MoreArgs = list(level = level),
                    SIMPLIFY = FALSE, USE.NAMES = FALSE)
      do.call(rbind, tmp)
    })
  }
  # return list of confidence intervals
  ci_list
}

get_ci_list.cov_fit_mediation <- function(object, level = 0.95,
                                          boot = NULL, ...) {
  # initializations
  alpha <- 1 - level
  # extract point estimates and standard errors
  if(is.null(boot)) {
    # combine point estimates
    estimates <- c(object$a, object$b, object$total, object$direct)
    # compute standard errors
    summary <- get_summary(object)
    se <- c(summary$a[1, 2], summary$b[1, 2], summary$total[1, 2],
            summary$direct[1, 2])
  } else {
    # compute means and standard errors from bootstrap replicates
    keep <- 1:4
    estimates <- colMeans(boot$t[, keep], na.rm = TRUE)
    se <- apply(boot$t[, keep], 2, sd, na.rm = TRUE)
  }
  # compute confidence intervals and combine into list
  list(a = confint_z(estimates[1], se[1], level = level),
       b = confint_z(estimates[2], se[2], level = level),
       total = confint_z(estimates[3], se[3], level = level),
       direct = confint_z(estimates[4], se[4], level = level))
}


## internal function to create names for confidence bounds
get_ci_names <- function(level = 0.95, alternative = "twosided") {
  if (alternative == "twosided") {
    alpha <- 1 - level
    paste(format(100 * c(alpha/2, 1 - alpha/2), trim = TRUE), "%")
  } else c("Lower", "Upper")
}


## internal function to compute confidence interval based on normal distribution
confint_z <- function(mean = 0, sd = 1, level = 0.95,
                      alternative = "twosided") {
  # initializations
  alpha <- 1 - level
  # compute confidence interval
  switch(alternative,
         twosided = qnorm(c(alpha/2, 1-alpha/2), mean = mean, sd = sd),
         less = c(-Inf, qnorm(level, mean = mean, sd = sd)),
         greater = c(qnorm(alpha, mean = mean, sd = sd), Inf))
}
