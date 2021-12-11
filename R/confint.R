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
#' @param parm  an integer, character or logical vector specifying the
#' coefficients for which to extract or compute confidence intervals, or
#' \code{NULL} to extract or compute confidence intervals for all coefficients.
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
#' # run standard bootstrap test
#' standard_boot <- test_mediation(BSG2014,
#'                                 x = "ValueDiversity",
#'                                 y = "TeamCommitment",
#'                                 m = "TaskConflict",
#'                                 robust = FALSE)
#' confint(standard_boot, type = "data")
#'
#' @keywords utilities

NULL


#' @rdname confint.test_mediation
#' @method confint boot_test_mediation
#' @export

## argument 'level' is ignored
confint.boot_test_mediation <- function(object, parm = NULL, level = NULL,
                                        type = c("boot", "data"), ...) {
  # number of hypothesized mediators
  nr_indirect <- get_nr_indirect(length(object$fit$x), length(object$fit$m),
                                 model = object$fit$model)
  # confidence interval of other effects
  type <- match.arg(type)
  if (type == "boot") {
    ci <- get_confint(object$fit, level = object$level, boot = object$reps)
  } else ci <- get_confint(object$fit, level = object$level)
  # combine with confidence interval of indirect effect
  if (nr_indirect == 1L) ci <- rbind(ci, "Indirect" = object$ci)
  else {
    ci_indirect <- object$ci
    rownames(ci_indirect) <- paste("Indirect", rownames(ci_indirect), sep = "_")
    ci <- rbind(ci, ci_indirect)
  }
  if (object$alternative != "twosided") colnames(ci) <- c("Lower", "Upper")
  # if requested, take subset of effects
  if (!is.null(parm)) ci <- ci[parm, , drop = FALSE]
  ci
}


#' @rdname confint.test_mediation
#' @method confint sobel_test_mediation
#' @export

confint.sobel_test_mediation <- function(object, parm = NULL, level = 0.95,
                                         ...) {
  # initializations
  level <- rep(as.numeric(level), length.out = 1)
  if(is.na(level) || level < 0 || level > 1) level <- formals()$level
  # confidence interval of indirect effect
  ci <- confint_z(object$fit$indirect, object$se, level = level,
                  alternative = object$alternative)
  # combine with confidence intervalse of other effects
  ci <- rbind(get_confint(object$fit, level = level), "Indirect" = ci)
  if(object$alternative != "twosided") colnames(ci) <- c("Lower", "Upper")
  # if requested, take subset of effects
  if(!is.null(parm)) ci <- ci[parm, , drop = FALSE]
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
  cn <- paste(format(100 * a, trim=TRUE), "%")
  # construct confidence interval
  ci <- array(NA_real_, dim = c(length(parm), 2L), dimnames = list(parm, cn))
  ci[] <- coef[parm] + se[parm] %o% q
  ci
}


## internal function to compute confidence intervals for estimated effects

get_confint <- function(object, parm, level = 0.95, ...) {
  UseMethod("get_confint")
}

get_confint.cov_fit_mediation <- function(object, parm = NULL, level = 0.95,
                                          boot = NULL, ...) {
  # initializations
  alpha <- 1 - level
  # extract point estimates and standard errors
  if(is.null(boot)) {
    # combine point estimates
    estimates <- c(object$a, object$b, object$direct, object$total)
    # compute standard errors
    summary <- get_summary(object)
    se <- c(summary$a[1, 2], summary$b[1, 2], summary$direct[1, 2],
            summary$total[1, 2])
  } else {
    # compute means and standard errors from bootstrap replicates
    keep <- c(3L, 5L:7L)
    estimates <- colMeans(boot$t[, keep], na.rm = TRUE)
    se <- apply(boot$t[, keep], 2, sd, na.rm = TRUE)
  }
  # compute confidence intervals and combine into one matrix
  ci <- rbind(confint_z(estimates[1], se[1], level=level),
              confint_z(estimates[2], se[2], level=level),
              confint_z(estimates[4], se[4], level=level),
              confint_z(estimates[3], se[3], level=level))
  rn <- get_effect_names(effects = object[c("a", "b", "total", "direct")])
  cn <- paste(format(100 * c(alpha/2, 1 - alpha/2), trim = TRUE), "%")
  dimnames(ci) <- list(rn, cn)
  # if requested, take subset of effects
  if(!is.null(parm)) ci <- ci[parm, , drop = FALSE]
  ci
}

get_confint.reg_fit_mediation <- function(object, parm = NULL, level = 0.95,
                                          boot = NULL, ...) {
  # initializations
  alpha <- 1 - level
  p_x <- length(object$x)
  p_m <- length(object$m)
  model <- object$model
  # row names to be used for confidence intervals
  keep <- c("a", "b", "d", "total", "direct")
  rn <- get_effect_names(effects = object[keep])
  # compute confidence intervals
  if (is.null(boot)) {
    # compute confidence intervals for a path
    if (p_m == 1L) {
      # single mediator
      confint_a <- confint(object$fit_mx, parm = 1L + seq_len(p_x),
                           level = level)
    } else {
      # multiple mediators
      if (model == "serial") {
        confint_a <- mapply(confint, object$fit_mx, parm = 1L + seq_len(p_m),
                            SIMPLIFY = FALSE, USE.NAMES = FALSE)
      } else {
        confint_a <- lapply(object$fit_mx, confint, parm = 1L + seq_len(p_x),
                            level = level)
      }
      confint_a <- do.call(rbind, confint_a)
    }
    # for serial multiple mediators, compute confidence intervals for d path
    if (model == "serial") {
      j_list <- lapply(seq_len(p_m-1L), function(j) 1L + seq_len(j))
      confint_d <- mapply(confint, object$fit_mx[-1L], parm = j_list,
                          SIMPLIFY = FALSE, USE.NAMES = FALSE)
      confint_d <- do.call(rbind, confint_d)
    } else confint_d <- NULL
    # compute confidence intervals for b path and direct effect
    confint_b <- confint(object$fit_ymx, parm = 1L + seq_len(p_m),
                         level = level)
    confint_direct <- confint(object$fit_ymx, parm = 1L + p_m + seq_len(p_x),
                              level = level)
    # compute confidence intervals for total effect
    if (is.null(object$fit_yx)) {
      # confidence interval not available
      confint_total <- matrix(NA_real_, nrow = p_x, ncol = 2L)
    } else {
      # extract confidence interval from regression model
      confint_total <- confint(object$fit_yx, parm = 1L + seq_len(p_x),
                               level = level)
    }
    # combine confidence intervals
    ci <- rbind(confint_a, confint_b, confint_d, confint_total, confint_direct)
    rownames(ci) <- rn
  } else {
    # get indices of columns of bootstrap replicates that that correspond to
    # the respective models
    p_covariates <- length(object$fit$covariates)
    index_list <- get_index_list(p_x, p_m, p_covariates, model = model)
    # keep indices for a path in model m ~ x + covariates
    if (p_m == 1L) keep_a <- index_list$fit_mx[1L + seq_len(p_x)]
    else if (model == "serial") {
      keep_a <- mapply("[", index_list$fit_mx, 1L + seq_len(p_m),
                       USE.NAMES = FALSE)
    } else keep_a <- sapply(index_list$fit_mx, "[", 1L + seq_len(p_x))
    # for serial multiple mediators, keep indices for d path
    if (model == "serial") {
      j_list <- lapply(seq_len(p_m-1L), function(j) 1L + seq_len(j))
      keep_d <- unlist(mapply("[", index_list$fit_mx[-1L], parm = j_list,
                              SIMPLIFY = FALSE, USE.NAMES = FALSE))
    } else keep_d <- NULL
    # keep indeces of b path and direct effect in model y ~ m + x + covariates
    keep_b <- index_list$fit_ymx[1L + seq_len(p_m)]
    keep_direct <- index_list$fit_ymx[1L + p_m + seq_len(p_x)]
    # index of total effect is stored separately in this list
    keep <- c(keep_a, keep_b, keep_d, index_list$total, keep_direct)
    # compute means and standard errors from bootstrap replicates
    estimates <- colMeans(boot$t[, keep], na.rm = TRUE)
    se <- apply(boot$t[, keep], 2L, sd, na.rm = TRUE)
    # compute confidence intervals and combine into one matrix
    ci <- mapply(confint_z, mean = estimates, sd = se,
                 MoreArgs = list(level = level),
                 SIMPLIFY = FALSE, USE.NAMES = FALSE)
    ci <- do.call(rbind, ci)
    # add row and column names
    cn <- paste(format(100 * c(alpha/2, 1 - alpha/2), trim = TRUE), "%")
    dimnames(ci) <- list(rn, cn)
  }
  # if requested, take subset of effects
  if(!is.null(parm)) ci <- ci[parm, , drop = FALSE]
  ci
}


# extract confidence interval from bootstrap results
# (argument 'parm' can be used for completeness; we only need the confidence
# interval for the indirect effect in the first column of the bootstrap results)
confint.boot <- function(object, parm = 1L, level = 0.95,
                         alternative = c("twosided", "less", "greater"),
                         type = c("bca", "perc"), ...) {
  # initializations
  alternative <- match.arg(alternative)
  type <- match.arg(type)
  component <- if(type == "perc") "percent" else type
  # extract confidence interval
  if(level == 0) {
    ci <- rep.int(mean(object$t[, parm], na.rm=TRUE), 2L)
  } else if(level == 1) {
    ci <- c(-Inf, Inf)
  } else if(alternative == "twosided") {
    ci <- boot.ci(object, conf=level, type=type, index=parm)[[component]][4:5]
  } else {
    alpha <- 1 - level
    ci <- boot.ci(object, conf=1-2*alpha, type=type, index=parm)[[component]][4:5]
    if(alternative == "less") ci[1] <- -Inf
    else ci[2] <- Inf
  }
  # return confidence interval
  ci
}

## internal function to compute confidence interval based on normal distribution
confint_z <- function(mean = 0, sd = 1, level = 0.95,
                      alternative = c("twosided", "less", "greater")) {
  # initializations
  alternative <- match.arg(alternative)
  # compute confidence interval
  alpha <- 1 - level
  switch(alternative,
         twosided = qnorm(c(alpha/2, 1-alpha/2), mean = mean, sd = sd),
         less = c(-Inf, qnorm(level, mean = mean, sd = sd)),
         greater = c(qnorm(alpha, mean = mean, sd = sd), Inf))
}
