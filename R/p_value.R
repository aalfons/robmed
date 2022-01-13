# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

#' p-Values from (robust) mediation analysis
#'
#' Compute or extract the p-values for effects in (robust) mediation
#' analysis.
#'
#' For bootstrap tests, the p-value of the indirect effect is computed as the
#' smallest significance level \eqn{\alpha}{alpha} for which the
#' \eqn{(1 - \alpha) * 100\%}{(1 - alpha) * 100\%} confidence interval obtained
#' from the bootstrapped distribution does not contain 0.
#'
#' This is a simple implementation, where each digit after the comma is
#' determined via a grid search.  Hence computation time can be long if
#' confidence intervals are computed via the bias-corrected and accelerated
#' method (\code{"bca"}).
#'
#' For Sobel tests, the p-value of the indirect effect is already stored in the
#' object returned by \code{\link{test_mediation}()} and is simply extracted.
#'
#' @param object  an object inheriting from class
#' \code{"\link{test_mediation}"} containing results
#' from (robust) mediation analysis.
#' @param parm  an integer, character or logical vector specifying the paths
#' for which to extract or compute p-values, or \code{NULL} to extract or
#' compute p-values for all coefficients.  In case of a character vector,
#' possible values are \code{"a"}, \code{"b"}, \code{"d"} (only serial
#' multiple mediator models), \code{"total"}, \code{"direct"}, and
#' \code{"indirect"}.
#' @param type  a character string specifying how to compute the p-values of
#' the effects other than the indirect effect(s).  Possible values are
#' \code{"boot"} (the default) to compute bootstrap p-values using the normal
#' approximation (i.e., to assume a normal distribution of the corresponding
#' effect with the standard deviation computed from the bootstrap replicates),
#' or \code{"data"} to compute p-values via statistical theory based on the
#' original data (e.g., based on a t-distribution if the coefficients are
#' estimated via regression).  Note that this is only relevant for mediation
#' analysis via a bootstrap test, where the p-value of the indirect effect is
#' always computed as described in \sQuote{Details}.
#' @param digits  an integer determining how many digits to compute for the
#' p-values of the indirect effects (see \sQuote{Details}).  The default is to
#' compute 4 digits after the comma.
#' @param \dots  for the generic function, additional arguments to be passed
#' down to methods.  For the methods, additional arguments are currently
#' ignored.
#'
#' @return A numeric vector containing the requested p-values.
#'
#' @note The behavior of argument \code{parm} has changed in version 0.10.0.
#'
#' @author Andreas Alfons
#'
#' @seealso
#' \code{\link{test_mediation}()}, \code{\link[=coef.test_mediation]{coef}()},
#' \code{\link[=confint.test_mediation]{confint}()}
#'
#' @examples
#' data("BSG2014")
#'
#' \donttest{
#' # BCa intervals are recommended, but take a while to run
#' test_bca <- test_mediation(BSG2014,
#'                            x = "ValueDiversity",
#'                            y = "TeamCommitment",
#'                            m = "TaskConflict",
#'                            type = "bca")
#' p_value(test_bca)
#' }
#'
#' @keywords utilities
#'
#' @export

p_value <- function(object, ...) UseMethod("p_value")


#' @rdname p_value
#' @method p_value boot_test_mediation
#' @export

p_value.boot_test_mediation <- function(object, parm = NULL,
                                        type = c("boot", "data"),
                                        digits = 4L, ...) {
  # p-values of other effects
  type <- match.arg(type)
  if (type == "boot") {
    p_value_list <- get_p_value_list(object$fit, boot = object$reps)
  } else p_value_list <- get_p_value_list(object$fit)
  # add placeholder for p-values of indirect effects: those are only computed
  # if selected, as they take longer to compute
  p_value_list$indirect <- numeric()
  # if requested, take subset of p-values
  if (!is.null(parm)) p_value_list <- p_value_list[check_parm(parm)]
  # compute p-values of indirect effect if they are selected
  if (!is.null(p_value_list$indirect)) {
    # number of hypothesized mediators
    nr_indirect <- get_nr_indirect(length(object$fit$x), length(object$fit$m),
                                   model = object$fit$model)
    contrast <- object$fit$contrast          # only implemented for regression fit
    have_contrast <- is.character(contrast)  # but this always works
    # preparations to modify bootstrap object if contrasts are requested
    bootstrap <- object$reps
    indices_indirect <- if (nr_indirect == 1L) 1L else seq_len(1L + nr_indirect)
    # if contrasts are requested, modify bootstrap object to contain only
    # indirect effects and contrasts
    if (have_contrast) {
      # list of all combinations of indices of the relevant indirect effects
      combinations <- combn(indices_indirect[-1L], 2, simplify = FALSE)
      # modify bootstrap object to add contrasts
      bootstrap$t0 <- c(bootstrap$t0[indices_indirect],
                        get_contrasts(bootstrap$t0, combinations,
                                      type = contrast))
      bootstrap$t <- cbind(bootstrap$t[, indices_indirect],
                           get_contrasts(bootstrap$t, combinations,
                                         type = contrast))
      # update indices of indirect effects
      indices_indirect <- seq_along(bootstrap$t0)
    }
    # compute p-value of indirect effects as the smallest significance
    # level where 0 is not in the confidence interval
    p_value_list$indirect <- sapply(indices_indirect, function(j) {
      p_value(bootstrap, parm = j, digits = digits,
              alternative = object$alternative,
              type = object$type)
    })
  }
  # convert list to vector
  p_values <- unlist(p_value_list, use.names = FALSE)
  # add names
  keep <- names(p_value_list)
  names(p_values) <- get_effect_names(effects = object[keep])
  # return p-values
  p_values
}


#' @rdname p_value
#' @method p_value sobel_test_mediation
#' @export

p_value.sobel_test_mediation <- function(object, parm = NULL, ...) {
  # p-values of other effects
  p_value_list <- get_p_value_list(object$fit)
  # add p-value of indirect effect
  p_value_list$indirect <- object$p_value
  # if requested, take subset of p-values
  if (!is.null(parm)) p_value_list <- p_value_list[check_parm(parm)]
  # convert list to vector
  p_values <- unlist(p_value_list, use.names = FALSE)
  # add names
  keep <- names(p_value_list)
  names(p_values) <- get_effect_names(effects = object$fit[keep])
  # return p-values
  p_values
}


## methods to extract p-values from model fits

#' @export
p_value.lm <- function(object, parm = NULL, ...) {
  # compute the usual summary and extract coefficient matrix
  summary <- summary(object)
  coef_mat <- coef(summary)
  coef_names <- rownames(coef_mat)
  # extract p-values
  p_values <- coef_mat[, 4L]
  if (is.null(names(p_values))) names(p_values) <- coef_names
  # check parameters to extract
  if (missing(parm) || is.null(parm)) parm <- coef_names
  else if (is.numeric(parm)) parm <- coef_names[parm]
  # extract p-values of selected parameters
  p_values[parm]
}

#' @export
p_value.lmrob <- function(object, parm = NULL, ...) {
  # compute the usual summary and extract coefficient matrix
  summary <- summary(object)
  coef_mat <- coef(summary)
  coef_names <- rownames(coef_mat)
  # extract p-values
  p_values <- coef_mat[, 4L]
  if (is.null(names(p_values))) names(p_values) <- coef_names
  # check parameters to extract
  if (missing(parm) || is.null(parm)) parm <- coef_names
  else if (is.numeric(parm)) parm <- coef_names[parm]
  # extract p-values of selected parameters
  p_values[parm]
}

#' @export
p_value.rq <- function(object, parm = NULL, ...) {
  # compute the usual summary and extract coefficient matrix
  summary <- summary(object, se = "iid")
  coef_mat <- coef(summary)
  coef_names <- rownames(coef_mat)
  # extract p-values
  p_values <- coef_mat[, 4L]
  if (is.null(names(p_values))) names(p_values) <- coef_names
  # check parameters to extract
  if (missing(parm) || is.null(parm)) parm <- coef_names
  else if (is.numeric(parm)) parm <- coef_names[parm]
  # extract p-values of selected parameters
  p_values[parm]
}


# extract p-value from bootstrap results
# (argument 'parm' can be used for completeness; we only need the p-value for
# the indirect effect in the first column of the bootstrap results)
p_value.boot <- function(object, parm = 1L, digits = 4L,
                         alternative = c("twosided", "less", "greater"),
                         type = c("bca", "perc"), ...) {
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
      ci <- confint(object, parm = parm, level = 1 - alpha,
                    alternative = alternative, type = type)
      # reject if 0 is not in the confidence interval
      reject <- prod(ci) > 0
    }
  }
  # return smallest significance level where 0 is not in the confidence interval
  alpha
}


## internal function to compute p-values for estimated effects

get_p_value_list <- function(object, ...) UseMethod("get_p_value_list")

get_p_value_list.cov_fit_mediation <- function(object, boot = NULL, ...) {
  # extract p-values
  if(is.null(boot)) {
    # compute summary of model fit
    summary <- get_summary(object)
    # extract p-values as list
    list(a = summary$a[1, 4], b = summary$b[1, 4], total = summary$total[1, 4],
         direct = summary$direct[1, 4])
  } else {
    # compute means, standard errors and z-statistics from bootstrap replicates
    keep <- c(3L, 5L:7L)
    estimates <- colMeans(boot$t[, keep], na.rm = TRUE)
    se <- apply(boot$t[, keep], 2, sd, na.rm = TRUE)
    z <- estimates / se
    # compute p-values and return as list
    p_values <- p_value_z(z)
    list(a = p_values[1], b = p_values[2],
         total = p_values[4], direct = p_values[3])
  }
}

get_p_value_list.reg_fit_mediation <- function(object, boot = NULL, ...) {
  # initializations
  p_x <- length(object$x)
  p_m <- length(object$m)
  model <- object$model
  # extract point estimates and standard errors
  if(is.null(boot)) {
    # compute p-values for a path
    if(p_m == 1L) p_value_a <- p_value(object$fit_mx, parm = 1L + seq_len(p_x))
    else if (model == "serial") {
      p_value_a <- mapply(p_value, object$fit_mx, parm = 1L + seq_len(p_m),
                          USE.NAMES = FALSE)
    } else {
      p_value_a <- sapply(object$fit_mx, p_value, parm = 1L + seq_len(p_x),
                          USE.NAMES = FALSE)
    }
    # for serial multiple mediator models, compute p-values for d path
    if (model == "serial") {
      j_list <- lapply(seq_len(p_m-1L), function(j) 1L + seq_len(j))
      p_value_d <- mapply(p_value, object$fit_mx[-1L], parm = j_list,
                          SIMPLIFY = FALSE, USE.NAMES = FALSE)
      p_value_d <- unlist(p_value_d, use.names = FALSE)
    }
    # compute p-values for b path and direct effect
    p_value_b <- p_value(object$fit_ymx, parm = 1L + seq_len(p_m))
    p_value_direct <- p_value(object$fit_ymx, parm = 1L + p_m + seq_len(p_x))
    # compute p-value for total effect
    if(is.null(object$fit_yx)) {
      # p-value not available
      p_value_total <- rep.int(NA_real_, p_x)
    } else {
      # extract p-value from regression model
      p_value_total <- p_value(object$fit_yx, parm = 1L + seq_len(p_x))
    }
    # combine p-values into list
    if (model == "serial") {
      p_value_list <- list(a = p_value_a, b = p_value_b, d = p_value_d,
                           total = p_value_total, direct = p_value_direct)
    } else {
      p_value_list <- list(a = p_value_a, b = p_value_b, total = p_value_total,
                           direct = p_value_direct)
    }
  } else {
    # get indices of columns of bootstrap replicates that that correspond to
    # the respective models
    p_covariates <- length(object$covariates)
    index_list <- get_index_list(p_x, p_m, p_covariates, model = model)
    # keep indices for a path in model m ~ x + covariates
    if(p_m == 1L) keep_a <- index_list$fit_mx[1L + seq_len(p_x)]
    else if (model == "serial") {
      keep_a <- mapply("[", index_list$fit_mx, 1L + seq_len(p_m))
    } else keep_a <- sapply(index_list$fit_mx, "[", 1L + seq_len(p_x))
    # for serial multiple mediators, keep indices for d path
    if (model == "serial") {
      j_list <- lapply(seq_len(p_m-1L), function(j) 1L + seq_len(j))
      keep_d <- unlist(mapply("[", index_list$fit_mx[-1L], parm = j_list,
                              SIMPLIFY = FALSE, USE.NAMES = FALSE))
    }
    # keep indeces of b path and direct effect in model y ~ m + x + covariates
    keep_b <- index_list$fit_ymx[1L + seq_len(p_m)]
    keep_direct <- index_list$fit_ymx[1L + p_m + seq_len(p_x)]
    # create list of indices to keep for each path
    # (index of total effect is stored separately in 'index_list')
    if (model == "serial") {
      keep_list <- list(a = keep_a, b = keep_b, d = keep_d,
                        total = index_list$total, direct = keep_direct)
    } else {
      keep_list <- list(a = keep_a, b = keep_b, total = index_list$total,
                     direct = keep_direct)
    }
    # construct list of p_values
    p_value_list <- lapply(keep_list, function(keep) {
      # compute means, standard errors and z-statistics from bootstrap replicates
      estimates <- colMeans(boot$t[, keep, drop = FALSE], na.rm = TRUE)
      se <- apply(boot$t[, keep, drop = FALSE], 2L, sd, na.rm = TRUE)
      z <- estimates / se
      # compute p-values
      p_value_z(z)
    })
  }
  # return list of p-values
  p_value_list
}


## internal function to compute p-value based on normal distribution
p_value_z <- function(z, alternative = c("twosided", "less", "greater")) {
  # initializations
  alternative <- match.arg(alternative)
  # compute p-value
  switch(alternative, twosided = 2 * pnorm(abs(z), lower.tail = FALSE),
         less = pnorm(z), greater = pnorm(z, lower.tail = FALSE))
}
