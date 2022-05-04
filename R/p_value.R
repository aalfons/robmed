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
#' boot <- test_mediation(BSG2014,
#'                        x = "ValueDiversity",
#'                        y = "TeamCommitment",
#'                        m = "TaskConflict",
#'                        type = "bca")
#' p_value(boot)
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
    # regression fits store the bootstrap replicates of the regression
    # coefficients, whereas covariance fits directly store the bootstrap
    # replicates of the effects in the mediation model
    if (inherits(object$fit, "reg_fit_mediation")) {
      # extract bootstrap replicates for the indirect effects
      boot_indirect <- extract_boot(object$fit, boot = object$reps)$indirect
      # compute p-value of indirect effects as the smallest significance
      # level where 0 is not in the confidence interval
      p_value_list$indirect <- boot_p_value(object$fit$indirect, boot_indirect,
                                            object = object$reps,
                                            digits = digits,
                                            alternative = object$alternative,
                                            type = object$type)
    } else if (inherits(object$fit, "cov_fit_mediation")) {
      # compute p-value of indirect effects as the smallest significance
      # level where 0 is not in the confidence interval
      p_value_list$indirect <- extract_p_value(parm = 5L,
                                               object = object$reps,
                                               digits = digits,
                                               alternative = object$alternative,
                                               type = object$type)
    } else stop("not implemented for this type of model fit")
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


## internal function to compute p-values for estimated effects other than the
## indirect effect

get_p_value_list <- function(object, ...) UseMethod("get_p_value_list")

get_p_value_list.reg_fit_mediation <- function(object, boot = NULL, ...) {
  # initializations
  have_serial <- object$model == "serial"
  # either extract p-values from regression models or compute bootstrap
  # p-values using the normal approximation
  if(is.null(boot)) {
    # further initializations
    x <- object$x
    p_x <- length(x)
    m <- object$m
    p_m <- length(m)
    have_yx <- !is.null(object$fit_yx)
    # compute p-values of effects for a path
    if (p_m == 1L) p_value_a <- p_value(object$fit_mx, parm = x)
    else {
      p_value_a <- sapply(object$fit_mx, p_value, parm = x, USE.NAMES = FALSE)
    }
    # compute p-values of effects for b path
    p_value_b <- p_value(object$fit_ymx, parm = m)
    # compute p-values of total effects
    if (have_yx) {
      # extract p-values from regression model
      p_value_total <- p_value(object$fit_yx, parm = x)
    } else {
      # p-values not available
      p_value_total <- rep.int(NA_real_, p_x)
    }
    # compute p-values of direct effects
    p_value_direct <- p_value(object$fit_ymx, parm = x)
    # construct list of p-values
    if (have_serial) {
      # compute p-values of effects for d path
      j_list <- lapply(seq_len(p_m-1L), function(j) 1L + seq_len(j))
      p_value_d <- mapply(p_value, object$fit_mx[-1L], parm = j_list,
                          SIMPLIFY = FALSE, USE.NAMES = FALSE)
      p_value_d <- unlist(p_value_d, use.names = FALSE)
      # construct list
      p_value_list <- list(a = p_value_a, b = p_value_b, d = p_value_d,
                           total = p_value_total, direct = p_value_direct)
    } else {
      p_value_list <- list(a = p_value_a, b = p_value_b, total = p_value_total,
                           direct = p_value_direct)
    }
  } else {
    # extract bootstrap replicates for effects other than the indirect effect
    keep <- c("a", "b", if (have_serial) "d", "total", "direct")
    boot_list <- extract_boot(object, boot = boot)[keep]
    # construct list of p_values for the different effects
    p_value_list <- lapply(boot_list, function(boot) {
      # compute means, standard errors, and z-statistics of bootstrap replicates
      estimates <- colMeans(boot, na.rm = TRUE)
      se <- apply(boot, 2L, sd, na.rm = TRUE)
      z <- estimates / se
      # compute p-values
      p_value_z(z)
    })
  }
  # return list of p-values
  p_value_list
}

get_p_value_list.cov_fit_mediation <- function(object, boot = NULL, ...) {
  # extract p-values
  if(is.null(boot)) {
    # compute summary of model fit
    summary <- get_summary(object)
    # extract p-values as list
    list(a = summary$a[1, 4], b = summary$b[1, 4], total = summary$total[1, 4],
         direct = summary$direct[1, 4])
  } else {
    # compute means, standard errors, and z-statistics from bootstrap replicates
    keep <- 1:4
    estimates <- colMeans(boot$t[, keep], na.rm = TRUE)
    se <- apply(boot$t[, keep], 2, sd, na.rm = TRUE)
    z <- estimates / se
    # compute p-values and return as list
    p_values <- p_value_z(z)
    list(a = p_values[1], b = p_values[2], total = p_values[3],
         direct = p_values[4])
  }
}


## internal function to compute p-value based on normal distribution
p_value_z <- function(z, alternative = "twosided") {
  switch(alternative, twosided = 2 * pnorm(abs(z), lower.tail = FALSE),
         less = pnorm(z), greater = pnorm(z, lower.tail = FALSE))
}
