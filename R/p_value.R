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
#' @param parm  an integer, character or logical vector specifying the
#' coefficients for which to extract or compute p-values, or
#' \code{NULL} to extract or compute p-values for all coefficients.
#' @param type  a character string specifying how to compute the p-values of
#' the effects other than the indirect effect(s).  Possible values are
#' \code{"boot"} (the default) to compute bootstrap p-values using the normal
#' approximation (i.e., to assume a normal distribution of the corresponding
#' effect with the standard deviation computed from the bootstrap replicates),
#' or \code{"data"} to compute p-values via statistical theory based on the
#' original data (e.g., based on a t-distribution the coefficients are
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
  # number of hypothesized mediators
  p_m <- length(object$fit$m)
  # p-values of other effects
  type <- match.arg(type)
  if(type == "boot") {
    p_values <- get_p_value(object$fit, boot = object$reps)
  } else p_values <- get_p_value(object$fit)
  # add temporary NA's for p-values of indirect effect,
  # as those take longer to compute
  if (p_m == 1) {
    # only one mediator
    indirect_names <- "ab"
    alpha <- NA_real_
  } else {
    # multiple mediators
    indirect_names <- paste("ab", rownames(object$ci), sep = "_")
    alpha <- rep.int(NA_real_, nrow(object$ci))
  }
  names(alpha) <- indirect_names
  p_values <- c(p_values, alpha)
  # if requested, take subset of effects
  if(!is.null(parm)) {
    p_values <- p_values[parm]
    # check which p-values of indirect effects need to be computed
    which_names <- grep("ab", names(p_values), value = TRUE)
    which_indices <- match(which_names, indirect_names, nomatch = integer())
  } else {
    which_names <- indirect_names
    which_indices <- seq_along(which_names)
  }
  # compute p-value of requested indirect effects as the smallest significance
  # level where 0 is not in the confidence interval
  if (length(which_names) > 0) {
    p_values[which_names] <- sapply(which_indices, function(j) {
      p_value(object$reps, parm = j, digits = digits,
              alternative = object$alternative,
              type = object$type)
    })
  }
  p_values
}

# p_value.boot_test_mediation <- function(object, digits = 4L, ...) {
#   # number of hypothesized mediators
#   p_m <- length(object$fit$m)
#   # compute p-value
#   if(p_m == 1L) {
#     # only one mediator
#     # set lower bound of significance level to 0
#     lower <- 0
#     # loop over the number of digits and determine the corresponding digit after
#     # the comma of the p-value
#     for (digit in seq_len(digits)) {
#       # set step size
#       step <- 1 / 10^digit
#       # reset the significance level to the lower bound as we continue from there
#       # with a smaller stepsize
#       alpha <- lower
#       # there is no rejection at the lower bound, so increase significance level
#       # until there is rejection
#       reject <- FALSE
#       while(!reject) {
#         # update lower bound and significance level
#         lower <- alpha
#         alpha <- alpha + step
#         # retest at current significance level and extract confidence interval
#         ci <- confint(object$reps, parm = 1L, level = 1 - alpha,
#                       alternative = object$alternative, type = object$type)
#         # reject if 0 is not in the confidence interval
#         reject <- prod(ci) > 0
#       }
#     }
#   } else {
#     # multiple mediators
#     rn <- rownames(object$ci)
#     alpha <- sapply(seq_along(rn), function(j) {
#       # set lower bound of significance level to 0
#       lower <- 0
#       # loop over the number of digits and determine the corresponding digit after
#       # the comma of the p-value
#       for (digit in seq_len(digits)) {
#         # set step size
#         step <- 1 / 10^digit
#         # reset the significance level to the lower bound as we continue from there
#         # with a smaller stepsize
#         alpha_j <- lower
#         # there is no rejection at the lower bound, so increase significance level
#         # until there is rejection
#         reject <- FALSE
#         while(!reject) {
#           # update lower bound and significance level
#           lower <- alpha_j
#           alpha_j <- alpha_j + step
#           # retest at current significance level and extract confidence interval
#           # ci <- retest(object, level = 1 - alpha)$ci[m, ]
#           ci <- confint(object$reps, parm = j, level = 1 - alpha_j,
#                         alternative = object$alternative, type = object$type)
#           # reject if 0 is not in the confidence interval
#           reject <- prod(ci) > 0
#         }
#       }
#       # return p-value for current indirect effect
#       alpha_j
#     })
#     names(alpha) <- rn
#   }
#   # return smallest significance level where 0 is not in the confidence interval
#   alpha
# }


#' @rdname p_value
#' @method p_value sobel_test_mediation
#' @export

p_value.sobel_test_mediation <- function(object, parm = NULL, ...) {
  # old behavior by default if new arguments are missing
  if (missing(parm)) {
    warning("default behavior will change in a future version, see the ",
            sQuote("Note"), " section of the help file")
    object$p_value
  } else {
    # combine p-value of indirect effect with that of other effects
    p_values <- c(get_p_value(object$fit), ab = object$p_value)
    # if requested, take subset of effects
    if(!is.null(parm)) p_values <- p_values[parm]
    p_values
  }
}

# p_value.sobel_test_mediation <- function(object, ...) object$p_value


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


## internal function to compute p-values for estimated effects

get_p_value <- function(object, parm, ...) {
  UseMethod("get_p_value")
}

get_p_value.cov_fit_mediation <- function(object, parm = NULL, boot = NULL,
                                          ...) {
  # extract p-values
  if(is.null(boot)) {
    # combine point estimates
    estimates <- c(object$a, object$b, object$direct, object$total)
    # compute standard errors
    summary <- get_summary(object)
    p_values <- c(summary$a[1, 4], summary$b[1, 4], summary$direct[1, 4],
                  summary$total[1, 4])
    names(p_values) <- c("a", "b", "Direct", "Total")
  } else {
    # compute means, standard errors and z-statistics from bootstrap replicates
    keep <- c(3L, 5L:7L)
    estimates <- colMeans(boot$t[, keep], na.rm = TRUE)
    se <- apply(boot$t[, keep], 2, sd, na.rm = TRUE)
    z <- estimates / se
    # compute p-values
    p_values <- p_value_z(z)
    names(p_values) <- c("a", "b", "Direct", "Total")
  }
  # if requested, take subset of effects
  if(!is.null(parm)) p_values <- p_values[parm]
  p_values
}

get_p_value.reg_fit_mediation <- function(object, parm = NULL, boot = NULL,
                                          ...) {
  # initializations
  p_m <- length(object$m)
  # extract point estimates and standard errors
  if(is.null(boot)) {
    # extract p-values from regression models
    if(p_m == 1L) p_value_mx <- p_value(object$fit_mx, parm = 2L)
    else p_value_mx <- sapply(object$fit_mx, p_value, parm = 2L)
    p_value_ymx <- p_value(object$fit_ymx, parm = 1L + seq_len(p_m + 1L))
    # compute p-valuel for total effect
    if(is_robust(object)) {
      # p-value not available
      p_value_yx <- NA_real_
    } else {
      # extract p-value from regression model
      p_value_yx <- p_value(object$fit_yx, parm = 2L)
    }
    # combine p-values
    p_values <- c(p_value_mx, p_value_ymx, p_value_yx)
    names(p_values) <- get_effect_names(object$m)
  } else {
    # get indices of columns of bootstrap replicates that that correspond to
    # the respective models
    p_covariates <- length(object$fit$covariates)
    index_list <- get_index_list(p_m, p_covariates)
    # the a path is the second coefficient in the model m ~ x + covariates
    if (p_m == 1) keep_mx <- index_list$fit_mx[2L]
    else keep_mx <- sapply(index_list$fit_mx, "[", 2L)
    # keep b and c coefficients of model y ~ m + x + covariates
    keep_ymx <- index_list$fit_ymx[1L + seq_len(p_m + 1)]
    # index of total effect is stored separately in this list
    keep <- c(keep_mx, keep_ymx, index_list$total)
    # compute means, standard errors and z-statistics from bootstrap replicates
    estimates <- colMeans(boot$t[, keep], na.rm = TRUE)
    se <- apply(boot$t[, keep], 2L, sd, na.rm = TRUE)
    z <- estimates / se
    # compute p-values
    p_values <- p_value_z(z)
    names(p_values) <- get_effect_names(object$m)
  }
  # if requested, take subset of effects
  if(!is.null(parm)) p_values <- p_values[parm]
  p_values
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
    while(!reject) {
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
