# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

#' Retest for mediation
#'
#' Re-perform a (fast and robust) bootstrap test or Sobel's test for the
#' indirect effect(s) based on results from (robust) mediation analysis.
#'
#' @param object  an object inheriting from class
#' \code{"\link{test_mediation}"} containing results from (robust) mediation
#' analysis.
#' @param alternative  a character string specifying the alternative hypothesis
#' in the test for the indirect effect.  Possible values are \code{"twosided"},
#' \code{"less"} or \code{"greater"}.
#' @param level  numeric; the confidence level of the confidence interval in
#' the bootstrap test.
#' @param type  a character string specifying the type of confidence interval
#' to be computed in the bootstrap test.  Possible values are \code{"bca"} for
#' the bias-corrected and accelerated bootstrap, or \code{"perc"} for the
#' percentile bootstrap.
#' @param contrast  a logical indicating whether to compute pairwise contrasts
#' of the indirect effects.  This can also be a character string, with
#' \code{"estimates"} for computing the differences of the indirect effects
#' (such that it is tested whether two indirect effects are equal), and
#' \code{"absolute"} for computing the differences of the absolute values of
#' the indirect effects (such that it is tested whether two indirect effects
#' are equal in magnitude).  This is only relevant for models with multiple
#' hypothesized mediators, which are currently only implemented for bootstrap
#' tests and estimation via regressions.
#' @param \dots  additional arguments to be passed down to methods.
#'
#' @return An object of the same class as \code{object} with updated test
#' results (see \code{\link{test_mediation}()}).
#'
#' @note From version 0.8.0 onwards, the behavior of this function changed.
#' For arguments that are not supplied, the corresponding values of
#' \code{object} are now used as defaults.
#'
#' @author Andreas Alfons
#'
#' @seealso \code{\link{test_mediation}()}
#'
#' @examples
#' data("BSG2014")
#'
#' # run fast and robust bootstrap test
#' test <- test_mediation(BSG2014,
#'                        x = "ValueDiversity",
#'                        y = "TeamCommitment",
#'                        m = "TaskConflict")
#' summary(test)
#'
#' # now compute 97.5% confidence interval
#' retest(test, level = 0.975)
#'
#' @keywords multivariate
#'
#' @export

retest <- function(object, ...) UseMethod("retest")


#' @rdname retest
#' @method retest boot_test_mediation
#' @export

retest.boot_test_mediation <- function(object, alternative, level,
                                       type, contrast, ...) {
  # initializations
  m <- object$fit$m
  p_m <- length(m)
  have_covariance <- inherits(object$fit, "cov_fit_mediation")
  # check alternative hypothesis
  if (missing(alternative)) alternative <- object$alternative
  else {
    alternative <- match.arg(alternative,
                             choices = c("twosided", "less", "greater"))
  }
  # check confidence level
  if (missing(level)) level <- object$level
  else {
    level <- rep(as.numeric(level), length.out = 1)
    if(is.na(level) || level < 0 || level > 1) level <- object$level
  }
  # check type of confidence intervals
  if (missing(type)) type <- object$type
  else type <- match.arg(type, choices = c("bca", "perc"))
  # check contrasts of indirect effect
  if (inherits(object$fit, "reg_fit_mediation")) {
    # regression model fit (multiple mediators and contrasts are supported)
    if (missing(contrast)) contrast <- object$fit$contrast
    else {
      if (p_m == 1L) contrast <- FALSE
      else if (is.logical(contrast)) {
        contrast <- isTRUE(contrast)
        if (contrast) contrast <- "estimates"
      } else {
        contrast <- match.arg(contrast, choices = c("estimates", "absolute"))
      }
    }
    have_contrast <- is.character(contrast)
    update_contrast <- contrast != object$fit$contrast
  } else if (inherits(object$fit, "cov_fit_mediation")) {
    # covariance model fit (only implemented for a simple mediation model)
    have_contrast <- FALSE
    update_contrast <- FALSE
  } else stop("not implemented for this type of model fit")
  # check for any new arguments
  update_alternative <- alternative != object$alternative
  update_level <- level != object$level
  update_type <- type != object$type
  update <- update_alternative || update_level || update_type || update_contrast
  # reperform test if necessary
  if (update) {
    # recompute confidence interval
    if(p_m == 1L) {
      # only one mediator
      ci <- confint(object$reps, parm = 1L, level = level,
                    alternative = alternative, type = type)
    } else {
      # multiple mediators
      indices_ab <- seq_len(1L + p_m)
      bootstrap <- object$reps
      ci <- lapply(indices_ab, function(j) {
        confint(bootstrap, parm = j, level = level,
                alternative = alternative, type = type)
      })
      ci <- do.call(rbind, ci)
      # if requested, compute contrasts of indirect effects
      if (have_contrast) {
        # list of all combinations of indices of the relevant indirect effects
        combinations <- combn(indices_ab[-1], 2, simplify = FALSE)
        # prepare "boot" object for the calculation of the confidence intervals
        contrast_bootstrap <- bootstrap
        contrast_bootstrap$t0 <- get_contrasts(bootstrap$t0, combinations,
                                               type = contrast)
        contrast_bootstrap$t <- get_contrasts(bootstrap$t, combinations,
                                              type = contrast)
        # compute confidence intervals of contrasts
        indices_contrasts <- seq_along(combinations)
        contrast_ci <- lapply(indices_contrasts, function(j) {
          confint(contrast_bootstrap, parm = j, level = level,
                  alternative = alternative, type = type)
        })
        contrast_ci <- do.call(rbind, contrast_ci)
        # combine confidence intervals of indirect effects and contrasts
        ci <- rbind(ci, contrast_ci)
        # add rownames to conficence intervals
        rownames(ci) <- rownames(object$ci)
        # compute estimates of the contrasts
        if (update_contrast) {
          # compute bootstrap estimates
          ab_boot <- c(object$ab[indices_ab],
                       colMeans(contrast_bootstrap$t, na.rm = TRUE))
          # compute estimates on original data
          ab_data <- object$fit$ab[indices_ab]
          contrasts_ab <- get_contrasts(ab_data, combinations, type = contrast)
          ab_data <- c(ab_data, contrasts_ab)
          # add names
          names(ab_boot) <- names(ab_data) <- names(object$ab)
        }
      } else if (update_contrast) {
        ## the new object doesn't have contrasts, but the old object did
        # add rownames to conficence intervals
        rownames(ci) <- rownames(object$ci[indices_ab])
        # removing contrasts from estimates
        ab_boot <- object$ab[indices_ab]
        ab_data <- object$fit$ab[indices_ab]
      } else {
        ## neither the old nor the new object had contrasts
        # add rownames to conficence intervals
        rownames(ci) <- rownames(object$ci)
      }
    }
    # modify object
    object$ci <- ci
    if (update_alternative) object$alternative <- alternative
    if (update_level) object$level <- level
    if (update_type) object$type <- type
    if (update_contrast) {
      object$ab <- ab_boot
      object$fit$ab <- ab_data
      object$fit$contrast <- contrast
    }
  } else warning("no new argument values; returning original object")
  # return modified object
  object
}


#' @rdname retest
#' @method retest sobel_test_mediation
#' @export

retest.sobel_test_mediation <- function(object, alternative, ...) {
  # initializations
  if (missing(alternative)) update <- FALSE
  else {
    alternative <- match.arg(alternative,
                             choices = c("twosided", "less", "greater"))
    update <- alternative != object$alternative
  }
  # reperform test if necessary
  if (update) {
    # recompute confidence interval and modify object
    object$p_value <- p_value_z(object$statistic, alternative = alternative)
    object$alternative <- alternative
  } else warning("no new argument values; returning original object")
  # return modified object
  object
}
