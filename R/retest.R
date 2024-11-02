# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

#' Retest for mediation
#'
#' Re-perform a test for the indirect effect(s) based on results from (robust)
#' mediation analysis.  This function is purely available for computational
#' convenience if the analysis was accidentally run with the wrong parameter
#' settings, as it avoids having to re-run the bootstrap procedure.  It must
#' not be abused for \eqn{p}{p}-hacking.
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
#' to be computed in the bootstrap test.  Possible values are \code{"perc"}
#' for the percentile bootstrap, or \code{"bca"} for the bias-corrected and
#' accelerated (BCa) bootstrap.
#' @param contrast  a logical indicating whether to compute pairwise contrasts
#' of the indirect effects.  This can also be a character string, with
#' \code{"estimates"} for computing the pairwise differences of the indirect
#' effects (such that it is tested whether two indirect effects are equal),
#' and \code{"absolute"} for computing the pairwise differences of the absolute
#' values of the indirect effects (such that it is tested whether two indirect
#' effects are equal in magnitude).  This is only relevant for models with
#' multiple indirect effects, which are currently only implemented for
#' bootstrap tests and estimation via regressions.
#' @param order  a character string specifying the order of approximation of
#' the standard error in Sobel's test.  Possible values are \code{"first"}
#' for a first-order approximation, and \code{"second"} for a second-order
#' approximation.
#' @param \dots  additional arguments to be passed down to methods.
#'
#' @return An object of the same class as \code{object} with updated test
#' results (see \code{\link{test_mediation}()}).
#'
#' @note
#' From version 0.9.0 onwards, the behavior of this function changed. For
#' arguments that are not supplied, the corresponding values of \code{object}
#' are now used as defaults.
#'
#' Since version 1.1.0, bias-corrected and accelerated (BCa) bootstrap
#' confidence intervals are no longer recommended, and the option
#' \code{type = "bca"} may disappear in the future.
#'
#' @author Andreas Alfons
#'
#' @seealso \code{\link{test_mediation}()}
#'
#' @examples
#' data("BSG2014")
#'
#' # run fast-and-robust bootstrap test
#' boot <- test_mediation(BSG2014,
#'                        x = "SharedLeadership",
#'                        y = "TeamPerformance",
#'                        m = c("ProceduralJustice",
#'                              "InteractionalJustice"),
#'                        covariates = c("AgeDiversity",
#'                                       "GenderDiversity"))
#' summary(boot)
#'
#' # now include comparison of indirect effects
#' retest(boot, contrast = "estimates")
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
  defaults <- list(alternative = c("twosided", "less", "greater"),
                   type = c("perc", "bca"),
                   contrast = c("estimates", "absolute"))
  fit <- object$fit
  # check alternative hypothesis
  if (missing(alternative)) alternative <- object$alternative
  else alternative <- match.arg(alternative, choices = defaults$alternative)
  # check confidence level
  if (missing(level)) level <- object$level
  else {
    level <- rep(as.numeric(level), length.out = 1L)
    if(is.na(level) || level < 0 || level > 1) {
      level <- object$level
      warning("confidence level must be between 0 and 1; not updating it")
    }
  }
  # check type of confidence intervals
  if (missing(type)) type <- object$type
  else {
    # check for a valid value
    type <- match.arg(type, choices = defaults$type)
    # check type if BCa confidence intervals can be computed
    if (type == "bca" && object$R < nrow(fit$data)) {
      type <- "perc"
      warning("cannot compute BCa confidence intervals as number of ",
              "bootstrap samples is smaller than number of observations; ",
              "using percentile confidence intervals")
    }
  }
  # check contrasts of indirect effect
  if (inherits(fit, "reg_fit_mediation")) {
    # further initializations
    model <- fit$model
    # regression model fit (multiple mediators and contrasts are supported)
    if (missing(contrast)) contrast <- fit$contrast
    else {
      if (model == "simple") contrast <- FALSE
      else if (is.logical(contrast)) {
        contrast <- isTRUE(contrast)
        if (contrast) contrast <- defaults$contrast[1L]
      } else contrast <- match.arg(contrast, choices = defaults$contrast)
    }
    update_contrast <- contrast != fit$contrast
  } else if (inherits(fit, "cov_fit_mediation")) {
    # covariance model fit (only implemented for a simple mediation model)
    update_contrast <- FALSE
  } else stop("not implemented for this type of model fit")
  # check for any new arguments
  update_alternative <- alternative != object$alternative
  update_level <- level != object$level
  update_type <- type != object$type
  update <- update_alternative || update_level || update_type || update_contrast
  # reperform test if necessary
  if (update) {
    # if contrasts are updated, first recompute the point estimates of the
    # indirect effects before recomputing bootstrap replicates
    if (update_contrast) {
      # extract estimates of the indirect effects
      indirect_data <- extract_effects(fit$x, fit$m, family = fit$family,
                                       model = model, contrast = contrast,
                                       fit_mx = fit$fit_mx,
                                       fit_ymx = fit$fit_ymx,
                                       fit_yx = fit$fit_yx)$indirect
      # update the object for the model fit
      fit$indirect <- fit$ab <- indirect_data
      fit$contrast <- contrast
      # recompute bootstrap estimates of the indirect effects
      boot_indirect <- extract_boot(fit, boot = object$reps)$indirect
      indirect_boot <- colMeans(boot_indirect, na.rm = TRUE)
      # update the object for the bootstrap test
      object$indirect <- object$ab <- indirect_boot
      object$fit <- fit
    } else if (inherits(fit, "reg_fit_mediation")) {
      # extract bootstrap replicates of the indirect effects
      boot_indirect <- extract_boot(fit, boot = object$reps)$indirect
    }
    # recompute confidence intervals of indirect effects with updated arguments
    if (inherits(fit, "reg_fit_mediation")) {
      # compute confidence intervals of indirect effects
      ci <- boot_ci(fit$indirect, boot_indirect, object = object$reps,
                    alternative = alternative, level = level, type = type)
    } else {
      ci <- extract_ci(parm = 5L, object = object$reps,
                       alternative = alternative,
                       level = level, type = type)
    }
    # update the object for the bootstrap test
    object$ci <- ci
    if (update_alternative) object$alternative <- alternative
    if (update_level) object$level <- level
    if (update_type) object$type <- type
  } else warning("no new argument values; returning original object")
  # return modified object
  object
}


#' @rdname retest
#' @method retest sobel_test_mediation
#' @export

retest.sobel_test_mediation <- function(object, alternative, order, ...) {
  # initializations
  defaults <- list(alternative = c("twosided", "less", "greater"),
                   order = c("first", "second"))
  if (missing(alternative)) alternative <- object$alternative
  else alternative <- match.arg(alternative, choices = defaults$alternative)
  if (missing(order)) order <- object$order
  else order <- match.arg(order, choices = defaults$order)
  # reperform test if necessary
  if (order != object$order) {
    # entire test needs to be re-run (standard error and test statistic change)
    object <- sobel_test_mediation(object$fit, alternative = alternative,
                                   order = order)
  } else if (alternative != object$alternative) {
    # only recompute p-value and modify object (test statistic is unchanged)
    object$p_value <- p_value_z(object$statistic, alternative = alternative)
    object$alternative <- alternative
  } else warning("no new argument values; returning original object")
  # return modified object
  object
}
