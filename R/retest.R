# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

#' Retest for mediation
#'
#' Reperform a (fast and robust) bootstrap test or Sobel's test for the
#' indirect effect(s) based on results from (robust) mediation analysis.
#'
#' @param object  an object inheriting from class
#' \code{"\link{test_mediation}"} containing results from (robust) mediation
#' analysis.
#' @param alternative  a character string specifying the alternative hypothesis
#' in the test for the indirect effect.  Possible values are \code{"twosided"}
#' (the default), \code{"less"} or \code{"greater"}.
#' @param level  numeric; the confidence level of the confidence interval in
#' the bootstrap test.  The default is to compute a 95\% confidence interval.
#' @param type  a character string specifying the type of confidence interval
#' to be computed in the bootstrap test.  Possible values are \code{"bca"} (the
#' default) for the bias-corrected and accelerated bootstrap, or \code{"perc"}
#' for the percentile bootstrap.
#' @param \dots  additional arguments to be passed down to methods.
#'
#' @return An object of the same class as \code{object} with updated test
#' results (see \code{\link{test_mediation}}).
#'
#' @author Andreas Alfons
#'
#' @seealso \code{\link{test_mediation}}
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

retest.boot_test_mediation <- function(object,
                                       alternative = c("twosided", "less", "greater"),
                                       level = 0.95, type = c("bca", "perc"),
                                       ...) {
  # initializations
  alternative <- match.arg(alternative)
  level <- rep(as.numeric(level), length.out = 1)
  if(is.na(level) || level < 0 || level > 1) level <- formals()$level
  type <- match.arg(type)
  # recompute confidence interval
  m <- object$fit$m
  p_m <- length(m)
  if(p_m == 1L) {
    # only one mediator
    ci <- confint(object$reps, parm = 1L, level = level,
                  alternative = alternative, type = type)
  } else {
    # multiple mediators
    ci <- lapply(seq_len(1L + p_m), function(j) {
      confint(object$reps, parm = j, level = level,
              alternative = alternative, type = type)
    })
    ci <- do.call(rbind, ci)
    rownames(ci) <- c("Total", m)
  }
  # modify object with updated confidence interval
  object$ci <- ci
  object$level <- level
  object$alternative <- alternative
  object$type <- type
  # return modified object
  object
}


#' @rdname retest
#' @method retest sobel_test_mediation
#' @export

retest.sobel_test_mediation <- function(object,
                                        alternative = c("twosided", "less", "greater"),
                                        ...) {
  # initializations
  alternative <- match.arg(alternative)
  # recompute confidence interval and modify object
  object$p_value <- p_value_z(object$statistic, alternative=alternative)
  object$alternative <- alternative
  # return modified object
  object
}
