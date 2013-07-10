# ------------------------------------
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ------------------------------------

#' @export
retest <- function(object, ...) UseMethod("retest")

#' @method retest bootMA
#' @export
retest.bootMA <- function(object, alpha = 0.05, 
                          alternative = c("twosided", "less", "greater"), 
                          ...) {
  # initializations
  alpha <- rep(as.numeric(alpha), length.out=1)
  if(is.na(alpha) || alpha < 0 || alpha > 1) alpha <- formals()$alpha
  alternative <- match.arg(alternative)
  # recompute confidence interval and modify object
  object$ci <- confint(object$reps, level=1-alpha, alternative=alternative)
  object$alpha <- alpha
  object$alternative <- alternative
  # return modified object
  object
}

#' @method retest sobelMA
#' @export
retest.sobelMA <- function(object, 
                           alternative = c("twosided", "less", "greater"), 
                           ...) {
  # initializations
  alternative <- match.arg(alternative)
  # recompute confidence interval and modify object
  object$pValue <- pvalZ(object$statistic, alternative=alternative)
  object$alternative <- alternative
  # return modified object
  object
}
