# ------------------------------------
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ------------------------------------

#' @export
retest <- function(object, ...) UseMethod("retest")

#' @method retest bootMA
#' @export
retest.bootMA <- function(object, 
                          alternative = c("twosided", "less", "greater"), 
                          level = 0.95, type = c("bca", "perc"), ...) {
  # initializations
  alternative <- match.arg(alternative)
  level <- rep(as.numeric(level), length.out=1)
  if(is.na(level) || level < 0 || level > 1) level <- formals()$level
  type <- match.arg(type)
  # recompute confidence interval and modify object
  object$ci <- confint(object$reps, level=level, alternative=alternative, 
                       type=type)
  object$level <- level
  object$alternative <- alternative
  object$type <- type
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
