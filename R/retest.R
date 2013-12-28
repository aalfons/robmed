# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

#' Retest for mediation
#' 
#' Reperform a (fast and robust) bootstrap test or Sobel's test for the 
#' indirect effect based on results from (robust) mediation analysis.
#' 
#' @param object  an object inheriting from class \code{"\link{testMediation}"} 
#' containing results from (robust) mediation analysis.
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
#' results (see \code{\link{testMediation}}).
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link{testMediation}}
#' 
#' @keywords multivariate
#' 
#' @export

retest <- function(object, ...) UseMethod("retest")


#' @rdname retest
#' @method retest bootTestMediation
#' @export

retest.bootTestMediation <- function(object, 
                                     alternative = c("twosided", "less", "greater"), 
                                     level = 0.95, type = c("bca", "perc"), 
                                     ...) {
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


#' @rdname retest
#' @method retest sobelTestMediation
#' @export

retest.sobelTestMediation <- function(object, 
                                      alternative = c("twosided", "less", "greater"), 
                                      ...) {
  # initializations
  alternative <- match.arg(alternative)
  # recompute confidence interval and modify object
  object$pValue <- pValueZ(object$statistic, alternative=alternative)
  object$alternative <- alternative
  # return modified object
  object
}
