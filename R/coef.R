# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

#' Coefficients in (robust) mediation analysis
#' 
#' Extract coefficients from regression models computed in (robust) mediation 
#' analysis.
#' 
#' @method coef testMA
#' 
#' @param object  an object of class \code{"bootMA"} or \code{"sobelMA"} 
#' containing results from (robust) mediation analysis, as returned by 
#' \code{\link{mediate}}.
#' @param parm  an integer, character or logical vector specifying the 
#' coefficients to be extracted, or \code{NULL} to extract all coefficients.
#' @param \dots  additional arguments are currently ignored.
#' 
#' @return A numeric vetor containing the requested coefficients.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link{mediate}}, \code{\link[=confint.bootMA]{confint}}
#' 
#' @keywords utilities
#' 
#' @export

coef.testMA <- function(object, parm = NULL, ...) {
  # extract effects other than indirect effect from mediation model fit
  coef <- coef(object$fit)
  # add coefficient of indirect effect
  coef <- c(coef, ab=object$ab)
  # if requested, take subset of effects
  if(!is.null(parm))  coef <- coef[parm]
  coef
}


## internal function to extract coefficients from mediation model fit
coef.fitMA <- function(object, parm = NULL, ...) {
  # extract effects
  coef <- c(object$a, object$b, object$c, object$cPrime)
  names(coef) <- c("a", "b", "c", "c'")
  # if requested, take subset of effects
  if(!is.null(parm))  coef <- coef[parm]
  coef
}
