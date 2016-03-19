# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

#' Coefficients in (robust) mediation analysis
#' 
#' Extract coefficients from models computed in (robust) mediation analysis.
#' 
#' @method coef testMediation
#' 
#' @param object  an object inheriting from class \code{"\link{testMediation}"} 
#' containing results from (robust) mediation analysis, or an object inheriting 
#' from class \code{"\link{fitMediation}"} containing a (robust) mediation 
#' model fit.
#' @param parm  an integer, character or logical vector specifying the 
#' coefficients to be extracted, or \code{NULL} to extract all coefficients.
#' @param \dots  additional arguments are currently ignored.
#' 
#' @return A numeric vetor containing the requested coefficients.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link{testMediation}}, \code{\link{fitMediation}}, 
#' \code{\link[=confint.testMediation]{confint}}
#' 
#' @keywords utilities
#' 
#' @export

coef.testMediation <- function(object, parm = NULL, ...) {
  # extract effects (including indirect effect)
  coef <- c(coef(object$fit), ab=object$ab)
  # if requested, take subset of effects
  if(!is.null(parm))  coef <- coef[parm]
  coef
}


#' @rdname coef.testMediation
#' @method coef fitMediation
#' @export

coef.fitMediation <- function(object, parm = NULL, ...) {
  # extract effects
  coef <- c(object$a, object$b, object$c, object$cPrime)
  names(coef) <- c("a", "b", "c", "c'")
  # if requested, take subset of effects
  if(!is.null(parm))  coef <- coef[parm]
  coef
}