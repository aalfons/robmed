# ------------------------------------
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ------------------------------------

## internal function to extract coefficients
coefMA <- function(object, parm = NULL, ...) {
  # extract effects
  a <- unname(coef(object$fitMX)[-1])
  bc <- unname(coef(object$fitYMX)[-1])
  cPrime <- unname(coef(object$fitYX)[-1])
  coef <- c(a, bc[1], cPrime, bc[2], object$ab)
  names(coef) <- c("a", "b", "c'", "c", "ab")
  # if requested, take subset of effects
  if(!is.null(parm))  coef <- coef[parm]
  # return effects
  coef
}


#' Coefficients in (robust) mediation analysis
#' 
#' Extract coefficients from regression models computed in (robust) mediation 
#' analysis.
#' 
#' @method coef bootMA
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

coef.bootMA <- coefMA


#' @rdname coef.bootMA
#' @method coef sobelMA
#' @export

coef.sobelMA <- coefMA
