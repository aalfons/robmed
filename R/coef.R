# ------------------------------------
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ------------------------------------

## internal function
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


#' @method coef bootMA
#' @export
coef.bootMA <- coefMA


#' @method coef sobelMA
#' @export
coef.sobelMA <- coefMA
