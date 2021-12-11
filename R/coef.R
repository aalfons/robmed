# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

#' Coefficients in (robust) mediation analysis
#'
#' Extract coefficients from models computed in (robust) mediation analysis.
#'
#' @method coef test_mediation
#'
#' @param object  an object inheriting from class \code{"\link{test_mediation}"}
#' containing results from (robust) mediation analysis, or an object inheriting
#' from class \code{"\link{fit_mediation}"} containing a (robust) mediation
#' model fit.
#' @param type  a character string specifying whether to extract the means
#' of the bootstrap distribution (\code{"boot"}; the default), or the
#' coefficient estimates based on the original data set (\code{"data"}).
#' @param parm  an integer, character or logical vector specifying the
#' coefficients to be extracted, or \code{NULL} to extract all coefficients.
#' @param \dots  additional arguments are currently ignored.
#'
#' @return A numeric vector containing the requested coefficients.
#'
#' @author Andreas Alfons
#'
#' @seealso \code{\link{test_mediation}()}, \code{\link{fit_mediation}()},
#' \code{\link[=confint.test_mediation]{confint}()}, \code{\link{p_value}()}
#'
#' @examples
#' data("BSG2014")
#'
#' # fit robust mediation model and extract coefficients
#' fit <- fit_mediation(BSG2014,
#'                      x = "ValueDiversity",
#'                      y = "TeamCommitment",
#'                      m = "TaskConflict")
#' coef(fit)
#'
#' # run fast-and-robust bootstrap test and extract coefficients
#' test <- test_mediation(fit)
#' coef(test, type = "data")  # from original sample
#' coef(test, type = "boot")  # means of bootstrap replicates
#'
#' @keywords utilities
#'
#' @export

coef.test_mediation <- function(object, parm = NULL, ...) {
  coef(object$fit, parm = parm, ...)
}

#' @rdname coef.test_mediation
#' @method coef boot_test_mediation
#' @export

coef.boot_test_mediation <- function(object, parm = NULL,
                                     type = c("boot", "data"),
                                     ...) {
  # initializations
  type <- match.arg(type)
  # extract effects
  if(type == "boot") {
    # TODO: 'parm' should only allow the following values:
    #       "a", "b", "d", "total", "direct", "indirect"
    #       Then all effects of such a path should be returned, as the names
    #       of the effects can become quite complex.  At least this way it is
    #       predictable for the user.  For backwards compatibility, "ab" should
    #       be allowed as a synonym for "indirect".
    # extract bootstrap estimates (also works if d path does not exist)
    keep <- c("a", "b", "d", "total", "direct", "indirect")
    coef <- unlist(object[keep], use.names = FALSE)
    names(coef) <- get_effect_names(effects = object[keep])
    # if requested, take subset of effects
    if (!is.null(parm))  coef <- coef[parm]
  } else coef <- coef(object$fit, parm = parm, ...)
  # return effects
  coef
}


#' @rdname coef.test_mediation
#' @method coef fit_mediation
#' @export

coef.fit_mediation <- function(object, parm = NULL, ...) {
  # TODO: 'parm' should only allow the following values:
  #       "a", "b", "d", "total", "direct", "indirect"
  #       Then all effects of such a path should be returned, as the names
  #       of the effects can become quite complex.  At least this way it is
  #       predictable for the user.  For backwards compatibility, "ab" should
  #       be allowed as a synonym for "indirect".
  # extract effect estimates (also works if d path does not exist)
  keep <- c("a", "b", "d", "total", "direct", "indirect")
  coef <- unlist(object[keep], use.names = FALSE)
  names(coef) <- get_effect_names(effects = object[keep])
  # if requested, take subset of effects
  if (!is.null(parm))  coef <- coef[parm]
  coef
}
