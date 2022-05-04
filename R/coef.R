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
#' @param parm  an integer, character or logical vector specifying the
#' paths for which to extract coefficients, or \code{NULL} to extract all
#' coefficients.  In case of a character vector, possible values are
#' \code{"a"}, \code{"b"}, \code{"d"} (only serial multiple mediator
#' models), \code{"total"}, \code{"direct"}, and \code{"indirect"}.
#' @param type  a character string specifying whether to extract the means
#' of the bootstrap distribution (\code{"boot"}; the default), or the
#' coefficient estimates based on the original data set (\code{"data"}).
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
#' boot <- test_mediation(fit)
#' coef(boot, type = "data")  # from original sample
#' coef(boot, type = "boot")  # means of bootstrap replicates
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
  fit <- object$fit
  # extract effects
  if(type == "boot") {
    # call workhorse function with list of effect estimates
    have_d <- inherits(fit, "reg_fit_mediation") && fit$model == "serial"
    keep <- c("a", "b", if (have_d) "d", "total", "direct", "indirect")
    coef_mediation(object[keep], parm = parm)
  } else coef(fit, parm = parm, ...)
}


#' @rdname coef.test_mediation
#' @method coef fit_mediation
#' @export

coef.fit_mediation <- function(object, parm = NULL, ...) {
  # call workhorse function with list of effect estimates
  have_d <- inherits(object, "reg_fit_mediation") && object$model == "serial"
  keep <- c("a", "b", if (have_d) "d", "total", "direct", "indirect")
  coef_mediation(object[keep], parm = parm)
}


# workhorse function to extract coefficients
coef_mediation <- function(coef_list, parm = NULL) {
  # if requested, take subset of effects
  if (!is.null(parm)) coef_list <- coef_list[check_parm(parm)]
  # convert effect estimates to vector and add names
  coef <- unlist(coef_list, use.names = FALSE)
  names(coef) <- get_effect_names(effects = coef_list)
  coef
}

# internal function to check argument 'parm' for backwards compatibility
check_parm <- function(parm = NULL) {
  # some checks if effects are selected via a character vector
  if (is.character(parm)) {
    # for backwards compatibility, check if 'ab' is used for the indirect effect
    which_ab <- which(parm == "ab")
    if (length(which_ab) > 0) {
      parm[which_ab] <- "indirect"
      warning("component 'ab' is deprecated, please use 'indirect' instead",
              call. = FALSE)
    }
    # allow for capitalized names
    parm[which(parm == "Total")] <- "total"
    parm[which(parm == "Direct")] <- "direct"
    parm[which(parm == "Indirect")] <- "indirect"
  }
  # return checked object
  parm
}
