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
#' coefficient estimates based on the full data set (\code{"data"}).
#' @param parm  an integer, character or logical vector specifying the
#' coefficients to be extracted, or \code{NULL} to extract all coefficients.
#' @param \dots  additional arguments are currently ignored.
#'
#' @return A numeric vector containing the requested coefficients.
#'
#' @author Andreas Alfons
#'
#' @seealso \code{\link{test_mediation}}, \code{\link{fit_mediation}},
#' \code{\link[=confint.test_mediation]{confint}}
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
#' # run fast and robust bootstrap test and extract coefficients
#' test <- test_mediation(fit)
#' coef(test, type = "data")  # from original sample
#' coef(test, type = "boot")  # means of bootstrap replicates
#'
#' @keywords utilities
#'
#' @export

coef.test_mediation <- function(object, parm = NULL, ...) {
  # extract effects (including indirect effect)
  ab <- object$ab
  names(ab) <- if(length(ab) == 1L) "ab" else paste("ab", names(ab), sep = "_")
  coef <- c(coef(object$fit), ab)
  # if requested, take subset of effects
  if(!is.null(parm))  coef <- coef[parm]
  coef
}

#' @rdname coef.test_mediation
#' @method coef boot_test_mediation
#' @export

coef.boot_test_mediation <- function(object, parm = NULL,
                                     type = c("boot", "data"),
                                     ...) {
  # initializations
  type <- match.arg(type)
  # extract effects (including indirect effect)
  if(type == "boot") {
    # number of mediators and covariates
    p_m <- length(object$fit$m)
    p_covariates <- length(object$fit$covariates)
    # get indices of columns of bootstrap replicates that that correspond to
    # the respective models
    index_list <- get_index_list(p_m, p_covariates)
    # the a path is the second coefficient in the model m ~ x + covariates
    if (p_m == 1L) keep_mx <- index_list$fit_mx[2L]
    else keep_mx <- sapply(index_list$fit_mx, "[", 2L)
    # keep b and c coefficients of model y ~ m + x + covariates
    keep_ymx <- index_list$fit_ymx[1L + seq_len(p_m + 1)]
    # index of c' is stored separately in this list
    keep <- c(keep_mx, keep_ymx, index_list$c_prime)
    # compute means of bootstrap replicates
    coef <- colMeans(object$reps$t[, keep], na.rm = TRUE)
    names(coef) <- get_effect_names(object$fit$m)
    ab <- object$ab
  } else {
    coef <- coef(object$fit)
    ab <- object$fit$a * object$fit$b
    if(length(ab) > 1L) ab <- c(Total = sum(ab), ab)
  }
  # set names and combine coefficients with indirect effect
  names(ab) <- if(length(ab) == 1L) "ab" else paste("ab", names(ab), sep = "_")
  coef <- c(coef, ab)
  # if requested, take subset of effects
  if(!is.null(parm))  coef <- coef[parm]
  coef
}


#' @rdname coef.test_mediation
#' @method coef fit_mediation
#' @export

coef.fit_mediation <- function(object, parm = NULL, ...) {
  # extract effects
  coef <- c(object$a, object$b, object$c, object$c_prime)
  names(coef) <- get_effect_names(object$m)
  # if requested, take subset of effects
  if(!is.null(parm))  coef <- coef[parm]
  coef
}


# utility function to get names of coefficients
get_effect_names <- function(m, sep = "_") {
  if(length(m) == 1L) c("a", "b", "c", "c'")
  else c(paste("a", m, sep = sep), paste("b", m, sep = sep), "c", "c'")
}
