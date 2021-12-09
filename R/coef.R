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
  x <- object$fit$x
  m <- object$fit$m
  nr_indirect <- length(x) * length(m)
  type <- match.arg(type)
  # extract effects (including indirect effect)
  if(type == "boot") {
    # construct vector of bootstrap estimates
    ab <- object$ab
    coef <- c(object$a, object$b, object$direct, object$total, ab)
    # add coefficient names
    if (nr_indirect == 1L) indirect_names <- "ab"
    else indirect_names <- paste("ab", names(ab), sep = "_")
    names(coef) <- c(get_effect_names(x, m), indirect_names)
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
  # TODO: 'parm' should only allow the values "a", "b", "d", "direct", "total", "indirect"
  #       and then return all effects of that type, as the names of the effects
  #       can become quite complex.  At least this way it is predictable for
  #       the user.  For backwards compatibility, "ab" should be allowed as
  #       a synonym for "indirect".
  # extract effects
  coef <- c(object$a, object$b, object$d, object$direct,
            object$total, object$indirect)
  names(coef) <- get_effect_names(object$a, object$b, object$d, object$direct,
                                  object$total, object$indirect)
  # if requested, take subset of effects
  if (!is.null(parm))  coef <- coef[parm]
  coef
}
# coef.fit_mediation <- function(object, parm = NULL, ...) {
#   # initializations
#   x <- object$x
#   m <- object$m
#   nr_indirect <- length(x) * length(m)
#   # extract effects
#   ab <- object$ab
#   coef <- c(object$a, object$b, object$d, object$direct, object$total, ab)
#   # add coefficient names
#   if (nr_indirect == 1L) indirect_names <- "ab"
#   else indirect_names <- paste("ab", names(ab), sep = "_")
#   names(coef) <- c(get_effect_names(x, m), indirect_names)
#   # if requested, take subset of effects
#   if (!is.null(parm))  coef <- coef[parm]
#   coef
# }


# utility function to get names of coefficients
get_effect_names <- function(a = NULL, b = NULL, d = NULL, direct = NULL,
                             total = NULL, indirect = NULL, sep = "_") {
  # names for effect(s) a
  if (is.null(a)) a_names <- NULL
  else a_names <- if (length(a) == 1L) "a" else paste("a", names(a), sep = sep)
  # names for effect(s) b
  if (is.null(b)) b_names <- NULL
  else b_names <- if (length(b) == 1L) "b" else paste("b", names(b), sep = sep)
  # names for effect(s) d
  if (is.null(d)) d_names <- NULL
  else d_names <- paste("d", names(d), sep = sep)
  # names for direct effect(s)
  if (is.null(direct)) direct_names <- NULL
  else if (length(direct) == 1L) direct_names <- "Direct"
  else direct_names <- paste("Direct", names(direct), sep = sep)
  # names for total effect(s)
  if (is.null(total)) total_names <- NULL
  else if (length(total) == 1L) total_names <- "Total"
  else total_names <- paste("Total", names(total), sep = sep)
  # names for indirect effect(s)
  if (is.null(indirect)) indirect_names <- NULL
  else if (length(indirect) == 1L) indirect_names <- "Indirect"
  else indirect_names <- paste("Indirect", names(indirect), sep = sep)
  # return requested names
  c(a_names, b_names, d_names, direct_names, total_names, indirect_names)
}
# get_effect_names <- function(x, m, sep = "_") {
#   # initializations
#   p_x <- length(x)
#   p_m <- length(m)
#   # construct names
#   if (p_x == 1L) {
#     if (p_m == 1L) c("a", "b", "Direct", "Total")
#     else {
#       c(paste("a", m, sep = sep), paste("b", m, sep = sep), "Direct", "Total")
#     }
#   } else {
#     if (p_m == 1L) {
#       c(paste("a", x, sep = sep), "b", paste("Direct", x, sep = sep),
#         paste("Total", x, sep = sep))
#     } else {
#       c(paste("a", sapply(m, paste, x, sep = "."), sep = sep),
#         paste("b", m, sep = sep),
#         paste("Direct", x, sep = sep),
#         paste("Total", x, sep = sep))
#     }
#   }
# }
