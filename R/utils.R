# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------


## check whether an object corresponds to a robust model fit

is_robust <- function(object) UseMethod("is_robust")

is_robust.reg_fit_mediation <- function(object) {
  robust <- object$robust
  if (is.character(robust)) TRUE else robust
}

is_robust.cov_fit_mediation <- function(object) object$robust
