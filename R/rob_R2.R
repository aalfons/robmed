# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

## compute robust R-squared from Renaud & Victoria-Feser (2010)
## (this is no longer necessary since summary method for "lmrob" objects now
## computes robust R-squared as well)

# # generic function
# rob_R2 <- function(object, ...) UseMethod("rob_R2")
#
# # method for "lmrob" objects
# rob_R2.lmrob <- function(object, corrected = TRUE, ...) {
#   # initializations
#   corrected <- isTRUE(corrected)
#   psi_control <- get_psi_control(object)
#   # check correction
#   if(corrected) {
#     # compute correction factor for given weight function and tuning
#     # parameters via numerical integration
#     integrand <- function(r, control) {
#       Mwgt(r, cc=control$tuning.psi, psi=control$psi) * dnorm(r)
#     }
#     E1 <- integrate(integrand, -Inf, Inf, control=psi_control)$value
#     integrand <- function(r, control) {
#       r * Mpsi(r, cc=control$tuning.psi, psi=control$psi) * dnorm(r)
#     }
#     E2 <- integrate(integrand, -Inf, Inf, control=psi_control)$value
#     a <- E1 / E2
#   } else a <- 1
#   # extract fitted values and residuals and compute response
#   fitted <- fitted(object)
#   residuals <- residuals(object)
#   y <- fitted + residuals
#   # extract weights
#   w <- weights(object, type="robustness")
#   # compute regression sum of weighted squares
#   SSR <- sum(w * (residuals)^2)
#   # compute robust R-squared
#   if(corrected) {
#     # compute total sum of weighted squares for fitted values
#     mean_fitted <- weighted.mean(fitted, w)
#     SSF <- sum(w * (fitted - mean_fitted)^2)
#     # compute R-squared
#     R2 <- SSF / (SSF + a*SSR)
#   } else {
#     # compute total sum of weighted squares
#     meanY <- weighted.mean(y, w)
#     SST <- sum(w * (y - meanY)^2)
#     # compute robust R-squared
#     R2 <- 1 - SSR / SST
#   }
#   # compute adjusted R-squared
#   n <- length(y)
#   adj_R2 <- 1 - (1 - R2) * (n-1) / object$df.residual  # we always use intercept
#   # return robust R-squared and correction factor
#   list(R2=R2, adj_R2=adj_R2, a=a)
# }
