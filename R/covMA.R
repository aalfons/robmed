# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

## fit a mediation model based on a scatter matrix
covMA <- function(x, y, m, data, robust = TRUE, control=covHuber.control(...), 
                  ...) {
  # initializations
  robust <- isTRUE(robust)
  # compute scatter matrix (Huber-type M estimator or MLE of covariance matrix)
  cov <- if(robust) covHuber(data, control=control) else covMLE(data)
  S <- cov$cov
  # compute coefficients of mediation model
  a <- S[m, x] / S[x, x]
  det <- S[x, x] * S[m, m] - S[m, x]^2
  b <- (-S[m, x] * S[y, x] + S[x, x] * S[y, m]) / det
  c <- (S[m, m] * S[y, x] - S[m, x] * S[y, m]) / det
  cPrime <- S[y, x] / S[x, x]
  # return results
  result <- list(a=a, b=b, c=c, cPrime=cPrime, robust=robust, cov=cov)
  class(result) <- c("covMA", "fitMA")
  result
}


## compute duplication matrix according to Magnus & Neudecker (1999, p.49)
duplicationMatrix <- function(p){
  D <- diag(p)
  index <- seq(p*(p+1)/2)
  D[lower.tri(D, diag=TRUE)] <- index
  D[upper.tri(D)] <- D[lower.tri(D)]
  outer(c(D), index, function(i, j) ifelse(i == j, 1, 0 ))
}
