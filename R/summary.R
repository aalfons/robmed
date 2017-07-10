# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

# ## summary of a mediation model fit based on a scatter matrix
# ## @export
# summary.covFitMediation <- function(object, ...) {
#   # extract variable names
#   cn <- names(object$data)
#   x <- cn[1]
#   y <- cn[2]
#   m <- cn[3]
#   # extract covariance matrix
#   S <- object$cov$cov[c(x, m, y), c(x, m, y)]
#   # extract number of observations
#   n <- nobs(object$cov)
#   # compute the inverse of the Fisher information matrix for the unique
#   # elements of the covariance matrix
#   D <- duplicationMatrix(3)
#   invS <- solve(S)
#   W <- t(D) %*% kronecker(invS, invS) %*% D / 2
#   OmegaSigma <- solve(W)
#   # parameters in mediation model
#   a <- object$a
#   b <- object$b
#   c <- object$c
#   sEpsilonMX <- S[m,m] - a^2 * S[x,x]
#   sEpsilonYMX <- S[y,y] - b^2*S[m,m] - c^2*S[x,x] - 2*b*c*S[m,x]
#   # apply the delta method to obtain the inverse of the Fisher information
#   # matrix for the coefficients in the mediation model, see Zu & Yuan (2010)
#   hDot <- matrix(c(-a/S[x,x], a*c/S[x,x], -a^2*c/sEpsilonMX-c/S[x,x], 1, a^2, c^2,
#                    1/S[x,x], (a*b-c)/sEpsilonMX, -b/S[x,x]-a*(a*b-c)/sEpsilonMX, 0, -2*a, 2*b*c,
#                    0, -a/sEpsilonMX, a^2/sEpsilonMX+1/S[x,x], 0, 0, -2*c,
#                    0, -b/sEpsilonMX, a*b/sEpsilonMX, 0, 1, b^2,
#                    0, 1/sEpsilonMX, -a/sEpsilonMX, 0, 0, -2*b,
#                    0, 0, 0, 0, 0, 1),
#                  nrow=6, ncol=6)
#   OmegaTheta <- hDot %*% OmegaSigma %*% t(hDot)
#   # total effect
#   cPrime <- object$cPrime
#   OmegaSigmaYX <- OmegaSigma[c(1, 3, 6), c(1, 3, 6)]
#   hDotYX <- matrix(c(-cPrime/S[x,x], 1, 0, 1/S[x,x], 0, -1, 0, 0, 1),
#                    nrow=3, ncol=3)
#   OmegaThetaYX <- hDotYX %*% OmegaSigmaYX %*% t(hDotYX)
#   # compute standard errors
#   se <- sqrt(c(diag(OmegaTheta)[1:3], OmegaThetaYX[1,1]) / n)
#   # perform asymptotic tests
#   coefficients <- c(a, b, c, cPrime)
#   z <- coefficients / se
#   pValue <- pValueZ(z)
#   coefficients <- cbind(coefficients, se, z, pValue)
#   tn <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
#   dimnames(coefficients) <- list(c(x, m, x, x), tn)
#   # residual standard error as list (for compatibility with regression method)
#   s <- list(value=sEpsilonYMX)
#   # return results
#   result <- list(a=coefficients[1, , drop=FALSE],
#                  b=coefficients[2, , drop=FALSE],
#                  c=coefficients[3, , drop=FALSE],
#                  cPrime=coefficients[4, , drop=FALSE],
#                  robust=object$robust, s=s, n=n, variables=cn)
#   class(result) <- "summaryFitMediation"
#   result
# }
#
# ## summary of a mediation model fit based on regression
# ## @export
# summary.regFitMediation <- function(object, ...) {
#   # initializations
#   robust <- object$robust
#   # compute summaries of regression models and extract t-tests for coefficients
#   if(!robust) {
#     tmp <- summary(object$fitYX)
#     cPrime <- tmp$coefficients[2, , drop=FALSE]
#   }
#   tmp <- summary(object$fitMX)
#   a <- tmp$coefficients[2, , drop=FALSE]
#   tmp <- summary(object$fitYMX)
#   b <- tmp$coefficients[2, , drop=FALSE]
#   c <- tmp$coefficients[3, , drop=FALSE]
#   if(robust) {
#     cPrime <- matrix(c(object$cPrime, rep.int(NA, ncol(c)-1)), nrow=1)
#     dimnames(cPrime) <- dimnames(c)
#   }
#   # construct return object
#   result <- list(a=a, b=b, c=c, cPrime=cPrime)
#   # add partial effects of control variables if they exist
#   p <- nrow(tmp$coefficients)
#   if(p > 3) result$covariates <- tmp$coefficients[4:p, , drop=FALSE]
#   # add residual standard error
#   result$robust <- robust
#   result$s <- list(value=tmp$sigma, df=tmp$df[2])
#   # add R-squared and F-test for nonrobust fit
#   if(robust) {
#     # robust F-test not yet implemented, only robust R-squared
#     result$FTest <- robR2(object$fitYMX)
#   } else {
#     statistic <- unname(tmp$fstatistic[1])
#     df <- unname(tmp$fstatistic[-1])
#     pValue <- pf(statistic, df[1], df[2], lower.tail=FALSE)
#     result$FTest <- list(R2=tmp$r.squared, adjR2=tmp$adj.r.squared,
#                          statistic=statistic, df=df, pValue=pValue)
#   }
#   # add number of observations and variable names
#   result$n <- nobs(object$fitYMX)
#   result$variables <- names(object$data)
#   ## add class and return results
#   class(result) <- "summaryFitMediation"
#   result
# }

## summary of a mediation model fit based on a scatter matrix
## (currently doesn't do anything)
#' @export
summary.covFitMediation <- function(object, ...) object

## summary of a mediation model fit based on regression
## (currently doesn't do anything)
#' @export
summary.regFitMediation <- function(object, ...) object

## summary of mediation analysis objects

# summary.testMediation <- function(object, ...) {
#   result <- list(object=object, summary=getSummary(object$fit))
#   if(object$fit$robust && inherits(object$fit, "regFitMediation") &&
#      inherits(object, "bootTestMediation")) {
#     # correct standard deviation of bootstrap estimates for estimation of
#     # parameters a, b and c, and perform two-sided t-test
#     cPrime <- result$summary$cPrime
#     n <- result$summary$n
#     cPrime[1, 2] <- sd(object$reps$t[, 2]) * (n-1)/(n-4)
#     cPrime[1, 3] <- cPrime[1, 1] / cPrime[1, 2]
#     cPrime[1, 4] <- 2*pt(abs(cPrime[1, 3]), df=n-4, lower.tail=FALSE)
#     result$summary$cPrime <- cPrime
#   }
#   class(result) <- "summaryTestMediation"
#   result
# }

#' Summary of results from (robust) mediation analysis
#'
#' Summarize results from (robust) mediation analysis for proper interpretation.
#'
#' @name summary.testMediation
#'
#' @param object  an object inheriting from class \code{"\link{testMediation}"}
#' containing results from (robust) mediation analysis.
#' @param other  a character string specifying how to summarize the effects
#' other than the indirect effect.  Possible values are \code{"boot"} (the
#' default) to compute significance tests using the normal approximation of the
#' bootstrap distribution (i.e., to assume a normal distribution of the
#' corresponding effect with the standard deviation computed from the bootstrap
#' replicates), or \code{"theory"} to compute significance tests via
#' statistical theory (e.g., t-tests if the coefficients are estimated via
#' regression).  Note that this is only relevant for mediation analysis via a
#' bootstrap test, where significance of the indirect effect is always assessed
#' via a percentile-based confidence interval due to the asymmetry of its
#' distribution.
#' @param \dots  additional arguments are currently ignored.
#'
#' @return An object of class \code{"summaryTestMediation"} with the following
#' components:
#' \item{object}{the \code{object} passed to the \code{summary} method, which
#' contains the results from testing the indirect effect.}
#' \item{summary}{an object containing all necessary information to summarize
#' the effects other than the indirect effect.}
#'
#' @author Andreas Alfons
#'
#' @seealso \code{\link{testMediation}}
#'
#' @keywords utilities

NULL


#' @rdname summary.testMediation
#' @method summary bootTestMediation
#' @export

summary.bootTestMediation <- function(object, other = c("boot", "theory"),
                                      ...) {
  # get significance of effects and summary of model fit
  # component 'boot' only exists for bootstrap test, otherwise NULL
  other <- match.arg(other)
  if(other == "boot") summary <- getSummary(object$fit, boot=object$reps)
  else summary <- getSummary(object$fit)
  # construct return object
  result <- list(object=object, summary=summary)
  class(result) <- "summaryTestMediation"
  result
}


#' @rdname summary.testMediation
#' @method summary sobelTestMediation
#' @export

summary.sobelTestMediation <- function(object, ...) {
  # get significance of effects and summary of model fit
  # component 'se' only exists for Sobel test, otherwise NULL
  summary <- getSummary(object$fit)
  # construct return object
  result <- list(object=object, summary=summary)
  class(result) <- "summaryTestMediation"
  result
}


## internal function to compute significance of effects

getSummary <- function(object, ...) UseMethod("getSummary")

getSummary.covFitMediation <- function(object, boot = NULL, ...) {
  # extract variable names
  x <- object$x
  y <- object$y
  m <- object$m
  # extract coefficients
  a <- object$a
  b <- object$b
  c <- object$c
  cPrime <- object$cPrime
  coefficients <- c(a, b, c, cPrime)
  # extract covariance matrix
  S <- object$cov$cov[c(x, m, y), c(x, m, y)]
  # extract number of observations
  n <- nobs(object$cov)
  # compute standard errors and z-statistics
  if(is.null(boot)) {
    # compute the inverse of the Fisher information matrix for the unique
    # elements of the covariance matrix
    D <- duplicationMatrix(3)
    invS <- solve(S)
    W <- t(D) %*% kronecker(invS, invS) %*% D / 2
    OmegaSigma <- solve(W)
    # parameters in mediation model
    # apply the delta method to obtain the inverse of the Fisher information
    # matrix for the coefficients in the mediation model, see Zu & Yuan (2010)
    sEpsilonMX <- S[m,m] - a^2 * S[x,x]
    hDot <- matrix(c(-a/S[x,x], a*c/S[x,x], -a^2*c/sEpsilonMX-c/S[x,x], 1, a^2, c^2,
                     1/S[x,x], (a*b-c)/sEpsilonMX, -b/S[x,x]-a*(a*b-c)/sEpsilonMX, 0, -2*a, 2*b*c,
                     0, -a/sEpsilonMX, a^2/sEpsilonMX+1/S[x,x], 0, 0, -2*c,
                     0, -b/sEpsilonMX, a*b/sEpsilonMX, 0, 1, b^2,
                     0, 1/sEpsilonMX, -a/sEpsilonMX, 0, 0, -2*b,
                     0, 0, 0, 0, 0, 1),
                   nrow=6, ncol=6)
    OmegaTheta <- hDot %*% OmegaSigma %*% t(hDot)
    # total effect
    OmegaSigmaYX <- OmegaSigma[c(1, 3, 6), c(1, 3, 6)]
    hDotYX <- matrix(c(-cPrime/S[x,x], 1, 0, 1/S[x,x], 0, -1, 0, 0, 1),
                     nrow=3, ncol=3)
    OmegaThetaYX <- hDotYX %*% OmegaSigmaYX %*% t(hDotYX)
    # compute standard errors and z-statistics
    means <- NULL
    se <- sqrt(c(diag(OmegaTheta)[1:3], OmegaThetaYX[1,1]) / n)
    z <- coefficients / se
    tn <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  } else {
    # compute means, standard errors and z-statistics from bootstrap replicates
    means <- colMeans(boot$t[, -1], na.rm=TRUE)
    se <- apply(boot$t[, -1], 2, sd, na.rm=TRUE)
    z <- means / se
    tn <- c("Data", "Boot", "Std. Error", "z value", "Pr(>|z|)")
  }
  # perform z-tests and combine results
  pValue <- pValueZ(z)
  coefficients <- cbind(coefficients, means, se, z, pValue)
  dimnames(coefficients) <- list(c(x, m, x, x), tn)
  # residual standard error as list (for compatibility with regression method)
  sEpsilonYMX <- S[y,y] - b^2*S[m,m] - c^2*S[x,x] - 2*b*c*S[m,x]
  s <- list(value=sEpsilonYMX)
  # return results
  result <- list(a=coefficients[1, , drop=FALSE],
                 b=coefficients[2, , drop=FALSE],
                 c=coefficients[3, , drop=FALSE],
                 cPrime=coefficients[4, , drop=FALSE],
                 robust=object$robust, s=s, n=n,
                 variables=c(x, y, m))
  class(result) <- "summaryFitMediation"
  result
}

getSummary.regFitMediation <- function(object, boot = NULL, ...) {
  # initializations
  x <- object$x
  y <- object$y
  m <- object$m
  covariates <- object$covariates
  robust <- object$robust
  haveBoot <- !is.null(boot)
  # compute summary of y ~ m + x + covariates
  summaryYMX <- summary(object$fitYMX)
  # extract number of observations
  n <- nobs(object$fitYMX)
  # perform tests for significance of effects
  if(haveBoot) {
    # extract coefficients
    coefficients <- c(coefficients(object), coef(object$fitYMX)[-(1:3)])
    # compute standard errors and z-statistics from bootstrap replicates
    means <- colMeans(boot$t[, -1], na.rm=TRUE)
    se <- apply(boot$t[, -1], 2, sd, na.rm=TRUE)
    z <- means / se
    # perform z-tests and combine results
    pValue <- pValueZ(z)
    coefficients <- cbind(coefficients, means, se, z, pValue)
    tn <- c("Data", "Boot", "Std. Error", "z value", "Pr(>|z|)")
    dimnames(coefficients) <- list(c(x, m, x, x, covariates), tn)
    # split up effect summaries
    a <- coefficients[1, , drop=FALSE]
    b <- coefficients[2, , drop=FALSE]
    c <- coefficients[3, , drop=FALSE]
    cPrime <- coefficients[4, , drop=FALSE]
  } else {
    # compute summaries of regression models and extract t-tests for coefficients
    tmp <- summary(object$fitMX)
    a <- tmp$coefficients[2, , drop=FALSE]
    b <- summaryYMX$coefficients[2, , drop=FALSE]
    c <- summaryYMX$coefficients[3, , drop=FALSE]
    if(robust) {
      # standard errors and t-test not available
      cPrime <- matrix(c(object$cPrime, rep.int(NA_real_, 3)), nrow=1)
      dimnames(cPrime) <- dimnames(c)
    } else {
      tmp <- summary(object$fitYX)
      cPrime <- tmp$coefficients[2, , drop=FALSE]
    }
  }
  # initialize return object
  result <- list(a=a, b=b, c=c, cPrime=cPrime)
  # add partial effects of control variables if they exist
  if(length(covariates) > 0) {
    if(haveBoot) result$covariates <- coefficients[-(1:4), , drop=FALSE]
    else result$covariates <- summaryYMX$coefficients[-(1:3), , drop=FALSE]
  }
  # add residual standard error
  result$robust <- robust
  result$s <- list(value=summaryYMX$sigma, df=summaryYMX$df[2])
  # add R-squared and F-test for nonrobust fit
  if(robust) {
    # compute robust R-squared
    result$FTest <- robR2(object$fitYMX)
    # TODO: implement robust F-test
  } else {
    # add R-squared and F-test for nonrobust fit
    statistic <- unname(summaryYMX$fstatistic[1])
    df <- unname(summaryYMX$fstatistic[-1])
    pValue <- pf(statistic, df[1], df[2], lower.tail=FALSE)
    result$FTest <- list(R2=summaryYMX$r.squared,
                         adjR2=summaryYMX$adj.r.squared,
                         statistic=statistic, df=df,
                         pValue=pValue)
  }
  # add number of observations and variable names
  result$n <- n
  result$variables <- names(object$data)
  ## add class and return results
  class(result) <- "summaryFitMediation"
  result
}


## compute duplication matrix according to Magnus & Neudecker (1999, p.49)
## (required for computing the Fischer information matrix of a mediation model
## fit based on a scatter matrix)
duplicationMatrix <- function(p){
  D <- diag(p)
  index <- seq(p*(p+1)/2)
  D[lower.tri(D, diag=TRUE)] <- index
  D[upper.tri(D)] <- D[lower.tri(D)]
  outer(c(D), index, function(i, j) ifelse(i == j, 1, 0 ))
}
