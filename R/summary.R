# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

## summary of a mediation model fit based on a scatter matrix
#' @S3method summary covFitMediation
summary.covFitMediation <- function(object, ...) {
  # extract variable names
  cn <- names(object$data)
  x <- cn[1]
  y <- cn[2]
  m <- cn[3]
  # extract covariance matrix
  S <- object$cov$cov[c(x, m, y), c(x, m, y)]
  # extract number of observations
  n <- nobs(object$cov)
  # compute the inverse of the Fisher information matrix for the unique 
  # elements of the covariance matrix
  D <- duplicationMatrix(3)
  invS <- solve(S)
  W <- t(D) %*% kronecker(invS, invS) %*% D / 2
  OmegaSigma <- solve(W)
  # parameters in mediation model
  a <- object$a
  b <- object$b
  c <- object$c
  sEpsilonMX <- S[m,m] - a^2 * S[x,x]
  sEpsilonYMX <- S[y,y] - b^2*S[m,m] - c^2*S[x,x] - 2*b*c*S[m,x]
  # apply the delta method to obtain the inverse of the Fisher information 
  # matrix for the coefficients in the mediation model, see Zu & Yuan (2010)
  hDot <- matrix(c(-a/S[x,x], a*c/S[x,x], -a^2*c/sEpsilonMX-c/S[x,x], 1, a^2, c^2, 
                   1/S[x,x], (a*b-c)/sEpsilonMX, -b/S[x,x]-a*(a*b-c)/sEpsilonMX, 0, -2*a, 2*b*c, 
                   0, -a/sEpsilonMX, a^2/sEpsilonMX+1/S[x,x], 0, 0, -2*c, 
                   0, -b/sEpsilonMX, a*b/sEpsilonMX, 0, 1, b^2, 
                   0, 1/sEpsilonMX, -a/sEpsilonMX, 0, 0, -2*b, 
                   0, 0, 0, 0, 0, 1), 
                 nrow=6, ncol=6)
  OmegaTheta <- hDot %*% OmegaSigma %*% t(hDot)
  # total effect
  cPrime <- object$cPrime
  OmegaSigmaYX <- OmegaSigma[c(1, 3, 6), c(1, 3, 6)]
  hDotYX <- matrix(c(-cPrime/S[x,x], 1, 0, 1/S[x,x], 0, -1, 0, 0, 1), 
                   nrow=3, ncol=3)
  OmegaThetaYX <- hDotYX %*% OmegaSigmaYX %*% t(hDotYX)
  # compute standard errors
  se <- sqrt(c(diag(OmegaTheta)[1:3], OmegaThetaYX[1,1]) / n)
  # perform asymptotic tests
  coefficients <- c(a, b, c, cPrime)
  z <- coefficients / se
  pValue <- pValueZ(z)
  coefficients <- cbind(coefficients, se, z, pValue)
  tn <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  dimnames(coefficients) <- list(c(x, m, x, x), tn)
  # residual standard error as list (for compatibility with regression method)
  s <- list(value=sEpsilonYMX)
  # return results
  result <- list(a=coefficients[1, , drop=FALSE], 
                 b=coefficients[2, , drop=FALSE], 
                 c=coefficients[3, , drop=FALSE], 
                 cPrime=coefficients[4, , drop=FALSE], 
                 robust=object$robust, s=s, n=n, variables=cn)
  class(result) <- "summaryFitMediation"
  result
}

## summary of a mediation model fit based on regression
#' @S3method summary regFitMediation
summary.regFitMediation <- function(object, ...) {
  # initializations
  robust <- object$robust
  # compute summaries of regression models and extract t-tests for coefficients
  tmp <- summary(object$fitYX)
  cPrime <- tmp$coefficients[2, , drop=FALSE]
  tmp <- summary(object$fitMX)
  a <- tmp$coefficients[2, , drop=FALSE]
  tmp <- summary(object$fitYMX)
  b <- tmp$coefficients[2, , drop=FALSE]
  c <- tmp$coefficients[3, , drop=FALSE]
  # construct return object
  result <- list(a=a, b=b, c=c, cPrime=cPrime)
  # add partial effects of control variables if they exist
  p <- nrow(tmp$coefficients)
  if(p > 3) result$covariates <- tmp$coefficients[4:p, , drop=FALSE]
  # add residual standard error
  result$robust <- robust
  result$s <- list(value=tmp$sigma, df=tmp$df[2])
  # add R-squared and F-test for nonrobust fit
  if(robust) {
    # robust F-test not yet implemented, only robust R-squared
    result$FTest <- robR2(object$fitYMX)
  } else {
    statistic <- unname(tmp$fstatistic[1])
    df <- unname(tmp$fstatistic[-1])
    pValue <- pf(statistic, df[1], df[2], lower.tail=FALSE)
    result$FTest <- list(R2=tmp$r.squared, adjR2=tmp$adj.r.squared, 
                         statistic=statistic, df=df, pValue=pValue)
  }
  # add number of observations and variable names
  result$n <- nobs(object$fitYMX)
  result$variables <- names(object$data)
  ## add class and return results
  class(result) <- "summaryFitMediation"
  result
}

## summary of mediation analysis objects
#' @S3method summary testMediation
summary.testMediation <- function(object, ...) {
  result <- list(object=object, summary=summary(object$fit))
  class(result) <- "summaryTestMediation"
  result
}
