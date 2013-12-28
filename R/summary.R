# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

## summary of a mediation model fit based on a scatter matrix
#' @S3method summary covFitMediation
summary.covFitMediation <- function(object, ...) {
  # extract covariance matrix
  S <- object$cov$cov
  # extract variable names
  cn <- names(object$data)
  x <- cn[1]
  y <- cn[2]
  m <- cn[3]
  # extract number of observations
  n <- nobs(object$cov)
  # compute the inverse of the Fisher information matrix for the unique 
  # elements of the covariance matrix
  D <- duplicationMatrix(3)
  invS <- solve(S)
  W <- t(D) %*% kronecker(invS, invS) %*% D / 2
  OmegaSigma <- solve(W)
  # apply the delta method to obtain the inverse of the Fisher information 
  # matrix for the coefficients in the mediation model
  d <- S[x,x]*S[m,m] - S[m,x]^2
  dsq <- d^2
  hDot <- matrix(c(-S[m,x]/S[x,x]^2, 0, 1/S[x,x], 0, 0, 0,   # a
                   (S[m,x]^2*S[y,m]+S[m,x]*S[y,x]*S[m,m])/dsq, -S[m,x]/d, 
                   (-S[y,x]*S[x,x]*S[m,m]-S[m,x]^2*S[y,x]+2*S[x,x]*S[y,m]*S[m,x])/dsq,
                   (S[m,x]*S[y,x]*S[x,x]-S[x,x]^2*S[y,m])/dsq, S[x,x]/d, 0,  # b
                   (-S[m,m]^2*S[y,x]+S[m,x]*S[y,m]*S[m,m])/dsq, S[m,m]/d, 
                   (-S[y,m]*S[x,x]*S[m,m]-S[m,x]^2*S[y,m]+2*S[m,m]*S[y,x]*S[m,x])/dsq, 
                   (-S[m,x]^2*S[y,x]+S[m,x]*S[y,m]*S[x,x])/dsq, -S[m,x]/d, 0,  # c
                   -S[y,x]/S[x,x]^2, 1/S[x,x], 0, 0, 0, 0),  # c'
                 6, 4)
  OmegaTheta <- t(hDot) %*% OmegaSigma %*% hDot
  # compute standard errors
  se <- sqrt(diag(OmegaTheta) / n)
  # perform asymptotic tests
  coefficients <- c(object$a, object$b, object$c, object$cPrime)
  z <- coefficients / se
  pValue <- pValueZ(z)
  coefficients <- cbind(coefficients, se, z, pValue)
  cn <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  dimnames(coefficients) <- list(c(x, m, x, x), cn)
  # compute residual standard error
  s <- S[y,y]-object$b^2*S[m,m]-object$c^2*S[x,x]-2*object$b*object$c*S[m,x]
  s <- list(value=s)
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
  # extract residual standard error
  s <- list(value=tmp$sigma, df=tmp$df[2])
  # construct return object
  result <- list(a=a, b=b, c=c, cPrime=cPrime, robust=robust, s=s)
  # add R-squared and F-test for nonrobust fit
  if(!robust) {
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
