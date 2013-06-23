# ------------------------------------
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ------------------------------------

#' (Robust) mediation analysis
#' 
#' Perform (robust) mediation analysis via a (fast and robust) bootstrap test 
#' or Sobel's test.
#' 
#' @aliases summary.bootMA summary.sobelMA
#' 
#' @param x  the explanatory variable.
#' @param y  the response variable.
#' @param m  the mediator variable.
#' @param method  the type of test to be performed.  Possible values are 
#' \code{"boot"} for the bootstrap or \code{"sobel"} for Sobel's test.
#' @param R  the number of bootstrap replicates.
#' @param alpha  the significance level of the test.
#' @param \dots  additional arguments to be passed to \code{\link[boot]{boot}}.
#' 
#' @return An object containing the necessary information for mediation 
#' analysis.  A \code{summary} method to produce summaries is available.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link[robustbase]{lmrob}}, \code{\link[stats]{lm}}, 
#' \code{\link[boot]{boot}}
#' 
#' @keywords multivariate
#' 
#' @export

mediate <- function(x, y, m, method = c("boot", "sobel"), 
                    alternative = c("twosided", "less", "greater"), 
                    R = 5000, alpha = 0.05, robust = TRUE, 
                    transform = FALSE, control, ...) {
  ## initializations
  # prepare data frame containing all variables with original names
  x <- substitute(x)
  y <- substitute(y)
  m <- substitute(m)
  data <- eval.parent(call("data.frame", x, y, m))
  # make sure that variables are numeric
  convert <- !sapply(data, is.numeric)
  data[convert] <- lapply(data[convert], as.numeric)
  # check if there are enough observations
  n <- nrow(data)
  if(n <= 3) stop("not enough observations")
  # check other arguments
  method <- match.arg(method)
  alternative <- match.arg(alternative)
  robust <- isTRUE(robust)
  transform <- isTRUE(transform)
  if(robust && missing(control)) {
    control <- if(transform) covHuber.control() else lmrob.control()
  }
  ## check if the robust transformation of Zu & Yuan (2010) should be applied
  if(robust && transform) {
    huberEst <- covHuber(data, control=control)
    data[] <- mapply("-", data, huberEst$mu, SIMPLIFY=FALSE, USE.NAMES=FALSE)
    data <- huberEst$weights * data
  }
  ## compute regression coefficients
  fYX <- as.formula(paste(y, "~", x))
  fMX <- as.formula(paste(m, "~", x))
  fYMX <- as.formula(paste(y, "~", m, "+", x))
  if(robust && !transform) {
    fitMX <- lmrob(fMX, data=data, control=control)
    fitYX <- lmrob(fYX, data=data, control=control)
    fitYMX <- lmrob(fYMX, data=data, control=control)
  } else {
    fitYX <- lm(fYX, data=data)
    fitMX <- lm(fMX, data=data)
    fitYMX <- lm(fYMX, data=data)
  }
  ## perform mediation analysis
  if(method == "boot") {
    ## bootstrap test
    if(robust && !transform) {
      # extract (square root of) robustness weights and combine data into matrix
      wM <- sqrt(weights(fitMX, type="robustness"))
      wY <- sqrt(weights(fitYMX, type="robustness"))
      z <- cbind(rep.int(1, n), as.matrix(data), wM, wY)
      # compute matrices for linear corrections
      psiControl <- getPsiControl(fitMX)  # the same for both model fits
      corrM <- correctionMatrix(z[, 1:2], weights=wM, 
                                residuals=residuals(fitMX), 
                                scale=fitMX$scale, 
                                psiControl=psiControl)
      coefM <- coef(fitMX)
      corrY <- correctionMatrix(z[, c(1, 4, 2)], weights=wY, 
                                  residuals=residuals(fitYMX), 
                                  scale=fitYMX$scale, 
                                  psiControl=psiControl)
      coefY <- coef(fitYMX)
      # perform fast and robust bootstrap
      bootstrap <- localBoot(z, function(z, i, corrM, coefM, corrY, coefY) {
        # extract bootstrap sample from the data
        zi <- z[i, , drop=FALSE]
        wMi <- zi[, "wM"]
        wYi <- zi[, "wY"]
        # check whether there are enough observations with nonzero weights
        if(sum(wMi > 0) <= 2 || sum(wYi > 0) <= 3) return(NA)
        # compute coefficients from weighted regression m ~ x
        wxi <- wMi * zi[, 1:2]
        wmi <- wMi * zi[, 4]
        coefMi <- solve(crossprod(wxi)) %*% crossprod(wxi, wmi)
        # compute coefficients from weighted regression y ~ m + x
        wmxi <- wYi * zi[, c(1, 4, 2)]
        wyi <- wYi * zi[, 3]
        coefYi <- solve(crossprod(wmxi)) %*% crossprod(wmxi, wyi)
        # compute corrected coefficients
        coefMi <- drop(coefM + corrM %*% (coefMi - coefM))
        coefYi <- drop(coefY + corrY %*% (coefYi - coefY))
        # compute indirect effect
        unname(coefMi[2]) * unname(coefYi[2])
      }, R=R, corrM=corrM, coefM=coefM, corrY=corrY, coefY=coefY, ...)
      R <- sum(!is.na(bootstrap$t))  # adjust number of replicates for NAs
    } else {
      z <- cbind(rep.int(1, n), as.matrix(data))
      bootstrap <- localBoot(z, function(z, i) {
        # extract bootstrap sample from the data
        zi <- z[i, , drop=FALSE]
        # compute coefficients from regression m ~ x
        xi <- zi[, 1:2]
        mi <- zi[, 4]
        coefMi <- drop(solve(crossprod(xi)) %*% crossprod(xi, mi))
        # compute coefficients from regression y ~ m + x
        mxi <- zi[, c(1, 4, 2)]
        yi <- zi[, 3]
        coefYi <- drop(solve(crossprod(mxi)) %*% crossprod(mxi, yi))
        # compute indirect effect
        unname(coefMi[2]) * unname(coefYi[2])
      }, R=R, ...)
    }
    # extract bootstrap replicates and confidence interval for indirect effect
    reps <- bootstrap$t
    if(alternative == "twosided") {
      ci <- boot.ci(bootstrap, conf=1-alpha, type="perc")$percent[4:5]
    } else {
      ci <- boot.ci(bootstrap, conf=1-2*alpha, type="perc")$percent[4:5]
      if(alternative == "less") ci[1] <- -Inf
      else ci[2] <- Inf
    }
    # construct return object
    result <- list(ab=mean(reps, na.rm=TRUE), ci=ci, reps=reps, R=R, 
                   alpha=alpha, alternative=alternative, robust=robust, 
                   transform=transform, fitYX=fitYX, fitMX=fitMX, 
                   fitYMX=fitYMX, data=data)
    class(result) <- "bootMA"
  } else {
    ## Sobel test
    a <- unname(coef(fitMX)[2])
    b <- unname(coef(fitYMX)[2])
    # compute standard errors
    summaryMX <- summary(fitMX)
    sa <- coef(summaryMX)[2, 2]
    summaryYMX <- summary(fitYMX)
    sb <- coef(summaryYMX)[2, 2]
    # compute test statistic and p-Value
    ab <- a * b
    se <- sqrt(b^2 * sa^2 + a^2 * sb^2)
    z <- ab / se
    pValue <- switch(alternative, twosided=2*pnorm(-abs(z)), less=pnorm(z), 
                     greater=pnorm(z, lower.tail=FALSE))
    # construct return item
    result <- list(ab=ab, se=se, statistic=z, pValue=pValue, 
                   alternative=alternative, robust=robust, transform=transform, 
                   fitYX=fitYX, fitMX=fitMX, fitYMX=fitYMX, data=data)
    class(result) <- "sobelMA"
  }
  ## return test results
  result
}

# wrapper function for boot() that ignores unused arguments, but allows 
# arguments for parallel computing to be passed down
localBoot <- function(..., sim, stype, L, m, ran.gen, mle) boot(...)

# psi function (derivative of the rho function) as used by lmrob()
psi <- function(x, control, derivative = 0) {
  Mpsi(x, cc=control$tuning.psi, psi=control$psi, deriv=derivative)
}

# get control arguments for psi function as used in a given model fit
getPsiControl <- function(object) object$control[c("tuning.psi", "psi")]

# compute matrix for linear correction
correctionMatrix <- function(X, weights, residuals, scale, psiControl) {
  tmp <- psi(residuals / scale, control=psiControl, derivative=1)
  solve(crossprod(X, tmp * X)) %*% crossprod(weights * X)
}
