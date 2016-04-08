# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

#' Confidence intervals for (robust) mediation analysis
#'
#' Extract or compute confidence intervals for coefficients from (robust)
#' mediation analysis.
#'
#' @name confint.testMediation
#'
#' @param object  an object inheriting from class \code{"\link{testMediation}"}
#' containing results from (robust) mediation analysis, or an object inheriting
#' from class \code{"\link{fitMediation}"} containing a (robust) mediation
#' model fit.
#' @param parm  an integer, character or logical vector specifying the
#' coefficients for which to extract or compute confidence intervals, or
#' \code{NULL} to extract or compute confidence intervals for all coefficients.
#' @param level  for the \code{"bootTestMediation"} method, this is ignored and
#' the confidence level of the bootstrap confidence interval for the indirect
#' effect is used.  For the other methods, the confidence level of the
#' confidence intervals to be computed.  The default is to compute 95\%
#' confidence intervals.
#' @param \dots  additional arguments are currently ignored.
#'
#' @return A numeric matrix containing the requested confidence intervals.
#'
#' @author Andreas Alfons
#'
#' @seealso \code{\link{testMediation}}, \code{\link{fitMediation}},
#' \code{\link[=coef.testMediation]{coef}}
#'
#' @keywords utilities

NULL


#' @rdname confint.testMediation
#' @method confint bootTestMediation
#' @export

## argument 'level' is ignored
confint.bootTestMediation <- function(object, parm = NULL, level = NULL, ...) {
  # combine confidence interval of indirect effect with those of other effects
  ci <- rbind(getConfint(object$fit, level=object$level, boot=object$reps),
              ab=object$ci)
  if(object$alternative != "twosided") colnames(ci) <- c("Lower", "Upper")
  # if requested, take subset of effects
  if(!is.null(parm)) ci <- ci[parm, , drop=FALSE]
  ci
}


#' @rdname confint.testMediation
#' @method confint sobelTestMediation
#' @export

confint.sobelTestMediation <- function(object, parm = NULL, level = 0.95, ...) {
  # initializations
  level <- rep(as.numeric(level), length.out=1)
  if(is.na(level) || level < 0 || level > 1) level <- formals()$level
  # confidence interval of indirect effect
  ci <- confintZ(object$ab, object$se, level=level,
                 alternative=object$alternative)
  # combine with confidence intervalse of other effects
  ci <- rbind(getConfint(object$fit, level=level, sAB=object$se), ab=ci)
  if(object$alternative != "twosided") colnames(ci) <- c("Lower", "Upper")
  # if requested, take subset of effects
  if(!is.null(parm)) ci <- ci[parm, , drop=FALSE]
  ci
}


# ## @rdname confint.testMediation
# ## @method confint regFitMediation
# ## @export
#
# confint.regFitMediation <- function(object, parm = NULL, level = 0.95, ...) {
#   # extract confidence intervals and combine into one matrix
#   ci <- rbind(confint(object$fitMX, parm=2, level=level),
#               confint(object$fitYMX, parm=2:3, level=level),
#               confint(object$fitYX, parm=2, level=level))
#   rownames(ci) <- c("a", "b", "c", "c'")
#   # if requested, take subset of effects
#   if(!is.null(parm)) ci <- ci[parm, , drop=FALSE]
#   ci
# }
#
#
# ## @rdname confint.testMediation
# ## @method confint covFitMediation
# ## @export
#
# confint.covFitMediation <- function(object, parm = NULL, level = 0.95, ...) {
#   # initializations
#   alpha <- 1 - level
#   # compute standard errors
#   summary <- summary(object)
#   # compute confidence intervals and combine into one matrix
#   ci <- rbind(confintZ(object$a, summary$a[1,2], level=level),
#               confintZ(object$b, summary$b[1,2], level=level),
#               confintZ(object$c, summary$c[1,2], level=level),
#               confintZ(object$cPrime, summary$cPrime[1,2], level=level))
#   cn <- paste(format(100 * c(alpha/2, 1-alpha/2), trim=TRUE), "%")
#   dimnames(ci) <- list(c("a", "b", "c", "c'"), cn)
#   # if requested, take subset of effects
#   if(!is.null(parm)) ci <- ci[parm, , drop=FALSE]
#   ci
# }


## internal function to compute confidence intervals for estimated effects

getConfint <- function(object, parm, level = 0.95, ...) UseMethod("getConfint")

getConfint.covFitMediation <- function(object, parm = NULL, level = 0.95,
                                       boot = NULL, ...) {
  # initializations
  alpha <- 1 - level
  # extract point estimates and standard erroers
  if(is.null(boot)) {
    # combine point estimates
    estimates <- c(object$a, object$b, object$c, object$cPrime)
    # compute standard errors
    summary <- getSummary(object)
    se <- c(summary$a[1,2], summary$b[1,2], summary$c[1,2], summary$cPrime[1,2])
  } else {
    # compute means and standard erroers from bootstrap replicates
    estimates <- colMeans(boot$t[, -1], na.rm=TRUE)
    se <- apply(boot$t[, -1], 2, sd, na.rm=TRUE)
  }
  # compute confidence intervals and combine into one matrix
  ci <- rbind(confintZ(estimates[1], se[1], level=level),
              confintZ(estimates[2], se[2], level=level),
              confintZ(estimates[3], se[3], level=level),
              confintZ(estimates[4], se[4], level=level))
  cn <- paste(format(100 * c(alpha/2, 1-alpha/2), trim=TRUE), "%")
  dimnames(ci) <- list(c("a", "b", "c", "c'"), cn)
  # if requested, take subset of effects
  if(!is.null(parm)) ci <- ci[parm, , drop=FALSE]
  ci
}

getConfint.regFitMediation <- function(object, parm = NULL, level = 0.95,
                                       boot = NULL, sAB = NULL, ...) {
  # initializations
  alpha <- 1 - level
  # extract point estimates and standard erroers
  if(is.null(boot)) {
    # extract confidence intervals from regression models
    confintMX <- confint(object$fitMX, parm=2, level=level)
    confintYMX <- confint(object$fitYMX, parm=2:3, level=level)
    # compute confidence interval for total effect
    if(object$robust) {
      # compute the variance estimate of c' = a*b + c assuming independence
      summaryYMX <- summary(object$fitYMX)
      sCPrime <- sqrt(sAB^2 + summaryYMX$coefficients[3, 2]^2)
      # compute the degrees of freedom
      # m ~ x + covariates: 2 + # covariates
      # y ~ m + x + covariates: 3 + # covariates
      # y ~ x + covariates: 2 + # covariates
      n <- nobs(object$fitYMX)
      df <- max(1, n - 7 - 3*length(object$covariates))
      # compute confidence interval
      confintYX <- object$cPrime + qt(c(alpha/2, 1-alpha/2), df=df) * sCPrime
    } else {
      # extract confidence interval from regression model
      confintYX <- confint(object$fitYX, parm=2, level=level)
    }
    ci <- rbind(confintMX, confintYMX, confintYX)
    rownames(ci) <- c("a", "b", "c", "c'")
  } else {
    # compute means and standard erroers from bootstrap replicates
    estimates <- colMeans(boot$t[, 2:5], na.rm=TRUE)
    se <- apply(boot$t[, 2:5], 2, sd, na.rm=TRUE)
    # compute confidence intervals and combine into one matrix
    ci <- rbind(confintZ(estimates[1], se[1], level=level),
                confintZ(estimates[2], se[2], level=level),
                confintZ(estimates[3], se[3], level=level),
                confintZ(estimates[4], se[4], level=level))
    # add row and column names
    cn <- paste(format(100 * c(alpha/2, 1-alpha/2), trim=TRUE), "%")
    dimnames(ci) <- list(c("a", "b", "c", "c'"), cn)
  }
  # if requested, take subset of effects
  if(!is.null(parm)) ci <- ci[parm, , drop=FALSE]
  ci
}


# extract confidence interval from bootstrap results
# (argument 'parm' can be used for completeness; we only need the confidence
# interval for the indirect effect in the first column of the bootstrap results)
confint.boot <- function(object, parm = 1, level = 0.95,
                         alternative = c("twosided", "less", "greater"),
                         type = c("bca", "perc"), ...) {
  # initializations
  alternative <- match.arg(alternative)
  type <- match.arg(type)
  component <- if(type == "perc") "percent" else type
  # extract confidence interval
  if(alternative == "twosided") {
    ci <- boot.ci(object, conf=level, type=type, index=parm)[[component]][4:5]
  } else {
    alpha <- 1 - level
    ci <- boot.ci(object, conf=1-2*alpha, type=type, index=parm)[[component]][4:5]
    if(alternative == "less") ci[1] <- -Inf
    else ci[2] <- Inf
  }
  # return confidence interval
  ci
}

## internal function to compute confidence interval based on normal distribution
confintZ <- function(mean = 0, sd = 1, level = 0.95,
                     alternative = c("twosided", "less", "greater")) {
  # initializations
  alternative <- match.arg(alternative)
  # compute confidence interval
  alpha <- 1 - level
  switch(alternative, twosided=qnorm(c(alpha/2, 1-alpha/2), mean=mean, sd=sd),
         less=c(-Inf, qnorm(level, mean=mean, sd=sd)),
         greater=c(qnorm(alpha, mean=mean, sd=sd), Inf))
}
