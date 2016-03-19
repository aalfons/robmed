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
  ci <- rbind(confint(object$fit, level=object$level), ab=object$ci)
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
  ci <- rbind(confint(object$fit, level=level), ab=ci)
  if(object$alternative != "twosided") colnames(ci) <- c("Lower", "Upper")
  # if requested, take subset of effects
  if(!is.null(parm)) ci <- ci[parm, , drop=FALSE]
  ci
}


#' @rdname confint.testMediation
#' @method confint regFitMediation
#' @export

confint.regFitMediation <- function(object, parm = NULL, level = 0.95, ...) {
  # extract confidence intervals and combine into one matrix
  ci <- rbind(confint(object$fitMX, parm=2, level=level), 
              confint(object$fitYMX, parm=2:3, level=level), 
              confint(object$fitYX, parm=2, level=level))
  rownames(ci) <- c("a", "b", "c", "c'")
  # if requested, take subset of effects
  if(!is.null(parm)) ci <- ci[parm, , drop=FALSE]
  ci
}


#' @rdname confint.testMediation
#' @method confint covFitMediation
#' @export

confint.covFitMediation <- function(object, parm = NULL, level = 0.95, ...) {
  # initializations
  alpha <- 1 - level
  # compute standard errors
  summary <- summary(object)
  # compute confidence intervals and combine into one matrix
  ci <- rbind(confintZ(object$a, summary$a[1,2], level=level), 
              confintZ(object$b, summary$b[1,2], level=level), 
              confintZ(object$c, summary$c[1,2], level=level), 
              confintZ(object$cPrime, summary$cPrime[1,2], level=level))
  cn <- paste(format(100 * c(alpha/2, 1-alpha/2), trim=TRUE), "%")
  dimnames(ci) <- list(c("a", "b", "c", "c'"), cn)
  # if requested, take subset of effects
  if(!is.null(parm)) ci <- ci[parm, , drop=FALSE]
  ci
}


## internal functions

# extract confidence interval from bootstrap results 
# (argument 'parm' is ignored)
confint.boot <- function(object, parm, level = 0.95, 
                         alternative = c("twosided", "less", "greater"), 
                         type = c("bca", "perc"), ...) {
  # initializations
  alternative <- match.arg(alternative)
  type <- match.arg(type)
  component <- if(type == "perc") "percent" else type
  # extract confidence interval
  if(alternative == "twosided") {
    ci <- boot.ci(object, conf=level, type=type)[[component]][4:5]
  } else {
    alpha <- 1 - level
    ci <- boot.ci(object, conf=1-2*alpha, type=type)[[component]][4:5]
    if(alternative == "less") ci[1] <- -Inf
    else ci[2] <- Inf
  }
  # return confidence interval
  ci
}

# compute confidence interval based on normal distribution
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