# ------------------------------------
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ------------------------------------

#' @method confint bootMA
#' @export

# confint.bootMA <- function(object, parm = NULL, level = NULL, 
#                            alternative = c("twosided", "less", "greater"), 
#                            ...) {
#   # initializations
#   if(!is.null(level)) {
#     level <- rep(as.numeric(level), length.out=1)
#     if(is.na(level) || level < 0 || level > 1) level <- NULL
#   }
#   if(missing(alternative)) alternative <- object$alternative
#   else alternative <- match.arg(alternative)
#   recompute <- !(is.null(level) && alternative == object$alternative)
#   if(is.null(level)) level <- 1 - object$alpha
#   # confidence interval for indirect effect
#   if(recompute) ci <- confint(object$reps, level=level, alternative=alternative)
#   else ci <- object$ci
#   ci <- t(ci)
#   rownames(ci) <- "ab"
#   # confidence interval for other effects
#   ciYX <- confint(object$fitYX, parm=2, level=level)
#   rownames(ciYX) <- "c'"
#   ciMX <- confint(object$fitMX, parm=2, level=level)
#   rownames(ciMX) <- "a"
#   ciYMX <- confint(object$fitYMX, parm=2:3, level=level)
#   rownames(ciYMX) <- c("b", "c")
#   # combine confidence intervals into one matrix
#   ci <- rbind(ciMX, ciYMX["b", , drop=FALSE], ciYX, 
#               ciYMX["c", , drop=FALSE], ci)
#   if(alternative != "twosided") colnames(ci) <- c("Lower", "Upper")
#   # if requested, take subset of effects
#   if(!is.null(parm)) ci <- ci[parm, , drop=FALSE]
#   ci
# }

## argument 'level' is ignored
confint.bootMA <- function(object, parm = NULL, level = NULL, ...) {
  # initializations
  level <- 1 - object$alpha
  # combine confidence interval of indirect effect with those of other effects
  ci <- rbind(confintEffects(object, level=level), ab=object$ci)
  if(object$alternative != "twosided") colnames(ci) <- c("Lower", "Upper")
  # if requested, take subset of effects
  if(!is.null(parm)) ci <- ci[parm, , drop=FALSE]
  ci
}


#' @method confint sobelMA
#' @export

confint.sobelMA <- function(object, parm = NULL, level = 0.95, ...) {
  # initializations
  level <- rep(as.numeric(level), length.out=1)
  if(is.na(level) || level < 0 || level > 1) level <- formals()$level
  # confidence interval of indirect effect
  ci <- confintZ(object$ab, object$se, level=level, 
                 alternative=object$alternative)
  # combine with confidence intervalse of other effects
  ci <- rbind(confintEffects(object, level=level), ab=ci)
  if(object$alternative != "twosided") colnames(ci) <- c("Lower", "Upper")
  # if requested, take subset of effects
  if(!is.null(parm)) ci <- ci[parm, , drop=FALSE]
  ci
}


## internal function to extract confidence interval from bootstrap results
## argument 'parm' is ignored
confint.boot <- function(object, parm, level = 0.95, 
                         alternative = c("twosided", "less", "greater"), 
                         ...) {
  # initializations
  alternative <- match.arg(alternative)
  # extract confidence interval
  if(alternative == "twosided") {
    ci <- boot.ci(object, conf=level, type="perc")$percent[4:5]
  } else {
    alpha <- 1-level
    ci <- boot.ci(object, conf=1-2*alpha, type="perc")$percent[4:5]
    if(alternative == "less") ci[1] <- -Inf
    else ci[2] <- Inf
  }
  # return confidence interval
  ci
}


## internal function for confidence interval based on normal distribution
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


## internal function to extract confidence intervals for effects other than 
## the indirect effect
confintEffects <- function(object, level = 0.95) {
  ciYX <- confint(object$fitYX, parm=2, level=level)
  rownames(ciYX) <- "c'"
  ciMX <- confint(object$fitMX, parm=2, level=level)
  rownames(ciMX) <- "a"
  ciYMX <- confint(object$fitYMX, parm=2:3, level=level)
  rownames(ciYMX) <- c("b", "c")
  # combine confidence intervals into one matrix
  rbind(ciMX, ciYMX["b", , drop=FALSE], ciYX, ciYMX["c", , drop=FALSE])
}
