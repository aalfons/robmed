# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------


# regression with skew-elliptical errors
#' @importFrom sn selm.fit
## @importFrom methods new
selm_fit <- function(x, y, intercept = TRUE, family = "ST",
                     fixed.param = list()) {
  # # selm.fit() requires matrix of predictors, not data frame
  # x <- as.matrix(x)
  # selm.fit() always requires constant for intercept, so argument is ignored
  n <- nrow(x)
  x <- cbind("(Intercept)" = rep.int(1, n), x)
  # fit the linear model
  control <- list(method = "MLE")
  fit <- selm.fit(x, y, family = family, fixed.param = fixed.param,
                  selm.control = control)
  # # return model fit as object of S4 class "selm"
  # new("selm", call = call("selm"), family = family, logL = fit$logL,
  #     method = control$method, param = fit$param, param.var = fit$param.var,
  #     size = fit$size, residuals.dp = fit$resid.dp,
  #     fitted.values.dp = fit$fitted.dp, control = list(), input = list(),
  #     opt.method = fit$opt.method)
  # return model fit as an S3 clone of S4 class "selm"
  out <- list(call = NULL, family = family, logL = fit$logL,
              method = control$method, param = fit$param,
              param.var = fit$param.var, size = fit$size,
              residuals.dp = fit$resid.dp, fitted.values.dp = fit$fitted.dp,
              control = list(), input = list(), opt.method = fit$opt.method)
  class(out) <- "selm"
  out
}


## function to convert 'family' argument to arguments as in package 'sn'
get_selm_args <- function(family = "student") {
  # convert argument
  if (family == "student") {
    family <- "ST"
    fixed.param <- list(alpha = 0)  # no skewness
  } else {
    family <- if (family == "skewnormal") "SN" else "ST"
    fixed.param <- list()
    # fixed.param <- if (family == "skewt") list(nu = 3) else list()
  }
  # return list of arguments
  list(family = family, fixed.param = fixed.param)
}


## function to extract regression coefficients from 'selm' objects
get_coef <- function(param) {
  coefficients <- param$dp
  keep <- setdiff(names(coefficients), c("omega", "alpha", "nu"))
  coefficients[keep]
}

## define a coef() method as the one from package 'sn' in some cases doesn't
## actually return the regression coefficients (and also returns additional
## parameters)
# ## @importFrom methods setMethod
# setMethod("coef", "selm", function(object, ...) {
#   get_coef(object@param)
# })
coef.selm <- function(object, ...) get_coef(object$param)


## functions to extract all parameters in direct or centered parametrization
## this is need to extract starting values for bootstrap replicates
get_dp <- function(object) object$param$dp
get_cp <- function(object) object$param$cp


## for regression with skew-elliptical errors: return list that contains
## coefficient matrix and information on parameters of error distribution
#' @importFrom methods new
#' @importFrom sn summary
get_summary.selm <- function(object, ...) {
  # construct S4 object and compute the usual summary
  objectS4 <- new("selm", call = call("selm"), family = object$family,
                  logL = object$logL, method = object$method,
                  param = object$param, param.var = object$param.var,
                  size = object$size, residuals.dp = object$residuals.dp,
                  fitted.values.dp = object$fitted.values.dp,
                  control = object$control, input = object$control,
                  opt.method = object$opt.method)
  summary <- summary(objectS4, param = "DP")
  # extract coefficients
  param_table <- summary@param.table
  # change column names to be consistent with other models
  colnames(param_table) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  # split up in regression coefficients and parameters of error distribution
  keep <- 1:object$size["p"]
  coefficients <- param_table[keep, , drop = FALSE]
  parameters <- param_table[-keep, , drop = FALSE]
  # return results
  result <- list(coefficients = coefficients, parameters = parameters)
  class(result) <- "summary_selm"
  result
}

#' @export
print.summary_selm <- function(x, digits = max(3, getOption("digits")-3),
                               signif.stars = getOption("show.signif.stars"),
                               signif.legend = signif.stars, ...) {
  # print coefficient matrix
  cat("Coefficients:\n")
  printCoefmat(x$coefficients, digits = digits, signif.stars = signif.stars,
               signif.legend = FALSE, ...)
  cat("\nParameters of error distribution:\n")
  printCoefmat(x$parameters, digits = digits, signif.stars = signif.stars,
               signif.legend = signif.legend, ...)
  # return object invisibly
  invisible(x)
}


## extract number of observations from median regression
#' @export
#' @import stats
nobs.selm <- function(object, ...) unname(object$size["n.obs"])


# there is no confint() method for median regression results
#' @importFrom methods new
#' @importFrom sn summary
#' @export
confint.selm <- function(object, parm = NULL, level = 0.95, ...) {
  # construct S4 object and compute the usual summary
  objectS4 <- new("selm", call = call("selm"), family = object$family,
                  logL = object$logL, method = object$method,
                  param = object$param, param.var = object$param.var,
                  size = object$size, residuals.dp = object$residuals.dp,
                  fitted.values.dp = object$fitted.values.dp,
                  control = object$control, input = object$control,
                  opt.method = object$opt.method)
  summary <- summary(objectS4, param = "DP")
  # extract coefficient matrix
  coef_mat <- summary@param.table
  keep <- 1:object$size["p"]
  coef_mat <- coef_mat[keep, , drop = FALSE]
  coef_names <- rownames(coef_mat)
  # extract point estimates and standard errors
  coef <- coef_mat[, 1L]
  se <- coef_mat[, 2L]
  # check parameters to extract
  if (missing(parm)) parm <- coef_names
  else if (is.numeric(parm)) parm <- coef_names[parm]
  # significance level and quantile of the t distribution
  a <- (1 - level) / 2
  a <- c(a, 1 - a)
  q <- qnorm(a)
  # column names for output should contain percentages
  cn <- paste(format(100 * a, trim=TRUE), "%")
  # construct confidence interval
  ci <- array(NA_real_, dim = c(length(parm), 2L), dimnames = list(parm, cn))
  ci[] <- coef[parm] + se[parm] %o% q
  ci
}
