# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------


## regression with skew-elliptical errors
#' @importFrom sn selm.fit
## @importFrom methods new
selm_fit <- function(x, y, intercept = TRUE, family = "ST",
                     fixed.param = list()) {
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
  # (called "lmse" to avoid conflicts with package 'sn')
  out <- list(call = NULL, family = family, logL = fit$logL,
              method = control$method, param = fit$param,
              param.var = fit$param.var, size = fit$size,
              residuals.dp = fit$resid.dp, fitted.values.dp = fit$fitted.dp,
              control = list(), input = list(), opt.method = fit$opt.method)
  class(out) <- "lmse"
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
  }
  # return list of arguments
  list(family = family, fixed.param = fixed.param)
}

## convert from class structure as in package 'sn'
get_family <- function(family, param) {
  if (is.null(family) && is.null(param)) "gaussian"
  else if (family == "SN") "skewnormal"
  else if (family == "ST") {
    fixed <- param$fixed
    alpha <- fixed$alpha
    if (length (fixed) == 0 || is.null(alpha) || alpha != 0) "skewt"
    else "student"
  }
}


## workhorse function to extract coefficients from lists returned by selm.fit()
## (this is needed for bootstrap replicates)
get_coef <- function(param, family) {
  # initializations
  family <- get_family(family, param)
  # extract coefficients
  if (family == "student") {
    coefficients <- param$dp
    remove <- c("omega", "alpha", "nu")
  } else {
    coefficients <- param$cp
    if (is.null(coefficients)) {
      # sometimes the centered parametrization of the skew-t distribution
      # doesn't exist (this is a hack to make bootstrap work; it may not be
      # fully justified, but it only affects the intercept)
      coefficients <- param[["pseudo-cp"]]
      names(coefficients) <- gsub("~", "", names(coefficients))
    }
    remove <- c("s.d.", "gamma1", "gamma2")
  }
  keep <- setdiff(names(coefficients), remove)
  coefficients[keep]
}

## method to extract regression coefficients from 'lmse' objects
#' @export
coef.lmse <- function(object, ...) get_coef(object$param, object$family)


## functions to extract all parameters in direct or centered parametrization
## (this is needed to extract starting values for bootstrap replicates)
get_dp <- function(object) object$param$dp
get_cp <- function(object) object$param$cp


## summary of regression with skew-elliptical errors: return list that contains
## coefficient matrix and information on parameters of error distribution
#' @importFrom methods new
#' @importFrom sn summary
get_summary.lmse <- function(object, ...) {
  # construct S4 object as in package 'sn'
  objectS4 <- new("selm", call = call("selm"), family = object$family,
                  logL = object$logL, method = object$method,
                  param = object$param, param.var = object$param.var,
                  size = object$size, residuals.dp = object$residuals.dp,
                  fitted.values.dp = object$fitted.values.dp,
                  control = object$control, input = object$control,
                  opt.method = object$opt.method)
  # use summary methods from package 'sn'
  family <- get_family(object$family, object$param)
  have_student <- family == "student"
  if (have_student) summary <- summary(objectS4, param.type = "DP")
  else summary <- summary(objectS4, param.type = "CP")
  # extract coefficients
  param_table <- summary@param.table
  # change column names to be consistent with other models
  colnames(param_table) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  # split up in regression coefficients and parameters of error distribution
  keep <- 1:object$size["p"]
  coefficients <- param_table[keep, , drop = FALSE]
  parameters <- param_table[-keep, 1:2, drop = FALSE]
  # parameters <- param_table[-keep, , drop = FALSE]
  # # for degrees of freedom, test if 1/nu > 0 rather than nu = 0
  # if (have_student) {
  #   which <- which(rownames(parameters) == "nu")
  #   nu <- parameters[which, 1]
  #   se <- parameters[which, 2]
  #   nu_inv <- 1 / nu
  #   se_inv <- nu_inv^2 * se
  #   z_inv <- nu_inv / se_inv
  #   p_inv <- p_value_z(z_inv, alternative = "greater")
  #   # define extra parameter matrix for degrees of freedom
  #   dn <- list("nu   ", c("Estimate", "Std. Error", "z value", "Pr(>z)"))
  #   nu <- matrix(c(nu, se, z_inv, p_inv), nrow = 1, ncol = 4, dimnames = dn)
  #   # return list instead of parameter matrix
  #   parameters <- list(other = parameters[-which,], nu = nu)
  # }
  # return results
  result <- list(coefficients = coefficients, parameters = parameters)
  class(result) <- "summary_lmse"
  result
}

#' @export
print.summary_lmse <- function(x, digits = max(3, getOption("digits")-3),
                               signif.stars = getOption("show.signif.stars"),
                               signif.legend = signif.stars, ...) {
  # print coefficient matrix
  cat("Coefficients:\n")
  printCoefmat(x$coefficients, digits = digits, signif.stars = signif.stars,
               signif.legend = signif.legend, ...)
  # printCoefmat(x$coefficients, digits = digits, signif.stars = signif.stars,
  #              signif.legend = FALSE, ...)
  # print information on other parameters
  cat("\nParameters of error distribution:\n")
  print(x$parameters, digits = digits, ...)
  # if (is.list(x$parameters)) {
  #   print(x$parameters$other, digits = digits, signif.stars = signif.stars,
  #         signif.legend = FALSE, ...)
  #   print(x$parameters$nu, digits = digits, signif.stars = signif.stars,
  #         signif.legend = signif.legend, ...)
  # } else {
  #   print(x$parameters, digits = digits, signif.stars = signif.stars,
  #         signif.legend = signif.legend, ...)
  # }
  # return object invisibly
  invisible(x)
}


## extract number of observations from regression with skew-elliptical errors
#' @export
#' @import stats
nobs.lmse <- function(object, ...) unname(object$size["n.obs"])


## extract confidence intervals from regression with skew-elliptical errors
#' @importFrom methods new
#' @importFrom sn confint
#' @export
confint.lmse <- function(object, parm = NULL, level = 0.95, ...) {
  # construct S4 object as in package 'sn'
  objectS4 <- new("selm", call = call("selm"), family = object$family,
                  logL = object$logL, method = object$method,
                  param = object$param, param.var = object$param.var,
                  size = object$size, residuals.dp = object$residuals.dp,
                  fitted.values.dp = object$fitted.values.dp,
                  control = object$control, input = object$control,
                  opt.method = object$opt.method)
  # use summary methods from package 'sn'
  family <- get_family(object$family, object$param)
  have_student <- family == "student"
  if (have_student) ci <- confint(objectS4, level = level, param.type = "DP")
  else ci <- confint(objectS4, level = level, param.type = "CP")
  # extract confidence intervals for regression coefficients
  keep <- 1:object$size["p"]
  ci <- ci[keep, , drop = FALSE]
  coef_names <- rownames(ci)
  # check parameters to extract
  if (missing(parm)) parm <- coef_names
  else if (is.numeric(parm)) parm <- coef_names[parm]
  # return confidence interval
  ci[parm, , drop = FALSE]
}


## function to extract log-likelihood (necessary to compute AIC/BIC)
#' @export
logLik.lmse <- function(object, ....) {
  logL <- object$logL
  attr(logL, "df") <- unname(object$size["n.param"])
  class(logL) <- "logLik"
  logL
}


## regression with selecting family of error distribution via BIC
## selm.fit() always requires constant for intercept, so argument is ignored
#' @importFrom sn selm.fit
## @importFrom methods new
lmselect_fit <- function(x, y, intercept = TRUE) {
  # selm.fit() always requires constant for intercept, so argument is ignored
  # fit model with normal errors
  fit_gaussian <- lm_fit(x, y)
  bic_gaussian <- BIC(fit_gaussian)
  # start with model for normal errors
  bic <- bic_gaussian
  fit <- fit_gaussian
  # fit model with skew-normal errors and check if it fits better
  fit_skewnormal <- selm_fit(x, y, family = "SN")
  bic_skewnormal <- BIC(fit_skewnormal)
  if (bic_skewnormal < bic) {
    bic <- bic_skewnormal
    fit <- fit_skewnormal
  }
  # fit model with t errors and check if it fits better
  fit_student <- selm_fit(x, y, family = "ST", fixed.param = list(alpha = 0))
  bic_student <- BIC(fit_student)
  if (bic_student < bic) {
    bic <- bic_student
    fit <- fit_student
  }
  # fit model with skew-t errors only if both skew-normal and t distribution
  # are an improvement to normal errors
  if (bic_skewnormal < bic_gaussian && bic_student < bic_gaussian) {
    fit_skewt <- selm_fit(x, y, family = "ST")
    bic_skewt <- BIC(fit_skewt)
    if (bic_skewt < bic) {
      bic <- bic_skewt
      fit <- fit_skewt
    }
  } else fit_skewt <- NULL
  # add starting values for bootstrap
  fit$start <- list(skewnormal = get_cp(fit_skewnormal),
                    student = get_dp(fit_student),
                    skewt = get_dp(fit_skewt))
  # return best model fit according to BIC
  fit
}

## regression with selecting family of error distribution via BIC
## (this is a more barebones version to be used within bootstrap)
lmselect_boot <- function(x, y, start = NULL, control = list(method = "MLE")) {
  # initializations
  d <- dim(x)
  n <- d[1]
  log_n <- log(n)
  p <- d[2]
  # fit model with normal errors
  coef_gaussian <- drop(solve(crossprod(x)) %*% crossprod(x, y))
  # compute BIC
  residuals <- y - x %*% coef_gaussian
  sigma2 <- sum(residuals^2) / n
  bic_gaussian <- n * (log(2 * pi) + 1 + log(sigma2)) + (p + 1) * log_n
  # start with model for normal errors
  bic <- bic_gaussian
  coefficients <- coef_gaussian
  # fit model with skew-normal errors and check if it fits better
  fit_skewnormal <- selm.fit(x, y, family = "SN", start = start$skewnormal,
                             selm.control = control)
  bic_skewnormal <- -2 * fit_skewnormal$logL + (p + 2) * log_n
  if (bic_skewnormal < bic) {
    bic <- bic_skewnormal
    coefficients <- get_coef(fit_skewnormal$param, family = "SN")
  }
  # fit model with t errors and check if it fits better
  fit_student <- selm.fit(x, y, family = "ST", start = start$student,
                          fixed.param = list(alpha = 0),
                          selm.control = control)
  bic_student <- -2 * fit_student$logL + (p + 2) * log_n
  if (bic_student < bic) {
    bic <- bic_student
    coefficients <- get_coef(fit_student$param, family = "ST")
  }
  # fit model with skew-t errors only if both skew-normal and t distribution
  # are an improvement to normal errors
  if (bic_skewnormal < bic_gaussian && bic_student < bic_gaussian) {
    fit_skewt <- selm.fit(x, y, family = "ST", start = start$skewt,
                          selm.control = control)
    bic_skewt <- -2 * fit_skewt$logL + (p + 3) * log_n
    if (bic_skewt < bic) {
      bic <- bic_skewt
      coefficients <- get_coef(fit_skewt$param, family = "ST")
    }
  }
  # return best model fit according to BIC
  unname(coefficients)
}
