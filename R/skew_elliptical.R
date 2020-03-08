# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------


# regression with skew-elliptical errors
#' @importFrom sn selm.fit
#' @importFrom methods new
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
  # return model fit as object of S4 class "selm"
  # setClass(fit) <- "selm"
  new("selm", call = call("selm"), family = family, logL = fit$logL,
      method = control$method, param = fit$param, param.var = fit$param.var,
      size = fit$size, residuals.dp = fit$resid.dp,
      fitted.values.dp = fit$fitted.dp, control = list(), input = list(),
      opt.method = fit$opt.method)
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
#' @importFrom methods setMethod
setMethod("coef", "selm", function(object, ...) {
  get_coef(object@param)
})
