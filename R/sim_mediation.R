# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------


## function to simulate data from a fitted mediation model

#' @export
sim_mediation <- function(object, n, ...) UseMethod("sim_mediation")

#' @export
# explanatory ... character string indicating whether to draw the explanatory
#                 variables from a (normal) distribution, or whether to
#                 bootstrap the explanatory variables.
# errors ........ character string indicating whether to draw the error terms
#                 from the fitted model distribution, or whether to bootstrap
#                 the error terms from the observed residuals.
sim_mediation.fit_mediation <- function(object, n = NULL,
                                        explanatory = c("sim", "boot"),
                                        errors = c("sim", "boot"), ...) {

  # TODO: variables with only a few unique values should be drawn from a
  #       multinomial distribution with the observed probabilities

  # initializations
  explanatory <- match.arg(explanatory)
  errors <- match.arg(errors)
  if (is.null(n)) n <- nrow(object$data)
  # extract relevant information
  x <- object$x
  y <- object$y
  m <- object$m
  p_m <- length(object$m)
  covariates <- object$covariates
  model <- object$model
  if (is.null(model)) model <- "simple"

  # simulate or bootstrap explanatory variables
  if (explanatory == "sim") X <- sim_explanatory(object, n = n)
  else X <- boot_explanatory(object, n = n)

  # simulate or bootstrap error terms
  if (errors == "sim") e <- sim_errors(object, n = n)
  else e <- boot_errors(object, n = n)

  # extract coefficients
  coef <- get_coefficients(object)

  # simulate mediators under the model
  if (model == "parallel") {
    # parallel mediators
    M <- sapply(seq_len(p_m), function(j) {
      coef_M <- coef$M[[j]]
      coef_M[1L] + X %*% coef_M[-1L] + e$M[[j]]
    })
  } else if (model == "serial") {
    # serial mediators
    M <- matrix(NA_real_, nrow = n, ncol = p_m)
    for (j in seq_len(p_m)) {
      coef_M <- coef$M[[j]]
      MX <- cbind(M[, seq_len(j-1)], X)
      M[, j] <- coef_M[1L] + MX %*% coef_M[-1L] + e$M[[j]]
    }
  } else {
    # single mediator
    coef_M <- coef$M
    M <- coef_M[1L] + X %*% coef_M[-1L] + e$M
  }
  # add column names to mediators
  colnames(M) <- m
  # compute the dependent variable under the model
  coef_Y <- coef$Y
  MX <- cbind(M, X)
  Y <- coef_Y[1L] + MX %*% coef_Y[-1L] + e$Y
  colnames(Y) <- y

  # return simulated data with variables in same order as 'data' component
  data.frame(X[, x, drop = FALSE], Y, M, X[, covariates, drop = FALSE])

}

#' @export
sim_mediation.test_mediation <- function(object, n = NULL, ...) {
  # call method for model fit
  sim_mediation(object$fit, n = n, ...)
}


## wrapper to conform to R convention regarding name and first argument
#' @export
rmediation <- function(n, object, ...) sim_mediation(object, n = n, ...)


## simulate explanatory variables based on a mediation model fit

sim_explanatory <- function(object, n) UseMethod("sim_explanatory")

sim_explanatory.reg_fit_mediation <- function(object, n) {
  # initializations
  data <- object$data
  n_data <- nrow(data)
  predictors <- c(object$x, object$covariates)
  robust <- object$robust
  family <- object$family
  # perform separate regressions with each explanatory variable as response
  # variable and no explanatory variables (just an intercept)
  null_matrix <- matrix(NA_real_, nrow = n_data, ncol = 0)
  if (robust == "MM") {
    # MM-estimator for robust regression
    fits <- lapply(predictors, function(y) {
      lmrob_fit(null_matrix, data[, y], control = object$control)
    })
  } else if (robust == "median") {
    # LAD-estimator for median regression
    fits <- lapply(predictors, function(y) {
      rq_fit(null_matrix, data[, y], tau = 0.5)
    })
  } else if (family == "gaussian") {
    # OLS estimator
    fits <- lapply(predictors, function(y) lm_fit(null_matrix, data[, y]))
  } else if (family == "select") {
    # select among normal, skew-normal, t and skew-t errors
    fits <- lapply(predictors, function(y) lmselect_fit(null_matrix, data[, y]))
  } else {
    # obtain parameters as required for package 'sn'
    selm_args <- get_selm_args(family)
    # perform regression with skew-elliptical errors
    fits <- lapply(predictors, function(y) {
      selm_fit(null_matrix, data[, y], family = selm_args$family,
               fixed.param = selm_args$fixed.param)
    })
  }
  # add names to list of regression fits
  names(fits) <- predictors
  # draw explanatory variables by drawing error terms and adding intercept
  X <- sapply(fits, function(fit) unname(coef(fit)) + sim_errors(fit, n = n))
  # return explanatory variables
  X
}

sim_explanatory.cov_fit_mediation <- function(object, n) {
  # initializations
  x <- object$x
  # extract means and standard deviations of the independent variable from
  # the already computed means and covariance matrix stored in the model fit
  center <- object$cov$center[x]
  scale <- sqrt(object$cov$cov[x, x])
  # draw explanatory variable from normal distribution
  X <- rnorm(n, mean = center, sd = scale)
  # return explanatory variable as matrix
  X <- as.matrix(X)
  colnames(X) <- x
  X
}


## simulate error terms based on a mediation model fit

sim_errors <- function(object, n) UseMethod("sim_errors")

sim_errors.reg_fit_mediation <- function(object, n) {
  # initializations
  p_m <- length(object$m)
  # draw errors for mediators
  if (p_m == 1L) e_M <- sim_errors(object$fit_mx, n = n)
  else e_M <- lapply(object$fit_mx, sim_errors, n = n)
  # draw errors for dependent variable
  e_Y <- sim_errors(object$fit_ymx, n = n)
  # return list of errors
  list(M = e_M, Y = e_Y)
}

sim_errors.lmrob <- function(object, n) {
  sigma <- object$scale           # residual scale
  rnorm(n, mean = 0, sd = sigma)  # return error terms
}

sim_errors.rq <- function(object, n) {
  if (object$tau != 0.5) stop("only implemented for median regression")
  sigma <- mad(residuals(object), center = 0)  # residual scale
  rnorm(n, mean = 0, sd = sigma)               # return error terms
}

sim_errors.lm <- function(object, n) {
  rss <- sum(residuals(object)^2)          # residual sum of squares
  sigma <- sqrt(rss / object$df.residual)  # residual scale
  rnorm(n, mean = 0, sd = sigma)           # return error terms
}

#' @importFrom sn rsn rst
sim_errors.lmse <- function(object, n) {
  # family specification as in package 'sn' (not 'robmed'): either "SN" or "ST"
  which <- object$family
  param <- object$param
  family <- get_family(which, param)
  # extract parameters of residual distribution
  dp <- get_dp(object)
  # draw error terms with fitted parameters in direct parametrization
  if (which == "SN") {
    e <- rsn(n, xi = 0, omega = dp["omega"], alpha = dp["alpha"])
  } else {
    e <- rst(n, xi = 0, omega = dp["omega"], alpha = dp["alpha"], nu = dp["nu"])
  }
  # remove attributes since rsn() and rst() store family and parameters
  attributes(e) <- NULL
  # transform to centered parametrization
  # (doesn't matter for the symmetric t-distribution, where 'mu0' is zero)
  if (family != "student") e <- e - unname(param$mu0)
  # return error terms
  e
}

sim_errors.cov_fit_mediation <- function(object, n) {
  # -----
  # Note: since Zu & Yuan (2010) use the maximum likelihood estimator of the
  # covariance matrix, we get the maximum likelihood estimates of the residual
  # scale.  For the nonrobust version, the resulting scales therefore differ
  # slightly from those of the OLS regression fit.
  # -----
  # initializations
  x <- object$x
  y <- object$y
  m <- object$m
  # extract variances of the variables
  cov <- object$cov$cov
  # extract effects
  a <- object$a
  b <- object$b
  direct <- object$direct
  # extract residual scales and draw errors for mediators
  sigma_M <- sqrt(cov[m, m] - a^2 * cov[x, x])
  e_M <- rnorm(n, mean = 0, sd = sigma_M)
  # extract residual scales and draw errors for dependent variable
  sigma_Y <- sqrt(cov[y, y] - b^2 * cov[m, m] - direct^2 * cov[x, x] -
                    b * direct * cov[m, x])
  e_Y <- rnorm(n, mean = 0, sd = sigma_Y)
  # return list of errors
  list(M = e_M, Y = e_Y)
}


## extract coefficients

get_coefficients <- function(object) UseMethod("get_coefficients")

get_coefficients.reg_fit_mediation <- function(object) {
  # initializations
  p_m <- length(object$m)
  # extract coefficients of models for mediators
  if (p_m == 1L) coef_M <- coef(object$fit_mx)
  else coef_M <- lapply(object$fit_mx, coef)
  # extract coefficients of model for dependent variable
  coef_Y <- coef(object$fit_ymx)
  # return list of coefficients
  list(M = coef_M, Y = coef_Y)
}

get_coefficients.cov_fit_mediation <- function(object) {
  # initializations
  x <- object$x
  y <- object$y
  m <- object$m
  # extract means of variables
  centers <- object$cov$center
  # extract coefficients of model for mediator
  a <- object$a
  coef_M <- c(centers[m] - a * centers[x], a)
  names(coef_M) <- c("(Intercept)", x)
  # extract coefficients of model for dependent variable
  b <- object$b
  direct <- object$direct
  coef_Y <- c(centers[y] - b * centers[m] - direct * centers[x], b, direct)
  names(coef_Y) <- c("(Intercept)", m, x)
  # return list of coefficients
  list(M = coef_M, Y = coef_Y)
}


## bootstrap explanatory variables from a mediation model fit
boot_explanatory <- function(object, n) {
  # initializations
  data <- as.matrix(object$data, rownames.force = FALSE)
  n_data <- nrow(data)
  predictors <- c(object$x, object$covariates)
  # draw bootstrap sample of explanatory variables
  i <- sample.int(n_data, size = n, replace = TRUE)
  data[i, predictors, drop = FALSE]
}


## bootstrap error terms from residuals of a mediation model fit

boot_errors <- function(object, n) UseMethod("boot_errors")

boot_errors.reg_fit_mediation <- function(object, n) {
  # initializations
  p_m <- length(object$m)
  # extract residuals and bootstrap errors for mediators
  if (p_m == 1L) {
    residuals_M <- unname(residuals(object$fit_mx))
    e_M <- sample(residuals_M, size = n, replace = TRUE)
  } else {
    e_M <- lapply(object$fit_mx, function(fit) {
      residuals_M <- unname(residuals(fit))
      sample(residuals_M, size = n, replace = TRUE)
    })
  }
  # extract residuals and bootstrap errors for dependent variable
  residuals_Y <- unname(residuals(object$fit_ymx))
  e_Y <- sample(residuals_Y, size = n, replace = TRUE)
  # return list of errors
  list(M = e_M, Y = e_Y)
}

boot_errors.cov_fit_mediation <- function(object, n) {
  # initializations
  x <- object$x
  y <- object$y
  m <- object$m
  data <- as.matrix(object$data, rownames.force = FALSE)
  # extract coefficients
  coef <- get_coefficients(object)
  # compute residuals and bootstrap errors for mediators
  residuals_M <- data[, m] - coef$M[1L] - coef$M[x] * data[, x]
  e_M <- sample(residuals_M, size = n, replace = TRUE)
  # compute residuals and bootstrap errors for dependent variable
  mx <- c(m, x)
  residuals_Y <- data[, y] - coef$Y[1L] - drop(data[, mx] %*% coef$Y[mx])
  e_Y <- sample(residuals_Y, size = n, replace = TRUE)
  # return list of errors
  list(M = e_M, Y = e_Y)
}
