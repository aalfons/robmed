# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------


## function to simulate data from a fitted mediation model

#' @export
sim_mediation <- function(object, n, ...) UseMethod("sim_mediation")

#'
sim_mediation.default <- function(object, n, model = "simple",
                                  family = "gaussian", ...) {
  # 'object' should be a list of effects? or the original data set?
}

#' @export
sim_mediation.fit_mediation <- function(object, n = NULL, ...) {
  # get summary of model fit to obtain standard deviations and standard errors
  summary <- get_summary(object)
  if (is.null(n)) n <- summary$n
  # temporarily return summary
  summary
}

#' @export
sim_mediation.test_mediation <- function(object, n = NULL, ...) {
  # call method for model fit
  sim_mediation(object$fit, n = n, ...)
}

#' @export
# explanatory ... character string indicating whether to draw the explanatory
#                 variables from a (normal) distribution, or whether to
#                 bootstrap the explanatory variables.
# errors ........ character string indicating whether to draw the error terms
#                 from the fitted model distribution, or whether to bootstrap
#                 the error terms from the observed residuals.
sim_mediation.boot_test_mediation <- function(object, n = NULL,
                                              explanatory = c("sim", "boot"),
                                              errors = c("sim", "boot"),
                                              ...) {

  # initializations
  explanatory <- match.arg(explanatory)
  errors <- match.arg(errors)
  # get summary of model fit to obtain standard deviations and standard errors
  # component 'boot' only exists for bootstrap test, otherwise NULL
  fit <- object$fit
  if (is.null(n)) n <- nobs(fit)
  # extract relevant information
  x <- fit$x
  y <- fit$y
  m <- fit$m
  p_m <- length(fit$m)
  covariates <- fit$covariates
  model <- fit$model
  if (is.null(model)) model <- "simple"

  # simulate or bootstrap explanatory variables
  if (explanatory == "sim") X <- sim_explanatory(fit, n = n)
  else stop("not implemented yet")

  # simulate or bootstrap error terms
  if (errors == "sim") e <- sim_errors(fit, n = n)
  else stop("not implemented yet")

  # extract coefficients
  coef <- get_coefficients(fit)

  # simulate mediators under the model
  if (model == "parallel") {
    # parallel mediators
    M <- sapply(seq_len(p_m), function(j) {
      coef_M <- coef$M[[j]]
      coef_M[1] + X %*% coef_M[-1] + e$M[[j]]
    })
  } else if (model == "serial") {
    # serial mediators
    M <- matrix(NA_real_, nrow = n, ncol = p_m)
    for (j in seq_len(p_m)) {
      coef_M <- coef$M[[j]]
      MX <- cbind(M[, seq_len(j-1)], X)
      M[, j] <- coef_M[1] + MX %*% coef_M[-1] + e$M[[j]]
    }
  } else {
    # single mediator
    coef_M <- coef$M
    M <- coef_M[1] + X %*% coef_M[-1] + e$M
  }
  # add column names to mediators
  colnames(M) <- m
  # compute the dependent variable under the model
  coef_Y <- coef$Y
  MX <- cbind(M, X)
  Y <- coef_Y[1] + MX %*% coef_Y[-1] + e$Y
  colnames(Y) <- y

  # return simulated data with variables in same order as 'data' component
  data.frame(X[, x, drop = FALSE], Y, M, X[, covariates, drop = FALSE])

}


## wrapper to conform to R convention regarding name and first argument
#' @export
rmediation <- function(n, object, ...) sim_mediation(object, n = n, ...)


## simulate explanatory variables based on a mediation model fit

sim_explanatory <- function(object, n) UseMethod("sim_explanatory")

sim_explanatory.reg_fit_mediation <- function(object, n) {
  # initializations
  data <- object$data
  predictors <- c(object$x, object$covariates)
  have_robust <- is_robust(fit)
  family <- object$family
  # extract relevant parameters and draw from respective distribution
  if (family == "gaussian") {
    if (have_robust) {
      # TODO: use the corresponding regression method with each variable as
      # response and a constant term as explanatory variable to determine the
      # centers and scales of the normal distribution
      centers <- sapply(data[, predictors, drop = FALSE], median)
      scales <- sapply(data[, predictors, drop = FALSE], mad)
    } else {
      # compute the mean and standard deviations of explanatory variables to
      # determine the parameters of the normal distributions
      centers <- sapply(data[, predictors, drop = FALSE], mean)
      scales <- sapply(data[, predictors, drop = FALSE], sd)
    }
    # draw explanatory variables from normal distributions
    X <- mapply(rnorm, n, mean = centers, sd = scales)
  } else {
    stop("not implemented yet")
    # TODO: use the corresponding regression method with each variable as
    # response and a constant term as explanatory variable to determine the
    # parameters to draw from the respective distribution
    # TODO: draw explanatory variables from skew-elliptical distributions
  }
  # return explanatory variables
  colnames(X) <- predictors
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
  m <- object$m
  p_m <- length(m)
  family <- object$family
  # draw error terms from the respective error distributions
  if (family == "gaussian") {
    # extract residual scales and draw errors for mediators
    if (p_m == 1L) {
      sigma_M <- get_scale(object$fit_mx)
      e_M <- rnorm(n, sd = sigma_M)
    } else {
      e_M <- lapply(object$fit_mx, function(fit) {
        sigma_M <- get_scale(fit)
        rnorm(n, sd = sigma_M)
      })
    }
    # extract residual scale and draw errors for dependent variable
    sigma_Y <- get_scale(object$fit_ymx)
    e_Y <- rnorm(n, sd = sigma_Y)
  } else {
    stop("not implemented yet")
  }
  # return list of errors
  list(M = e_M, Y = e_Y)
}

sim_errors.cov_fit_mediation <- function(object, n) {
  stop("not implemented yet")
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
  stop("not implemented yet")
}


## extract residual scale

get_scale <- function(object) UseMethod("get_scale")

get_scale.lmrob <- function(object) object$scale

get_scale.lm <- function(object) {
  rss <- sum(residuals(object)^2)  # residual sum of squares
  sqrt(rss / object$df.residual)   # residual scale
}

get_scale.rq <- function(object) {
  if (object$tau != 0.5) stop("only implemented for median regression")
  mad(residuals(object), center = 0)  # MAD with median residual set to 0
}
