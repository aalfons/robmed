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
# fitted ... character string indicating whether to simulate the fitted values
#            based on drawing the explanatory variables from a distribution and
#            the observed coefficient estimates, or whether to bootstrap the
#            fitted values from the model estimates.
# ('fitted' is not the best name for this, as this refers to the fitted values
# under the true simulated model, and not to the fitted values of an estimated
# model.)
# errors ... character string indicating whether to simulate the error terms
#            from the fitted model distribution, or whether to bootstrap the
#            error terms from the observed residuals.
sim_mediation.boot_test_mediation <- function(object, n = NULL,
                                              type = c("boot", "data"),
                                              explanatory = c("sim", "boot"),
                                              errors = c("sim", "boot"),
                                              ...) {

  # initializations
  type <- match.arg(type)
  explanatory <- match.arg(explanatory)
  errors <- match.arg(errors)
  # get summary of model fit to obtain standard deviations and standard errors
  # component 'boot' only exists for bootstrap test, otherwise NULL
  fit <- object$fit
  if (type == "boot") summary <- get_summary(fit, boot = object$reps)
  else summary <- get_summary(fit)
  if (is.null(n)) n <- summary$n
  # extract relevant information
  x <- fit$x
  covariates <- fit$covariates

  # simulate or bootstrap explanatory variables
  if (explanatory == "sim") X <- sim_explanatory(fit, n = n)
  else stop("not implemented yet")

  # simulate or bootstrap error terms
  if (errors == "sim") e <- sim_errors(summary, n = n)
  else stop("not implemented yet")

  # which type of coefficients to be used as true values
  which <- switch(type, boot = "Boot", data = "Estimate")
  # extract coefficients
  if (inherits(fit, "reg_fit_mediation")) {
    # additional information
    y <- fit$y
    m <- fit$m
    p_m <- length(m)
    model <- fit$model
    # compute the mediators under the model
    if (model == "parallel") {
      # parallel mediators
      M <- sapply(seq_len(p_m), function(j) {
        coef_M <- summary$fit_mx[[j]]$coefficients[, which]
        coef_M[1] + X %*% coef_M[-1] + e$M[, j]
      })
    } else if (model == "serial") {
      # serial mediators
      M <- matrix(NA_real_, nrow = n, ncol = p_m)
      for (j in seq_len(p_m)) {
        coef_M <- summary$fit_mx[[j]]$coefficients[, which]
        MX <- cbind(M[, seq_len(j-1)], X)
        M[, j] <- coef_M[1] + MX %*% coef_M[-1] + e$M[, j]
      }
    } else {
      # single mediator
      coef_M <- summary$fit_mx$coefficients[, which]
      M <- coef_M[1] + X %*% coef_M[-1] + e$M
    }
    # add column names to mediators
    colnames(M) <- m
    # compute the dependent variable under the model
    coef_Y <- summary$fit_ymx$coefficients[, which]
    MX <- cbind(M, X)
    Y <- coef_Y[1] + MX %*% coef_Y[-1] + e$Y
    colnames(Y) <- y
  } else if (inherits(fit, "cov_fit_mediation")) {
    stop("not implemented yet")
  } else stop("not implemented for this type of model fit")

  # return simulated data with variables in same order as 'data' component
  data.frame(X[, x, drop = FALSE], Y, M, X[, covariates, drop = FALSE])

}


## wrapper to conform to R convention regarding name and first argument
#' @export
rmediation <- function(n, object, ...) sim_mediation(object, n = n, ...)



## utility functions

# simulate explanatory variables based on a mediation model fit

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


# simulate error terms based on a mediation model fit

sim_errors <- function(object, n) UseMethod("sim_errors")

sim_errors.summary_reg_fit_mediation <- function(object, n) {
  # initializations
  m <- object$m
  p_m <- length(m)
  robust <- object$robust
  family <- object$family
  # draw error terms from the respective error distributions
  if (robust == "median") {
    stop("not implemented yet")
  } else if (family == "gaussian") {
    # extract residual scales for mediators
    if (p_m == 1L) sigma_M <- object$fit_mx$s$value
    else sigma_M <- sapply(object$fit_mx, function(fit) fit$s$value)
    # draw errors for mediators
    e_M <- mapply(rnorm, n, sd = sigma_M)
    colnames(e_M) <- m
    # extract residual scale and draw errors for dependent variable
    sigma_Y <- object$fit_ymx$s$value
    e_Y <- rnorm(n, sd = sigma_Y)
  } else {
    stop("not implemented yet")
  }
  # return errors
  list(M = e_M, Y = e_Y)
}

sim_errors.summary_cov_fit_mediation <- function(object, n) {
  stop("not implemented yet")
}
