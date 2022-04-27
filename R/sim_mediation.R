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
                                        errors = c("sim", "boot"),
                                        num_discrete = 10, ...) {

  # initializations
  if (is.null(n)) n <- nrow(object$data)
  else {
    n <- rep(as.integer(n), length.out = 1L)
    if (is.na(n) || n <= 0L) stop("'n' must be a positive integer")
  }
  explanatory <- match.arg(explanatory)
  errors <- match.arg(errors)
  # extract relevant information
  x <- object$x
  y <- object$y
  m <- object$m
  p_m <- length(object$m)
  covariates <- object$covariates
  model <- object$model
  if (is.null(model)) model <- "simple"

  # simulate or bootstrap explanatory variables
  if (explanatory == "sim") {
    # explanatory variables are simulated independently from each other and
    # therefore there are no correlations between them
    X <- sim_explanatory(object, n = n, num_discrete = num_discrete)
  } else X <- boot_explanatory(object, n = n)

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

sim_explanatory <- function(object, n, ...) UseMethod("sim_explanatory")

sim_explanatory.reg_fit_mediation <- function(object, n, num_discrete = 10,
                                              ...) {

  # initializations
  num_discrete <- rep(as.integer(num_discrete), length.out = 1L)
  if (is.na(num_discrete) || num_discrete <= 0L) {
    num_discrete <- formals()$num_discrete
  }
  # extract relevant information
  data <- object$data
  n_data <- nrow(data)
  predictors <- c(object$x, object$covariates)
  robust <- object$robust
  family <- object$family
  # check which explanatory variables are discrete
  is_discrete <- sapply(data[, predictors, drop = FALSE], function(x) {
    length(unique(x)) <= num_discrete
  })
  discrete <- predictors[is_discrete]
  continuous <- predictors[!is_discrete]

  # for continuous explanatory variables, perform separate regressions of each
  # variable on just a constant term (intercept-only model) and draw from the
  # fitted model
  if (length(continuous) == 0L) X_continuous <- NULL
  else {
    null_matrix <- matrix(NA_real_, nrow = n_data, ncol = 0)
    if (robust == "MM") {
      # MM-estimator for robust regression
      fits <- lapply(continuous, function(y) {
        lmrob_fit(null_matrix, data[, y], control = object$control)
      })
    } else if (robust == "median") {
      # LAD-estimator for median regression
      fits <- lapply(continuous, function(y) {
        rq_fit(null_matrix, data[, y], tau = 0.5)
      })
    } else if (family == "gaussian") {
      # OLS estimator
      fits <- lapply(continuous, function(y) lm_fit(null_matrix, data[, y]))
    } else if (family == "select") {
      # select among normal, skew-normal, t and skew-t errors
      fits <- lapply(continuous, function(y) lmselect_fit(null_matrix, data[, y]))
    } else {
      # obtain parameters as required for package 'sn'
      selm_args <- get_selm_args(family)
      # perform regression with skew-elliptical errors
      fits <- lapply(continuous, function(y) {
        selm_fit(null_matrix, data[, y], family = selm_args$family,
                 fixed.param = selm_args$fixed.param)
      })
    }
    # add names to list of regression fits
    names(fits) <- continuous
    # draw explanatory variables by drawing error terms and adding intercept
    X_continuous <- sapply(fits, function(fit) {
      unname(coef(fit)) + sim_errors(fit, n = n)
    })
  }

  # draw discrete explanatory variables from multinomial distributions
  if (length(discrete) == 0L) X_discrete <- NULL
  else {
    X_discrete <- sapply(discrete, function(x) {
      # Note: tabulate(as.factor(.)) can give a different number of elements from
      # unique(.) as the latter seems to be more sensitive to numerical precision
      # -----
      # # obtain observed frequencies
      # frequencies <- tabulate(as.factor(data[, x]))
      # -----
      # obtain unique values
      unique <- sort(unique(data[, x]))
      # draw from unique values based on observed frequencies
      if (length(unique) == 1L) rep.int(unique, n)
      else {
        frequencies <- sapply(unique, function(value) sum(data[, x] == value))
        sample(unique, size = n, replace = TRUE, prob = frequencies)
      }
    })
  }

  # return explanatory variables
  X <- cbind(X_continuous, X_discrete)
  X[, predictors, drop = FALSE]

}

sim_explanatory.cov_fit_mediation <- function(object, n, ...) {
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
  # draw errors for mediators
  e_M <- sim_errors(object$fit_mx, n = n)
  # draw errors for dependent variable
  e_Y <- sim_errors(object$fit_ymx, n = n)
  # return list of errors
  list(M = e_M, Y = e_Y)
}

# -----
# FIXME: I think the code below only works when the function is exported
# -----
# sim_errors.list <- function(object, ...) lapply(object, sim_errors, ...)
# -----
# dirty hack:
sim_errors.list <- function(object, n) {
  lapply(object, function(x, n) {
    # FIXME: this throws error when MM-estimator doesn't converge
    #        The reason is that the (unconverged) object has class "lmrob.S"
    #        instead of class "lmrob".
    if (inherits(x, "lmrob")) sim_errors.lmrob(x, n)
    else if (inherits(x, "rq")) sim_errors.rq(x, n)
    else if (inherits(x, "lm")) sim_errors.lm(x, n)
    else if (inherits(x, "lmse")) sim_errors.lmse(x, n)
    else stop("not implemented yet")
  }, n = n)
}
# -----

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
