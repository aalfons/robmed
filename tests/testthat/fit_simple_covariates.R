context("model fit: simple mediation, covariates")


## load package
library("robmed", quietly = TRUE)

## control parameters
n <- 250            # number of observations
a <- c <- 0.2       # true effects
b <- 0              # true effect
seed <- 20190201    # seed for the random number generator

## set seed for reproducibility
set.seed(seed)

## generate data
X <- rnorm(n)
M <- a * X + rnorm(n)
Y <- b * M + c * X + rnorm(n)
C1 <- rnorm(n)
C2 <- rnorm(n)
test_data <- data.frame(X, Y, M, C1, C2)

## control parameters for methods
max_iterations <- 500  # for MM-regression estimator

## fit mediation models
fit_list <- list(
  robust = {
    set.seed(seed)
    fit_mediation(test_data, x = "X", y = "Y", m = "M",
                  covariates = c("C1", "C2"), method = "regression",
                  robust = TRUE, max_iterations = max_iterations)
  },
  median = {
    fit_mediation(test_data, x = "X", y = "Y", m = "M",
                  covariates = c("C1", "C2"), method = "regression",
                  robust = "median")
  },
  OLS = {
    fit_mediation(test_data, x = "X", y = "Y", m = "M",
                  covariates = c("C1", "C2"), method = "regression",
                  robust = FALSE, family = "gaussian")
  },
  student = {
    fit_mediation(test_data, x = "X", y = "Y", m = "M",
                  covariates = c("C1", "C2"), method = "regression",
                  robust = FALSE, family = "student")
  },
  select = {
    fit_mediation(test_data, x = "X", y = "Y", m = "M",
                  covariates = c("C1", "C2"), method = "regression",
                  robust = FALSE, family = "select")
  }
)

## compute summaries
summary_list <- lapply(fit_list, summary)

## relevant information
classes <- c(robust = "lmrob", median = "rq", OLS = "lm",
             student = "lmse", select = "lm")

## correct values
coef_names <- c("a", "b", "Total", "Direct", "Indirect")


## common tests for all model fits

# loop over methods
methods <- names(fit_list)
for (method in methods) {

  # extract information for current method
  fit <- fit_list[[method]]
  summary <- summary_list[[method]]

  ## correct values
  class <- classes[method]
  family <- if (method %in% c("student", "select")) method else "gaussian"


  ## run tests

  test_that("output has correct structure", {

    # regression fit
    expect_s3_class(fit, "reg_fit_mediation")
    expect_s3_class(fit, "fit_mediation")
    # regressions m ~ x and y ~ m + x
    expect_s3_class(fit$fit_mx, class)
    expect_s3_class(fit$fit_ymx, class)
    # regression y ~ x is only computed for nonrobust methods
    if (method %in% c("robust", "median")) {
      expect_null(fit$fit_yx)
    } else {
      expect_s3_class(fit$fit_yx, class)
    }

  })

  test_that("arguments are correctly passed", {

    # variable names
    expect_identical(fit$x, "X")
    expect_identical(fit$y, "Y")
    expect_identical(fit$m, "M")
    expect_identical(fit$covariates, c("C1", "C2"))
    # robust or nonrobust fit
    if (method == "robust") {
      expect_identical(fit$robust, "MM")
      expect_equal(fit$control, reg_control(max_iterations = max_iterations))
    } else if (method == "median") {
      expect_identical(fit$robust, "median")
      expect_null(fit$control)
    } else {
      expect_false(fit$robust)
      expect_null(fit$control)
    }
    # assumed error distribution
    expect_identical(fit$family, family)
    # mediation model
    expect_identical(fit$model, "simple")
    # no contrasts
    expect_false(fit$contrast)

  })

  test_that("dimensions are correct", {

    # effects are scalars
    expect_length(fit$a, 1L)
    expect_length(fit$b, 1L)
    expect_length(fit$direct, 1L)
    expect_length(fit$total, 1L)
    expect_length(fit$indirect, 1L)
    # individual regressions
    expect_length(coef(fit$fit_mx), 4L)
    expect_length(coef(fit$fit_ymx), 5L)
    if (!(method %in% c("robust", "median"))) {
      expect_length(coef(fit$fit_yx), 4L)
    }
    # dimensions of data
    expect_identical(dim(fit$data), c(as.integer(n), 5L))

  })

  test_that("values of coefficients are correct", {

    expect_equivalent(fit$a, coef(fit$fit_mx)["X"])
    expect_equivalent(fit$b, coef(fit$fit_ymx)["M"])
    expect_equivalent(fit$direct, coef(fit$fit_ymx)["X"])
    # total effect
    if (method %in% c("robust", "median")) {
      expect_equivalent(fit$total, fit$indirect + fit$direct)
    } else if (method == "OLS") {
      expect_equivalent(fit$total, coef(fit$fit_yx)["X"])
      expect_equivalent(fit$total, fit$indirect + fit$direct)
    } else {
      expect_equivalent(fit$total, coef(fit$fit_yx)["X"])
    }

  })

  test_that("output of coef() method has correct attributes", {

    coefficients <- coef(fit)
    expect_length(coefficients, 5L)
    expect_named(coefficients, coef_names)

  })

  test_that("coef() method returns correct values of coefficients", {

    expect_equivalent(coef(fit, parm = "a"), fit$a)
    expect_equivalent(coef(fit, parm = "b"), fit$b)
    expect_equivalent(coef(fit, parm = "Total"), fit$total)
    expect_equivalent(coef(fit, parm = "Direct"), fit$direct)
    expect_equivalent(coef(fit, parm = "Indirect"), fit$indirect)

  })

  test_that("summary returns original object", {
    expect_identical(fit, summary)
  })

}


## covariance fits only implemented for simple mediation without covariates

# argument values
cov_args <- list(winsorized = list(robust = TRUE),
                 ML = list(robust = FALSE))

# loop over methods
cov_methods <- names(cov_args)
for (method in cov_methods) {

  # argument values
  robust <- cov_args[[method]]$robust


  # run tests

  test_that("covariates not implemented", {

    # run regression fit
    set.seed(seed)
    reg_fit <- fit_mediation(test_data, x = "X", y = "Y", m = "M",
                             covariates = c("C1", "C2"),
                             method = "regression",
                             robust = robust)

    # try to run covariance fit with covariates (should give warning)
    set.seed(seed)
    expect_warning(
      cov_fit <- fit_mediation(test_data, x = "X", y = "Y", m = "M",
                               covariates = c("C1", "C2"),
                               method = "covariance",
                               robust = robust)
    )

    # these should be the same
    expect_equal(cov_fit, reg_fit)

  })

}
