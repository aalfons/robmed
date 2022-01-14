context("model fit: simple mediation, no covariates")


## load package
library("robmed", quietly = TRUE)

## control parameters
n <- 250            # number of observations
a <- c <- 0.2       # true effects
b <- 0              # true effect
seed <- 20150601    # seed for the random number generator

## set seed for reproducibility
set.seed(seed)

## generate data
X <- rnorm(n)
M <- a * X + rnorm(n)
Y <- b * M + c * X + rnorm(n)
test_data <- data.frame(X, Y, M)

## control parameters for methods
efficiency <- 0.95  # for MM-regression estimator
prob <- 0.9         # for winsorized covariance matrix

## fit mediation models
fit_list <- list(
  robust = {
    set.seed(seed)
    fit_mediation(test_data, x = "X", y = "Y", m = "M",
                  method = "regression", robust = TRUE,
                  efficiency = efficiency)
  },
  median = {
    fit_mediation(test_data, x = "X", y = "Y", m = "M",
                  method = "regression", robust = "median")
  },
  OLS = {
    fit_mediation(test_data, x = "X", y = "Y", m = "M",
                  method = "regression", robust = FALSE,
                  family = "gaussian")
  },
  student = {
    fit_mediation(test_data, x = "X", y = "Y", m = "M",
                  method = "regression", robust = FALSE,
                  family = "student")
  },
  select = {
    fit_mediation(test_data, x = "X", y = "Y", m = "M",
                  method = "regression", robust = FALSE,
                  family = "select")
  },
  winsorized = {
    fit_mediation(test_data, x = "X", y = "Y", m = "M",
                  method = "covariance", robust = TRUE,
                  prob = prob)
  },
  ML = {
    fit_mediation(test_data, x = "X", y = "Y", m = "M",
                  method = "covariance", robust = FALSE)
  }
)

## compute summaries
summary_list <- lapply(fit_list, summary)

## correct values
coef_names <- c("a", "b", "Total", "Direct", "Indirect")


## common tests for all model fits

# loop over methods
methods <- names(fit_list)
for (method in methods) {

  # extract information for current method
  fit <- fit_list[[method]]
  summary <- summary_list[[method]]


  ## run tests

  test_that("output has correct structure", {

    # model fit
    expect_s3_class(fit, "fit_mediation")

  })

  test_that("arguments are correctly passed", {

    # variable names
    expect_identical(fit$x, "X")
    expect_identical(fit$y, "Y")
    expect_identical(fit$m, "M")
    expect_identical(fit$covariates, character())
    # robust or nonrobust fit
    if (method == "robust") {
      expect_identical(fit$robust, "MM")
      expect_equal(fit$control, reg_control(efficiency = efficiency))
    } else if (method == "median") {
      expect_identical(fit$robust, "median")
      expect_null(fit$control)
    } else if (method == "winsorized") {
      expect_true(fit$robust)
      expect_equal(fit$control, cov_control(prob = prob))
    } else {
      expect_false(fit$robust)
      expect_null(fit$control)
    }

  })

  test_that("dimensions are correct", {

    # effects are scalars
    expect_length(fit$a, 1L)
    expect_length(fit$b, 1L)
    expect_length(fit$direct, 1L)
    expect_length(fit$total, 1L)
    expect_length(fit$indirect, 1L)
    # dimensions of data
    expect_identical(dim(fit$data), c(as.integer(n), 3L))

  })

  test_that("values of coefficients are correct", {

    expect_null(fit$d)
    expect_equivalent(fit$indirect, fit$a * fit$b)

  })

  test_that("output of coef() method has correct attributes", {

    coefficients <- coef(fit)
    expect_length(coefficients, 5L)
    expect_named(coefficients, coef_names)

  })

  test_that("coef() method returns correct values of coefficients", {

    expect_equivalent(coef(fit, parm = "a"), fit$a)
    expect_equivalent(coef(fit, parm = "b"), fit$b)
    expect_null(coef(fit, parm = "d"))
    expect_equivalent(coef(fit, parm = "total"), fit$total)
    expect_equivalent(coef(fit, parm = "direct"), fit$direct)
    expect_equivalent(coef(fit, parm = "indirect"), fit$indirect)

  })

  test_that("summary returns original object", {
    expect_identical(fit, summary)
  })

}



## additional tests for regression fits

# relevant information
reg_methods <- c("robust", "median", "OLS", "student", "select")
classes <- c(robust = "lmrob", median = "rq", OLS = "lm",
             student = "lmse", select = "lm")

# loop over methods
for (method in reg_methods) {

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
    # regressions m ~ x and y ~ m + x
    expect_s3_class(fit$fit_mx, class)
    expect_s3_class(fit$fit_ymx, class)
    # regression y ~ x is only computed for nonrobust methods
    if (method %in% c("robust", "median")) {
      expect_null(fit$fit_yx)
    } else {
      expect_s3_class(fit$fit_yx, class)
    }
    # no covariance matrix
    expect_null(fit[["cov"]])  # can't use $ because of component 'covariates'

  })

  test_that("arguments are correctly passed", {

    # assumed error distribution
    expect_identical(fit$family, family)
    # mediation model
    expect_identical(fit$model, "simple")
    # no contrasts
    expect_false(fit$contrast)

  })

  test_that("dimensions are correct", {

    # individual regressions
    expect_length(coef(fit$fit_mx), 2L)
    expect_length(coef(fit$fit_ymx), 3L)
    if (!(method %in% c("robust", "median"))) {
      expect_length(coef(fit$fit_yx), 2L)
    }

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

}


## additional tests for covariance fits

# relevant information
cov_methods <- c("winsorized", "ML")
classes <- c(winsorized = "cov_Huber", ML = "cov_ML")

# loop over methods
for (method in cov_methods) {

  # extract information for current method
  fit <- fit_list[[method]]
  summary <- summary_list[[method]]


  ## run tests

  test_that("output has correct structure", {

    # covariance fit
    expect_s3_class(fit, "cov_fit_mediation")
    # covariance matrix
    expect_s3_class(fit$cov, classes[method])
    # no regression fits
    expect_null(fit$fit_mx)
    expect_null(fit$fit_ymx)
    expect_null(fit$fit_yx)

  })

  test_that("dimensions are correct", {

    # center and covariance matrix
    expect_length(fit$cov$center, 3L)
    expect_identical(dim(fit$cov$cov), c(3L, 3L))

  })

  test_that("values of coefficients are correct", {

    expect_equivalent(fit$total, fit$indirect + fit$direct)

  })

}
