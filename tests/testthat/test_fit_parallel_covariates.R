context("model fit: parallel multiple mediators, covariates")


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
M1 <- a * X + rnorm(n)
M2 <- rnorm(n)
Y <- b * M1 + c * X + rnorm(n)
C1 <- rnorm(n)
C2 <- rnorm(n)
test_data <- data.frame(X, Y, M1, M2, C1, C2)

## control parameters for methods
max_iterations <- 500  # for MM-regression estimator

## fit mediation models
fit_list <- list(
  robust = {
    set.seed(seed)
    fit_mediation(test_data, x = "X", y = "Y", m = c("M1", "M2"),
                  covariates = c("C1", "C2"), model = "parallel",
                  method = "regression", robust = TRUE,
                  max_iterations = max_iterations)
  },
  median = {
    fit_mediation(test_data, x = "X", y = "Y", m = c("M1", "M2"),
                  covariates = c("C1", "C2"), model = "parallel",
                  method = "regression", robust = "median")
  },
  OLS = {
    fit_mediation(test_data, x = "X", y = "Y", m = c("M1", "M2"),
                  covariates = c("C1", "C2"), model = "parallel",
                  method = "regression", robust = FALSE,
                  family = "gaussian")
  },
  student = {
    fit_mediation(test_data, x = "X", y = "Y", m = c("M1", "M2"),
                  covariates = c("C1", "C2"), model = "parallel",
                  method = "regression", robust = FALSE,
                  family = "student")
  },
  select = {
    fit_mediation(test_data, x = "X", y = "Y", m = c("M1", "M2"),
                  covariates = c("C1", "C2"), model = "parallel",
                  method = "regression", robust = FALSE,
                  family = "select")
  }
)

## compute summaries
summary_list <- lapply(fit_list, summary)

# relevant information
classes <- c(robust = "lmrob", median = "rq", OLS = "lm",
             student = "lmse", select = "lm")

## correct values
coef_names <- c("a_M1", "a_M2", "b_M1", "b_M2", "Total", "Direct",
                "Indirect_Total", "Indirect_M1", "Indirect_M2")


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
    expect_type(fit$fit_mx, "list")
    expect_length(fit$fit_mx, 2L)
    expect_named(fit$fit_mx, c("M1", "M2"))
    expect_s3_class(fit$fit_mx$M1, class)
    expect_s3_class(fit$fit_mx$M2, class)
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
    expect_identical(fit$m, c("M1", "M2"))
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
    expect_identical(fit$model, "parallel")
    # no contrasts
    expect_false(fit$contrast)

  })

  test_that("dimensions are correct", {

    # effect estimates
    expect_length(fit$a, 2L)
    expect_named(fit$a, c("M1", "M2"))
    expect_length(fit$b, 2L)
    expect_named(fit$b, c("M1", "M2"))
    expect_length(fit$direct, 1L)
    expect_length(fit$total, 1L)
    expect_length(fit$indirect, 3L)
    expect_named(fit$indirect, c("Total", "M1", "M2"))
    # individual regressions
    expect_length(coef(fit$fit_mx$M1), 4L)
    expect_length(coef(fit$fit_mx$M2), 4L)
    expect_length(coef(fit$fit_ymx), 6L)
    if (!(method %in% c("robust", "median"))) {
      expect_length(coef(fit$fit_yx), 4L)
    }
    # dimensions of data
    expect_identical(dim(fit$data), c(as.integer(n), 6L))

  })

  test_that("values of coefficients are correct", {

    # extract correct values
    a <- c(M1 = unname(coef(fit$fit_mx$M1)["X"]),
           M2 = unname(coef(fit$fit_mx$M2)["X"]))
    b <- coef(fit$fit_ymx)[c("M1", "M2")]
    direct <- coef(fit$fit_ymx)["X"]
    indirect <- a * b
    sum_indirect <- sum(indirect)
    indirect <- c(Total = sum_indirect, indirect)
    # compare with stored values
    expect_equivalent(fit$a, a)
    expect_equivalent(fit$b, b)
    expect_null(fit$d)
    expect_equivalent(fit$direct, direct)
    expect_equivalent(fit$indirect, indirect)
    # total effect
    if (method %in% c("robust", "median")) {
      expect_equivalent(fit$total, sum_indirect + direct)
    } else if (method == "OLS") {
      expect_equivalent(fit$total, coef(fit$fit_yx)["X"])
      expect_equivalent(fit$total, sum_indirect + direct)
    } else {
      expect_equivalent(fit$total, coef(fit$fit_yx)["X"])
    }

  })

  test_that("output of coef() method has correct attributes", {

    coefficients <- coef(fit)
    expect_length(coefficients, 9L)
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
    reg_fit <- fit_mediation(test_data, x = "X", y = "Y", m = c("M1", "M2"),
                             covariates = c("C1", "C2"), model = "parallel",
                             method = "regression", robust = robust)

    # try to run covariance fit (should give warning)
    set.seed(seed)
    expect_warning(
      cov_fit <- fit_mediation(test_data, x = "X", y = "Y", m = c("M1", "M2"),
                               covariates = c("C1", "C2"), model = "parallel",
                               method = "covariance", robust = robust)
    )

    # these should be the same
    expect_equal(cov_fit, reg_fit)

  })

}
