context("regression fit: simple mediation, no covariates")


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

## fit mediation models
fit_list <- list(
  robust = {
    set.seed(seed)
    fit_mediation(test_data, x = "X", y = "Y", m = "M",
                  method = "regression", robust = TRUE,
                  efficiency = 0.95)
  },
  median = {
    set.seed(seed)
    fit_mediation(test_data, x = "X", y = "Y", m = "M",
                  method = "regression", robust = "median")
  },
  standard = {
    set.seed(seed)
    fit_mediation(test_data, x = "X", y = "Y", m = "M",
                  method = "regression", robust = FALSE,
                  family = "gaussian")
  },
  student = {
    set.seed(seed)
    fit_mediation(test_data, x = "X", y = "Y", m = "M",
                  method = "regression", robust = FALSE,
                  family = "student")
  },
  select = {
    set.seed(seed)
    fit_mediation(test_data, x = "X", y = "Y", m = "M",
                  method = "regression", robust = FALSE,
                  family = "select")
  }
)

## compute summaries
summary_list <- lapply(fit_list, summary)

## other relevant information
methods <- names(fit_list)
classes <- c(robust = "lmrob", median = "rq", standard = "lm",
             student = "lmse", select = "lm")
families <- c(robust = "gaussian", median = "gaussian", standard = "gaussian",
              student = "student", select = "select")

## stuff needed to check correctness
coef_names <- c("a", "b", "Total", "Direct", "Indirect")


## loop over methods
for (method in methods) {

  ## extract information for current method
  fit <- fit_list[[method]]
  summary <- summary_list[[method]]

  ## correct values
  class <- classes[method]
  family <- unname(families[method])


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
    expect_identical(fit$covariates, character())
    # robust or nonrobust fit
    if (method == "robust") {
      expect_identical(fit$robust, "MM")
      expect_equal(fit$control, reg_control(efficiency = 0.95))
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
    expect_length(coef(fit$fit_mx), 2L)
    expect_length(coef(fit$fit_ymx), 3L)
    if (!(method %in% c("robust", "median"))) {
      expect_length(coef(fit$fit_yx), 2L)
    }
    # dimensions of data
    expect_identical(dim(fit$data), c(as.integer(n), 3L))

  })

  test_that("values of coefficients are correct", {

    expect_equivalent(fit$a, coef(fit$fit_mx)["X"])
    expect_equivalent(fit$b, coef(fit$fit_ymx)["M"])
    expect_equivalent(fit$direct, coef(fit$fit_ymx)["X"])
    expect_equivalent(fit$indirect, fit$a * fit$b)
    # total effect
    if (method %in% c("robust", "median")) {
      expect_equivalent(fit$total, fit$indirect + fit$direct)
    } else if (method == "standard") {
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
