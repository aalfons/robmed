context("model fit: multiple independent variables to be mediated, no covariates")


## load package
library("robmed", quietly = TRUE)

## control parameters
n <- 250               # number of observations
a <- c <- c(0.2, 0.2)  # true effects
b <- 0                 # true effect
seed <- 20150601       # seed for the random number generator

## set seed for reproducibility
set.seed(seed)

## generate data
X1 <- rnorm(n)
X2 <- rnorm(n)
M <- a[1] * X1 + a[2] * X2 + rnorm(n)
Y <- b * M + c[1] * X1 + c[2] * X2 + rnorm(n)
test_data <- data.frame(X1, X2, Y, M)

## control parameters for methods
x <- c("X1", "X2")           # independent variables
y <- "Y"                     # dependent variable
m <- "M"                     # mediator variable
covariates <- character()    # control variables
efficiency <- 0.95           # for MM-regression
algorithm <- "fn"            # for median regression

## fit mediation models
fit_list <- list(
  robust = {
    set.seed(seed)
    fit_mediation(test_data, x = x, y = y, m = m, covariates = covariates,
                  method = "regression", robust = TRUE, efficiency = efficiency)
  },
  median = {
    fit_mediation(test_data, x = x, y = y, m = m, covariates = covariates,
                  method = "regression", robust = "median",
                  algorithm = algorithm)
  },
  skewnormal = {
    fit_mediation(test_data, x = x, y = y, m = m, covariates = covariates,
                  method = "regression", robust = FALSE, family = "skewnormal")
  },
  select = {
    fit_mediation(test_data, x = x, y = y, m = m, covariates = covariates,
                  method = "regression", robust = FALSE, family = "select")
  },
  OLS = {
    fit_mediation(test_data, x = x, y = y, m = m, covariates = covariates,
                  method = "regression", robust = FALSE, family = "gaussian")
  }
)

## compute summaries
summary_list <- lapply(fit_list, summary)

## correct values
effect_names <- c("a_X1", "a_X2", "b", "Total_X1", "Total_X2", "Direct_X1",
                  "Direct_X2", "Indirect_X1", "Indirect_X2")
classes <- c(robust = "lmrob", median = "rq", skewnormal = "lmse",
             select = "lm", OLS = "lm")


## run tests

# loop over methods
methods <- names(fit_list)
for (method in methods) {

  # extract information for current method
  fit <- fit_list[[method]]
  summary <- summary_list[[method]]

  # correct values
  class <- classes[method]
  family <- if (method %in% c("skewnormal", "select")) method else "gaussian"


  # run tests

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
    expect_identical(fit$x, x)
    expect_identical(fit$y, y)
    expect_identical(fit$m, m)
    expect_identical(fit$covariates, covariates)
    # robust or nonrobust fit
    if (method == "robust") {
      expect_identical(fit$robust, "MM")
      expect_equal(fit$control, MM_reg_control(efficiency = efficiency))
    } else if (method == "median") {
      expect_identical(fit$robust, "median")
      expect_equal(fit$control, median_reg_control(algorithm = algorithm))
    } else {
      expect_false(fit$robust)
      expect_null(fit$control)
    }
    # assumed error distribution
    expect_identical(fit$family, family)
    # mediation model
    expect_identical(fit$model, "multiple")
    # no contrasts
    expect_false(fit$contrast)

  })

  test_that("dimensions are correct", {

    # effect estimates
    expect_length(fit$a, 2L)
    expect_named(fit$a, x)
    expect_length(fit$b, 1L)
    expect_named(fit$b, NULL)
    expect_length(fit$total, 2L)
    expect_named(fit$total, x)
    expect_length(fit$direct, 2L)
    expect_named(fit$direct, x)
    expect_length(fit$indirect, 2L)
    expect_named(fit$indirect, x)
    # individual regressions
    expect_length(coef(fit$fit_mx), 3L)
    expect_length(coef(fit$fit_ymx), 4L)
    if (!(method %in% c("robust", "median"))) {
      expect_length(coef(fit$fit_yx), 3L)
    }
    # dimensions of data
    expect_identical(dim(fit$data), c(as.integer(n), 4L))

  })

  test_that("values of coefficients are correct", {

    # extract correct values
    a <- coef(fit$fit_mx)[x]
    b <- unname(coef(fit$fit_ymx)[m])
    direct <- coef(fit$fit_ymx)[x]
    indirect <- a * b
    # compare with stored values
    expect_equivalent(fit$a, a)
    expect_equivalent(fit$b, b)
    expect_null(fit[["d"]])
    expect_equivalent(fit$direct, direct)
    expect_equivalent(fit$indirect, indirect)
    # total effect
    if (method %in% c("robust", "median")) {
      expect_equivalent(fit$total, indirect + direct)
    } else if (method == "OLS") {
      expect_equivalent(fit$total, coef(fit$fit_yx)[x])
      expect_equivalent(fit$total, indirect + direct)
    } else {
      expect_equivalent(fit$total, coef(fit$fit_yx)[x])
    }

  })

  test_that("output of coef() method has correct attributes", {

    coefficients <- coef(fit)
    expect_length(coefficients, 9L)
    expect_named(coefficients, effect_names)

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


  # tests for weight_plot(), which is only implemented for ROBMED
  if (method == "robust") {

    # loop over settings for weight plot
    for (setting in c("default", m, y)) {

      # obtain setup object for weight plot
      if (setting == "default") weight <- setup_weight_plot(fit)
      else weight <- setup_weight_plot(fit, outcome = setting)

      # run tests
      test_that("object returned by setup_weight_plot() has correct structure", {

        # check data frame for weight percentages to be plotted
        expect_s3_class(weight$data, "data.frame")

        # additional checks
        names <- c("Outcome", "Tail", "Weights", "Threshold", "Percentage")
        if (setting == "default") {
          # check dimensions and column names
          expect_identical(ncol(weight$data), 5L)
          expect_named(weight$data, names)
          # check if variables are passed correctly
          expect_identical(weight$outcome, c(m, y))
        } else {
          # check dimensions and column names
          expect_identical(ncol(weight$data), 4L)
          expect_named(weight$data, names[-1])
          # check if variables are passed correctly
          expect_identical(weight$outcome, setting)
        }

      })

    }

  } else {

    # run tests
    test_that("setup_weight_plot() gives error", {

      # not implemented
      expect_error(setup_weight_plot(fit))

    })

  }


  # tests for setup_ellipse_plot(), which is only implemented for some methods
  if (method %in% c("robust", "OLS", "winsorized", "ML")) {

    # loop over settings for ellipse plot
    for (setting in c("mx", "ym", "partial")) {

      # obtain setup object for ellipse plot
      if (setting == "mx") ellipse <- setup_ellipse_plot(fit)
      else if (setting == "ym") {
        ellipse <- setup_ellipse_plot(fit, horizontal = m, vertical = y,
                                      partial = FALSE)
      } else {
        ellipse <- setup_ellipse_plot(fit, horizontal = m, vertical = y,
                                      partial = TRUE)
      }

      # run tests
      test_that("object returned by setup_ellipse_plot() has correct structure", {

        # check data frame for data to be plotted
        expect_s3_class(ellipse$data, "data.frame")
        # check dimensions and column names
        if (method %in% c("robust", "winsorized")) {
          expect_identical(dim(ellipse$data), c(as.integer(n), 3L))
          expect_named(ellipse$data, c("x", "y", "Weight"))
        } else {
          expect_identical(dim(ellipse$data), c(as.integer(n), 2L))
          expect_named(ellipse$data, c("x", "y"))
        }

        # check data frame for ellipse
        expect_s3_class(ellipse$ellipse, "data.frame")
        # check dimensions
        expect_identical(ncol(ellipse$ellipse), 2L)
        expect_gt(nrow(ellipse$ellipse), 0L)
        # check column names
        expect_named(ellipse$ellipse, c("x", "y"))

        # check line to be plotted
        if (setting == "partial") {
          # check data frame for line representing the coefficient
          expect_s3_class(ellipse$line, "data.frame")
          # check dimensions
          expect_identical(dim(ellipse$line), c(1L, 2L))
          # check column names
          expect_named(ellipse$line, c("intercept", "slope"))
          # check if intercept is 0 for partial residuals
          expect_identical(ellipse$line$intercept, 0)
        } else {
          # no line to be plotted
          expect_null(ellipse$line)
        }

        # check if variables are passed correctly
        if (setting == "mx") {
          expect_identical(ellipse$horizontal, x[1])
          expect_identical(ellipse$vertical, m)
        } else {
          expect_identical(ellipse$horizontal, m)
          expect_identical(ellipse$vertical, y)
        }

        # check logical for partial residuals on the vertical axis
        if (setting == "partial") {
          expect_true(ellipse$partial)
        } else {
          expect_false(ellipse$partial)
        }

        # check logical for robust method
        if (method %in% c("robust", "winsorized")) {
          expect_true(ellipse$robust)
        } else {
          expect_false(ellipse$robust)
        }

        # check logical for multiple methods
        expect_false(ellipse$have_methods)

      })

    }

  } else {

    # run tests
    test_that("setup_ellipse_plot() gives error", {

      # not implemented
      expect_error(setup_ellipse_plot(fit))

    })

  }


  # tests for ci_plot() and density_plot()
  test_that("setup_ci_plot()  and setup_density_plot() give error", {

    # not meaningful
    expect_error(setup_ci_plot(fit))
    expect_error(setup_density_plot(fit))

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
    reg_fit <- fit_mediation(test_data, x = x, y = y, m = m,
                             covariates = covariates,
                             method = "regression",
                             robust = robust)

    # try to run covariance fit (should give warning)
    set.seed(seed)
    expect_warning(
      cov_fit <- fit_mediation(test_data, x = x, y = y, m = m,
                               covariates = covariates,
                               method = "covariance",
                               robust = robust)
    )

    # hack in case ATLAS is used as BLAS: signs of eigenvectors may be arbitrary
    if (robust) {
      attr(reg_fit$fit_mx$cov, "eigen") <- NULL
      attr(reg_fit$fit_ymx$cov, "eigen") <- NULL
      attr(cov_fit$fit_mx$cov, "eigen") <- NULL
      attr(cov_fit$fit_ymx$cov, "eigen") <- NULL
    }

    # these should be the same
    expect_equal(cov_fit, reg_fit)

  })

}
