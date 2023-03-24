context("model fit: serial multiple mediators, no covariates")


## load package
library("robmed", quietly = TRUE)

## control parameters
n <- 250            # number of observations
a <- c(0.2, 0.2)    # true effects
c <- d <- 0.2       # true effects
b <- c(0, 0.2)      # true effect
seed <- 20150601    # seed for the random number generator

## set seed for reproducibility
set.seed(seed)

## generate data
X <- rnorm(n)
M1 <- a[1] * X + rnorm(n)
M2 <- d * M1 + a[2] * X + rnorm(n)
Y <- b[1] * M1 + b[2] * M2 + c * X + rnorm(n)
test_data <- data.frame(X, Y, M1, M2)

## control parameters for methods
x <- "X"                     # independent variable
y <- "Y"                     # dependent variable
m <- c("M1", "M2")           # serial mediator variables
covariates <- character()    # control variables
efficiency <- 0.95           # for MM-regression estimator

## fit mediation models
fit_list <- list(
  robust = {
    set.seed(seed)
    fit_mediation(test_data, x = x, y = y, m = m, covariates = covariates,
                  model = "serial", method = "regression", robust = TRUE,
                  efficiency = efficiency)
  },
  median = {
    fit_mediation(test_data, x = x, y = y, m = m, covariates = covariates,
                  model = "serial", method = "regression", robust = "median")
  },
  skewnormal = {
    fit_mediation(test_data, x = x, y = y, m = m, covariates = covariates,
                  model = "serial", method = "regression", robust = FALSE,
                  family = "skewnormal")
  },
  select = {
    fit_mediation(test_data, x = x, y = y, m = m, covariates = covariates,
                  model = "serial", method = "regression", robust = FALSE,
                  family = "select")
  },
  OLS = {
    fit_mediation(test_data, x = x, y = y, m = m, covariates = covariates,
                  model = "serial", method = "regression", robust = FALSE,
                  family = "gaussian")
  }
)

## compute summaries
summary_list <- lapply(fit_list, summary)

## correct values
effect_names <- c("a_M1", "a_M2", "b_M1", "b_M2", "d", "Total", "Direct",
                  "Indirect_Total", "Indirect_M1", "Indirect_M2",
                  "Indirect_M1->M2")
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
    expect_type(fit$fit_mx, "list")
    expect_length(fit$fit_mx, 2L)
    expect_named(fit$fit_mx, m)
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
    expect_identical(fit$x, x)
    expect_identical(fit$y, y)
    expect_identical(fit$m, m)
    expect_identical(fit$covariates, covariates)
    # robust or nonrobust fit
    if (method == "robust") {
      expect_identical(fit$robust, "MM")
      expect_equal(fit$control, reg_control(efficiency = efficiency))
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
    expect_identical(fit$model, "serial")
    # no contrasts
    expect_false(fit$contrast)

  })

  test_that("dimensions are correct", {

    # effect estimates
    expect_length(fit$a, 2L)
    expect_named(fit$a, m)
    expect_length(fit$b, 2L)
    expect_named(fit$b, m)
    expect_length(fit[["d"]], 1L)
    expect_named(fit[["d"]], NULL)
    expect_length(fit$total, 1L)
    expect_named(fit$total, NULL)
    expect_length(fit$direct, 1L)
    expect_named(fit$direct, NULL)
    expect_length(fit$indirect, 4L)
    expect_named(fit$indirect, c("Total", m, paste(m, collapse = "->")))
    # individual regressions
    expect_length(coef(fit$fit_mx$M1), 2L)
    expect_length(coef(fit$fit_mx$M2), 3L)
    expect_length(coef(fit$fit_ymx), 4L)
    if (!(method %in% c("robust", "median"))) {
      expect_length(coef(fit$fit_yx), 2L)
    }
    # dimensions of data
    expect_identical(dim(fit$data), c(as.integer(n), 4L))

  })

  test_that("values of coefficients are correct", {

    # extract correct values
    a <- c(M1 = unname(coef(fit$fit_mx$M1)[x]),
           M2 = unname(coef(fit$fit_mx$M2)[x]))
    b <- coef(fit$fit_ymx)[m]
    d <- unname(coef(fit$fit_mx$M2)[m[1]])
    direct <- coef(fit$fit_ymx)[x]
    indirect <- c(a * b, "M1->M2" = unname(a[1]) * d * unname(b[2]))
    sum_indirect <- sum(indirect)
    indirect <- c(Total = sum_indirect, indirect)
    # compare with stored values
    expect_equivalent(fit$a, a)
    expect_equivalent(fit$b, b)
    expect_equivalent(fit[["d"]], d)
    expect_equivalent(fit$direct, direct)
    expect_equivalent(fit$indirect, indirect)
    # total effect
    if (method %in% c("robust", "median")) {
      expect_equivalent(fit$total, sum_indirect + direct)
    } else if (method == "OLS") {
      expect_equivalent(fit$total, coef(fit$fit_yx)[x])
      expect_equivalent(fit$total, sum_indirect + direct)
    } else {
      expect_equivalent(fit$total, coef(fit$fit_yx)[x])
    }

  })

  test_that("output of coef() method has correct attributes", {

    coefficients <- coef(fit)
    expect_length(coefficients, 11L)
    expect_named(coefficients, effect_names)

  })

  test_that("coef() method returns correct values of coefficients", {

    expect_equivalent(coef(fit, parm = "a"), fit$a)
    expect_equivalent(coef(fit, parm = "b"), fit$b)
    expect_equivalent(coef(fit, parm = "d"), fit$d)
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
        ellipse <- setup_ellipse_plot(fit, horizontal = m[1], vertical = y,
                                      partial = FALSE)
      } else {
        ellipse <- setup_ellipse_plot(fit, horizontal = m[1], vertical = y,
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
        if (setting == "ym") {
          # no line to be plotted
          expect_null(ellipse$line)
        } else {
          # check data frame for line representing the coefficient
          expect_s3_class(ellipse$line, "data.frame")
          # check dimensions
          expect_identical(dim(ellipse$line), c(1L, 2L))
          # check column names
          expect_named(ellipse$line, c("intercept", "slope"))
          if (setting == "partial") {
            # check if intercept is 0 for partial residuals
            expect_identical(ellipse$line$intercept, 0)
          }
        }

        # check if variables are passed correctly
        if (setting == "mx") {
          expect_identical(ellipse$horizontal, x)
          expect_identical(ellipse$vertical, m[1])
        } else {
          expect_identical(ellipse$horizontal, m[1])
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
                             covariates = covariates, model = "serial",
                             method = "regression", robust = robust)

    # try to run covariance fit (should give warning)
    set.seed(seed)
    expect_warning(
      cov_fit <- fit_mediation(test_data, x = x, y = y, m = m,
                               covariates = covariates, model = "serial",
                               method = "covariance", robust = robust)
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
