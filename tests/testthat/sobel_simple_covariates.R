context("Sobel test: simple mediation, covariates")


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
ctrl <- reg_control(max_iterations = 500)  # for MM-regression estimator

## perform Sobel tests
sobel_list <- list(
  robust = {
    set.seed(seed)
    test_mediation(test_data, x = "X", y = "Y", m = "M",
                   covariates = c("C1", "C2"), test = "sobel",
                   method = "regression", robust = TRUE, control = ctrl)
  },
  median = {
    test_mediation(test_data, x = "X", y = "Y", m = "M",
                   covariates = c("C1", "C2"), test = "sobel",
                   method = "regression", robust = "median")
  },
  OLS = {
    test_mediation(test_data, x = "X", y = "Y", m = "M",
                   covariates = c("C1", "C2"), test = "sobel",
                   method = "regression", robust = FALSE, family = "gaussian")
  },
  student = {
    test_mediation(test_data, x = "X", y = "Y", m = "M",
                   covariates = c("C1", "C2"), test = "sobel",
                   method = "regression", robust = FALSE, family = "student")
  },
  select = {
    test_mediation(test_data, x = "X", y = "Y", m = "M",
                   covariates = c("C1", "C2"), test = "sobel",
                   method = "regression", robust = FALSE, family = "select")
  }
)

## compute summaries
summary_list <- lapply(sobel_list, summary)

## relevant information
classes <- c(robust = "summary_lmrob", median = "summary_rq",
             OLS = "summary_lm", student = "summary_lmse",
             select = "summary_lm")
level <- 0.9

## correct values
coef_names <- c("a", "b", "Total", "Direct", "Indirect")


## common tests for all model fits

# loop over methods
methods <- names(sobel_list)
for (method in methods) {

  # extract information for current method
  sobel <- sobel_list[[method]]
  summary <- summary_list[[method]]

  ## correct values
  class <- classes[method]
  family <- if (method %in% c("student", "select")) method else "gaussian"
  if (method == "student") {
    mx_names <- c("(Intercept.DP)", "X", "C1", "C2")
    ymx_names <- c("(Intercept.DP)", "M", "X", "C1", "C2")
  } else {
    mx_names <- c("(Intercept)", "X", "C1", "C2")
    ymx_names <- c("(Intercept)", "M", "X", "C1", "C2")
  }


  ## run tests

  test_that("output has correct structure", {

    # Sobel test
    expect_s3_class(sobel, "sobel_test_mediation")
    expect_s3_class(sobel, "test_mediation")
    # regression fit
    expect_s3_class(sobel$fit, "reg_fit_mediation")
    expect_s3_class(sobel$fit, "fit_mediation")
    # standard error
    expect_is(sobel$se, "numeric")
    # test statistic
    expect_is(sobel$statistic, "numeric")
    # p-value
    expect_is(sobel$p_value, "numeric")

  })

  test_that("arguments are correctly passed", {

    # alternative hypothesis
    expect_identical(sobel$alternative, "twosided")
    # order of normal approximation
    expect_identical(sobel$order, "first")
    # variable names
    expect_identical(sobel$fit$x, "X")
    expect_identical(sobel$fit$y, "Y")
    expect_identical(sobel$fit$m, "M")
    expect_identical(sobel$fit$covariates, c("C1", "C2"))
    # robust or nonrobust fit and test
    if (method == "robust") {
      expect_identical(sobel$fit$robust, "MM")
      expect_equal(sobel$fit$control, ctrl)
    } else if (method == "median") {
      expect_identical(sobel$fit$robust, "median")
      expect_null(sobel$fit$control)
    } else {
      expect_false(sobel$fit$robust)
      expect_null(sobel$fit$control)
    }
    # assumed error distribution
    expect_identical(sobel$fit$family, family)
    # mediation model
    expect_identical(sobel$fit$model, "simple")
    # no contrasts
    expect_false(sobel$fit$contrast)

  })

  test_that("dimensions are correct", {

    # standard error
    expect_length(sobel$se, 1L)
    # test statistic
    expect_length(sobel$statistic, 1L)
    # p-value
    expect_length(sobel$p_value, 1L)

  })

  test_that("coef() method returns correct values of coefficients", {

    expect_identical(coef(sobel), coef(sobel$fit))

  })

  test_that("output of confint() method has correct attributes", {

    ci_sobel <- confint(sobel, level = level)
    expect_equal(dim(ci_sobel), c(5L, 2L))
    expect_equal(rownames(ci_sobel), coef_names)
    expect_equal(colnames(ci_sobel), c("5 %", "95 %"))

  })

  test_that("confint() method returns correct values of confidence intervals", {

    ci_90 <- confint(sobel, parm = "Indirect", level = level)
    if (method == "student") {
      # Sobel test not meaningful as higher moments may not exist
      expect_true(is.na(ci_90["Indirect", 1]))
      expect_true(is.na(ci_90["Indirect", 2]))
    } else {
      # default CI should be wider
      ci_default <- confint(sobel, parm = "Indirect")  # should be 95%
      expect_lt(ci_default["Indirect", 1], ci_90["Indirect", 1])
      expect_gt(ci_default["Indirect", 2], ci_90["Indirect", 2])
    }

  })

  test_that("summary has correct structure", {

    # summary
    expect_s3_class(summary, "summary_test_mediation")
    # original output of test for indirect effect
    expect_identical(summary$object, sobel)
    # summary of the model fit
    expect_s3_class(summary$summary, "summary_reg_fit_mediation")
    expect_s3_class(summary$summary, "summary_fit_mediation")
    # no diagnostic plot
    expect_null(summary$plot)

  })

  # model summaries of regressions m ~ x and y ~ m + x
  fit_names <- c("fit_mx", "fit_ymx")
  # correct degrees of freedom in the numerator for F-test
  df1 <- c(3, 4)

  # loop over model summaries
  for (i in 1:length(fit_names)) {

    # extract model summary
    fit_name <- fit_names[i]
    model_summary <- summary$summary[[fit_name]]

    # perform tests
    test_that("model summary has correct structure", {

      # model summary
      expect_s3_class(model_summary, class)
      # regression standard error, R-squared, and F-test
      if (method %in% c("robust", "OLS", "select")) {
        # regression standard error
        expect_type(model_summary$s, "list")
        expect_named(model_summary$s, c("value", "df"))
        # R-squared
        expect_type(model_summary$R2, "list")
        expect_named(model_summary$R2, c("R2", "adj_R2"))
        # F-test
        expect_type(model_summary$F_test, "list")
        expect_named(model_summary$F_test, c("statistic", "df", "p_value"))
        # degrees of freedom of F-test
        df_test <- model_summary$F_test$df
        if (method == "robust") {
          expect_identical(df_test[1], df1[i])
          expect_identical(df_test[2], Inf)
        } else {
          expect_identical(df_test[1], as.integer(df1[i]))
          expect_identical(df_test[2], model_summary$s$df)
        }
      } else {
        # regression standard error
        expect_null(model_summary$s)
        # R-squared
        expect_null(model_summary$R2)
        # F-test
        expect_null(model_summary$F_test)
      }
      # additional information for robust Sobel test
      if (method == "robust") {
        fit <- sobel$fit[[fit_name]]
        # information on convergence
        expect_type(model_summary$algorithm, "list")
        expect_named(model_summary$algorithm, c("converged", "method"))
        expect_identical(model_summary$algorithm$converged,
                         fit$converged)
        expect_identical(model_summary$algorithm$method,
                         fit$control$method)
        # information on outliers
        expect_type(model_summary$outliers, "list")
        expect_named(model_summary$outliers,
                     c("indices", "weights", "threshold"))
        expect_type(model_summary$outliers$indices, "integer")
        expect_identical(model_summary$outliers$weights,
                         weights(fit, type = "robustness"))
        expect_identical(model_summary$outliers$threshold,
                         summary(fit)$control$eps.outlier)
      }

    })

  }

  test_that("attributes are correctly passed through summary", {

    # robustness
    if (method == "robust") {
      expect_identical(summary$summary$robust, "MM")
    } else if (method == "median") {
      expect_identical(summary$summary$robust, "median")
    } else {
      expect_false(summary$summary$robust)
    }
    # number of observations
    expect_identical(summary$summary$n, as.integer(n))
    # variable names
    expect_identical(summary$summary$x, "X")
    expect_identical(summary$summary$y, "Y")
    expect_identical(summary$summary$m, "M")
    expect_identical(summary$summary$covariates, c("C1", "C2"))

  })

  test_that("effect summaries have correct names", {

    # regression m ~ x
    expect_identical(dim(summary$summary$fit_mx$coefficients),
                     c(4L, 4L))
    expect_identical(rownames(summary$summary$fit_mx$coefficients),
                     mx_names)
    expect_identical(colnames(summary$summary$fit_mx$coefficients)[1],
                     "Estimate")
    # regression y ~ m + x
    expect_identical(dim(summary$summary$fit_ymx$coefficients),
                     c(5L, 4L))
    expect_identical(rownames(summary$summary$fit_ymx$coefficient),
                     ymx_names)
    expect_identical(colnames(summary$summary$fit_ymx$coefficient)[1],
                     "Estimate")
    # total effect
    expect_identical(dim(summary$summary$total), c(1L, 4L))
    expect_identical(rownames(summary$summary$total), "X")
    expect_identical(colnames(summary$summary$total)[1], "Estimate")
    # direct effect
    expect_identical(dim(summary$summary$direct), c(1L, 4L))
    expect_identical(rownames(summary$summary$direct), "X")
    expect_identical(colnames(summary$summary$direct)[1], "Estimate")

  })

  test_that("effect summaries contain correct coefficient values", {

    expect_identical(summary$summary$fit_mx$coefficients["X", "Estimate"],
                     sobel$fit$a)
    expect_identical(summary$summary$fit_ymx$coefficients["M", "Estimate"],
                     sobel$fit$b)
    expect_identical(summary$summary$total["X", "Estimate"], sobel$fit$total)
    expect_identical(summary$summary$direct["X", "Estimate"], sobel$fit$direct)

  })

  test_that("output of p_value() method has correct attributes", {

    p_data <- p_value(sobel)
    expect_length(p_data, 5L)
    expect_named(p_data, coef_names)
    expect_equivalent(p_data["Indirect"], sobel$p_value)

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

    # use regression fit
    set.seed(seed)
    reg_sobel <- test_mediation(test_data, x = "X", y = "Y", m = "M",
                                covariates = c("C1", "C2"), test = "sobel",
                                method = "regression", robust = robust)

    # try to use covariance fit with covariates (should give warning)
    set.seed(seed)
    expect_warning(
      cov_sobel <- test_mediation(test_data, x = "X", y = "Y", m = "M",
                                  covariates = c("C1", "C2"), test = "sobel",
                                  method = "covariance", robust = robust)
    )

    # these should be the same
    expect_equal(cov_sobel, reg_sobel)

  })

}
