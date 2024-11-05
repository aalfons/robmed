context("Sobel test: simple mediation, no covariates")


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
x <- "X"                                    # independent variable
y <- "Y"                                    # dependent variable
m <- "M"                                    # mediator variable
covariates <- character()                   # control variables
robust_ctrl <- MM_reg_control(efficiency = 0.95)     # for MM-regression
median_ctrl <- median_reg_control(algorithm = "fn")  # for median regression
cov_ctrl <- cov_control(prob = 0.9)         # for winsorized covariance matrix

## perform Sobel tests
sobel_list <- list(
  robust = {
    set.seed(seed)
    test_mediation(test_data, x = x, y = y, m = m, covariates = covariates,
                   test = "sobel", method = "regression", robust = TRUE,
                   control = robust_ctrl)
  },
  median = {
    test_mediation(test_data, x = x, y = y, m = m, covariates = covariates,
                   test = "sobel", method = "regression", robust = "median",
                   control = median_ctrl)
  },
  skewnormal = {
    test_mediation(test_data, x = x, y = y, m = m, covariates = covariates,
                   test = "sobel", method = "regression", robust = FALSE,
                   family = "skewnormal")
  },
  select = {
    test_mediation(test_data, x = x, y = y, m = m, covariates = covariates,
                   test = "sobel", method = "regression", robust = FALSE,
                   family = "select")
  },
  OLS = {
    test_mediation(test_data, x = x, y = y, m = m, covariates = covariates,
                   test = "sobel", method = "regression", robust = FALSE,
                   family = "gaussian")
  },
  winsorized = {
    test_mediation(test_data, x = x, y = y, m = m, covariates = covariates,
                   test = "sobel", method = "covariance", robust = TRUE,
                   control = cov_ctrl)
  },
  ML = {
    test_mediation(test_data, x = x, y = y, m = m, covariates = covariates,
                   test = "sobel", method = "covariance", robust = FALSE)
  }
)

## compute summaries
summary_list <- lapply(sobel_list, summary)

## relevant information
level <- 0.9

## correct values
effect_names <- c("a", "b", "Total", "Direct", "Indirect")


## common tests for all model fits

# loop over methods
methods <- names(sobel_list)
for (method in methods) {

  # extract information for current method
  sobel <- sobel_list[[method]]
  summary <- summary_list[[method]]


  # run tests

  test_that("output has correct structure", {

    # Sobel test
    expect_s3_class(sobel, "sobel_test_mediation")
    expect_s3_class(sobel, "test_mediation")
    # model fit
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
    expect_identical(sobel$fit$x, x)
    expect_identical(sobel$fit$y, y)
    expect_identical(sobel$fit$m, m)
    expect_identical(sobel$fit$covariates, covariates)
    # robust or nonrobust fit and test
    if (method == "robust") {
      expect_identical(sobel$fit$robust, "MM")
      expect_equal(sobel$fit$control, robust_ctrl)
    } else if (method == "median") {
      expect_identical(sobel$fit$robust, "median")
      expect_equal(sobel$fit$control, median_ctrl)
    } else if (method == "winsorized") {
      expect_true(sobel$fit$robust)
      expect_equal(sobel$fit$control, cov_ctrl)
    } else {
      expect_false(sobel$fit$robust)
      expect_null(sobel$fit$control)
    }

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

    ci <- confint(sobel, level = level)
    expect_equal(dim(ci), c(5L, 2L))
    expect_equal(rownames(ci), effect_names)
    expect_equal(colnames(ci), c("5 %", "95 %"))

  })

  test_that("confint() method returns correct values of confidence intervals", {

    # mean of confidence intervals should equal point estimate
    if (method %in% c("robust", "median")) {
      keep <- c("a", "b", "direct", "indirect")
    } else keep <- c("a", "b", "total", "direct", "indirect")
    ci <- confint(sobel, parm = keep)
    coef <- coef(sobel, parm = keep)
    expect_equivalent(rowMeans(ci), coef)
    # default (95%) CI should be wider thn 90% CI
    ci_default <- confint(sobel, parm = "indirect")
    ci_90 <- confint(sobel, parm = "indirect", level = level)
    expect_lt(ci_default["Indirect", 1], ci_90["Indirect", 1])
    expect_gt(ci_default["Indirect", 2], ci_90["Indirect", 2])

  })

  test_that("output of p_value() method has correct attributes", {

    p_val <- p_value(sobel)
    expect_length(p_val, 5L)
    expect_named(p_val, effect_names)
    expect_equivalent(p_val["Indirect"], sobel$p_value)

  })

  test_that("summary has correct structure", {

    # summary
    expect_s3_class(summary, "summary_test_mediation")
    # original output of test for indirect effect
    expect_identical(summary$object, sobel)
    # summary of the model fit
    expect_s3_class(summary$summary, "summary_fit_mediation")
    # no diagnostic plot
    expect_null(summary$plot)

  })

  test_that("attributes are correctly passed through summary", {

    # robustness
    if (method == "robust") {
      expect_identical(summary$summary$robust, "MM")
    } else if (method == "median") {
      expect_identical(summary$summary$robust, "median")
    } else if (method == "winsorized") {
      expect_true(summary$summary$robust)
    } else {
      expect_false(summary$summary$robust)
    }
    # number of observations
    expect_identical(summary$summary$n, as.integer(n))
    # variable names
    expect_identical(summary$summary$x, x)
    expect_identical(summary$summary$y, y)
    expect_identical(summary$summary$m, m)
    expect_identical(summary$summary$covariates, covariates)

  })

  test_that("effect summaries have correct names", {

    # total effect
    expect_identical(dim(summary$summary$total), c(1L, 4L))
    expect_identical(rownames(summary$summary$total), x)
    expect_identical(colnames(summary$summary$total)[1], "Estimate")
    # direct effect
    expect_identical(dim(summary$summary$direct), c(1L, 4L))
    expect_identical(rownames(summary$summary$direct), x)
    expect_identical(colnames(summary$summary$direct)[1], "Estimate")

  })

  test_that("effect summaries contain correct coefficient values", {

    expect_identical(summary$summary$total[x, "Estimate"], sobel$fit$total)
    expect_identical(summary$summary$direct[x, "Estimate"], sobel$fit$direct)

  })


  # loop over settings for retest()
  for (setting in c("order", "less", "greater")) {

    # use retest() to change one of the arguments
    if (setting == "order") resobel <- retest(sobel, order = "second")
    else resobel <- retest(sobel, alternative = setting)

    test_that("output of retest() has correct structure", {

      # bootstrap test
      expect_identical(class(resobel), class(sobel))
      # regression fit
      expect_identical(resobel$fit, sobel$fit)
      # bootstrap replicates
      expect_identical(resobel$reps, sobel$reps)

    })

    test_that("arguments of retest() are correctly passed", {

      # alternative hypothesis and order of approximation
      if (setting == "order") {
        expect_identical(resobel$alternative, sobel$alternative)
        expect_identical(resobel$order, "second")
      } else {
        expect_identical(resobel$alternative, setting)
        expect_identical(resobel$order, sobel$order)
      }

    })

    test_that("retest() yields correct values", {

      # standard error, test statistic, and p-value
      if (setting == "order") {
        expect_gt(resobel$se, sobel$se)
        expect_lt(abs(resobel$statistic), abs(sobel$statistic))
        expect_gt(resobel$p_value, sobel$p_value)
      } else {
        # standard error and test statistic
        expect_equal(resobel$se, sobel$se)
        expect_equal(resobel$statistic, sobel$statistic)
        # p-value
        is_half <- (setting == "less" && sobel$statistic < 0) ||
          (setting == "greater" && sobel$statistic > 0)
        if (is_half) {
          expect_equal(resobel$p_value, sobel$p_value / 2)
        } else {
          expect_equal(resobel$p_value, 1 - sobel$p_value / 2)
        }
      }

    })

  }


  # tests for weight_plot(), which is only implemented for ROBMED
  test_that("setup_weight_plot() behaves as expected", {

    if (method == "robust") {
      # plot is inherited from model fit
      expect_identical(setup_weight_plot(sobel), setup_weight_plot(sobel$fit))
    } else {
      # plot is not implemented
      expect_error(setup_weight_plot(sobel))
    }

  })


  # tests for ellipse_plot(), which is only implemented for some methods
  test_that("setup_ellipse_plot() behaves as expected", {

    if (method %in% c("robust", "OLS", "winsorized", "ML")) {
      # plot is inherited from model fit
      expect_identical(setup_ellipse_plot(sobel), setup_ellipse_plot(sobel$fit))
    } else {
      # plot is not implemented
      expect_error(setup_ellipse_plot(sobel))
    }

  })


  # tests for density_plot()

  # obtain setup object for ellipse plot
  density <- setup_density_plot(sobel, level = level)

  # run tests
  test_that("object returned by setup_density_plot() has correct structure", {

    # check data frame for confidence interval
    expect_s3_class(density$density, "data.frame")
    # check dimensions
    expect_identical(ncol(density$density), 2L)
    expect_gt(nrow(density$density), 0L)
    # check column names
    expect_named(density$density, c("Indirect", "Density"))
    # check data frame confidence interval
    expect_s3_class(density$ci, "data.frame")
    # check dimensions
    expect_identical(dim(density$ci), c(1L, 3L))
    # check column names
    expect_named(density$ci, c("Estimate", "Lower", "Upper"))
    # check type of test
    expect_identical(density$test, "sobel")
    # check confidence level
    expect_identical(density$level, level)
    # check logical for multiple effects
    expect_false(density$have_effect)
    # check logical for multiple methods
    expect_false(density$have_methods)

  })

}


## additional tests for regression fits

# relevant information
reg_methods <- c("robust", "median", "OLS", "skewnormal", "select")
model_summary_classes <- c(robust = "summary_lmrob", median = "summary_rq",
                           skewnormal = "summary_lmse", select = "summary_lm",
                           OLS = "summary_lm")

# loop over methods
for (method in intersect(methods, reg_methods)) {

  # extract information for current method
  sobel <- sobel_list[[method]]
  summary <- summary_list[[method]]

  # correct values
  model_summary_class <- model_summary_classes[method]
  family <- if (method %in% c("skewnormal", "select")) method else "gaussian"
  intercept_name <- if (method == "skewnormal") "(Intercept.CP)" else "(Intercept)"


  # run tests

  test_that("output has correct structure", {

    # regression fit
    expect_s3_class(sobel$fit, "reg_fit_mediation")

  })

  test_that("arguments are correctly passed", {

    # assumed error distribution
    expect_identical(sobel$fit$family, family)
    # mediation model
    expect_identical(sobel$fit$model, "simple")
    # no contrasts
    expect_false(sobel$fit$contrast)

  })

  test_that("summary has correct structure", {

    # summary of the model fit
    expect_s3_class(summary$summary, "summary_reg_fit_mediation")

  })

  # loop over model summaries
  for (fit_name in c("fit_mx", "fit_ymx")) {

    # extract model summary
    model_summary <- summary$summary[[fit_name]]

    # correct values
    if (fit_name == "fit_mx") {
      coef_names <- c(intercept_name, x)
      n_coef <- 2L
      df1 <- 1
    } else {
      coef_names <- c(intercept_name, m, x)
      n_coef <- 3L
      df1 <- 2
    }

    # run tests
    test_that("model summary has correct structure", {

      # model summary
      expect_s3_class(model_summary, model_summary_class)
      # tests on coefficients
      expect_identical(dim(model_summary$coefficients), c(n_coef, 4L))
      expect_identical(rownames(model_summary$coefficients), coef_names)
      expect_identical(colnames(model_summary$coefficients)[1], "Estimate")
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
          expect_identical(df_test[1], df1)
          expect_identical(df_test[2], Inf)
        } else {
          expect_identical(df_test[1], as.integer(df1))
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

  test_that("effect summaries contain correct coefficient values", {

    expect_identical(summary$summary$fit_mx$coefficients[x, "Estimate"],
                     sobel$fit$a)
    expect_identical(summary$summary$fit_ymx$coefficients[m, "Estimate"],
                     sobel$fit$b)

  })

}


## additional tests for covariance fits

# relevant information
cov_methods <- c("winsorized", "ML")

# loop over methods
for (method in intersect(methods, cov_methods)) {

  # extract information for current method
  sobel <- sobel_list[[method]]
  summary <- summary_list[[method]]


  # run tests

  test_that("output has correct structure", {

    # covariance fit
    expect_s3_class(sobel$fit, "cov_fit_mediation")

  })

  test_that("summary has correct structure", {

    # summary of the model fit
    expect_s3_class(summary$summary, "summary_cov_fit_mediation")
    # no regression model summaries
    expect_null(summary$summary$fit_mx)
    expect_null(summary$summary$fit_ymx)

  })

  test_that("effect summaries have correct names", {

    # a path
    expect_identical(dim(summary$summary$a), c(1L, 4L))
    expect_identical(rownames(summary$summary$a), x)
    expect_identical(colnames(summary$summary$a)[1], "Estimate")
    # b path
    expect_identical(dim(summary$summary$b), c(1L, 4L))
    expect_identical(rownames(summary$summary$b), m)
    expect_identical(colnames(summary$summary$b)[1], "Estimate")

  })

  test_that("effect summaries contain correct coefficient values", {

    expect_identical(summary$summary$a[x, "Estimate"], sobel$fit$a)
    expect_identical(summary$summary$b[m, "Estimate"], sobel$fit$b)

  })

}
