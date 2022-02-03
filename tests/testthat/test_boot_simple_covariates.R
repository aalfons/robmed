context("bootstrap test: simple mediation, covariates")


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

## arguments for bootstrap tests
x <- "X"                                   # independent variable
y <- "Y"                                   # dependent variable
m <- "M"                                   # mediator variable
covariates <- c("C1", "C2")                # control variables
R <- 100                                   # number of bootstrap samples
level <- 0.9                               # confidence level
ci_type <- "perc"                          # type of confidence intervals
ctrl <- reg_control(max_iterations = 500)  # for MM-regression estimator

## perform bootstrap tests
boot_list <- list(
  robust = {
    set.seed(seed)
    test_mediation(test_data, x = x, y = y, m = m, covariates = covariates,
                   test = "boot", R = R, level = level, type = ci_type,
                   method = "regression", robust = TRUE, control = ctrl)
  },
  median = {
    set.seed(seed)
    test_mediation(test_data, x = x, y = y, m = m, covariates = covariates,
                   test = "boot", R = R, level = level, type = ci_type,
                   method = "regression", robust = "median")
  },
  OLS = {
    set.seed(seed)
    test_mediation(test_data, x = x, y = y, m = m, covariates = covariates,
                   test = "boot", R = R, level = level, type = ci_type,
                   method = "regression", robust = FALSE, family = "gaussian")
    # },
    #   student = {
    #     set.seed(seed)
    #     suppressWarnings(
    #       test_mediation(test_data, x = x, y = y, m = m, covariates = covariates,
    #                      test = "boot", R = R, level = level, type = ci_type,
    #                      method = "regression", robust = FALSE, family = "student")
    #     )
    #   },
    #   select = {
    #     set.seed(seed)
    #     suppressWarnings(
    #       test_mediation(test_data, x = x, y = y, m = m, covariates = covariates,
    #                      test = "boot", R = R, level = level, type = ci_type,
    #                      method = "regression", robust = FALSE, family = "select")
    #     )
  }
)

## correct values
effect_names <- c("a", "b", "Total", "Direct", "Indirect")
model_summary_classes <- c(robust = "summary_lmrob", median = "summary_rq",
                           OLS = "summary_lm", student = "summary_lmse",
                           select = "summary_lm")


## run tests

# loop over methods
methods <- names(boot_list)
for (method in methods) {

  # extract information for current method
  boot <- boot_list[[method]]

  # correct values
  family <- if (method %in% c("student", "select")) method else "gaussian"
  model_summary_class <- model_summary_classes[method]
  intercept_name <- if (method == "student") "(Intercept.DP)" else "(Intercept)"


  # run tests

  test_that("output has correct structure", {

    # bootstrap test
    expect_s3_class(boot, "boot_test_mediation")
    expect_s3_class(boot, "test_mediation")
    # model fit
    expect_s3_class(boot$fit, "fit_mediation")
    # bootstrap replicates
    expect_s3_class(boot$reps, "boot")

  })

  test_that("arguments are correctly passed", {

    # alternative hypothesis
    expect_identical(boot$alternative, "twosided")
    # number of bootstrap replicates
    # (doesn't hold for some methods if there are computational issues)
    expect_identical(boot$R, as.integer(R))
    # confidence level
    expect_identical(boot$level, level)
    # type of confidence intervals
    expect_identical(boot$type, ci_type)
    # variable names
    expect_identical(boot$fit$x, x)
    expect_identical(boot$fit$y, y)
    expect_identical(boot$fit$m, m)
    expect_identical(boot$fit$covariates, covariates)
    # robust or nonrobust fit and test
    if (method == "robust") {
      expect_identical(boot$fit$robust, "MM")
      expect_equal(boot$fit$control, ctrl)
    } else if (method == "median") {
      expect_identical(boot$fit$robust, "median")
      expect_null(boot$fit$control)
    } else {
      expect_false(boot$fit$robust)
      expect_null(boot$fit$control)
    }
    # assumed error distribution
    expect_identical(boot$fit$family, family)
    # mediation model
    expect_identical(boot$fit$model, "simple")
    # no contrasts
    expect_false(boot$fit$contrast)

  })

  test_that("dimensions are correct", {

    # effects are scalars
    expect_length(boot$indirect, 1L)
    # only one indirect effect, so only one confidence interval
    expect_length(boot$ci, 2L)
    expect_named(boot$ci, c("Lower", "Upper"))
    # dimensions of bootstrap replicates
    expect_identical(dim(boot$reps$t), c(as.integer(R), 11L))

  })

  test_that("values of coefficients are correct", {

    expect_equivalent(boot$a, mean(boot$reps$t[, 3]))
    expect_equivalent(boot$b, mean(boot$reps$t[, 7]))
    expect_null(boot[["d"]])
    expect_equivalent(boot$direct, mean(boot$reps$t[, 8]))
    expect_equivalent(boot$indirect, mean(boot$reps$t[, 1]))
    expect_equivalent(boot$reps$t[, 1], boot$reps$t[, 3] * boot$reps$t[, 7])
    # total effect
    expect_equivalent(boot$total, mean(boot$reps$t[, 11]))
    if (method %in% c("robust", "median", "OLS")) {
      expect_equivalent(boot$reps$t[, 11], rowSums(boot$reps$t[, c(1, 8)]))
    }

  })


  # loop over types of estimation (bootstrap or on original data)
  for (type in c("boot", "data")) {

    # run tests

    test_that("output of coef() method has correct attributes", {

      # extract coefficients
      coef <- coef(boot, type = type)
      # tests
      expect_length(coef, 5L)
      expect_named(coef, effect_names)

    })

    test_that("coef() method returns correct values of coefficients", {

      # object containing the true values
      true <- if (type == "boot") boot else boot$fit
      # tests
      expect_equivalent(coef(boot, parm = "a", type = type),
                        true$a)
      expect_equivalent(coef(boot, parm = "b", type = type),
                        true$b)
      expect_null(coef(boot, parm = "d", type = type))
      expect_equivalent(coef(boot, parm = "total", type = type),
                        true$total)
      expect_equivalent(coef(boot, parm = "direct", type = type),
                        true$direct)
      expect_equivalent(coef(boot, parm = "indirect", type = type),
                        true$indirect)

    })

    test_that("output of confint() method has correct attributes", {

      # extract confidence intervals
      ci <- confint(boot, type = type)
      # tests
      expect_equal(dim(ci), c(5L, 2L))
      expect_equal(rownames(ci), effect_names)
      expect_equal(colnames(ci), c("5 %", "95 %"))

    })

    test_that("confint() method returns correct values of confidence intervals", {

      expect_equivalent(confint(boot, parm = "indirect", type = type), boot$ci)

    })

    test_that("output of p_value() method has correct attributes", {

      # extract p-value
      digits <- 3
      p_val <- p_value(boot, type = type, digits = digits)
      # tests
      expect_length(p_val, 5L)
      expect_named(p_val, effect_names)
      expect_equal(p_val["Indirect"], round(p_val["Indirect"], digits = digits))

    })


    # compute summary
    summary <- summary(boot, type = type)

    # correct values
    if (type == "boot") {
      n_col <- 5L
      col_names <- c("Data", "Boot")
      keep <- 1:2
    } else {
      n_col <- 4L
      col_names <- "Estimate"
      keep <- 1
    }


    # run tests

    test_that("summary has correct structure", {

      # summary
      expect_s3_class(summary, "summary_test_mediation")
      # original output of test for indirect effect
      expect_identical(summary$object, boot)
      # summary of the model fit
      expect_s3_class(summary$summary, "summary_reg_fit_mediation")
      expect_s3_class(summary$summary, "summary_fit_mediation")
      # diagnostic plot only for ROBMED
      if (method == "robust") {
        expect_s3_class(summary$plot, "gg_weight_plot")
      } else {
        expect_null(summary$plot)
      }

    })

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
      expect_identical(summary$summary$x, x)
      expect_identical(summary$summary$y, y)
      expect_identical(summary$summary$m, m)
      expect_identical(summary$summary$covariates, covariates)

    })


    # loop over response variables in regression model fits
    for (response in c(m, y)) {

      # extract information on current regression model fit
      if (response == "Y") {
        model_summary <- summary$summary$fit_ymx
        fit <- boot$fit$fit_ymx
      } else {
        model_summary <- summary$summary$fit_mx
        fit <- boot$fit$fit_mx
      }

      # correct values
      if (response == "Y") {
        coef_names <- c(intercept_name, m, x, covariates)
        n_coef <- 5L
        df1 <- 4
      } else {
        coef_names <- c(intercept_name, x, covariates)
        n_coef <- 4L
        df1 <- 3
      }

      # run tests
      test_that("model summary has correct structure", {

        # model summary
        expect_s3_class(model_summary, model_summary_class)
        # tests on coefficients
        expect_identical(dim(model_summary$coefficients), c(n_coef, n_col))
        expect_identical(rownames(model_summary$coefficients), coef_names)
        expect_identical(colnames(model_summary$coefficients)[keep], col_names)
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
        # additional information for robust bootstrap test
        if (method == "robust") {
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


    test_that("effect summaries have correct names", {

      # total effect
      expect_identical(dim(summary$summary$total), c(1L, n_col))
      expect_identical(rownames(summary$summary$total), x)
      expect_identical(colnames(summary$summary$total)[keep], col_names)
      # direct effect
      expect_identical(dim(summary$summary$direct), c(1L, n_col))
      expect_identical(rownames(summary$summary$direct), x)
      expect_identical(colnames(summary$summary$direct)[keep], col_names)

    })

    test_that("effect summaries contain correct coefficient values", {

      # effects computed on original sample
      expect_equivalent(summary$summary$fit_mx$coefficients[x, 1], boot$fit$a)
      expect_identical(summary$summary$fit_ymx$coefficients[m, 1], boot$fit$b)
      expect_identical(summary$summary$total[x, 1], boot$fit$total)
      expect_identical(summary$summary$direct[x, 1], boot$fit$direct)

      # bootstrapped effects
      if (type == "boot") {
        expect_equivalent(summary$summary$fit_mx$coefficients[x, "Boot"],
                          boot$a)
        expect_equivalent(summary$summary$fit_ymx$coefficients[m, "Boot"],
                          boot$b)
        expect_equal(summary$summary$total[x, "Boot"], boot$total)
        expect_equal(summary$summary$direct[x, "Boot"], boot$direct)
      }

    })

  }

}


# ## covariance fits only implemented for simple mediation without covariates
#
# # argument values
# cov_args <- list(winsorized = list(robust = TRUE),
#                  ML = list(robust = FALSE))
#
# # loop over methods
# cov_methods <- names(cov_args)
# for (method in cov_methods) {
#
#   # argument values
#   robust <- cov_args[[method]]$robust
#
#
#   # run tests
#
#   test_that("covariates not implemented", {
#
#     # use regression fit
#     set.seed(seed)
#     reg_boot <- test_mediation(test_data, x = x, y = y, m = m,
#                                covariates = covariates, test = "boot",
#                                method = "regression", robust = robust)
#
#     # try to use covariance fit (should give warning)
#     set.seed(seed)
#     expect_warning(
#       cov_boot <- test_mediation(test_data, x = x, y = y, m = m,
#                                  covariates = covariates, test = "boot",
#                                  method = "covariance", robust = robust)
#     )
#
#     # these should be the same
#     expect_equal(cov_boot, reg_boot)
#
#   })
#
# }
