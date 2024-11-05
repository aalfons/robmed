context("bootstrap test: parallel multiple mediators, covariates")


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

## arguments for bootstrap tests
x <- "X"                                 # independent variable
y <- "Y"                                 # dependent variable
m <- c("M1", "M2")                       # parallel mediators
covariates <- c("C1", "C2")              # control variables
R <- 100                                 # number of bootstrap samples
level <- c(0.9, 0.95)                    # confidence level
ci_type <- "perc"                        # type of confidence intervals
robust_ctrl <- MM_reg_control(max_iterations = 500)  # for MM-regression
median_ctrl <- median_reg_control(algorithm = "fn")  # for median regression

## perform bootstrap tests
boot_list <- list(
  robust = {
    set.seed(seed)
    test_mediation(test_data, x = x, y = y, m = m, covariates = covariates,
                   model = "parallel", test = "boot", R = R, level = level[1],
                   type = ci_type, method = "regression", robust = TRUE,
                   control = robust_ctrl)
  },
  median = {
    set.seed(seed)
    test_mediation(test_data, x = x, y = y, m = m, covariates = covariates,
                   model = "parallel", test = "boot", R = R, level = level[1],
                   type = ci_type, method = "regression", robust = "median",
                   control = median_ctrl)
  },
  # skewnormal = {
  #   set.seed(seed)
  #   suppressWarnings(
  #     test_mediation(test_data, x = x, y = y, m = m, covariates = covariates,
  #                    model = "parallel", test = "boot", R = R, level = level[1],
  #                    type = ci_type,  method = "regression", robust = FALSE,
  #                    family = "skewnormal")
  #   )
  # },
  # select = {
  #   set.seed(seed)
  #   suppressWarnings(
  #     test_mediation(test_data, x = x, y = y, m = m, covariates = covariates,
  #                    model = "parallel", test = "boot", R = R, level = level[1],
  #                    type = ci_type, method = "regression", robust = FALSE,
  #                    family = "select")
  #   )
  # },
  OLS = {
    set.seed(seed)
    test_mediation(test_data, x = x, y = y, m = m, covariates = covariates,
                   model = "parallel", test = "boot", R = R, level = level[1],
                   type = ci_type, method = "regression", robust = FALSE,
                   family = "gaussian")
  }
)

## correct values
effect_names <- c("a_M1", "a_M2", "b_M1", "b_M2", "Total", "Direct",
                  "Indirect_Total", "Indirect_M1", "Indirect_M2")
model_summary_classes <- c(robust = "summary_lmrob", median = "summary_rq",
                           skewnormal = "summary_lmse", select = "summary_lm",
                           OLS = "summary_lm")


## run tests

# loop over methods
methods <- names(boot_list)
for (method in methods) {

  # extract information for current method
  boot <- boot_list[[method]]

  # correct values
  family <- if (method %in% c("skewnormal", "select")) method else "gaussian"
  model_summary_class <- model_summary_classes[method]
  intercept_name <- if (method == "skewnormal") "(Intercept.CP)" else "(Intercept)"


  # run tests

  test_that("output has correct structure", {

    # bootstrap test
    expect_s3_class(boot, "boot_test_mediation")
    expect_s3_class(boot, "test_mediation")
    # regression fit
    expect_s3_class(boot$fit, "reg_fit_mediation")
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
    expect_identical(boot$level, level[1])
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
      expect_equal(boot$fit$control, robust_ctrl)
    } else if (method == "median") {
      expect_identical(boot$fit$robust, "median")
      expect_equal(boot$fit$control, median_ctrl)
    } else {
      expect_false(boot$fit$robust)
      expect_null(boot$fit$control)
    }
    # assumed error distribution
    expect_identical(boot$fit$family, family)
    # mediation model
    expect_identical(boot$fit$model, "parallel")
    # no contrasts
    expect_false(boot$fit$contrast)

  })

  test_that("dimensions are correct", {

    # bootstrap effect estimates
    expect_length(boot$a, 2L)
    expect_named(boot$a, m)
    expect_length(boot$b, 2L)
    expect_named(boot$b, m)
    expect_length(boot$direct, 1L)
    expect_named(boot$direct, NULL)
    expect_length(boot$total, 1L)
    expect_named(boot$total, NULL)
    expect_length(boot$indirect, 3L)
    expect_named(boot$indirect, c("Total", m))
    # bootstrap confidence interval
    expect_identical(dim(boot$ci), c(3L, 2L))
    expect_identical(rownames(boot$ci), names(boot$indirect))
    expect_identical(colnames(boot$ci), c("Lower", "Upper"))
    # dimensions of bootstrap replicates
    n_coef <- if (method %in% c("robust", "median", "OLS")) 14L else 18L
    expect_identical(dim(boot$reps$t), c(as.integer(R), n_coef))

  })

  test_that("values of coefficients are correct", {

    # extract bootstrap replicates
    boot_a <- boot$reps$t[, c(2, 6)]
    boot_b <- boot$reps$t[, 10:11]
    boot_direct <- boot$reps$t[, 12]
    boot_indirect <- boot_a * boot_b
    boot_indirect_total <- rowSums(boot_indirect)
    boot_indirect <- cbind(boot_indirect_total, boot_indirect)
    if (method %in% c("robust", "median", "OLS")) {
      boot_total <- boot_indirect_total + boot_direct
    } else boot_total <- boot$reps$t[, 16]
    # check bootstrap estimates
    expect_equivalent(boot$a, colMeans(boot_a))
    expect_equivalent(boot$b, colMeans(boot_b))
    expect_null(boot[["d"]])
    expect_equivalent(boot$total, mean(boot_total))
    expect_equivalent(boot$direct, mean(boot_direct))
    expect_equivalent(boot$indirect, colMeans(boot_indirect))

  })


  # loop over types of estimation (bootstrap or on original data)
  for (type in c("boot", "data")) {

    # run tests

    test_that("output of coef() method has correct attributes", {

      # extract coefficients
      coef <- coef(boot, type = type)
      # tests
      expect_length(coef, 9L)
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
      expect_equal(dim(ci), c(9L, 2L))
      expect_equal(rownames(ci), effect_names)
      expect_equal(colnames(ci), c("5 %", "95 %"))

    })

    test_that("confint() method returns correct values of confidence intervals", {

      # confidence interval of indirect effect
      ci_indirect <- confint(boot, parm = "indirect", type = type)
      expect_equivalent(ci_indirect, boot$ci)
      # extract confidence intervals for other effects
      if (type == "boot") keep <- c("a", "b", "total", "direct")
      else keep <- c("a", "b", "direct")
      ci_other <- confint(boot, parm = keep, type = type)
      coef_other <- coef(boot, parm = keep, type = type)
      expect_equivalent(rowMeans(ci_other), coef_other)

    })

    test_that("output of p_value() method has correct attributes", {

      # extract p-value
      digits <- 3
      p_val <- p_value(boot, type = type, digits = digits)
      # check dimensions
      expect_length(p_val, 9L)
      expect_named(p_val, effect_names)
      # check number of digits
      which <- grep("Indirect", names(p_val))
      expect_equal(p_val[which], round(p_val[which], digits = digits))

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
      # list of objects for summaries of regressions m ~ x
      expect_type(summary$summary$fit_mx, "list")

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
        model_summary <- summary$summary$fit_mx[[response]]
        fit <- boot$fit$fit_mx[[response]]
      }

      # correct values
      if (response == "Y") {
        coef_names <- c(intercept_name, m, x, covariates)
        n_coef <- 6L
        df1 <- 5
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
      a_data <- c(summary$summary$fit_mx$M1$coefficients[x, 1],
                  summary$summary$fit_mx$M2$coefficients[x, 1])
      expect_equivalent(a_data, boot$fit$a)
      expect_identical(summary$summary$fit_ymx$coefficients[m, 1], boot$fit$b)
      expect_identical(summary$summary$total[x, 1], boot$fit$total)
      expect_identical(summary$summary$direct[x, 1], boot$fit$direct)

      # bootstrapped effects
      if (type == "boot") {
        a_boot <- c(summary$summary$fit_mx$M1$coefficients[x, "Boot"],
                    summary$summary$fit_mx$M2$coefficients[x, "Boot"])
        expect_equivalent(a_boot, boot$a)
        expect_equivalent(summary$summary$fit_ymx$coefficients[m, "Boot"],
                          boot$b)
        expect_equal(summary$summary$total[x, "Boot"], boot$total)
        expect_equal(summary$summary$direct[x, "Boot"], boot$direct)
      }

    })

  }


  # loop over settings for retest()
  for (setting in c("wider", "less", "greater")) {

    # use retest() to change one of the arguments
    if (setting == "wider") reboot <- retest(boot, level = level[2])
    else reboot <- retest(boot, alternative = setting, level = level[2])

    test_that("output of retest() has correct structure", {

      # bootstrap test
      expect_identical(class(reboot), class(boot))
      # regression fit
      expect_identical(reboot$fit, boot$fit)
      # bootstrap replicates
      expect_identical(reboot$reps, boot$reps)

    })

    test_that("arguments of retest() are correctly passed", {

      # alternative hypothesis
      if (setting == "wider") {
        expect_identical(reboot$alternative, boot$alternative)
      } else {
        expect_identical(reboot$alternative, setting)
      }
      # confidence level
      expect_identical(reboot$level, level[2])
      # type of confidence intervals
      expect_identical(reboot$type, boot$type)

    })

    test_that("retest() yields correct values", {

      # indirect effect
      expect_identical(reboot$indirect, boot$indirect)
      # confidence interval
      expect_identical(dim(reboot$ci), dim(boot$ci))
      expect_identical(dimnames(reboot$ci), dimnames(boot$ci))
      if (setting == "wider") {
        expect_true(all(reboot$ci[, 1] < boot$ci[, 1]))
        expect_true(all(reboot$ci[, 2] > boot$ci[, 2]))
        expect_equal(colnames(confint(reboot)), c("2.5 %", "97.5 %"))
      } else {
        if (setting == "less") {
          expect_equivalent(reboot$ci[, 1], rep.int(-Inf, 3))
          expect_equal(reboot$ci[, 2], boot$ci[, 2])
        } else {
          expect_equal(reboot$ci[, 1], boot$ci[, 1])
          expect_equivalent(reboot$ci[, 2], rep.int(Inf, 3))
        }
        expect_identical(colnames(confint(reboot)), c("Lower", "Upper"))
      }

    })

  }


  # tests for weight_plot(), which is only implemented for ROBMED
  test_that("setup_weight_plot() behaves as expected", {

    if (method == "robust") {
      # plot is inherited from model fit
      expect_identical(setup_weight_plot(boot), setup_weight_plot(boot$fit))
    } else {
      # plot is not implemented
      expect_error(setup_weight_plot(boot))
    }

  })


  # tests for ellipse_plot(), which is only implemented for some methods
  test_that("setup_ellipse_plot() behaves as expected", {

    if (method %in% c("robust", "OLS")) {
      # plot is inherited from model fit
      expect_identical(setup_ellipse_plot(boot), setup_ellipse_plot(boot$fit))
    } else {
      # plot is not implemented
      expect_error(setup_ellipse_plot(boot))
    }

  })


  # loop over settings for ci_plot()
  for (setting in c("default", "p_value")) {

    # obtain setup object for CI plot
    if (setting == "default") ci <- setup_ci_plot(boot)
    else ci <- setup_ci_plot(boot, p_value = TRUE)

    # run tests
    test_that("object returned by setup_ci_plot() has correct structure", {

      # check data frame for confidence interval
      expect_s3_class(ci$ci, "data.frame")
      # check dimensions and column names
      if (setting == "default") {
        expect_identical(dim(ci$ci), c(4L, 4L))
        expect_named(ci$ci, c("Effect", "Estimate", "Lower", "Upper"))
      } else {
        expect_identical(dim(ci$ci), c(4L, 5L))
        expect_named(ci$ci, c("Label", "Effect", "Estimate", "Lower", "Upper"))
      }

      # check data frame for p-value
      if (setting == "default") {
        # p-value not plotted
        expect_null(ci$p_value)
      } else {
        # check data frame
        expect_s3_class(ci$p_value, "data.frame")
        # check dimensions and column names
        expect_identical(dim(ci$p_value), c(4L, 3L))
        expect_named(ci$p_value, c("Label", "Effect", "Value"))
        # check that labels are correct
        labels <- c("Confidence interval", "p-Value")
        expect_identical(ci$ci$Label,
                         factor(rep.int(labels[1], 4), levels = labels))
        expect_identical(ci$p_value$Label,
                         factor(rep.int(labels[2], 4), levels = labels))
      }

      # check that direct effect and indirect effect are plotted by default
      names <- c("Direct", "Indirect_Total", "Indirect_M1", "Indirect_M2")
      factor <- factor(names, levels = names)
      expect_identical(ci$ci$Effect, factor)
      # check confidence level
      expect_identical(ci$level, level[1])
      # check logical for multiple methods
      expect_false(ci$have_methods)

    })

  }


  # tests for density_plot()

  # obtain setup object for ellipse plot
  density <- setup_density_plot(boot)

  # run tests
  test_that("object returned by setup_density_plot() has correct structure", {

    # check data frame for confidence interval
    expect_s3_class(density$density, "data.frame")
    # check dimensions
    expect_identical(ncol(density$density), 3L)
    expect_gt(nrow(density$density), 0L)
    # check column names
    expect_named(density$density, c("Effect", "Indirect", "Density"))
    # check data frame confidence interval
    expect_s3_class(density$ci, "data.frame")
    # check dimensions
    expect_identical(dim(density$ci), c(3L, 4L))
    # check column names
    expect_named(density$ci, c("Effect", "Estimate", "Lower", "Upper"))
    # check type of test
    expect_identical(density$test, "boot")
    # check confidence level
    expect_identical(density$level, level[1])
    # check logical for multiple effects
    expect_true(density$have_effect)
    # check logical for multiple methods
    expect_false(density$have_methods)

  })

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
#                                covariates = covariates, model = "parallel",
#                                test = "boot", method = "regression",
#                                robust = robust)
#
#     # try to use covariance fit (should give warning)
#     set.seed(seed)
#     expect_warning(
#       cov_boot <- test_mediation(test_data, x = x, y = y, m = m,
#                                  covariates = covariates, model = "parallel",
#                                  test = "boot", method = "covariance",
#                                  robust = robust)
#     )
#
#     # these should be the same
#     expect_equal(cov_boot, reg_boot)
#
#   })
#
# }
