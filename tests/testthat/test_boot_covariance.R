context("bootstrap test: covariance fit")


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

## arguments for bootstrap tests
x <- "X"                         # independent variable
y <- "Y"                         # dependent variable
m <- "M"                         # mediator variable
covariates <- character()        # control variables
R <- 100                         # number of bootstrap samples
level <- c(0.9, 0.95)            # confidence level
ci_type <- "perc"                # type of confidence intervals
ctrl <- cov_control(prob = 0.9)  # for winsorized covariance matrix

## perform bootstrap tests
boot_list <- list(
  winsorized = {
    set.seed(seed)
    test_mediation(test_data, x = x, y = y, m = m, test = "boot", R = R,
                   level = level[1], type = ci_type, method = "covariance",
                   robust = TRUE, control = ctrl)
  },
  ML = {
    set.seed(seed)
    test_mediation(test_data, x = x, y = y, m = m, test = "boot", R = R,
                   level = level[1], type = ci_type, method = "covariance",
                   robust = FALSE)
  }
)

## correct values
effect_names <- c("a", "b", "Total", "Direct", "Indirect")


## run tests

# loop over methods
methods <- names(boot_list)
for (method in methods) {

  # extract information for current method
  boot <- boot_list[[method]]

  # compute summaries
  summary_list <- list(boot = summary(boot, type = "boot"),
                       data = summary(boot, type = "data"))


  ## run tests

  test_that("output has correct structure", {

    # bootstrap test
    expect_s3_class(boot, "boot_test_mediation")
    expect_s3_class(boot, "test_mediation")
    # model fit
    expect_s3_class(boot$fit, "cov_fit_mediation")
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
    if (method == "winsorized") {
      expect_true(boot$fit$robust)
      expect_equal(boot$fit$control, ctrl)
    } else {
      expect_false(boot$fit$robust)
      expect_null(boot$fit$control)
    }

  })

  test_that("dimensions are correct", {

    # effects are scalars
    expect_length(boot$indirect, 1L)
    # only one indirect effect, so only one confidence interval
    expect_length(boot$ci, 2L)
    expect_named(boot$ci, c("Lower", "Upper"))
    # dimensions of bootstrap replicates
    expect_identical(dim(boot$reps$t), c(as.integer(R), 5L))

  })

  test_that("values of coefficients are correct", {

    expect_equivalent(boot$a, mean(boot$reps$t[, 1]))
    expect_equivalent(boot$b, mean(boot$reps$t[, 2]))
    expect_null(boot[["d"]])
    expect_equivalent(boot$direct, mean(boot$reps$t[, 4]))
    expect_equivalent(boot$indirect, mean(boot$reps$t[, 5]))
    expect_equivalent(boot$reps$t[, 5], boot$reps$t[, 1] * boot$reps$t[, 2])
    # total effect
    expect_equivalent(boot$total, mean(boot$reps$t[, 3]))
    expect_equivalent(boot$reps$t[, 3], rowSums(boot$reps$t[, 4:5]))

  })


  # loop over types of estimation (bootstrap or on original data)
  types <- names(summary_list)
  for (type in types) {

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

      # confidence interval of indirect effect
      ci_indirect <- confint(boot, parm = "indirect", type = type)
      expect_equivalent(ci_indirect, boot$ci)
      # extract confidence intervals for other effects
      keep <- c("a", "b", "total", "direct")
      ci_other <- confint(boot, parm = keep, type = type)
      coef_other <- coef(boot, parm = keep, type = type)
      expect_equivalent(rowMeans(ci_other), coef_other)

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


    # extract current summary
    summary <- summary_list[[type]]

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
      expect_s3_class(summary$summary, "summary_cov_fit_mediation")
      expect_s3_class(summary$summary, "summary_fit_mediation")
      # no regression model summaries
      expect_null(summary$summary$fit_mx)
      expect_null(summary$summary$fit_ymx)

    })

    test_that("attributes are correctly passed through summary", {

      # robustness
      if (method == "winsorized") {
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

      # a path
      expect_identical(dim(summary$summary$a), c(1L, n_col))
      expect_identical(rownames(summary$summary$a), x)
      expect_identical(colnames(summary$summary$a)[keep], col_names)
      # b path
      expect_identical(dim(summary$summary$b), c(1L, n_col))
      expect_identical(rownames(summary$summary$b), m)
      expect_identical(colnames(summary$summary$b)[keep], col_names)
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
      expect_identical(summary$summary$a[x, 1], boot$fit$a)
      expect_identical(summary$summary$b[m, 1], boot$fit$b)
      expect_identical(summary$summary$total[x, 1], boot$fit$total)
      expect_identical(summary$summary$direct[x, 1], boot$fit$direct)

      # bootstrapped effects
      if (type == "boot") {
        expect_equal(summary$summary$a[x, "Boot"], boot$a)
        expect_equal(summary$summary$b[m, "Boot"], boot$b)
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
      expect_identical(length(reboot$ci), length(boot$ci))
      expect_identical(names(reboot$ci), names(boot$ci))
      if (setting == "wider") {
        expect_lt(reboot$ci[1], boot$ci[1])
        expect_gt(reboot$ci[2], boot$ci[2])
        expect_equal(colnames(confint(reboot)), c("2.5 %", "97.5 %"))
      } else {
        if (setting == "less") {
          expect_equivalent(reboot$ci[1], -Inf)
          expect_equal(reboot$ci[2], boot$ci[2])
        } else {
          expect_equal(reboot$ci[1], boot$ci[1])
          expect_equivalent(reboot$ci[2], Inf)
        }
        expect_identical(colnames(confint(reboot)), c("Lower", "Upper"))
      }

    })

  }


  # tests for weight_plot(), which is only implemented for ROBMED
  test_that("setup_weight_plot() gives error", {

    # plot is not implemented
    expect_error(setup_weight_plot(boot))

  })


  # tests for ellipse_plot()
  test_that("setup_ellipse_plot() behaves as expected", {

    # plot is inherited from model fit
    expect_identical(setup_ellipse_plot(boot), setup_ellipse_plot(boot$fit))

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
        expect_identical(dim(ci$ci), c(2L, 4L))
        expect_named(ci$ci, c("Effect", "Estimate", "Lower", "Upper"))
      } else {
        expect_identical(dim(ci$ci), c(2L, 5L))
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
        expect_identical(dim(ci$p_value), c(2L, 3L))
        expect_named(ci$p_value, c("Label", "Effect", "Value"))
        # check that labels are correct
        labels <- c("Confidence interval", "p-Value")
        expect_identical(ci$ci$Label,
                         factor(rep.int(labels[1], 2), levels = labels))
        expect_identical(ci$p_value$Label,
                         factor(rep.int(labels[2], 2), levels = labels))
      }

      # check that direct effect and indirect effect are plotted by default
      names <- c("Direct", "Indirect")
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
    expect_identical(density$test, "boot")
    # check confidence level
    expect_identical(density$level, level[1])
    # check logical for multiple effects
    expect_false(density$have_effect)
    # check logical for multiple methods
    expect_false(density$have_methods)

  })

}
