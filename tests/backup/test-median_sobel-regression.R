context("Sobel test: median regression")


## load package
library("robmed", quietly = TRUE)

## control parameters
n <- 250          # number of observations
a <- c <- 0.2     # true effects
b <- 0            # true effect
R <- 1000         # number of bootstrap samples
seed <- 20190201  # seed for the random number generator

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

## run bootstrap test
sobel <- test_mediation(test_data, x = "X", y = "Y", m = "M1", test = "sobel",
                        method = "regression", robust = "median")

## compute summary
summary_sobel <- summary(sobel)

## retest with different parameters
sobel_less <- retest(sobel, alternative = "less")
sobel_greater <- retest(sobel, alternative = "greater")
sobel_second <- retest(sobel, order = "second")

## create data for plotting
level <- 0.9
ci <- setup_ci_plot(sobel, level = level)
ci_p <- setup_ci_plot(sobel, level = level, p_value = TRUE)
density <- setup_density_plot(sobel, level = level)

## stuff needed to check correctness
coef_names <- c("a", "b", "Direct", "Total", "ab")
mx_names <- c("(Intercept)", "X")
ymx_names <- c("(Intercept)", "M1", "X")


## run tests

test_that("output has correct structure", {

  # bootstrap test
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
  # variable names
  expect_identical(sobel$fit$x, "X")
  expect_identical(sobel$fit$y, "Y")
  expect_identical(sobel$fit$m, "M1")
  expect_identical(sobel$fit$covariates, character())
  # robust fit and test
  expect_identical(sobel$fit$robust, "median")
  expect_identical(sobel$fit$family, "gaussian")
  expect_null(sobel$fit$control)

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

  ci_default <- confint(sobel, parm = "ab")  # should be 95%
  ci_90 <- confint(sobel, parm = "ab", level = level)
  # default CI should be wider
  expect_lt(ci_default["ab", 1], ci_90["ab", 1])
  expect_gt(ci_default["ab", 2], ci_90["ab", 2])

})

test_that("summary has correct structure", {

  # summary
  expect_s3_class(summary_sobel, "summary_test_mediation")
  # original output of test for indirect effect
  expect_identical(summary_sobel$object, sobel)
  # summary of the model fit
  expect_s3_class(summary_sobel$summary, "summary_reg_fit_mediation")
  expect_s3_class(summary_sobel$summary, "summary_fit_mediation")
  # summary for model m ~ x
  expect_s3_class(summary_sobel$summary$fit_mx, "summary_rq")
  # summary for model y ~ m + x
  expect_s3_class(summary_sobel$summary$fit_ymx, "summary_rq")
  # regression standard error for model y ~ m + x
  expect_null(summary_sobel$summary$fit_ymx$s)
  # R-squared for model y ~ m + x
  expect_null(summary_sobel$summary$fit_ymx$R2)
  # F-test for model y ~ m + x
  expect_null(summary_sobel$summary$fit_ymx$F_test)
  # no plot is created
  expect_null(summary_sobel$plot)

})

test_that("attributes are correctly passed through summary", {

  # robustness
  expect_identical(summary_sobel$summary$robust, "median")
  # number of observations
  expect_identical(summary_sobel$summary$n, as.integer(n))
  # variable names
  expect_identical(summary_sobel$summary$x, "X")
  expect_identical(summary_sobel$summary$y, "Y")
  expect_identical(summary_sobel$summary$m, "M1")
  expect_identical(summary_sobel$summary$covariates, character())

})

test_that("effect summaries have correct names", {

  # a path
  expect_identical(dim(summary_sobel$summary$fit_mx$coefficients),
                   c(2L, 4L))
  expect_identical(rownames(summary_sobel$summary$fit_mx$coefficients),
                   mx_names)
  expect_identical(colnames(summary_sobel$summary$fit_mx$coefficients)[1],
                   "Estimate")
  # b path
  expect_identical(dim(summary_sobel$summary$fit_ymx$coefficients),
                   c(3L, 4L))
  expect_identical(rownames(summary_sobel$summary$fit_ymx$coefficient),
                   ymx_names)
  expect_identical(colnames(summary_sobel$summary$fit_ymx$coefficient)[1],
                   "Estimate")
  # direct effect
  expect_identical(dim(summary_sobel$summary$direct),
                   c(1L, 4L))
  expect_identical(rownames(summary_sobel$summary$direct),
                   "X")
  expect_identical(colnames(summary_sobel$summary$direct)[1],
                   "Estimate")
  # total effect
  expect_identical(dim(summary_sobel$summary$total),
                   c(1L, 4L))
  expect_identical(rownames(summary_sobel$summary$total),
                   "X")
  expect_identical(colnames(summary_sobel$summary$total)[1],
                   "Estimate")

})

test_that("effect summaries contain correct coefficient values", {

  expect_identical(summary_sobel$summary$fit_mx$coefficients["X", "Estimate"],
                   sobel$fit$a)
  expect_identical(summary_sobel$summary$fit_ymx$coefficients["M1", "Estimate"],
                   sobel$fit$b)
  expect_identical(summary_sobel$summary$direct["X", "Estimate"],
                   sobel$fit$direct)
  expect_identical(summary_sobel$summary$total["X", "Estimate"],
                   sobel$fit$total)

})

test_that("output of retest() has correct structure", {

  # Sobel test
  expect_identical(class(sobel_less), class(sobel))
  expect_identical(class(sobel_greater), class(sobel))
  expect_identical(class(sobel_second), class(sobel))
  # mediation model fit
  expect_identical(sobel_less$fit, sobel$fit)
  expect_identical(sobel_greater$fit, sobel$fit)
  expect_identical(sobel_second$fit, sobel$fit)

})

test_that("arguments of retest() are correctly passed", {

  # standard error
  expect_identical(sobel_less$se, sobel$se)
  expect_identical(sobel_greater$se, sobel$se)
  expect_true(sobel_second$se > sobel$se)
  # test statistic
  expect_identical(sobel_less$statistic, sobel$statistic)
  expect_identical(sobel_greater$statistic, sobel$statistic)
  expect_true(abs(sobel_second$statistic) < abs(sobel$statistic))
  # p-value
  expect_equal(sobel_less$p_value, 1-sobel$p_value/2)
  expect_equal(sobel_greater$p_value, sobel$p_value/2)
  expect_true(sobel_second$p_value > sobel$p_value)
  # alternative hypothesis
  expect_identical(sobel_less$alternative, "less")
  expect_identical(sobel_greater$alternative, "greater")
  expect_identical(sobel_second$alternative, sobel$alternative)

})

test_that("output of p_value() method has correct attributes", {

  p_data <- p_value(sobel, parm = NULL)
  expect_length(p_data, 5L)
  expect_named(p_data, coef_names)
  expect_equivalent(p_data["ab"], sobel$p_value)

})

test_that("objects returned by setup_xxx_plot() have correct structure", {

  ## ci plot without p-value
  # check data frame for confidence interval
  expect_s3_class(ci$ci, "data.frame")
  # check dimensions
  expect_identical(dim(ci$ci), c(2L, 4L))
  # check column names
  expect_named(ci$ci, c("Effect", "Estimate", "Lower", "Upper"))
  # check that direct effect and indirect effect are plotted by default
  effect_names <- c("Direct", "ab")
  effect_factor <- factor(effect_names, levels = effect_names)
  expect_identical(ci$ci$Effect, effect_factor)
  # check confidence level
  expect_identical(ci$level, level)
  # check logical for multiple methods
  expect_false(ci$have_methods)

  ## ci plot with p-value
  # check data frame for confidence interval and p-value
  expect_s3_class(ci_p$ci, "data.frame")
  expect_s3_class(ci_p$p_value, "data.frame")
  # check dimensions
  expect_identical(dim(ci_p$ci), c(2L, 5L))
  expect_identical(dim(ci_p$p_value), c(2L, 3L))
  # check column names
  expect_named(ci_p$ci, c("Label", "Effect", "Estimate", "Lower", "Upper"))
  expect_named(ci_p$p_value, c("Label", "Effect", "Value"))
  # check that labels are correct
  label_names <- c("Confidence interval", "p-Value")
  expect_identical(ci_p$ci$Label,
                   factor(rep.int(label_names[1], 2), levels = label_names))
  expect_identical(ci_p$p_value$Label,
                   factor(rep.int(label_names[2], 2), levels = label_names))
  # check that direct effect and indirect effect are plotted by default
  effect_names <- c("Direct", "ab")
  effect_factor <- factor(effect_names, levels = effect_names)
  expect_identical(ci_p$ci$Effect, effect_factor)
  expect_identical(ci_p$p_value$Effect, effect_factor)
  # check confidence level
  expect_identical(ci$level, level)
  # check logical for multiple methods
  expect_false(ci$have_methods)

  ## density plot
  # check data frame for confidence interval
  expect_s3_class(density$density, "data.frame")
  # check dimensions
  expect_identical(ncol(density$density), 2L)
  expect_gt(nrow(density$density), 0L)
  # check column names
  expect_named(density$density, c("ab", "Density"))
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

  ## ellipse plot and weight plot
  expect_error(setup_ellipse_plot(sobel))
  expect_error(setup_weight_plot(sobel))

})


# run mediation analysis through formula interface with data argument
sobel_f1 <- test_mediation(Y ~ m(M1) + X, data = test_data, test = "sobel",
                           method = "regression", robust = "median")
# run mediation analysis through formula interface without data argument
sobel_f2 <- test_mediation(Y ~ m(M1) + X, test = "sobel",
                           method = "regression", robust = "median")
# define mediator outside formula
med <- m(M1)
sobel_f3 <- test_mediation(Y ~ med + X, data = test_data, test = "sobel",
                           method = "regression", robust = "median")


test_that("formula interface works correctly", {

  # check that results are the same as with default method
  expect_equal(sobel_f1, sobel)
  expect_equal(sobel_f2, sobel)
  expect_equal(sobel_f3, sobel)

})
