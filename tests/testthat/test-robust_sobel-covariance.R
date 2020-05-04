context("robust Sobel test: covariance")


## load package
library("robmed", quietly = TRUE)

## control parameters
n <- 250          # number of observations
a <- c <- 0.2     # true effects
b <- 0            # true effect
R <- 1000         # number of bootstrap samples
seed <- 20150601  # seed for the random number generator

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
ctrl <- cov_control(prob = 0.9)
sobel <- test_mediation(test_data, x = "X", y = "Y", m = "M1", test = "sobel",
                        method = "covariance", robust = TRUE, control = ctrl)

## compute summary
summary_sobel <- summary(sobel)

## retest with different parameters
sobel_less <- retest(sobel, alternative = "less")
sobel_greater <- retest(sobel, alternative = "greater")

## create data for plotting
level <- 0.9
ci <- setup_ci_plot(sobel, level = level)
density <- setup_density_plot(sobel, level = level)
ellipse <- setup_ellipse_plot(sobel)

## stuff needed to check correctness
coef_names <- c("a", "b", "Direct", "Total", "ab")


## run tests

test_that("output has correct structure", {

  # bootstrap test
  expect_s3_class(sobel, "sobel_test_mediation")
  expect_s3_class(sobel, "test_mediation")
  # regression fit
  expect_s3_class(sobel$fit, "cov_fit_mediation")
  expect_s3_class(sobel$fit, "fit_mediation")
  # indirect effect
  expect_is(sobel$ab, "numeric")
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
  expect_true(sobel$fit$robust)
  expect_equal(sobel$fit$control, ctrl)

})

test_that("dimensions are correct", {

  # indirect effect
  expect_length(sobel$ab, 1L)
  # standard error
  expect_length(sobel$se, 1L)
  # test statistic
  expect_length(sobel$statistic, 1L)
  # p-value
  expect_length(sobel$p_value, 1L)

})

test_that("values of coefficients are correct", {

  expect_equivalent(sobel$ab, sobel$fit$a * sobel$fit$b)

})

test_that("output of coef() method has correct attributes", {

  coef_sobel <- coef(sobel)
  expect_length(coef_sobel, 5L)
  expect_named(coef_sobel, coef_names)

})

test_that("coef() method returns correct values of coefficients", {

  expect_equivalent(coef(sobel, parm = "ab"), sobel$ab)

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
  expect_s3_class(summary_sobel$summary, "summary_cov_fit_mediation")
  expect_s3_class(summary_sobel$summary, "summary_fit_mediation")
  # regression standard error for model y ~ m + x
  expect_null(summary_sobel$summary$s)
  # R-squared for model y ~ m + x
  expect_null(summary_sobel$summary$R2)
  # F-test for model y ~ m + x
  expect_null(summary_sobel$summary$F_test)

})

test_that("attributes are correctly passed through summary", {

  # robustness
  expect_true(summary_sobel$summary$robust)
  # number of observations
  expect_identical(summary_sobel$summary$n, as.integer(n))
  # variable names
  expect_identical(summary_sobel$summary$x, "X")
  expect_identical(summary_sobel$summary$y, "Y")
  expect_identical(summary_sobel$summary$m, "M1")
  expect_null(summary_sobel$summary$covariates)

})

test_that("effect summaries have correct names", {

  # a path
  expect_identical(dim(summary_sobel$summary$a), c(1L, 4L))
  expect_identical(rownames(summary_sobel$summary$a), "X")
  expect_identical(colnames(summary_sobel$summary$a)[1], "Estimate")
  # b path
  expect_identical(dim(summary_sobel$summary$b), c(1L, 4L))
  expect_identical(rownames(summary_sobel$summary$b), "M1")
  expect_identical(colnames(summary_sobel$summary$b)[1], "Estimate")
  # direct effect
  expect_identical(dim(summary_sobel$summary$direct), c(1L, 4L))
  expect_identical(rownames(summary_sobel$summary$direct), "X")
  expect_identical(colnames(summary_sobel$summary$direct)[1], "Estimate")
  # total effect
  expect_identical(dim(summary_sobel$summary$total), c(1L, 4L))
  expect_identical(rownames(summary_sobel$summary$total), "X")
  expect_identical(colnames(summary_sobel$summary$total)[1], "Estimate")
  # no model fits
  expect_null(summary_sobel$summary$fit_mx)
  expect_null(summary_sobel$summary$fit_ymx)

})

test_that("effect summaries contain correct coefficient values", {

  expect_identical(summary_sobel$summary$a["X", "Estimate"], sobel$fit$a)
  expect_identical(summary_sobel$summary$b["M1", "Estimate"], sobel$fit$b)
  expect_identical(summary_sobel$summary$direct["X", "Estimate"], sobel$fit$direct)
  expect_identical(summary_sobel$summary$total["X", "Estimate"], sobel$fit$total)

})

test_that("output of retest() has correct structure", {

  # Sobel test
  expect_identical(class(sobel_less), class(sobel))
  expect_identical(class(sobel_greater), class(sobel))
  # regression fit
  expect_identical(sobel_less$fit, sobel$fit)
  expect_identical(sobel_greater$fit, sobel$fit)

})

test_that("arguments of retest() are correctly passed", {

  # alternative hypothesis
  expect_identical(sobel_less$alternative, "less")
  expect_identical(sobel_greater$alternative, "greater")
  # indirect effect
  expect_identical(sobel_less$ab, sobel$ab)
  expect_identical(sobel_greater$ab, sobel$ab)
  # standard error
  expect_identical(sobel_less$se, sobel$se)
  expect_identical(sobel_greater$se, sobel$se)
  # test statistic
  expect_identical(sobel_less$statistic, sobel$statistic)
  expect_identical(sobel_greater$statistic, sobel$statistic)
  # p-value
  expect_equal(sobel_less$p_value, 1-sobel$p_value/2)
  expect_equal(sobel_greater$p_value, sobel$p_value/2)

})

test_that("output of p_value() method has correct attributes", {

  p_data <- p_value(sobel, parm = NULL)
  expect_length(p_data, 5L)
  expect_named(p_data, coef_names)
  expect_equivalent(p_data["ab"], sobel$p_value)

})

test_that("objects returned by setup_xxx_plot() have correct structure", {

  ## ci plot
  # check data frame for confidence interval
  expect_s3_class(ci$ci, "data.frame")
  # check dimensions
  expect_identical(dim(ci$ci), c(2L, 4L))
  # check column names
  column_names <- c("Effect", "Estimate", "Lower", "Upper")
  expect_named(ci$ci, column_names)
  # check that direct effect and indirect effect are plotted by default
  effect_names <- c("Direct", "ab")
  expect_identical(ci$ci$Effect, factor(effect_names, levels = effect_names))
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
  column_names <- c("ab", "Density")
  expect_named(density$density, column_names)
  # check data frame confidence interval
  expect_s3_class(density$ci, "data.frame")
  # check dimensions
  expect_identical(dim(density$ci), c(1L, 3L))
  # check column names
  column_names <- c("Estimate", "Lower", "Upper")
  expect_named(density$ci, column_names)
  # check type of test
  expect_identical(density$test, "sobel")
  # check confidence level
  expect_identical(density$level, level)
  # check logical for multiple effects
  expect_false(density$have_effect)
  # check logical for multiple methods
  expect_false(density$have_methods)

  ## ellipse_plot
  expect_identical(ellipse, setup_ellipse_plot(sobel$fit))

})


# run mediation analysis through formula interface with data argument
sobel_f1 <- test_mediation(Y ~ m(M1) + X, data = test_data, test = "sobel",
                           method = "covariance", robust = TRUE, control = ctrl)
# run mediation analysis through formula interface without data argument
sobel_f2 <- test_mediation(Y ~ m(M1) + X, test = "sobel",
                           method = "covariance", robust = TRUE, control = ctrl)
# define mediator outside formula
med <- m(M1)
sobel_f3 <- test_mediation(Y ~ med + X, data = test_data, test = "sobel",
                           method = "covariance", robust = TRUE, control = ctrl)


test_that("formula interface works correctly", {

  # check that results are the same as with default method
  expect_equal(sobel_f1, sobel)
  expect_equal(sobel_f2, sobel)
  expect_equal(sobel_f3, sobel)

})
