context("standard covariance fit")


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
M1 <- a * X + rnorm(n)
M2 <- rnorm(n)
Y <- b * M1 + c * X + rnorm(n)
C1 <- rnorm(n)
C2 <- rnorm(n)
test_data <- data.frame(X, Y, M1, M2, C1, C2)

## fit mediation model and compute summary
foo <- fit_mediation(test_data, x = "X", y = "Y", m = "M1",
                     method = "covariance", robust = FALSE)
bar <- summary(foo)

## create data for plotting
ellipse_mx <- setup_ellipse_plot(foo)
ellipse_ym <- setup_ellipse_plot(foo, horizontal = "M1", vertical = "Y",
                                 partial = FALSE)
ellipse_partial <- setup_ellipse_plot(foo, horizontal = "M1", vertical = "Y",
                                      partial = TRUE)


## run tests

test_that("output has correct structure", {

  # covariance fit
  expect_s3_class(foo, "cov_fit_mediation")
  expect_s3_class(foo, "fit_mediation")
  # covariance matrix
  expect_s3_class(foo$cov, "cov_ML")
  expect_null(foo$control)

})

test_that("arguments are correctly passed", {

  # variable names
  expect_identical(foo$x, "X")
  expect_identical(foo$y, "Y")
  expect_identical(foo$m, "M1")
  expect_null(foo$fit$covariates)
  # standard fit
  expect_false(foo$robust)
  expect_null(foo$control)

})

test_that("dimensions are correct", {

  # effects are scalars
  expect_length(foo$a, 1L)
  expect_length(foo$b, 1L)
  expect_length(foo$direct, 1L)
  expect_length(foo$total, 1L)
  # dimensions of data
  expect_identical(dim(foo$data), c(as.integer(n), 3L))

})

test_that("values of coefficients are correct", {

  expect_equivalent(foo$total, foo$a * foo$b + foo$direct)

})

test_that("output of coef() method has correct attributes", {

  coefficients <- coef(foo)
  expect_length(coefficients, 4L)
  expect_named(coefficients, c("a", "b", "Direct", "Total"))

})

test_that("coef() method returns correct values of coefficients", {

  expect_equivalent(coef(foo, parm = "a"), foo$a)
  expect_equivalent(coef(foo, parm = "b"), foo$b)
  expect_equivalent(coef(foo, parm = "Direct"), foo$direct)
  expect_equivalent(coef(foo, parm = "Total"), foo$total)

})

test_that("summary returns original object", {
  expect_identical(foo, bar)
})

test_that("object returned by setup_ellipse_plot() has correct structure", {

  # check data frame for data to be plotted
  expect_s3_class(ellipse_mx$data, "data.frame")
  expect_s3_class(ellipse_ym$data, "data.frame")
  expect_s3_class(ellipse_partial$data, "data.frame")
  # check dimensions
  expect_identical(dim(ellipse_mx$data), c(as.integer(n), 2L))
  expect_identical(dim(ellipse_ym$data), c(as.integer(n), 2L))
  expect_identical(dim(ellipse_partial$data), c(as.integer(n), 2L))
  # check column names
  column_names <- c("x", "y")
  expect_named(ellipse_mx$data, column_names)
  expect_named(ellipse_ym$data, column_names)
  expect_named(ellipse_partial$data, column_names)

  # check data frame for ellipse
  expect_s3_class(ellipse_mx$ellipse, "data.frame")
  expect_s3_class(ellipse_ym$ellipse, "data.frame")
  expect_s3_class(ellipse_partial$ellipse, "data.frame")
  # check dimensions
  expect_identical(ncol(ellipse_mx$ellipse), 2L)
  expect_gt(nrow(ellipse_mx$ellipse), 0L)
  expect_identical(ncol(ellipse_ym$ellipse), 2L)
  expect_gt(nrow(ellipse_ym$ellipse), 0L)
  expect_identical(ncol(ellipse_partial$ellipse), 2L)
  expect_gt(nrow(ellipse_partial$ellipse), 0L)
  # check column names
  column_names <- c("x", "y")
  expect_named(ellipse_mx$ellipse, column_names)
  expect_named(ellipse_ym$ellipse, column_names)
  expect_named(ellipse_partial$ellipse, column_names)

  # check data frame for line representing the coefficient
  expect_s3_class(ellipse_mx$line, "data.frame")
  expect_null(ellipse_ym$line)
  expect_s3_class(ellipse_partial$line, "data.frame")
  # check dimensions
  expect_identical(dim(ellipse_mx$line), c(1L, 2L))
  expect_identical(dim(ellipse_partial$line), c(1L, 2L))
  # check column names
  column_names <- c("intercept", "slope")
  expect_named(ellipse_mx$line, column_names)
  expect_named(ellipse_partial$line, column_names)
  # check if intercept is 0 for partial residuals
  expect_identical(ellipse_partial$line$intercept, 0)

  # check if variables are passed correctly
  expect_identical(ellipse_mx$horizontal, "X")
  expect_identical(ellipse_mx$vertical, "M1")
  expect_identical(ellipse_ym$horizontal, "M1")
  expect_identical(ellipse_ym$vertical, "Y")
  expect_identical(ellipse_partial$horizontal, "M1")
  expect_identical(ellipse_partial$vertical, "Y")

  # check logical for partial residuals on the vertical axis
  expect_false(ellipse_mx$partial)
  expect_false(ellipse_ym$partial)
  expect_true(ellipse_partial$partial)

  # check logical for robust method
  expect_false(ellipse_mx$robust)
  expect_false(ellipse_ym$robust)
  expect_false(ellipse_partial$robust)

  # check logical for multiple methods
  expect_false(ellipse_mx$have_methods)
  expect_false(ellipse_ym$have_methods)
  expect_false(ellipse_partial$have_methods)

})


## only implemented for simple mediation without covariates

test_that("covariates not implemented", {

  # run regression fit
  set.seed(seed)
  suppressWarnings(
    reg_fit <- fit_mediation(test_data, x = "X", y = "Y", m = "M1",
                             covariates = c("C1", "C2"), method = "regression",
                             robust = FALSE)
  )

  # try to run with covariates (should give warning)
  set.seed(seed)
  expect_warning(
    cov_fit <- fit_mediation(test_data, x = "X", y = "Y", m = "M1",
                             covariates = c("C1", "C2"), method = "covariance",
                             robust = FALSE)
  )

  # these should be the same
  expect_equal(cov_fit, reg_fit)

})

test_that("multiple mediators not implemented", {

  # run regression fit
  set.seed(seed)
  suppressWarnings(
    reg_fit <- fit_mediation(test_data, x = "X", y = "Y", m = c("M1", "M2"),
                             method = "regression", robust = FALSE)
  )

  # try to run with multiple mediators (should give warning)
  set.seed(seed)
  expect_warning(
    cov_fit <- fit_mediation(test_data, x = "X", y = "Y", m = c("M1", "M2"),
                             method = "covariance", robust = FALSE)
  )

  # these should be the same
  expect_equal(cov_fit, reg_fit)

})


# fit mediation model through formula interface with data argument
fit_f1 <- fit_mediation(Y ~ m(M1) + X, data = test_data,
                        method = "covariance", robust = FALSE)
# fit mediation model through formula interface without data argument
fit_f2 <- fit_mediation(Y ~ m(M1) + X,
                        method = "covariance", robust = FALSE)
# define mediator outside formula
med <- m(M1)
fit_f3 <- fit_mediation(Y ~ med + X, data = test_data,
                        method = "covariance", robust = FALSE)


test_that("formula interface works correctly", {

  # check that results are the same as with default method
  expect_equal(fit_f1, foo)
  expect_equal(fit_f2, foo)
  expect_equal(fit_f3, foo)

})
