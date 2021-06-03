context("robust regression fit: one X, one M, no covariates")


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

## fit mediation model and compute summary
set.seed(seed)
foo <- fit_mediation(test_data, x = "X", y = "Y", m = "M",
                     method = "regression", robust = TRUE,
                     efficiency = 0.95)
bar <- summary(foo)

## create data for plotting
ellipse_mx <- setup_ellipse_plot(foo)
ellipse_ym <- setup_ellipse_plot(foo, horizontal = "M", vertical = "Y",
                                 partial = FALSE)
ellipse_partial <- setup_ellipse_plot(foo, horizontal = "M", vertical = "Y",
                                      partial = TRUE)
weight_default <- setup_weight_plot(foo)
weight_y <- setup_weight_plot(foo, outcome = "Y")


## run tests

test_that("output has correct structure", {

  # regression fit
  expect_s3_class(foo, "reg_fit_mediation")
  expect_s3_class(foo, "fit_mediation")
  # individual regressions
  expect_s3_class(foo$fit_mx, "lmrob")
  expect_s3_class(foo$fit_ymx, "lmrob")
  expect_null(foo$fit_yx)

})

test_that("arguments are correctly passed", {

  # variable names
  expect_identical(foo$x, "X")
  expect_identical(foo$y, "Y")
  expect_identical(foo$m, "M")
  expect_identical(foo$covariates, character())
  # robust fit
  expect_identical(foo$robust, "MM")
  expect_identical(foo$family, "gaussian")
  expect_equal(foo$control, reg_control(efficiency = 0.95))
  # no contrasts
  expect_false(foo$contrast)

})

test_that("dimensions are correct", {

  # effects are scalars
  expect_length(foo$a, 1L)
  expect_length(foo$b, 1L)
  expect_length(foo$direct, 1L)
  expect_length(foo$total, 1L)
  expect_length(foo$ab, 1L)
  # individual regressions
  expect_length(coef(foo$fit_mx), 2L)
  expect_length(coef(foo$fit_ymx), 3L)
  # dimensions of data
  expect_identical(dim(foo$data), c(as.integer(n), 3L))

})

test_that("values of coefficients are correct", {

  expect_equivalent(foo$a, coef(foo$fit_mx)["X"])
  expect_equivalent(foo$b, coef(foo$fit_ymx)["M"])
  expect_equivalent(foo$direct, coef(foo$fit_ymx)["X"])
  expect_equivalent(foo$total, foo$a * foo$b + foo$direct)
  expect_equivalent(foo$ab, foo$a * foo$b)

})

test_that("output of coef() method has correct attributes", {

  coefficients <- coef(foo)
  expect_length(coefficients, 5L)
  expect_named(coefficients, c("a", "b", "Direct", "Total", "ab"))

})

test_that("coef() method returns correct values of coefficients", {

  expect_equivalent(coef(foo, parm = "a"), foo$a)
  expect_equivalent(coef(foo, parm = "b"), foo$b)
  expect_equivalent(coef(foo, parm = "Direct"), foo$direct)
  expect_equivalent(coef(foo, parm = "Total"), foo$total)
  expect_equivalent(coef(foo, parm = "ab"), foo$ab)

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
  expect_identical(dim(ellipse_mx$data), c(as.integer(n), 3L))
  expect_identical(dim(ellipse_ym$data), c(as.integer(n), 3L))
  expect_identical(dim(ellipse_partial$data), c(as.integer(n), 3L))
  # check column names
  column_names <- c("x", "y", "Weight")
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
  expect_identical(ellipse_mx$vertical, "M")
  expect_identical(ellipse_ym$horizontal, "M")
  expect_identical(ellipse_ym$vertical, "Y")
  expect_identical(ellipse_partial$horizontal, "M")
  expect_identical(ellipse_partial$vertical, "Y")

  # check logical for partial residuals on the vertical axis
  expect_false(ellipse_mx$partial)
  expect_false(ellipse_ym$partial)
  expect_true(ellipse_partial$partial)

  # check logical for robust method
  expect_true(ellipse_mx$robust)
  expect_true(ellipse_ym$robust)
  expect_true(ellipse_partial$robust)

  # check logical for multiple methods
  expect_false(ellipse_mx$have_methods)
  expect_false(ellipse_ym$have_methods)
  expect_false(ellipse_partial$have_methods)

})

test_that("object returned by setup_weight_plot() has correct structure", {

  # check data frame for weight percentages to be plotted
  expect_s3_class(weight_default$data, "data.frame")
  expect_s3_class(weight_y$data, "data.frame")
  # check dimensions
  expect_identical(ncol(weight_default$data), 5L)
  expect_identical(ncol(weight_y$data), 4L)
  # check column names
  column_names <- c("Outcome", "Tail", "Weights", "Threshold", "Percentage")
  expect_named(weight_default$data, column_names)
  expect_named(weight_y$data, column_names[-1])

  # check if variables are passed correctly
  expect_identical(weight_default$outcome, c("M", "Y"))
  expect_identical(weight_y$outcome, "Y")

})


# fit mediation model through formula interface with data argument
set.seed(seed)
fit_f1 <- fit_mediation(Y ~ m(M) + X, data = test_data,
                        method = "regression", robust = TRUE,
                        efficiency = 0.95)
# fit mediation model through formula interface without data argument
set.seed(seed)
fit_f2 <- fit_mediation(Y ~ m(M) + X,
                        method = "regression", robust = TRUE,
                        efficiency = 0.95)
# define mediator outside formula
med <- m(M)
set.seed(seed)
fit_f3 <- fit_mediation(Y ~ med + X, data = test_data,
                        method = "regression", robust = TRUE,
                        efficiency = 0.95)


test_that("formula interface works correctly", {

  # check that results are the same as with default method
  expect_equal(fit_f1, foo)
  expect_equal(fit_f2, foo)
  expect_equal(fit_f3, foo)

})
