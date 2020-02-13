context("standard regression fit: single mediator, covariates")


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
C1 <- rnorm(n)
C2 <- rnorm(n)
test_data <- data.frame(X, Y, M, C1, C2)

## fit mediation model and compute summary
foo <- fit_mediation(test_data, x = "X", y = "Y", m = "M",
                     covariates = c("C1", "C2"), method = "regression",
                     robust = FALSE)
bar <- summary(foo)

## create data for plotting
ellipse_default <- setup_ellipse_plot(foo)
ellipse_partial <- setup_ellipse_plot(foo, horizontal = "M", vertical = "Y",
                                      partial = TRUE)


## run tests

test_that("output has correct structure", {

  # regression fit
  expect_s3_class(foo, "reg_fit_mediation")
  expect_s3_class(foo, "fit_mediation")
  # individual regressions
  expect_s3_class(foo$fit_mx, "lm")
  expect_s3_class(foo$fit_ymx, "lm")
  expect_s3_class(foo$fit_yx, "lm")

})

test_that("arguments are correctly passed", {

  # variable names
  expect_identical(foo$x, "X")
  expect_identical(foo$y, "Y")
  expect_identical(foo$m, "M")
  expect_identical(foo$covariates, c("C1", "C2"))
  # robust fit
  expect_false(foo$robust)
  expect_null(foo$control)

})

test_that("dimensions are correct", {

  # effects are scalars
  expect_length(foo$a, 1L)
  expect_length(foo$b, 1L)
  expect_length(foo$direct, 1L)
  expect_length(foo$total, 1L)
  # individual regressions
  expect_length(coef(foo$fit_mx), 4L)
  expect_length(coef(foo$fit_ymx), 5L)
  expect_length(coef(foo$fit_yx), 4L)
  # dimensions of data
  expect_identical(dim(foo$data), c(as.integer(n), 5L))

})

test_that("values of coefficients are correct", {

  expect_equivalent(foo$a, coef(foo$fit_mx)["X"])
  expect_equivalent(foo$b, coef(foo$fit_ymx)["M"])
  expect_equivalent(foo$direct, coef(foo$fit_ymx)["X"])
  expect_equivalent(foo$total, coef(foo$fit_yx)["X"])
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
  expect_s3_class(ellipse_default$data, "data.frame")
  expect_s3_class(ellipse_partial$data, "data.frame")
  # check dimensions
  expect_identical(dim(ellipse_default$data), c(as.integer(n), 2L))
  expect_identical(dim(ellipse_partial$data), c(as.integer(n), 2L))
  # check column names
  column_names <- c("x", "y")
  expect_named(ellipse_default$data, column_names)
  expect_named(ellipse_partial$data, column_names)

  # check data frame for ellipse
  expect_s3_class(ellipse_default$ellipse, "data.frame")
  expect_s3_class(ellipse_partial$ellipse, "data.frame")
  # check dimensions
  expect_identical(ncol(ellipse_default$ellipse), 2L)
  expect_gt(nrow(ellipse_default$ellipse), 0L)
  expect_identical(ncol(ellipse_partial$ellipse), 2L)
  expect_gt(nrow(ellipse_partial$ellipse), 0L)
  # check column names
  column_names <- c("x", "y")
  expect_named(ellipse_default$ellipse, column_names)
  expect_named(ellipse_partial$ellipse, column_names)

  # check data frame for line representing the coefficient
  expect_null(ellipse_default$line)
  expect_s3_class(ellipse_partial$line, "data.frame")
  # check dimensions
  expect_identical(dim(ellipse_partial$line), c(1L, 2L))
  # check column names
  column_names <- c("intercept", "slope")
  expect_named(ellipse_partial$line, column_names)
  # check if intercept is 0 for partial residuals
  expect_identical(ellipse_partial$line$intercept, 0)

  # check if variables are passed correctly
  expect_identical(ellipse_default$horizontal, "X")
  expect_identical(ellipse_default$vertical, "M")
  expect_identical(ellipse_partial$horizontal, "M")
  expect_identical(ellipse_partial$vertical, "Y")

  # check logical for partial residuals on the vertical axis
  expect_false(ellipse_default$partial)
  expect_true(ellipse_partial$partial)

  # check logical for robust method
  expect_false(ellipse_default$robust)
  expect_false(ellipse_partial$robust)

  # check logical for multiple methods
  expect_false(ellipse_default$have_methods)
  expect_false(ellipse_partial$have_methods)

})


# fit mediation model through formula interface with data argument
fit_f1 <- fit_mediation(Y ~ m(M) + X + covariates(C1, C2), data = test_data,
                        method = "regression", robust = FALSE)
# fit mediation model through formula interface without data argument
fit_f2 <- fit_mediation(Y ~ m(M) + X + covariates(C1, C2),
                        method = "regression", robust = FALSE)
# define mediator and covariates outside formula
med <- m(M)
cov <- covariates(C1, C2)
fit_f3 <- fit_mediation(Y ~ med + X + cov, data = test_data,
                        method = "regression", robust = FALSE)


test_that("formula interface works correctly", {

  # check that results are the same as with default method
  expect_equal(fit_f1, foo)
  expect_equal(fit_f2, foo)
  expect_equal(fit_f3, foo)

})
