context("robust regression fit: two X, one M, no covariates")


## load package
library("robmed", quietly = TRUE)

## control parameters
n <- 250               # number of observations
a <- c <- c(0.2, 0.2)  # true effects
b <- 0                 # true effect
seed <- 20150601       # seed for the random number generator

## set seed for reproducibility
set.seed(seed)

## generate data
X1 <- rnorm(n)
X2 <- rnorm(n)
M <- a[1] * X1 + a[2] * X2 + rnorm(n)
Y <- b * M + c[1] * X1 + c[2] * X2 + rnorm(n)
C1 <- rnorm(n)
C2 <- rnorm(n)
test_data <- data.frame(X1, X2, Y, M, C1, C2)

## fit mediation model and compute summary
set.seed(seed)
foo <- fit_mediation(test_data, x = c("X1", "X2"), y = "Y", m = "M",
                     covariates = c("C1", "C2"), method = "regression",
                     robust = TRUE, efficiency = 0.95)
bar <- summary(foo)

## create data for plotting
ellipse_default <- setup_ellipse_plot(foo)
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
  expect_identical(foo$x, c("X1", "X2"))
  expect_identical(foo$y, "Y")
  expect_identical(foo$m, "M")
  expect_identical(foo$covariates, c("C1", "C2"))
  # robust fit
  expect_identical(foo$robust, "MM")
  expect_identical(foo$family, "gaussian")
  expect_equal(foo$control, reg_control(efficiency = 0.95))
  # no contrasts
  expect_false(foo$contrast)

})

test_that("dimensions are correct", {

  # effects are scalars
  expect_length(foo$a, 2L)
  expect_length(foo$b, 1L)
  expect_length(foo$direct, 2L)
  expect_length(foo$total, 2L)
  expect_length(foo$ab, 3L)
  # individual regressions
  expect_length(coef(foo$fit_mx), 5L)
  expect_length(coef(foo$fit_ymx), 6L)
  # dimensions of data
  expect_identical(dim(foo$data), c(as.integer(n), 6L))

})

test_that("values of coefficients are correct", {

  # extract correct values
  a <- coef(foo$fit_mx)[c("X1", "X2")]
  b <- unname(coef(foo$fit_ymx)["M"])
  direct <- coef(foo$fit_ymx)[c("X1", "X2")]
  ab <- a * b
  sum_ab <- sum(ab)
  # compare with stored values
  expect_equal(foo$a, a)
  expect_equal(foo$b, b)
  expect_equal(foo$direct, direct)
  expect_equal(foo$total, ab + direct)
  expect_equal(foo$ab, c(Total = sum_ab, ab))

})

test_that("output of coef() method has correct attributes", {

  coefficients <- coef(foo)
  expect_length(coefficients, 10L)
  expect_named(coefficients,
               c("a_X1", "a_X2", "b", "Direct_X1", "Direct_X2",
                 "Total_X1", "Total_X2", "ab_Total", "ab_X1", "ab_X2"))

})

test_that("coef() method returns correct values of coefficients", {

  expect_equivalent(coef(foo, parm = "a_X1"), foo$a["X1"])
  expect_equivalent(coef(foo, parm = "a_X2"), foo$a["X2"])
  expect_equivalent(coef(foo, parm = "b"), foo$b)
  expect_equivalent(coef(foo, parm = "Direct_X1"), foo$direct["X1"])
  expect_equivalent(coef(foo, parm = "Direct_X2"), foo$direct["X2"])
  expect_equivalent(coef(foo, parm = "Total_X1"), foo$total["X1"])
  expect_equivalent(coef(foo, parm = "Total_X2"), foo$total["X2"])
  expect_equivalent(coef(foo, parm = "ab_Total"), foo$ab["Total"])
  expect_equivalent(coef(foo, parm = "ab_X1"), foo$ab["X1"])
  expect_equivalent(coef(foo, parm = "ab_X2"), foo$ab["X2"])

})

test_that("summary returns original object", {
  expect_identical(foo, bar)
})

test_that("object returned by setup_ellipse_plot() has correct structure", {

  # check data frame for data to be plotted
  expect_s3_class(ellipse_default$data, "data.frame")
  expect_s3_class(ellipse_partial$data, "data.frame")
  # check dimensions
  expect_identical(dim(ellipse_default$data), c(as.integer(n), 3L))
  expect_identical(dim(ellipse_partial$data), c(as.integer(n), 3L))
  # check column names
  column_names <- c("x", "y", "Weight")
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
  expect_identical(ellipse_default$horizontal, "X1")
  expect_identical(ellipse_default$vertical, "M")
  expect_identical(ellipse_partial$horizontal, "M")
  expect_identical(ellipse_partial$vertical, "Y")

  # check logical for partial residuals on the vertical axis
  expect_false(ellipse_default$partial)
  expect_true(ellipse_partial$partial)

  # check logical for robust method
  expect_true(ellipse_default$robust)
  expect_true(ellipse_partial$robust)

  # check logical for multiple methods
  expect_false(ellipse_default$have_methods)
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
fit_f1 <- fit_mediation(Y ~ m(M) + X1 + X2 + covariates(C1, C2),
                        data = test_data, method = "regression",
                        robust = TRUE, efficiency = 0.95)
# fit mediation model through formula interface without data argument
set.seed(seed)
fit_f2 <- fit_mediation(Y ~ m(M) + X1 + X2 + covariates(C1, C2),
                        method = "regression", robust = TRUE,
                        efficiency = 0.95)
# define mediators and covariates outside formula
med <- m(M)
cov <- covariates(C1, C2)
set.seed(seed)
fit_f3 <- fit_mediation(Y ~ med + X1 + X2 + cov, data = test_data,
                        method = "regression", robust = TRUE,
                        efficiency = 0.95)


test_that("formula interface works correctly", {

  # check that results are the same as with default method
  expect_equal(fit_f1, foo)
  expect_equal(fit_f2, foo)
  expect_equal(fit_f3, foo)

})
