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
  expect_length(foo$c, 1L)
  expect_length(foo$c_prime, 1L)
  # dimensions of data
  expect_identical(dim(foo$data), c(as.integer(n), 3L))

})

test_that("values of coefficients are correct", {

  expect_equivalent(foo$c_prime, foo$a * foo$b + foo$c)

})

test_that("output of coef() method has correct attributes", {

  coefficients <- coef(foo)
  expect_length(coefficients, 4L)
  expect_named(coefficients, c("a", "b", "c", "c'"))

})

test_that("coef() method returns correct values of coefficients", {

  expect_equivalent(coef(foo, parm = "a"), foo$a)
  expect_equivalent(coef(foo, parm = "b"), foo$b)
  expect_equivalent(coef(foo, parm = "c"), foo$c)
  expect_equivalent(coef(foo, parm = "c'"), foo$c_prime)

})

test_that("summary returns original object", {
  expect_identical(foo, bar)
})


## only implemented for simple mediation without covariates

test_that("covariates not implemented", {

  # run regression fit
  set.seed(seed)
  reg_fit <- fit_mediation(test_data, x = "X", y = "Y", m = "M1",
                           covariates = c("C1", "C2"), method = "regression",
                           robust = FALSE)

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
  reg_fit <- fit_mediation(test_data, x = "X", y = "Y", m = c("M1", "M2"),
                           method = "regression", robust = FALSE)

  # try to run with multiple mediators (should give warning)
  set.seed(seed)
  expect_warning(
    cov_fit <- fit_mediation(test_data, x = "X", y = "Y", m = c("M1", "M2"),
                             method = "covariance", robust = FALSE)
  )

  # these should be the same
  expect_equal(cov_fit, reg_fit)

})
