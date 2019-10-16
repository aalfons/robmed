context("median regression fit: single mediator, covariates")


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

## fit mediation model and compute summary
foo <- fit_mediation(test_data, x = "X", y = "Y", m = "M",
                     covariates = c("C1", "C2"), method = "regression",
                     robust = TRUE, median = TRUE)
bar <- summary(foo)


## run tests

test_that("output has correct structure", {

  # regression fit
  expect_s3_class(foo, "reg_fit_mediation")
  expect_s3_class(foo, "fit_mediation")
  # individual regressions
  expect_s3_class(foo$fit_mx, "rq")
  expect_s3_class(foo$fit_ymx, "rq")
  expect_null(foo$fit$fit_yx)

})

test_that("arguments are correctly passed", {

  # variable names
  expect_identical(foo$x, "X")
  expect_identical(foo$y, "Y")
  expect_identical(foo$m, "M")
  expect_identical(foo$covariates, c("C1", "C2"))
  # robust fit
  expect_true(foo$robust)
  expect_true(foo$median)
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
  # dimensions of data
  expect_identical(dim(foo$data), c(as.integer(n), 5L))

})

test_that("values of coefficients are correct", {

  expect_equivalent(foo$a, coef(foo$fit_mx)["X"])
  expect_equivalent(foo$b, coef(foo$fit_ymx)["M"])
  expect_equivalent(foo$direct, coef(foo$fit_ymx)["X"])
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
