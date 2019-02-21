context("standard regression fit: multiple mediators, no covariates")


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
test_data <- data.frame(X, Y, M1, M2)

## fit mediation model and compute summary
foo <- fit_mediation(test_data, x = "X", y = "Y", m = c("M1", "M2"),
                     method = "regression", robust = FALSE)
bar <- summary(foo)


## run tests

test_that("output has correct structure", {

  # regression fit
  expect_s3_class(foo, "reg_fit_mediation")
  expect_s3_class(foo, "fit_mediation")
  # individual regressions
  expect_type(foo$fit_mx, "list")
  expect_length(foo$fit_mx, 2L)
  expect_named(foo$fit_mx, c("M1", "M2"))
  expect_s3_class(foo$fit_mx$M1, "lm")
  expect_s3_class(foo$fit_mx$M2, "lm")
  expect_s3_class(foo$fit_ymx, "lm")
  expect_s3_class(foo$fit_yx, "lm")

})

test_that("arguments are correctly passed", {

  # variable names
  expect_identical(foo$x, "X")
  expect_identical(foo$y, "Y")
  expect_identical(foo$m, c("M1", "M2"))
  expect_identical(foo$covariates, character())
  # robust fit
  expect_false(foo$robust)
  expect_false(foo$median)
  expect_null(foo$control)

})

test_that("dimensions are correct", {

  # effects are scalars
  expect_length(foo$a, 2L)
  expect_length(foo$b, 2L)
  expect_length(foo$c, 1L)
  expect_length(foo$c_prime, 1L)
  # individual regressions
  expect_length(coef(foo$fit_mx$M1), 2L)
  expect_length(coef(foo$fit_mx$M2), 2L)
  expect_length(coef(foo$fit_ymx), 4L)
  expect_length(coef(foo$fit_yx), 2L)
  # dimensions of data
  expect_identical(dim(foo$data), c(as.integer(n), 4L))

})

test_that("values of coefficients are correct", {

  a <- c(coef(foo$fit_mx$M1)["X"], coef(foo$fit_mx$M2)["X"])
  expect_equivalent(foo$a, a)
  expect_named(foo$a, c("M1", "M2"))
  expect_equal(foo$b, coef(foo$fit_ymx)[c("M1", "M2")])  # also checks names
  expect_equivalent(foo$c, coef(foo$fit_ymx)["X"])
  expect_equivalent(foo$c_prime, coef(foo$fit_yx)["X"])
  expect_equivalent(foo$c_prime, sum(foo$a * foo$b) + foo$c)

})

test_that("output of coef() method has correct attributes", {

  coefficients <- coef(foo)
  expect_length(coefficients, 6L)
  expect_named(coefficients, c("a_M1", "a_M2", "b_M1", "b_M2", "c", "c'"))

})

test_that("coef() method returns correct values of coefficients", {

  expect_equivalent(coef(foo, parm = "a_M1"), foo$a["M1"])
  expect_equivalent(coef(foo, parm = "a_M2"), foo$a["M2"])
  expect_equivalent(coef(foo, parm = "b_M1"), foo$b["M1"])
  expect_equivalent(coef(foo, parm = "b_M2"), foo$b["M2"])
  expect_equivalent(coef(foo, parm = "c"), foo$c)
  expect_equivalent(coef(foo, parm = "c'"), foo$c_prime)

})

test_that("summary returns original object", {
  expect_identical(foo, bar)
})
