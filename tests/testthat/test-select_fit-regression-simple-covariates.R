context("regression fit with skewed errors: single mediator, covariates")


## load package
library("robmed", quietly = TRUE)
library("sn", quietly = TRUE)

## control parameters
n <- 250            # number of observations
a <- c <- 0.2       # true effects
b <- 0              # true effect
seed <- 20190201    # seed for the random number generator

## set seed for reproducibility
set.seed(seed)

## generate data
X <- rnorm(n)
M <- a * X + rsn(n, alpha = 10)
Y <- b * M + c * X + rsn(n, alpha = 10)
C1 <- rnorm(n)
C2 <- rnorm(n)
test_data <- data.frame(X, Y, M, C1, C2)

## fit mediation model and compute summary
foo <- fit_mediation(test_data, x = "X", y = "Y", m = "M",
                     covariates = c("C1", "C2"), method = "regression",
                     robust = FALSE, family = "select")
bar <- summary(foo)


## run tests

test_that("output has correct structure", {

  # regression fit
  expect_s3_class(foo, "reg_fit_mediation")
  expect_s3_class(foo, "fit_mediation")
  # individual regressions
  expect_s3_class(foo$fit_mx, "lmse")
  expect_s3_class(foo$fit_ymx, "lmse")
  expect_s3_class(foo$fit_yx, "lmse")

})

test_that("arguments are correctly passed", {

  # variable names
  expect_identical(foo$x, "X")
  expect_identical(foo$y, "Y")
  expect_identical(foo$m, "M")
  expect_identical(foo$covariates, c("C1", "C2"))
  # robust fit
  expect_false(foo$robust)
  expect_identical(foo$family, "select")
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
  expect_equivalent(foo$total, coef(foo$fit_yx)["X"])

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

  # not meaningful for median regression
  expect_error(setup_ellipse_plot(foo))

})


# fit mediation model through formula interface with data argument
fit_f1 <- fit_mediation(Y ~ m(M) + X + covariates(C1, C2), data = test_data,
                        method = "regression", robust = FALSE,
                        family = "select")
# fit mediation model through formula interface without data argument
fit_f2 <- fit_mediation(Y ~ m(M) + X + covariates(C1, C2),
                        method = "regression", robust = FALSE,
                        family = "select")
# define mediator and covariates outside formula
med <- m(M)
cov <- covariates(C1, C2)
fit_f3 <- fit_mediation(Y ~ med + X + cov, data = test_data,
                        method = "regression", robust = FALSE,
                        family = "select")


test_that("formula interface works correctly", {

  # check that results are the same as with default method
  expect_equal(fit_f1, foo)
  expect_equal(fit_f2, foo)
  expect_equal(fit_f3, foo)

})
