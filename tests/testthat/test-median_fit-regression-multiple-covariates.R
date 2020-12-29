context("median regression fit: multiple mediators, covariates")


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
M1 <- a * X + rnorm(n)
M2 <- rnorm(n)
Y <- b * M1 + c * X + rnorm(n)
C1 <- rnorm(n)
C2 <- rnorm(n)
test_data <- data.frame(X, Y, M1, M2, C1, C2)

## fit mediation model and compute summary
foo <- fit_mediation(test_data, x = "X", y = "Y", m = c("M1", "M2"),
                     covariates = c("C1", "C2"), method = "regression",
                     robust = "median")
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
  expect_s3_class(foo$fit_mx$M1, "rq")
  expect_s3_class(foo$fit_mx$M2, "rq")
  expect_s3_class(foo$fit_ymx, "rq")
  expect_null(foo$fit_yx)

})

test_that("arguments are correctly passed", {

  # variable names
  expect_identical(foo$x, "X")
  expect_identical(foo$y, "Y")
  expect_identical(foo$m, c("M1", "M2"))
  expect_identical(foo$covariates, c("C1", "C2"))
  # robust fit
  expect_identical(foo$robust, "median")
  expect_identical(foo$family, "gaussian")
  expect_null(foo$control)

})

test_that("dimensions are correct", {

  # effects are scalars
  expect_length(foo$a, 2L)
  expect_length(foo$b, 2L)
  expect_length(foo$direct, 1L)
  expect_length(foo$total, 1L)
  expect_length(foo$ab, 3L)
  # individual regressions
  expect_length(coef(foo$fit_mx$M1), 4L)
  expect_length(coef(foo$fit_mx$M2), 4L)
  expect_length(coef(foo$fit_ymx), 6L)
  # dimensions of data
  expect_identical(dim(foo$data), c(as.integer(n), 6L))

})

test_that("values of coefficients are correct", {

  # extract correct values
  a <- c(M1 = unname(coef(foo$fit_mx$M1)["X"]),
         M2 = unname(coef(foo$fit_mx$M2)["X"]))
  b <- coef(foo$fit_ymx)[c("M1", "M2")]
  direct <- coef(foo$fit_ymx)["X"]
  ab <- a * b
  sum_ab <- sum(ab)
  ab <- c(Total = sum_ab, ab)
  # compare with stored values
  expect_equal(foo$a, a)
  expect_equal(foo$b, b)
  expect_equivalent(foo$direct, direct)
  expect_equivalent(foo$total, sum_ab + direct)
  expect_equal(foo$ab, ab)

})

test_that("output of coef() method has correct attributes", {

  coefficients <- coef(foo)
  expect_length(coefficients, 9L)
  expect_named(coefficients,
               c("a_M1", "a_M2", "b_M1", "b_M2", "Direct", "Total",
                 "ab_Total", "ab_M1", "ab_M2"))

})

test_that("coef() method returns correct values of coefficients", {

  expect_equivalent(coef(foo, parm = "a_M1"), foo$a["M1"])
  expect_equivalent(coef(foo, parm = "a_M2"), foo$a["M2"])
  expect_equivalent(coef(foo, parm = "b_M1"), foo$b["M1"])
  expect_equivalent(coef(foo, parm = "b_M2"), foo$b["M2"])
  expect_equivalent(coef(foo, parm = "Direct"), foo$direct)
  expect_equivalent(coef(foo, parm = "Total"), foo$total)
  expect_equivalent(coef(foo, parm = "ab_Total"), foo$ab["Total"])
  expect_equivalent(coef(foo, parm = "ab_M1"), foo$ab["M1"])
  expect_equivalent(coef(foo, parm = "ab_M2"), foo$ab["M2"])

})

test_that("summary returns original object", {
  expect_identical(foo, bar)
})

test_that("object returned by setup_xxx_plot() has correct structure", {

  # not meaningful for median regression
  expect_error(setup_ellipse_plot(foo))
  expect_error(setup_weight_plot(foo))

})


# fit mediation model through formula interface with data argument
fit_f1 <- fit_mediation(Y ~ m(M1, M2) + X + covariates(C1, C2),
                        data = test_data, method = "regression",
                        robust = "median")
# fit mediation model through formula interface without data argument
fit_f2 <- fit_mediation(Y ~ m(M1, M2) + X + covariates(C1, C2),
                        method = "regression", robust = "median")
# define mediators and covariates outside formula
med <- m(M1, M2)
cov <- covariates(C1, C2)
fit_f3 <- fit_mediation(Y ~ med + X + cov, data = test_data,
                        method = "regression", robust = "median")


test_that("formula interface works correctly", {

  # check that results are the same as with default method
  expect_equal(fit_f1, foo)
  expect_equal(fit_f2, foo)
  expect_equal(fit_f3, foo)

})
