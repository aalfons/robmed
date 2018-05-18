context("robust bootstrap test: covariance")


## load package
library("robmed", quietly = TRUE)

## control parameters
n <- 250             # number of observations
a <- c <- 0.2        # true effects
b <- 0               # true effect
level <- 0.99        # confidence level
today <- 20180517    # seed for the random number generator
R <- 1000            # number of bootstrap samples

## set seed for reproducibility
set.seed(today)

## generate data
X <- rnorm(n)
M1 <- a * X + rnorm(n)
M2 <- rnorm(n)
Y <- b * M1 + c * X + rnorm(n)
C1 <- rnorm(n)
C2 <- rnorm(n)
test_data <- data.frame(X, Y, M1, M2, C1, C2)

## run bootstrap test and compute summary
foo <- test_mediation(test_data, x = "X", y = "Y", m = "M1",
                      test = "boot", R = R, level = level, type = "bca",
                      method = "covariance", robust = TRUE)
bar <- summary(foo)


## run tests

test_that("output has correct structure", {

  # bootstrap test
  expect_s3_class(foo, "boot_test_mediation")
  expect_s3_class(foo, "test_mediation")
  # regression fit
  expect_s3_class(foo$fit, "cov_fit_mediation")
  expect_s3_class(foo$fit, "fit_mediation")
  # individual robust regressions
  expect_s3_class(foo$fit$cov, "cov_Huber")
  # bootstrap replicates
  expect_s3_class(foo$reps, "boot")

})

test_that("arguments are correctly passed", {

  # number of bootstrap replicates
  expect_identical(foo$R, as.integer(R))  # doesn't hold for too many outliers
  # confidence level
  expect_identical(foo$level, level)
  # type of confidence intervals
  expect_identical(foo$type, "bca")
  # variable names
  expect_identical(foo$fit$x, "X")
  expect_identical(foo$fit$y, "Y")
  expect_identical(foo$fit$m, "M1")
  expect_null(foo$fit$covariates)
  # robust fit and test
  expect_true(foo$fit$robust)

})

test_that("dimensions are correct", {

  # effects are scalars
  expect_length(foo$ab, 1L)
  expect_length(foo$fit$a, 1L)
  expect_length(foo$fit$b, 1L)
  expect_length(foo$fit$c, 1L)
  expect_length(foo$fit$c_prime, 1L)
  # only one indirect effect, so only one confidence interval
  expect_length(foo$ci, 2L)
  # dimensions of bootstrap replicates
  d_boot <- dim(foo$reps$t)
  expect_identical(d_boot, c(as.integer(R), 5L))
  # dimensions of data
  d_data <- dim(foo$fit$data)
  expect_identical(d_data, c(as.integer(n), 3L))

})

test_that("values of coefficients are correct", {

  expect_equivalent(foo$ab, mean(foo$reps$t[, 1]))
  expect_equivalent(foo$fit$c_prime, foo$fit$a * foo$fit$b + foo$fit$c)

})

test_that("output of coef() method has correct attributes", {

  coef_boot <- coef(foo, type = "boot")
  coef_data <- coef(foo, type = "data")
  coef_names <- c("a", "b", "c", "c'", "ab")
  # bootstrapped effects
  expect_length(coef_boot, 5L)
  expect_named(coef_boot, coef_names)
  # effects computed on original sample
  expect_length(coef_data, 5L)
  expect_named(coef_data, coef_names)

})

test_that("coef() method returns correct values of coefficients", {

  # bootstrapped effects
  expect_equivalent(coef(foo, parm = "a", type = "boot"), mean(foo$reps$t[, 2]))
  expect_equivalent(coef(foo, parm = "b", type = "boot"), mean(foo$reps$t[, 3]))
  expect_equivalent(coef(foo, parm = "c", type = "boot"), mean(foo$reps$t[, 4]))
  expect_equivalent(coef(foo, parm = "c'", type = "boot"), mean(foo$reps$t[, 5]))
  expect_equivalent(coef(foo, parm = "ab", type = "boot"), foo$ab)

  # effects computed on original sample
  expect_equivalent(coef(foo, parm = "a", type = "data"), foo$fit$a)
  expect_equivalent(coef(foo, parm = "b", type = "data"), foo$fit$b)
  expect_equivalent(coef(foo, parm = "c", type = "data"), foo$fit$c)
  expect_equivalent(coef(foo, parm = "c'", type = "data"), foo$fit$c_prime)
  expect_equivalent(coef(foo, parm = "ab", type = "data"), foo$fit$a * foo$fit$b)

})

test_that("summary has correct structure", {

  # summary
  expect_s3_class(bar, "summary_test_mediation")
  # original output of test for indirect effect
  expect_identical(bar$object, foo)
  # summary of the model fit
  expect_s3_class(bar$summary, "summary_fit_mediation")
  # regression standard error for model y ~ m + x
  expect_type(bar$summary$s, "list")
  expect_named(bar$summary$s, "value")

})

test_that("attributes are correctly passed through summary", {
  # robustness
  expect_true(bar$summary$robust)
  # number of observations
  expect_identical(bar$summary$n, as.integer(n))
  # variables
  expect_identical(bar$summary$variables, c("X", "Y", "M1"))
})

test_that("effect summaries have correct names", {

  # a path
  expect_identical(dim(bar$summary$a), c(1L, 5L))
  expect_identical(rownames(bar$summary$a), "X")
  expect_identical(colnames(bar$summary$a)[1:2], c("Data", "Boot"))
  # b path
  expect_identical(dim(bar$summary$b), c(1L, 5L))
  expect_identical(rownames(bar$summary$b), "M1")
  expect_identical(colnames(bar$summary$b)[1:2], c("Data", "Boot"))
  # c path
  expect_identical(dim(bar$summary$c), c(1L, 5L))
  expect_identical(rownames(bar$summary$c), "X")
  expect_identical(colnames(bar$summary$c)[1:2], c("Data", "Boot"))
  # c' path
  expect_identical(dim(bar$summary$c_prime), c(1L, 5L))
  expect_identical(rownames(bar$summary$c_prime), "X")
  expect_identical(colnames(bar$summary$c_prime)[1:2], c("Data", "Boot"))

})

test_that("effect summaries contain correct coefficient values", {

  # effects computed on original sample
  expect_identical(bar$summary$a["X", "Data"], foo$fit$a)
  expect_identical(bar$summary$b["M1", "Data"], foo$fit$b)
  expect_identical(bar$summary$c["X", "Data"], foo$fit$c)
  expect_identical(bar$summary$c_prime["X", "Data"], foo$fit$c_prime)

  # bootstrapped effects
  expect_identical(bar$summary$a["X", "Boot"], mean(foo$reps$t[, 2]))
  expect_identical(bar$summary$b["M1", "Boot"], mean(foo$reps$t[, 3]))
  expect_identical(bar$summary$c["X", "Boot"], mean(foo$reps$t[, 4]))
  expect_identical(bar$summary$c_prime["X", "Boot"], mean(foo$reps$t[, 5]))

})


## only implemented for simple mediation without covariates

test_that("covariates not implemented", {

  # run test with regression method
  set.seed(today)
  test_reg <- test_mediation(test_data, x = "X", y = "Y", m = "M1",
                             covariates = c("C1", "C2"), test = "boot",
                             R = R, level = level, type = "bca",
                             method = "regression", robust = TRUE)

  # try to run test with covariates (should give warning)
  set.seed(today)
  expect_warning(
    test_cov <- test_mediation(test_data, x = "X", y = "Y", m = "M1",
                               covariates = c("C1", "C2"), test = "boot",
                               R = R, level = level, type = "bca",
                               method = "covariance", robust = TRUE)
  )

  # these should be identical
  expect_equal(test_cov, test_reg)

})


test_that("multiple mediators not implemented", {

  # run test with regression method
  set.seed(today)
  test_reg <- test_mediation(test_data, x = "X", y = "Y", m = c("M1", "M2"),
                             test = "boot", R = R, level = level, type = "bca",
                             method = "regression", robust = TRUE)

  # try to run test with covariates (should give warning)
  set.seed(today)
  expect_warning(
    test_cov <- test_mediation(test_data, x = "X", y = "Y", m = c("M1", "M2"),
                               test = "boot", R = R, level = level, type = "bca",
                               method = "covariance", robust = TRUE)
  )

  # these should be identical
  expect_equal(test_cov, test_reg)

})
