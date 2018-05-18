context("standard bootstrap test: regression, single mediator, no covariates")


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
M <- a * X + rnorm(n)
Y <- b * M + c * X + rnorm(n)
test_data <- data.frame(X, Y, M)

## run bootstrap test and compute summary
foo <- test_mediation(test_data, x = "X", y = "Y", m = "M",
                      test = "boot", R = R, level = level, type = "bca",
                      method = "regression", robust = FALSE)
bar <- summary(foo)


## run tests

test_that("output has correct structure", {

  # bootstrap test
  expect_s3_class(foo, "boot_test_mediation")
  expect_s3_class(foo, "test_mediation")
  # regression fit
  expect_s3_class(foo$fit, "reg_fit_mediation")
  expect_s3_class(foo$fit, "fit_mediation")
  # individual robust regressions
  expect_s3_class(foo$fit$fit_mx, "lm")
  expect_s3_class(foo$fit$fit_ymx, "lm")
  expect_s3_class(foo$fit$fit_yx, "lm")
  # bootstrap replicates
  expect_s3_class(foo$reps, "boot")

})

test_that("arguments are correctly passed", {

  # number of bootstrap replicates
  expect_identical(foo$R, as.integer(R))
  # confidence level
  expect_identical(foo$level, level)
  # type of confidence intervals
  expect_identical(foo$type, "bca")
  # variable names
  expect_identical(foo$fit$x, "X")
  expect_identical(foo$fit$y, "Y")
  expect_identical(foo$fit$m, "M")
  expect_identical(foo$fit$covariates, character())
  # nonrobust fit and test
  expect_false(foo$fit$robust)

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
  expect_equivalent(foo$fit$a, coef(foo$fit$fit_mx)["X"])
  expect_equivalent(foo$fit$b, coef(foo$fit$fit_ymx)["M"])
  expect_equivalent(foo$fit$c, coef(foo$fit$fit_ymx)["X"])
  expect_equivalent(foo$fit$c_prime, coef(foo$fit$fit_yx)["X"])

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
  expect_named(bar$summary$s, c("value", "df"))
  # R-squared for model y ~ m + x
  expect_type(bar$summary$R2, "list")
  expect_named(bar$summary$R2, c("R2", "adj_R2"))
  # F-test for model y ~ m + x
  expect_type(bar$summary$F_test, "list")
  expect_named(bar$summary$F_test, c("statistic", "df", "p_value"))
  df_test <- bar$summary$F_test$df
  expect_identical(df_test[1], 2L)
  expect_identical(df_test[2], bar$summary$s$df)

})

test_that("attributes are correctly passed through summary", {
  # robustness
  expect_false(bar$summary$robust)
  # number of observations
  expect_identical(bar$summary$n, as.integer(n))
  # variable names
  expect_identical(bar$summary$x, "X")
  expect_identical(bar$summary$y, "Y")
  expect_identical(bar$summary$m, "M")
  expect_identical(bar$summary$covariates, character())
})

test_that("effect summaries have correct names", {

  # a path
  expect_identical(dim(bar$summary$a), c(1L, 5L))
  expect_identical(rownames(bar$summary$a), "X")
  expect_identical(colnames(bar$summary$a)[1:2], c("Data", "Boot"))
  # b path
  expect_identical(dim(bar$summary$b), c(1L, 5L))
  expect_identical(rownames(bar$summary$b), "M")
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
  expect_identical(bar$summary$b["M", "Data"], foo$fit$b)
  expect_identical(bar$summary$c["X", "Data"], foo$fit$c)
  expect_identical(bar$summary$c_prime["X", "Data"], foo$fit$c_prime)

  # bootstrapped effects
  expect_identical(bar$summary$a["X", "Boot"], mean(foo$reps$t[, 2]))
  expect_identical(bar$summary$b["M", "Boot"], mean(foo$reps$t[, 3]))
  expect_identical(bar$summary$c["X", "Boot"], mean(foo$reps$t[, 4]))
  expect_identical(bar$summary$c_prime["X", "Boot"], mean(foo$reps$t[, 5]))

})
