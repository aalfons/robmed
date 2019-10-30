context("standard bootstrap test: covariance")


## load package
library("robmed", quietly = TRUE)

## control parameters
n <- 250          # number of observations
a <- c <- 0.2     # true effects
b <- 0            # true effect
R <- 1000         # number of bootstrap samples
seed <- 20150601  # seed for the random number generator

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

## run bootstrap test and compute summary
boot <- test_mediation(test_data, x = "X", y = "Y", m = "M1", test = "boot",
                       R = R, level = 0.9, type = "bca", method = "covariance",
                       robust = FALSE)

## compute summary
summary_boot <- summary(boot, other = "boot")
summary_theory <- summary(boot, other = "theory")

## create data for plotting
dot <- fortify(boot, method = "dot")
density <- fortify(boot, method = "density")

## stuff needed to check correctness
coef_names <- c("a", "b", "Direct", "Total", "ab")


## run tests

test_that("output has correct structure", {

  # bootstrap test
  expect_s3_class(boot, "boot_test_mediation")
  expect_s3_class(boot, "test_mediation")
  # regression fit
  expect_s3_class(boot$fit, "cov_fit_mediation")
  expect_s3_class(boot$fit, "fit_mediation")
  # bootstrap replicates
  expect_s3_class(boot$reps, "boot")

})

test_that("arguments are correctly passed", {

  # alternative hypothesis
  expect_identical(boot$alternative, "twosided")
  # number of bootstrap replicates
  expect_identical(boot$R, as.integer(R))
  # confidence level
  expect_identical(boot$level, 0.9)
  # type of confidence intervals
  expect_identical(boot$type, "bca")
  # variable names
  expect_identical(boot$fit$x, "X")
  expect_identical(boot$fit$y, "Y")
  expect_identical(boot$fit$m, "M1")
  expect_null(boot$fit$covariates)
  # nonrobust fit and test
  expect_false(boot$fit$robust)
  expect_null(boot$fit$control)

})

test_that("dimensions are correct", {

  # effects are scalars
  expect_length(boot$ab, 1L)
  # only one indirect effect, so only one confidence interval
  expect_length(boot$ci, 2L)
  # dimensions of bootstrap replicates
  d_boot <- dim(boot$reps$t)
  expect_identical(d_boot, c(as.integer(R), 7L))

})

test_that("values of coefficients are correct", {

  expect_equivalent(boot$ab, mean(boot$reps$t[, 1]))

})

test_that("output of coef() method has correct attributes", {

  coef_boot <- coef(boot, type = "boot")
  coef_data <- coef(boot, type = "data")
  # bootstrapped effects
  expect_length(coef_boot, 5L)
  expect_named(coef_boot, coef_names)
  # effects computed on original sample
  expect_length(coef_data, 5L)
  expect_named(coef_data, coef_names)

})

test_that("coef() method returns correct values of coefficients", {

  # bootstrapped effects
  expect_equivalent(coef(boot, parm = "a", type = "boot"), mean(boot$reps$t[, 3]))
  expect_equivalent(coef(boot, parm = "b", type = "boot"), mean(boot$reps$t[, 5]))
  expect_equivalent(coef(boot, parm = "Direct", type = "boot"), mean(boot$reps$t[, 6]))
  expect_equivalent(coef(boot, parm = "Total", type = "boot"), mean(boot$reps$t[, 7]))
  expect_equivalent(coef(boot, parm = "ab", type = "boot"), boot$ab)

  # effects computed on original sample
  expect_equivalent(coef(boot, parm = "a", type = "data"), boot$fit$a)
  expect_equivalent(coef(boot, parm = "b", type = "data"), boot$fit$b)
  expect_equivalent(coef(boot, parm = "Direct", type = "data"), boot$fit$direct)
  expect_equivalent(coef(boot, parm = "Total", type = "data"), boot$fit$total)
  expect_equivalent(coef(boot, parm = "ab", type = "data"), boot$fit$a * boot$fit$b)

})

test_that("output of confint() method has correct attributes", {

  ci_boot <- confint(boot, other = "boot")
  ci_theory <- confint(boot, other = "theory")
  # bootstrapped confidence intervals
  expect_equal(dim(ci_boot), c(5L, 2L))
  expect_equal(rownames(ci_boot), coef_names)
  expect_equal(colnames(ci_boot), c("5 %", "95 %"))
  # confidence intervals based on theory (except for indirect effect)
  expect_equal(dim(ci_theory), c(5L, 2L))
  expect_equal(rownames(ci_theory), coef_names)
  expect_equal(colnames(ci_theory), c("5 %", "95 %"))

})

test_that("confint() method returns correct values of confidence intervals", {

  # bootstrapped confidence intervals
  expect_equivalent(confint(boot, parm = "ab", other = "boot"), boot$ci)
  # confidence intervals based on theory (except for indirect effect)
  expect_equivalent(confint(boot, parm = "ab", other = "theory"), boot$ci)

})

test_that("summary has correct structure", {

  # summary
  expect_s3_class(summary_boot, "summary_test_mediation")
  expect_s3_class(summary_theory, "summary_test_mediation")
  # original output of test for indirect effect
  expect_identical(summary_boot$object, boot)
  expect_identical(summary_theory$object, boot)
  # summary of the model fit
  expect_s3_class(summary_boot$summary, "summary_cov_fit_mediation")
  expect_s3_class(summary_boot$summary, "summary_fit_mediation")
  expect_s3_class(summary_theory$summary, "summary_cov_fit_mediation")
  expect_s3_class(summary_theory$summary, "summary_fit_mediation")
  # regression standard error for model y ~ m + x
  expect_null(summary_boot$summary$s)
  expect_null(summary_theory$summary$s)
  # R-squared for model y ~ m + x
  expect_null(summary_boot$summary$R2)
  expect_null(summary_theory$summary$R2)
  # F-test for model y ~ m + x
  expect_null(summary_boot$summary$F_test)
  expect_null(summary_theory$summary$F_test)

})

test_that("attributes are correctly passed through summary", {

  # robustness
  expect_false(summary_boot$summary$robust)
  expect_false(summary_theory$summary$robust)
  # number of observations
  expect_identical(summary_boot$summary$n, as.integer(n))
  expect_identical(summary_theory$summary$n, as.integer(n))
  # variable names
  expect_identical(summary_boot$summary$x, "X")
  expect_identical(summary_boot$summary$y, "Y")
  expect_identical(summary_boot$summary$m, "M1")
  expect_null(summary_boot$summary$covariates)
  expect_identical(summary_theory$summary$x, "X")
  expect_identical(summary_theory$summary$y, "Y")
  expect_identical(summary_theory$summary$m, "M1")
  expect_null(summary_theory$summary$covariates)

})

test_that("effect summaries have correct names", {

  # a path
  expect_identical(dim(summary_boot$summary$a), c(1L, 5L))
  expect_identical(rownames(summary_boot$summary$a), "X")
  expect_identical(colnames(summary_boot$summary$a)[1:2], c("Data", "Boot"))
  expect_identical(dim(summary_theory$summary$a), c(1L, 4L))
  expect_identical(rownames(summary_theory$summary$a), "X")
  expect_identical(colnames(summary_theory$summary$a)[1], "Estimate")
  # b path
  expect_identical(dim(summary_boot$summary$b), c(1L, 5L))
  expect_identical(rownames(summary_boot$summary$b), "M1")
  expect_identical(colnames(summary_boot$summary$b)[1:2], c("Data", "Boot"))
  expect_identical(dim(summary_theory$summary$b), c(1L, 4L))
  expect_identical(rownames(summary_theory$summary$b), "M1")
  expect_identical(colnames(summary_theory$summary$b)[1], "Estimate")
  # direct effect
  expect_identical(dim(summary_boot$summary$direct), c(1L, 5L))
  expect_identical(rownames(summary_boot$summary$direct), "X")
  expect_identical(colnames(summary_boot$summary$direct)[1:2], c("Data", "Boot"))
  expect_identical(dim(summary_theory$summary$direct), c(1L, 4L))
  expect_identical(rownames(summary_theory$summary$direct), "X")
  expect_identical(colnames(summary_theory$summary$direct)[1], "Estimate")
  # total effect
  expect_identical(dim(summary_boot$summary$total), c(1L, 5L))
  expect_identical(rownames(summary_boot$summary$total), "X")
  expect_identical(colnames(summary_boot$summary$total)[1:2], c("Data", "Boot"))
  expect_identical(dim(summary_theory$summary$total), c(1L, 4L))
  expect_identical(rownames(summary_theory$summary$total), "X")
  expect_identical(colnames(summary_theory$summary$total)[1], "Estimate")
  # no model fits
  expect_null(summary_boot$summary$fit_mx)
  expect_null(summary_boot$summary$fit_ymx)
  expect_null(summary_theory$summary$fit_mx)
  expect_null(summary_theory$summary$fit_ymx)

})

test_that("effect summaries contain correct coefficient values", {

  # effects computed on original sample
  expect_identical(summary_boot$summary$a["X", "Data"], boot$fit$a)
  expect_identical(summary_boot$summary$b["M1", "Data"], boot$fit$b)
  expect_identical(summary_boot$summary$direct["X", "Data"], boot$fit$direct)
  expect_identical(summary_boot$summary$total["X", "Data"], boot$fit$total)
  expect_identical(summary_theory$summary$a["X", "Estimate"], boot$fit$a)
  expect_identical(summary_theory$summary$b["M1", "Estimate"], boot$fit$b)
  expect_identical(summary_theory$summary$direct["X", "Estimate"], boot$fit$direct)
  expect_identical(summary_theory$summary$total["X", "Estimate"], boot$fit$total)

  # bootstrapped effects
  expect_equal(summary_boot$summary$a["X", "Boot"], mean(boot$reps$t[, 3]))
  expect_equal(summary_boot$summary$b["M1", "Boot"], mean(boot$reps$t[, 5]))
  expect_equal(summary_boot$summary$direct["X", "Boot"], mean(boot$reps$t[, 6]))
  expect_equal(summary_boot$summary$total["X", "Boot"], mean(boot$reps$t[, 7]))

})

test_that("data returned by fortify() has correct structure", {

  ## dot plot
  # check dimensions
  expect_s3_class(dot, "data.frame")
  expect_identical(dim(dot), c(2L, 4L))
  # check column names
  column_names <- c("Effect", "Point", "Lower", "Upper")
  expect_named(dot, column_names)
  # check that direct effect and indirect effect are plotted by default
  effect_names <- c("Direct", "ab")
  expect_identical(dot$Effect, factor(effect_names, levels = effect_names))

  ## density plot
  # check dimensions
  expect_s3_class(density, "data.frame")
  expect_identical(ncol(density), 2L)
  expect_gt(nrow(density), 0L)
  # check column names
  column_names <- c("ab", "Density")
  expect_named(density, column_names)

})

test_that("data returned by fortify() has correct attributes", {

  ## dot plot
  # check aesthetic mapping
  mapping <- aes_string(x = "Effect", y = "Point",
                        ymin = "Lower", ymax = "Upper")
  expect_equal(attr(dot, "mapping"), mapping)
  # check default geom()
  expect_identical(attr(dot, "geom"), geom_pointrange)
  # check facets
  expect_null(attr(dot, "facets"))
  # check that method is stored correctly
  expect_identical(attr(dot, "method"), "dot")

  ## density plot
  # check aesthetic mapping
  mapping <- aes_string(x = "ab", y = "Density")
  expect_equal(attr(density, "mapping"), mapping)
  # check default geom()
  expect_equal(attr(density, "geom"), function(..., stat) {
    geom_density(..., stat = "identity")
  })
  # check facets
  expect_null(attr(density, "facets"))
  # check title
  expect_identical(attr(density, "main"), "Bootstrap distribution")
  # check confidence interval
  ci <- attr(density, "ci")
  expect_s3_class(ci, "data.frame")
  expect_identical(dim(ci), c(1L, 4L))
  expect_named(ci, c("ab", "Density", "Lower", "Upper"))
  # check that method is stored correctly
  expect_identical(attr(density, "method"), "density")

})


# ## only implemented for simple mediation without covariates
#
# test_that("covariates not implemented", {
#
#   # run test with regression method
#   set.seed(seed)
#   suppressWarnings(
#     test_reg <- test_mediation(test_data, x = "X", y = "Y", m = "M1",
#                                covariates = c("C1", "C2"), test = "boot",
#                                R = 500, level = 0.9, type = "perc",
#                                method = "regression", robust = FALSE)
#   )
#
#   # try to run test with covariates (should give warning)
#   set.seed(seed)
#   expect_warning(
#     test_cov <- test_mediation(test_data, x = "X", y = "Y", m = "M1",
#                                covariates = c("C1", "C2"), test = "boot",
#                                R = 500, level = 0.9, type = "perc",
#                                method = "covariance", robust = FALSE)
#   )
#
#   # these should be the same
#   expect_equal(test_cov, test_reg)
#
# })
#
#
# test_that("multiple mediators not implemented", {
#
#   # run test with regression method
#   set.seed(seed)
#   suppressWarnings(
#     test_reg <- test_mediation(test_data, x = "X", y = "Y", m = c("M1", "M2"),
#                                test = "boot", R = 500, level = 0.9,
#                                type = "perc", method = "regression",
#                                robust = FALSE)
#   )
#
#   # try to run test with multiple mediators (should give warning)
#   set.seed(seed)
#   expect_warning(
#     test_cov <- test_mediation(test_data, x = "X", y = "Y", m = c("M1", "M2"),
#                                test = "boot", R = 500, level = 0.9,
#                                type = "perc", method = "covariance",
#                                robust = FALSE)
#   )
#
#   # these should be the same
#   expect_equal(test_cov, test_reg)
#
# })
