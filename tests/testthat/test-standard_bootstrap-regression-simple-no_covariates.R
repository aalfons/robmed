context("standard bootstrap test: regression, single mediator, no covariates")


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
M <- a * X + rnorm(n)
Y <- b * M + c * X + rnorm(n)
test_data <- data.frame(X, Y, M)

## run bootstrap test
boot <- test_mediation(test_data, x = "X", y = "Y", m = "M", test = "boot",
                       R = R, level = 0.9, type = "bca", method = "regression",
                       robust = FALSE)

## compute summary
summary_boot <- summary(boot, other = "boot")
summary_theory <- summary(boot, other = "theory")

## create data for plotting
dot <- fortify(boot, method = "dot")
density <- fortify(boot, method = "density")

## stuff needed to check correctness
coef_names <- c("a", "b", "c", "c'", "ab")
mx_names <- c("(Intercept)", "X")
ymx_names <- c("(Intercept)", "M", "X")


## run tests

test_that("output has correct structure", {

  # bootstrap test
  expect_s3_class(boot, "boot_test_mediation")
  expect_s3_class(boot, "test_mediation")
  # regression fit
  expect_s3_class(boot$fit, "reg_fit_mediation")
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
  expect_identical(boot$fit$m, "M")
  expect_identical(boot$fit$covariates, character())
  # nonrobust fit and test
  expect_false(boot$fit$robust)
  expect_false(boot$fit$median)
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
  expect_equivalent(coef(boot, parm = "c", type = "boot"), mean(boot$reps$t[, 6]))
  expect_equivalent(coef(boot, parm = "c'", type = "boot"), mean(boot$reps$t[, 7]))
  expect_equivalent(coef(boot, parm = "ab", type = "boot"), boot$ab)

  # effects computed on original sample
  expect_equivalent(coef(boot, parm = "a", type = "data"), boot$fit$a)
  expect_equivalent(coef(boot, parm = "b", type = "data"), boot$fit$b)
  expect_equivalent(coef(boot, parm = "c", type = "data"), boot$fit$c)
  expect_equivalent(coef(boot, parm = "c'", type = "data"), boot$fit$c_prime)
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
  expect_s3_class(summary_boot$summary, "summary_reg_fit_mediation")
  expect_s3_class(summary_boot$summary, "summary_fit_mediation")
  expect_s3_class(summary_theory$summary, "summary_reg_fit_mediation")
  expect_s3_class(summary_theory$summary, "summary_fit_mediation")
  # summary for model m ~ x
  expect_s3_class(summary_boot$summary$fit_mx, "summary_lm")
  # summary for model y ~ m + x
  expect_s3_class(summary_boot$summary$fit_ymx, "summary_lm")
  # regression standard error for model y ~ m + x
  expect_type(summary_boot$summary$fit_ymx$s, "list")
  expect_named(summary_boot$summary$fit_ymx$s, c("value", "df"))
  expect_type(summary_theory$summary$fit_ymx$s, "list")
  expect_named(summary_theory$summary$fit_ymx$s, c("value", "df"))
  # R-squared for model y ~ m + x
  expect_type(summary_boot$summary$fit_ymx$R2, "list")
  expect_named(summary_boot$summary$fit_ymx$R2, c("R2", "adj_R2"))
  expect_type(summary_theory$summary$fit_ymx$R2, "list")
  expect_named(summary_theory$summary$fit_ymx$R2, c("R2", "adj_R2"))
  # F-test for model y ~ m + x
  expect_type(summary_boot$summary$fit_ymx$F_test, "list")
  expect_named(summary_boot$summary$fit_ymx$F_test, c("statistic", "df", "p_value"))
  df_test_boot <- summary_boot$summary$fit_ymx$F_test$df
  expect_identical(df_test_boot[1], 2L)
  expect_identical(df_test_boot[2], summary_boot$summary$fit_ymx$s$df)
  expect_type(summary_theory$summary$fit_ymx$F_test, "list")
  expect_named(summary_theory$summary$fit_ymx$F_test, c("statistic", "df", "p_value"))
  df_test_theory <- summary_theory$summary$fit_ymx$F_test$df
  expect_identical(df_test_theory[1], 2L)
  expect_identical(df_test_theory[2], summary_theory$summary$fit_ymx$s$df)

})

test_that("attributes are correctly passed through summary", {

  # robustness
  expect_false(summary_boot$summary$robust)
  expect_false(summary_boot$summary$median)
  expect_false(summary_theory$summary$robust)
  expect_false(summary_theory$summary$median)
  # number of observations
  expect_identical(summary_boot$summary$n, as.integer(n))
  expect_identical(summary_theory$summary$n, as.integer(n))
  # variable names
  expect_identical(summary_boot$summary$x, "X")
  expect_identical(summary_boot$summary$y, "Y")
  expect_identical(summary_boot$summary$m, "M")
  expect_identical(summary_boot$summary$covariates, character())
  expect_identical(summary_theory$summary$x, "X")
  expect_identical(summary_theory$summary$y, "Y")
  expect_identical(summary_theory$summary$m, "M")
  expect_identical(summary_theory$summary$covariates, character())

})

test_that("effect summaries have correct names", {

  # a path
  expect_identical(dim(summary_boot$summary$fit_mx$coefficients), c(2L, 5L))
  expect_identical(rownames(summary_boot$summary$fit_mx$coefficients), mx_names)
  expect_identical(colnames(summary_boot$summary$fit_mx$coefficients)[1:2], c("Data", "Boot"))
  expect_identical(dim(summary_theory$summary$fit_mx$coefficients), c(2L, 4L))
  expect_identical(rownames(summary_theory$summary$fit_mx$coefficients), mx_names)
  expect_identical(colnames(summary_theory$summary$fit_mx$coefficients)[1], "Estimate")
  # b path
  expect_identical(dim(summary_boot$summary$fit_ymx$coefficients), c(3L, 5L))
  expect_identical(rownames(summary_boot$summary$fit_ymx$coefficient), ymx_names)
  expect_identical(colnames(summary_boot$summary$fit_ymx$coefficient)[1:2], c("Data", "Boot"))
  expect_identical(dim(summary_theory$summary$fit_ymx$coefficient), c(3L, 4L))
  expect_identical(rownames(summary_theory$summary$fit_ymx$coefficient), ymx_names)
  expect_identical(colnames(summary_theory$summary$fit_ymx$coefficient)[1], "Estimate")
  # c path
  expect_identical(dim(summary_boot$summary$c), c(1L, 5L))
  expect_identical(rownames(summary_boot$summary$c), "X")
  expect_identical(colnames(summary_boot$summary$c)[1:2], c("Data", "Boot"))
  expect_identical(dim(summary_theory$summary$c), c(1L, 4L))
  expect_identical(rownames(summary_theory$summary$c), "X")
  expect_identical(colnames(summary_theory$summary$c)[1], "Estimate")
  # c' path
  expect_identical(dim(summary_boot$summary$c_prime), c(1L, 5L))
  expect_identical(rownames(summary_boot$summary$c_prime), "X")
  expect_identical(colnames(summary_boot$summary$c_prime)[1:2], c("Data", "Boot"))
  expect_identical(dim(summary_theory$summary$c_prime), c(1L, 4L))
  expect_identical(rownames(summary_theory$summary$c_prime), "X")
  expect_identical(colnames(summary_theory$summary$c_prime)[1], "Estimate")

})

test_that("effect summaries contain correct coefficient values", {

  # effects computed on original sample
  expect_equivalent(summary_boot$summary$fit_mx$coefficients[2, "Data"], boot$fit$a)
  expect_identical(summary_boot$summary$fit_ymx$coefficients[2, "Data"], boot$fit$b)
  expect_identical(summary_boot$summary$c["X", "Data"], boot$fit$c)
  expect_identical(summary_boot$summary$c_prime["X", "Data"], boot$fit$c_prime)
  expect_equivalent(summary_theory$summary$fit_mx$coefficients[2, "Estimate"], boot$fit$a)
  expect_identical(summary_theory$summary$fit_ymx$coefficients[2, "Estimate"], boot$fit$b)
  expect_identical(summary_theory$summary$c["X", "Estimate"], boot$fit$c)
  expect_identical(summary_theory$summary$c_prime["X", "Estimate"], boot$fit$c_prime)

  # bootstrapped effects
  expect_equivalent(summary_boot$summary$fit_mx$coefficients[2, "Boot"], mean(boot$reps$t[, 3]))
  expect_equivalent(summary_boot$summary$fit_ymx$coefficients[2, "Boot"], mean(boot$reps$t[, 5]))
  expect_equal(summary_boot$summary$c["X", "Boot"], mean(boot$reps$t[, 6]))
  expect_equal(summary_boot$summary$c_prime["X", "Boot"], mean(boot$reps$t[, 7]))

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
  effect_names <- c("c", "ab")
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
