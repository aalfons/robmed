context("robust bootstrap test: regression, multiple mediators, no covariates")


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
test_data <- data.frame(X, Y, M1, M2)

## run bootstrap test
ctrl <- reg_control(efficiency = 0.95)
boot <- test_mediation(test_data, x = "X", y = "Y", m = c("M1", "M2"),
                       test = "boot", R = R, level = 0.9, type = "bca",
                       method = "regression", robust = TRUE, control = ctrl)

## compute summary
summary_boot <- summary(boot, other = "boot")
summary_theory <- summary(boot, other = "theory")

## create data for plotting
dot <- fortify(boot, method = "dot")
density <- fortify(boot, method = "density")

## stuff needed to check correctness
indirect_names <- c("Total", "M1", "M2")
ab_names <- paste("ab", indirect_names, sep = "_")
coef_names <- c("a_M1", "a_M2", "b_M1", "b_M2", "c", "c'", ab_names)


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
  expect_identical(boot$R, as.integer(R))  # doesn't hold for too many outliers
  # confidence level
  expect_identical(boot$level, 0.9)
  # type of confidence intervals
  expect_identical(boot$type, "bca")
  # variable names
  expect_identical(boot$fit$x, "X")
  expect_identical(boot$fit$y, "Y")
  expect_identical(boot$fit$m, c("M1", "M2"))
  expect_identical(boot$fit$covariates, character())
  # robust fit and test
  expect_true(boot$fit$robust)
  expect_equal(boot$fit$control, ctrl)

})

test_that("dimensions are correct", {

  # effects are scalars
  expect_length(boot$ab, 3L)
  expect_named(boot$ab, indirect_names)
  # multiple indirect effects and confidence intervals
  expect_identical(dim(boot$ci), c(3L, 2L))
  expect_identical(rownames(boot$ci), indirect_names)
  # dimensions of bootstrap replicates
  d_boot <- dim(boot$reps$t)
  expect_identical(d_boot, c(as.integer(R), 9L))

})

test_that("values of coefficients are correct", {

  expect_equivalent(boot$ab, colMeans(boot$reps$t[, 1:3]))

})

test_that("output of coef() method has correct attributes", {

  coef_boot <- coef(boot, type = "boot")
  coef_data <- coef(boot, type = "data")
  # bootstrapped effects
  expect_length(coef_boot, 9L)
  expect_named(coef_boot, coef_names)
  # effects computed on original sample
  expect_length(coef_data, 9L)
  expect_named(coef_data, coef_names)

})

test_that("coef() method returns correct values of coefficients", {

  # bootstrapped effects
  expect_equivalent(coef(boot, parm = c("a_M1", "a_M2"), type = "boot"),
                    colMeans(boot$reps$t[, 4:5]))
  expect_equivalent(coef(boot, parm = c("b_M1", "b_M2"), type = "boot"),
                    colMeans(boot$reps$t[, 6:7]))
  expect_equivalent(coef(boot, parm = "c", type = "boot"), mean(boot$reps$t[, 8]))
  expect_equivalent(coef(boot, parm = "c'", type = "boot"), mean(boot$reps$t[, 9]))
  expect_equivalent(coef(boot, parm = ab_names, type = "boot"), boot$ab)

  # effects computed on original sample
  expect_equivalent(coef(boot, parm = c("a_M1", "a_M2"), type = "data"), boot$fit$a)
  expect_equivalent(coef(boot, parm = c("b_M1", "b_M2"), type = "data"), boot$fit$b)
  expect_equivalent(coef(boot, parm = "c", type = "data"), boot$fit$c)
  expect_equivalent(coef(boot, parm = "c'", type = "data"), boot$fit$c_prime)
  ab_data <- boot$fit$a * boot$fit$b
  expect_equivalent(coef(boot, parm = ab_names, type = "data"),
                    c(sum(ab_data), ab_data))

})

test_that("output of confint() method has correct attributes", {

  ci_boot <- confint(boot, other = "boot")
  ci_theory <- confint(boot, other = "theory")
  # bootstrapped confidence intervals
  expect_equal(dim(ci_boot), c(9L, 2L))
  expect_equal(rownames(ci_boot), coef_names)
  expect_equal(colnames(ci_boot), c("5 %", "95 %"))
  # confidence intervals based on theory (except for indirect effect)
  expect_equal(dim(ci_theory), c(9L, 2L))
  expect_equal(rownames(ci_theory), coef_names)
  expect_equal(colnames(ci_theory), c("5 %", "95 %"))

})

test_that("confint() method returns correct values of confidence intervals", {

  # bootstrapped confidence intervals
  expect_equivalent(confint(boot, parm = ab_names, other = "boot"), boot$ci)
  # confidence intervals based on theory (except for indirect effect)
  expect_equivalent(confint(boot, parm = ab_names, other = "theory"), boot$ci)

})

test_that("summary has correct structure", {

  # summary
  expect_s3_class(summary_boot, "summary_test_mediation")
  expect_s3_class(summary_theory, "summary_test_mediation")
  # original output of test for indirect effect
  expect_identical(summary_boot$object, boot)
  expect_identical(summary_theory$object, boot)
  # summary of the model fit
  expect_s3_class(summary_boot$summary, "summary_fit_mediation")
  expect_s3_class(summary_theory$summary, "summary_fit_mediation")
  # regression standard error for model y ~ m + x
  expect_type(summary_boot$summary$s, "list")
  expect_named(summary_boot$summary$s, c("value", "df"))
  expect_type(summary_theory$summary$s, "list")
  expect_named(summary_theory$summary$s, c("value", "df"))
  # R-squared for model y ~ m + x
  expect_type(summary_boot$summary$R2, "list")
  expect_named(summary_boot$summary$R2, c("R2", "adj_R2"))
  expect_type(summary_theory$summary$R2, "list")
  expect_named(summary_theory$summary$R2, c("R2", "adj_R2"))
  # F-test for model y ~ m + x
  expect_type(summary_boot$summary$F_test, "list")
  expect_named(summary_boot$summary$F_test, c("statistic", "df", "p_value"))
  df_test_boot <- summary_boot$summary$F_test$df
  expect_identical(df_test_boot[1], 3)
  expect_identical(df_test_boot[2], Inf)
  expect_type(summary_theory$summary$F_test, "list")
  expect_named(summary_theory$summary$F_test, c("statistic", "df", "p_value"))
  df_test_theory <- summary_theory$summary$F_test$df
  expect_identical(df_test_theory[1], 3)
  expect_identical(df_test_theory[2], Inf)

})

test_that("attributes are correctly passed through summary", {

  # robustness
  expect_true(summary_boot$summary$robust)
  expect_true(summary_theory$summary$robust)
  # number of observations
  expect_identical(summary_boot$summary$n, as.integer(n))
  expect_identical(summary_theory$summary$n, as.integer(n))
  # variable names
  expect_identical(summary_boot$summary$x, "X")
  expect_identical(summary_boot$summary$y, "Y")
  expect_identical(summary_boot$summary$m, c("M1", "M2"))
  expect_identical(summary_boot$summary$covariates, character())
  expect_identical(summary_theory$summary$x, "X")
  expect_identical(summary_theory$summary$y, "Y")
  expect_identical(summary_theory$summary$m, c("M1", "M2"))
  expect_identical(summary_theory$summary$covariates, character())

})

test_that("effect summaries have correct names", {

  # a path
  expect_identical(dim(summary_boot$summary$a), c(2L, 5L))
  expect_identical(rownames(summary_boot$summary$a), c("M1~X", "M2~X"))
  expect_identical(colnames(summary_boot$summary$a)[1:2], c("Data", "Boot"))
  expect_identical(dim(summary_theory$summary$a), c(2L, 4L))
  expect_identical(rownames(summary_theory$summary$a), c("M1~X", "M2~X"))
  expect_identical(colnames(summary_theory$summary$a)[1], "Estimate")
  # b path
  expect_identical(dim(summary_boot$summary$b), c(2L, 5L))
  expect_identical(rownames(summary_boot$summary$b), c("M1", "M2"))
  expect_identical(colnames(summary_boot$summary$b)[1:2], c("Data", "Boot"))
  expect_identical(dim(summary_theory$summary$b), c(2L, 4L))
  expect_identical(rownames(summary_theory$summary$b), c("M1", "M2"))
  expect_identical(colnames(summary_theory$summary$b)[1], "Estimate")
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
  # covariates
  expect_null(summary_boot$summary$covariate_effects)
  expect_null(summary_theory$summary$covariate_effects)

})

test_that("effect summaries contain correct coefficient values", {

  # effects computed on original sample
  expect_equivalent(summary_boot$summary$a[, "Data"], boot$fit$a)
  expect_identical(summary_boot$summary$b[, "Data"], boot$fit$b)
  expect_identical(summary_boot$summary$c["X", "Data"], boot$fit$c)
  expect_identical(summary_boot$summary$c_prime["X", "Data"], boot$fit$c_prime)
  expect_equivalent(summary_theory$summary$a[, "Estimate"], boot$fit$a)
  expect_identical(summary_theory$summary$b[, "Estimate"], boot$fit$b)
  expect_identical(summary_theory$summary$c["X", "Estimate"], boot$fit$c)
  expect_identical(summary_theory$summary$c_prime["X", "Estimate"], boot$fit$c_prime)

  # bootstrapped effects
  expect_equivalent(summary_boot$summary$a[, "Boot"], colMeans(boot$reps$t[, 4:5]))
  expect_equivalent(summary_boot$summary$b[, "Boot"], colMeans(boot$reps$t[, 6:7]))
  expect_equal(summary_boot$summary$c["X", "Boot"], mean(boot$reps$t[, 8]))
  expect_equal(summary_boot$summary$c_prime["X", "Boot"], mean(boot$reps$t[, 9]))

})

test_that("data returned by fortify() has correct structure", {

  ## dot plot
  # check dimensions
  expect_s3_class(dot, "data.frame")
  expect_identical(dim(dot), c(4L, 4L))
  # check column names
  column_names <- c("Effect", "Point", "Lower", "Upper")
  expect_named(dot, column_names)
  # check that direct effect and indirect effect are plotted by default
  effect_names <- c("c", ab_names)
  expect_identical(dot$Effect, factor(effect_names, levels = effect_names))

  ## density plot
  # check dimensions
  expect_s3_class(density, "data.frame")
  expect_identical(ncol(density), 3L)
  expect_gt(nrow(density), 0L)
  # check column names
  column_names <- c("ab", "Density", "Effect")
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
  expect_equal(attr(density, "facets"), ~Effect)
  # check title
  expect_identical(attr(density, "main"), "Bootstrap distribution")
  # check confidence interval
  ci <- attr(density, "ci")
  expect_s3_class(ci, "data.frame")
  expect_identical(dim(ci), c(3L, 5L))
  expect_named(ci, c("ab", "Density", "Lower", "Upper", "Effect"))
  # check that method is stored correctly
  expect_identical(attr(density, "method"), "density")

})
