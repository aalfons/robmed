context("standard bootstrap test: regression, multiple mediators, covariates")


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

## run bootstrap test
set.seed(seed)
level <- 0.9
boot <- test_mediation(test_data, x = "X", y = "Y", m = c("M1", "M2"),
                       covariates = c("C1", "C2"), test = "boot", R = R,
                       level = level, type = "bca", method = "regression",
                       robust = FALSE)

## compute summary
summary_boot <- summary(boot, type = "boot")
summary_data <- summary(boot, type = "data")

## create data for plotting
ci <- setup_ci_plot(boot)
density <- setup_density_plot(boot)
ellipse <- setup_ellipse_plot(boot)
# deprecated:
dot_deprecated <- suppressWarnings(fortify(boot, method = "dot"))
density_deprecated <- suppressWarnings(fortify(boot, method = "density"))

## stuff needed to check correctness
indirect_names <- c("Total", "M1", "M2")
ab_names <- paste("ab", indirect_names, sep = "_")
coef_names <- c("a_M1", "a_M2", "b_M1", "b_M2", "Direct", "Total", ab_names)
mx_names <- c("(Intercept)", "X", "C1", "C2")
ymx_names <- c("(Intercept)", "M1", "M2", "X", "C1", "C2")


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
  expect_identical(boot$level, level)
  # type of confidence intervals
  expect_identical(boot$type, "bca")
  # variable names
  expect_identical(boot$fit$x, "X")
  expect_identical(boot$fit$y, "Y")
  expect_identical(boot$fit$m, c("M1", "M2"))
  expect_identical(boot$fit$covariates, c("C1", "C2"))
  # nonrobust fit and test
  expect_false(boot$fit$robust)
  expect_null(boot$fit$control)

})

test_that("dimensions are correct", {

  # multiple indirect effects and confidence intervals
  expect_length(boot$ab, 3L)
  expect_named(boot$ab, indirect_names)
  expect_identical(dim(boot$ci), c(3L, 2L))
  expect_identical(rownames(boot$ci), indirect_names)
  # dimensions of bootstrap replicates
  d_boot <- dim(boot$reps$t)
  expect_identical(d_boot, c(as.integer(R), 18L))

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
                    colMeans(boot$reps$t[, c(5, 9)]))
  expect_equivalent(coef(boot, parm = c("b_M1", "b_M2"), type = "boot"),
                    colMeans(boot$reps$t[, 13:14]))
  expect_equivalent(coef(boot, parm = "Direct", type = "boot"), mean(boot$reps$t[, 15]))
  expect_equivalent(coef(boot, parm = "Total", type = "boot"), mean(boot$reps$t[, 18]))
  expect_equivalent(coef(boot, parm = ab_names, type = "boot"), boot$ab)

  # effects computed on original sample
  expect_equivalent(coef(boot, parm = c("a_M1", "a_M2"), type = "data"), boot$fit$a)
  expect_equivalent(coef(boot, parm = c("b_M1", "b_M2"), type = "data"), boot$fit$b)
  expect_equivalent(coef(boot, parm = "Direct", type = "data"), boot$fit$direct)
  expect_equivalent(coef(boot, parm = "Total", type = "data"), boot$fit$total)
  ab_data <- boot$fit$a * boot$fit$b
  expect_equivalent(coef(boot, parm = ab_names, type = "data"),
                    c(sum(ab_data), ab_data))

})

test_that("output of confint() method has correct attributes", {

  ci_boot <- confint(boot, type = "boot")
  ci_data <- confint(boot, type = "data")
  # bootstrapped confidence intervals
  expect_equal(dim(ci_boot), c(9L, 2L))
  expect_equal(rownames(ci_boot), coef_names)
  expect_equal(colnames(ci_boot), c("5 %", "95 %"))
  # confidence intervals based on theory (except for indirect effect)
  expect_equal(dim(ci_data), c(9L, 2L))
  expect_equal(rownames(ci_data), coef_names)
  expect_equal(colnames(ci_data), c("5 %", "95 %"))

})

test_that("confint() method returns correct values of confidence intervals", {

  # bootstrapped confidence intervals
  expect_equivalent(confint(boot, parm = ab_names, type = "boot"), boot$ci)
  # confidence intervals based on theory (except for indirect effect)
  expect_equivalent(confint(boot, parm = ab_names, type = "data"), boot$ci)

})

test_that("summary has correct structure", {

  # summary
  expect_s3_class(summary_boot, "summary_test_mediation")
  expect_s3_class(summary_data, "summary_test_mediation")
  # original output of test for indirect effect
  expect_identical(summary_boot$object, boot)
  expect_identical(summary_data$object, boot)
  # summary of the model fit
  expect_s3_class(summary_boot$summary, "summary_reg_fit_mediation")
  expect_s3_class(summary_boot$summary, "summary_fit_mediation")
  expect_s3_class(summary_data$summary, "summary_reg_fit_mediation")
  expect_s3_class(summary_data$summary, "summary_fit_mediation")
  # number of models m ~ x
  expect_type(summary_boot$summary$fit_mx, "list")
  expect_length(summary_boot$summary$fit_mx, 2L)
  expect_s3_class(summary_boot$summary$fit_mx[[1]], "summary_lm")
  expect_s3_class(summary_boot$summary$fit_mx[[2]], "summary_lm")
  # summary for model y ~ m + x
  expect_s3_class(summary_boot$summary$fit_ymx, "summary_lm")
  # regression standard error for model y ~ m + x
  expect_type(summary_boot$summary$fit_ymx$s, "list")
  expect_named(summary_boot$summary$fit_ymx$s, c("value", "df"))
  expect_type(summary_data$summary$fit_ymx$s, "list")
  expect_named(summary_data$summary$fit_ymx$s, c("value", "df"))
  # R-squared for model y ~ m + x
  expect_type(summary_boot$summary$fit_ymx$R2, "list")
  expect_named(summary_boot$summary$fit_ymx$R2, c("R2", "adj_R2"))
  expect_type(summary_data$summary$fit_ymx$R2, "list")
  expect_named(summary_data$summary$fit_ymx$R2, c("R2", "adj_R2"))
  # F-test for model y ~ m + x
  expect_type(summary_boot$summary$fit_ymx$F_test, "list")
  expect_named(summary_boot$summary$fit_ymx$F_test, c("statistic", "df", "p_value"))
  df_test_boot <- summary_boot$summary$fit_ymx$F_test$df
  expect_identical(df_test_boot[1], 5L)
  expect_identical(df_test_boot[2], summary_boot$summary$fit_ymx$s$df)
  expect_type(summary_data$summary$fit_ymx$F_test, "list")
  expect_named(summary_data$summary$fit_ymx$F_test, c("statistic", "df", "p_value"))
  df_test_data <- summary_data$summary$fit_ymx$F_test$df
  expect_identical(df_test_data[1], 5L)
  expect_identical(df_test_data[2], summary_data$summary$fit_ymx$s$df)

})

test_that("attributes are correctly passed through summary", {

  # robustness
  expect_false(summary_boot$summary$robust)
  expect_false(summary_data$summary$robust)
  # number of observations
  expect_identical(summary_boot$summary$n, as.integer(n))
  expect_identical(summary_data$summary$n, as.integer(n))
  # variable names
  expect_identical(summary_boot$summary$x, "X")
  expect_identical(summary_boot$summary$y, "Y")
  expect_identical(summary_boot$summary$m, c("M1", "M2"))
  expect_identical(summary_boot$summary$covariates, c("C1", "C2"))
  expect_identical(summary_data$summary$x, "X")
  expect_identical(summary_data$summary$y, "Y")
  expect_identical(summary_data$summary$m, c("M1", "M2"))
  expect_identical(summary_data$summary$covariates, c("C1", "C2"))

})

test_that("effect summaries have correct names", {

  # a path
  expect_identical(dim(summary_boot$summary$fit_mx[[1]]$coefficients), c(4L, 5L))
  expect_identical(rownames(summary_boot$summary$fit_mx[[1]]$coefficients), mx_names)
  expect_identical(colnames(summary_boot$summary$fit_mx[[1]]$coefficients)[1:2], c("Data", "Boot"))
  expect_identical(dim(summary_boot$summary$fit_mx[[2]]$coefficients), c(4L, 5L))
  expect_identical(rownames(summary_boot$summary$fit_mx[[2]]$coefficients), mx_names)
  expect_identical(colnames(summary_boot$summary$fit_mx[[2]]$coefficients)[1:2], c("Data", "Boot"))
  expect_identical(dim(summary_data$summary$fit_mx[[1]]$coefficients), c(4L, 4L))
  expect_identical(rownames(summary_data$summary$fit_mx[[1]]$coefficients), mx_names)
  expect_identical(colnames(summary_data$summary$fit_mx[[1]]$coefficients)[1], "Estimate")
  expect_identical(dim(summary_data$summary$fit_mx[[2]]$coefficients), c(4L, 4L))
  expect_identical(rownames(summary_data$summary$fit_mx[[2]]$coefficients), mx_names)
  expect_identical(colnames(summary_data$summary$fit_mx[[2]]$coefficients)[1], "Estimate")
  # b path
  expect_identical(dim(summary_boot$summary$fit_ymx$coefficients), c(6L, 5L))
  expect_identical(rownames(summary_boot$summary$fit_ymx$coefficient), ymx_names)
  expect_identical(colnames(summary_boot$summary$fit_ymx$coefficient)[1:2], c("Data", "Boot"))
  expect_identical(dim(summary_data$summary$fit_ymx$coefficient), c(6L, 4L))
  expect_identical(rownames(summary_data$summary$fit_ymx$coefficient), ymx_names)
  expect_identical(colnames(summary_data$summary$fit_ymx$coefficient)[1], "Estimate")
  # direct effect
  expect_identical(dim(summary_boot$summary$direct), c(1L, 5L))
  expect_identical(rownames(summary_boot$summary$direct), "X")
  expect_identical(colnames(summary_boot$summary$direct)[1:2], c("Data", "Boot"))
  expect_identical(dim(summary_data$summary$direct), c(1L, 4L))
  expect_identical(rownames(summary_data$summary$direct), "X")
  expect_identical(colnames(summary_data$summary$direct)[1], "Estimate")
  # total effect
  expect_identical(dim(summary_boot$summary$total), c(1L, 5L))
  expect_identical(rownames(summary_boot$summary$total), "X")
  expect_identical(colnames(summary_boot$summary$total)[1:2], c("Data", "Boot"))
  expect_identical(dim(summary_data$summary$total), c(1L, 4L))
  expect_identical(rownames(summary_data$summary$total), "X")
  expect_identical(colnames(summary_data$summary$total)[1], "Estimate")

})

test_that("effect summaries contain correct coefficient values", {

  # effects computed on original sample
  expect_equivalent(summary_boot$summary$fit_mx[[1]]$coefficients[2, "Data"], boot$fit$a[1])
  expect_equivalent(summary_boot$summary$fit_mx[[2]]$coefficients[2, "Data"], boot$fit$a[2])
  expect_identical(summary_boot$summary$fit_ymx$coefficients[2:3, "Data"], boot$fit$b)
  expect_identical(summary_boot$summary$direct["X", "Data"], boot$fit$direct)
  expect_identical(summary_boot$summary$total["X", "Data"], boot$fit$total)
  expect_equivalent(summary_data$summary$fit_mx[[1]]$coefficients[2, "Estimate"], boot$fit$a[1])
  expect_equivalent(summary_data$summary$fit_mx[[2]]$coefficients[2, "Estimate"], boot$fit$a[2])
  expect_identical(summary_data$summary$fit_ymx$coefficients[2:3, "Estimate"], boot$fit$b)
  expect_identical(summary_data$summary$direct["X", "Estimate"], boot$fit$direct)
  expect_identical(summary_data$summary$total["X", "Estimate"], boot$fit$total)

  # bootstrapped effects
  expect_equivalent(summary_boot$summary$fit_mx[[1]]$coefficients[2, "Boot"], mean(boot$reps$t[, 5]))
  expect_equivalent(summary_boot$summary$fit_mx[[2]]$coefficients[2, "Boot"], mean(boot$reps$t[, 9]))
  expect_equivalent(summary_boot$summary$fit_ymx$coefficients[2:3, "Boot"], colMeans(boot$reps$t[, 13:14]))
  expect_equal(summary_boot$summary$direct["X", "Boot"], mean(boot$reps$t[, 15]))
  expect_equal(summary_boot$summary$total["X", "Boot"], mean(boot$reps$t[, 18]))

})

test_that("objects returned by setup_xxx_plot() have correct structure", {

  ## ci plot
  # check data frame for confidence interval
  expect_s3_class(ci$ci, "data.frame")
  # check dimensions
  expect_identical(dim(ci$ci), c(4L, 4L))
  # check column names
  column_names <- c("Effect", "Estimate", "Lower", "Upper")
  expect_named(ci$ci, column_names)
  # check that direct effect and indirect effect are plotted by default
  effect_names <- c("Direct", ab_names)
  expect_identical(ci$ci$Effect, factor(effect_names, levels = effect_names))
  # check confidence level
  expect_identical(ci$level, level)
  # check logical for multiple methods
  expect_false(ci$have_methods)

  ## density plot
  # check data frame for confidence interval
  expect_s3_class(density$density, "data.frame")
  # check dimensions
  expect_identical(ncol(density$density), 3L)
  expect_gt(nrow(density$density), 0L)
  # check column names
  column_names <- c("Effect", "ab", "Density")
  expect_named(density$density, column_names)
  # check data frame confidence interval
  expect_s3_class(density$ci, "data.frame")
  # check dimensions
  expect_identical(dim(density$ci), c(3L, 4L))
  # check column names
  column_names <- c("Effect", "Estimate", "Lower", "Upper")
  expect_named(density$ci, column_names)
  # check type of test
  expect_identical(density$test, "boot")
  # check confidence level
  expect_identical(density$level, level)
  # check logical for multiple effects
  expect_true(density$have_effect)
  # check logical for multiple methods
  expect_false(density$have_methods)

  ## ellipse_plot
  expect_identical(ellipse, setup_ellipse_plot(boot$fit))

})


## deprecated:

test_that("data returned by fortify() has correct structure", {

  ## dot plot
  # check dimensions
  expect_s3_class(dot_deprecated, "data.frame")
  expect_identical(dim(dot_deprecated), c(4L, 4L))
  # check column names
  column_names <- c("Effect", "Point", "Lower", "Upper")
  expect_named(dot_deprecated, column_names)
  # check that direct effect and indirect effect are plotted by default
  effect_names <- c("Direct", ab_names)
  expect_identical(dot_deprecated$Effect, factor(effect_names, levels = effect_names))

  ## density plot
  # check dimensions
  expect_s3_class(density_deprecated, "data.frame")
  expect_identical(ncol(density_deprecated), 3L)
  expect_gt(nrow(density_deprecated), 0L)
  # check column names
  column_names <- c("ab", "Density", "Effect")
  expect_named(density_deprecated, column_names)

})

test_that("data returned by fortify() has correct attributes", {

  ## dot plot
  # check aesthetic mapping
  mapping <- aes_string(x = "Effect", y = "Point",
                        ymin = "Lower", ymax = "Upper")
  expect_equal(attr(dot_deprecated, "mapping"), mapping)
  # check default geom()
  expect_identical(attr(dot_deprecated, "geom"), geom_pointrange)
  # check facets
  expect_null(attr(dot_deprecated, "facets"))
  # check that method is stored correctly
  expect_identical(attr(dot_deprecated, "method"), "dot")

  ## density plot
  # check aesthetic mapping
  mapping <- aes_string(x = "ab", y = "Density")
  expect_equal(attr(density_deprecated, "mapping"), mapping)
  # check default geom()
  expect_equal(attr(density_deprecated, "geom"), robmed:::geom_densityline)
  # check facets
  expect_equal(attr(density_deprecated, "facets"), ~Effect)
  # check title
  expect_identical(attr(density_deprecated, "main"), "Bootstrap distribution")
  # check confidence interval
  ci <- attr(density_deprecated, "ci")
  expect_s3_class(ci, "data.frame")
  expect_identical(dim(ci), c(3L, 5L))
  expect_named(ci, c("ab", "Density", "Lower", "Upper", "Effect"))
  # check that method is stored correctly
  expect_identical(attr(density_deprecated, "method"), "density")

})


# run mediation analysis through formula interface with data argument
set.seed(seed)
boot_f1 <- test_mediation(Y ~ m(M1, M2) + X + covariates(C1, C2),
                          data = test_data, test = "boot", R = R,
                          level = 0.9, type = "bca", method = "regression",
                          robust = FALSE)
# run mediation analysis through formula interface without data argument
set.seed(seed)
boot_f2 <- test_mediation(Y ~ m(M1, M2) + X + covariates(C1, C2),
                          test = "boot", R = R, level = 0.9, type = "bca",
                          method = "regression", robust = FALSE)
# define mediators and covariates outside formula
med <- m(M1, M2)
cov <- covariates(C1, C2)
set.seed(seed)
boot_f3 <- test_mediation(Y ~ med + X + cov, data = test_data,
                          test = "boot", R = R, level = 0.9, type = "bca",
                          method = "regression", robust = FALSE)


test_that("formula interface works correctly", {

  # check that results are the same as with default method
  expect_equal(boot_f1, boot)
  expect_equal(boot_f2, boot)
  expect_equal(boot_f3, boot)

})
