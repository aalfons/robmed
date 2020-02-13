context("bootstrap test: median regression, single mediator, covariates")


## load package
library("robmed", quietly = TRUE)

## control parameters
n <- 250          # number of observations
a <- c <- 0.2     # true effects
b <- 0            # true effect
R <- 1000         # number of bootstrap samples
seed <- 20190201  # seed for the random number generator

## set seed for reproducibility
set.seed(seed)

## generate data
X <- rnorm(n)
M <- a * X + rnorm(n)
Y <- b * M + c * X + rnorm(n)
C1 <- rnorm(n)
C2 <- rnorm(n)
test_data <- data.frame(X, Y, M, C1, C2)

## run bootstrap test
set.seed(seed)
level <- 0.9
boot <- test_mediation(test_data, x = "X", y = "Y", m = "M",
                       covariates = c("C1", "C2"), test = "boot", R = R,
                       level = level, type = "bca", method = "regression",
                       robust = "median")

## compute summary
summary_boot <- summary(boot, type = "boot")
summary_data <- summary(boot, type = "data")

## create data for plotting
ci <- setup_ci_plot(boot)
density <- setup_density_plot(boot)
# deprecated:
dot_deprecated <- suppressWarnings(fortify(boot, method = "dot"))
density_deprecated <- suppressWarnings(fortify(boot, method = "density"))

## stuff needed to check correctness
coef_names <- c("a", "b", "Direct", "Total", "ab")
mx_names <- c("(Intercept)", "X", "C1", "C2")
ymx_names <- c("(Intercept)", "M", "X", "C1", "C2")


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
  expect_identical(boot$level, level)
  # type of confidence intervals
  expect_identical(boot$type, "bca")
  # variable names
  expect_identical(boot$fit$x, "X")
  expect_identical(boot$fit$y, "Y")
  expect_identical(boot$fit$m, "M")
  expect_identical(boot$fit$covariates, c("C1", "C2"))
  # robust fit and test
  expect_identical(boot$fit$robust, "median")
  expect_null(boot$fit$control)

})

test_that("dimensions are correct", {

  # effects are scalars
  expect_length(boot$ab, 1L)
  # only one indirect effect, so only one confidence interval
  expect_length(boot$ci, 2L)
  # dimensions of bootstrap replicates
  d_boot <- dim(boot$reps$t)
  expect_identical(d_boot, c(as.integer(R), 11L))

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
  expect_equivalent(coef(boot, parm = "b", type = "boot"), mean(boot$reps$t[, 7]))
  expect_equivalent(coef(boot, parm = "Direct", type = "boot"), mean(boot$reps$t[, 8]))
  expect_equivalent(coef(boot, parm = "Total", type = "boot"), mean(boot$reps$t[, 11]))
  expect_equivalent(coef(boot, parm = "ab", type = "boot"), boot$ab)

  # effects computed on original sample
  expect_equivalent(coef(boot, parm = "a", type = "data"), boot$fit$a)
  expect_equivalent(coef(boot, parm = "b", type = "data"), boot$fit$b)
  expect_equivalent(coef(boot, parm = "Direct", type = "data"), boot$fit$direct)
  expect_equivalent(coef(boot, parm = "Total", type = "data"), boot$fit$total)
  expect_equivalent(coef(boot, parm = "ab", type = "data"), boot$fit$a * boot$fit$b)

})

test_that("output of confint() method has correct attributes", {

  ci_boot <- confint(boot, type = "boot")
  ci_data <- confint(boot, type = "data")
  # bootstrapped confidence intervals
  expect_equal(dim(ci_boot), c(5L, 2L))
  expect_equal(rownames(ci_boot), coef_names)
  expect_equal(colnames(ci_boot), c("5 %", "95 %"))
  # confidence intervals based on theory (except for indirect effect)
  expect_equal(dim(ci_data), c(5L, 2L))
  expect_equal(rownames(ci_data), coef_names)
  expect_equal(colnames(ci_data), c("5 %", "95 %"))

})

test_that("confint() method returns correct values of confidence intervals", {

  # bootstrapped confidence intervals
  expect_equivalent(confint(boot, parm = "ab", type = "boot"), boot$ci)
  # confidence intervals based on theory (except for indirect effect)
  expect_equivalent(confint(boot, parm = "ab", type = "data"), boot$ci)

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
  # summary for model m ~ x
  expect_s3_class(summary_boot$summary$fit_mx, "summary_rq")
  # summary for model y ~ m + x
  expect_s3_class(summary_boot$summary$fit_ymx, "summary_rq")
  # regression standard error for model y ~ m + x
  expect_null(summary_boot$summary$fit_ymx$s)
  expect_null(summary_data$summary$fit_ymx$s)
  # R-squared for model y ~ m + x
  expect_null(summary_boot$summary$fit_ymx$R2)
  expect_null(summary_data$summary$fit_ymx$R2)
  # F-test for model y ~ m + x
  expect_null(summary_boot$summary$fit_ymx$F_test)
  expect_null(summary_data$summary$fit_ymx$F_test)

})

test_that("attributes are correctly passed through summary", {

  # robustness
  expect_identical(summary_boot$summary$robust, "median")
  expect_identical(summary_data$summary$robust, "median")
  # number of observations
  expect_identical(summary_boot$summary$n, as.integer(n))
  expect_identical(summary_data$summary$n, as.integer(n))
  # variable names
  expect_identical(summary_boot$summary$x, "X")
  expect_identical(summary_boot$summary$y, "Y")
  expect_identical(summary_boot$summary$m, "M")
  expect_identical(summary_boot$summary$covariates, c("C1", "C2"))
  expect_identical(summary_data$summary$x, "X")
  expect_identical(summary_data$summary$y, "Y")
  expect_identical(summary_data$summary$m, "M")
  expect_identical(summary_data$summary$covariates, c("C1", "C2"))

})

test_that("effect summaries have correct names", {

  # a path
  expect_identical(dim(summary_boot$summary$fit_mx$coefficients), c(4L, 5L))
  expect_identical(rownames(summary_boot$summary$fit_mx$coefficients), mx_names)
  expect_identical(colnames(summary_boot$summary$fit_mx$coefficients)[1:2], c("Data", "Boot"))
  expect_identical(dim(summary_data$summary$fit_mx$coefficients), c(4L, 4L))
  expect_identical(rownames(summary_data$summary$fit_mx$coefficients), mx_names)
  expect_identical(colnames(summary_data$summary$fit_mx$coefficients)[1], "Estimate")
  # b path
  expect_identical(dim(summary_boot$summary$fit_ymx$coefficients), c(5L, 5L))
  expect_identical(rownames(summary_boot$summary$fit_ymx$coefficient), ymx_names)
  expect_identical(colnames(summary_boot$summary$fit_ymx$coefficient)[1:2], c("Data", "Boot"))
  expect_identical(dim(summary_data$summary$fit_ymx$coefficient), c(5L, 4L))
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
  expect_equivalent(summary_boot$summary$fit_mx$coefficients[2, "Data"], boot$fit$a)
  expect_identical(summary_boot$summary$fit_ymx$coefficients[2, "Data"], boot$fit$b)
  expect_identical(summary_boot$summary$direct["X", "Data"], boot$fit$direct)
  expect_identical(summary_boot$summary$total["X", "Data"], boot$fit$total)
  expect_equivalent(summary_data$summary$fit_mx$coefficients[2, "Estimate"], boot$fit$a)
  expect_identical(summary_data$summary$fit_ymx$coefficients[2, "Estimate"], boot$fit$b)
  expect_identical(summary_data$summary$direct["X", "Estimate"], boot$fit$direct)
  expect_identical(summary_data$summary$total["X", "Estimate"], boot$fit$total)

  # bootstrapped effects
  expect_equivalent(summary_boot$summary$fit_mx$coefficients[2, "Boot"], mean(boot$reps$t[, 3]))
  expect_equivalent(summary_boot$summary$fit_ymx$coefficients[2, "Boot"], mean(boot$reps$t[, 7]))
  expect_equal(summary_boot$summary$direct["X", "Boot"], mean(boot$reps$t[, 8]))
  expect_equal(summary_boot$summary$total["X", "Boot"], mean(boot$reps$t[, 11]))

})

test_that("objects returned by setup_xxx_plot() have correct structure", {

  ## ci plot
  # check data frame for confidence interval
  expect_s3_class(ci$ci, "data.frame")
  # check dimensions
  expect_identical(dim(ci$ci), c(2L, 4L))
  # check column names
  column_names <- c("Effect", "Estimate", "Lower", "Upper")
  expect_named(ci$ci, column_names)
  # check that direct effect and indirect effect are plotted by default
  effect_names <- c("Direct", "ab")
  expect_identical(ci$ci$Effect, factor(effect_names, levels = effect_names))
  # check confidence level
  expect_identical(ci$level, level)
  # check logical for multiple methods
  expect_false(ci$have_methods)

  ## density plot
  # check data frame for confidence interval
  expect_s3_class(density$density, "data.frame")
  # check dimensions
  expect_identical(ncol(density$density), 2L)
  expect_gt(nrow(density$density), 0L)
  # check column names
  column_names <- c("ab", "Density")
  expect_named(density$density, column_names)
  # check data frame confidence interval
  expect_s3_class(density$ci, "data.frame")
  # check dimensions
  expect_identical(dim(density$ci), c(1L, 3L))
  # check column names
  column_names <- c("Estimate", "Lower", "Upper")
  expect_named(density$ci, column_names)
  # check type of test
  expect_identical(density$test, "boot")
  # check confidence level
  expect_identical(density$level, level)
  # check logical for multiple effects
  expect_false(density$have_effect)
  # check logical for multiple methods
  expect_false(density$have_methods)

  ## ellipse_plot
  expect_error(setup_ellipse_plot(boot))

})


## deprecated:

test_that("data returned by fortify() has correct structure", {

  ## dot plot
  # check dimensions
  expect_s3_class(dot_deprecated, "data.frame")
  expect_identical(dim(dot_deprecated), c(2L, 4L))
  # check column names
  column_names <- c("Effect", "Point", "Lower", "Upper")
  expect_named(dot_deprecated, column_names)
  # check that direct effect and indirect effect are plotted by default
  effect_names <- c("Direct", "ab")
  expect_identical(dot_deprecated$Effect, factor(effect_names, levels = effect_names))

  ## density plot
  # check dimensions
  expect_s3_class(density_deprecated, "data.frame")
  expect_identical(ncol(density_deprecated), 2L)
  expect_gt(nrow(density_deprecated), 0L)
  # check column names
  column_names <- c("ab", "Density")
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
  expect_null(attr(density_deprecated, "facets"))
  # check title
  expect_identical(attr(density_deprecated, "main"), "Bootstrap distribution")
  # check confidence interval
  ci <- attr(density_deprecated, "ci")
  expect_s3_class(ci, "data.frame")
  expect_identical(dim(ci), c(1L, 4L))
  expect_named(ci, c("ab", "Density", "Lower", "Upper"))
  # check that method is stored correctly
  expect_identical(attr(density_deprecated, "method"), "density")

})


# run mediation analysis through formula interface with data argument
set.seed(seed)
boot_f1 <- test_mediation(Y ~ m(M) + X + covariates(C1, C2), data = test_data,
                          test = "boot", R = R, level = 0.9, type = "bca",
                          method = "regression", robust = "median")
# run mediation analysis through formula interface without data argument
set.seed(seed)
boot_f2 <- test_mediation(Y ~ m(M) + X + covariates(C1, C2),
                          test = "boot", R = R, level = 0.9, type = "bca",
                          method = "regression", robust = "median")
# define mediator and covariates outside formula
med <- m(M)
cov <- covariates(C1, C2)
set.seed(seed)
boot_f3 <- test_mediation(Y ~ med + X + cov, data = test_data,
                          test = "boot", R = R, level = 0.9, type = "bca",
                          method = "regression", robust = "median")


test_that("formula interface works correctly", {

  # check that results are the same as with default method
  expect_equal(boot_f1, boot)
  expect_equal(boot_f2, boot)
  expect_equal(boot_f3, boot)

})
