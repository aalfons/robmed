context("standard bootstrap test: regression, single mediator, covariates")


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
C1 <- rnorm(n)
C2 <- rnorm(n)
test_data <- data.frame(X, Y, M, C1, C2)

## run bootstrap test
set.seed(seed)
level <- c(0.9, 0.95)
boot <- test_mediation(test_data, x = "X", y = "Y", m = "M",
                       covariates = c("C1", "C2"), test = "boot", R = R,
                       level = level[1], type = "bca", method = "regression",
                       robust = FALSE)

## compute summary
summary_boot <- summary(boot, type = "boot")
summary_data <- summary(boot, type = "data")

## retest with different parameters
boot_less <- retest(boot, alternative = "less", level = level[2])
boot_greater <- retest(boot, alternative = "greater", level = level[2])
boot_perc <- retest(boot, type = "perc", level = level[2])

## create data for plotting
ci <- setup_ci_plot(boot)
ci_perc <- setup_ci_plot(boot_perc, p_value = TRUE)
density <- setup_density_plot(boot)
ellipse <- setup_ellipse_plot(boot)

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
  expect_identical(boot$R, as.integer(R))
  # confidence level
  expect_identical(boot$level, level[1])
  # type of confidence intervals
  expect_identical(boot$type, "bca")
  # variable names
  expect_identical(boot$fit$x, "X")
  expect_identical(boot$fit$y, "Y")
  expect_identical(boot$fit$m, "M")
  expect_identical(boot$fit$covariates, c("C1", "C2"))
  # nonrobust fit and test
  expect_false(boot$fit$robust)
  expect_identical(boot$fit$family, "gaussian")
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
  expect_equivalent(coef(boot, parm = "a", type = "boot"),
                    mean(boot$reps$t[, 3]))
  expect_equivalent(coef(boot, parm = "b", type = "boot"),
                    mean(boot$reps$t[, 7]))
  expect_equivalent(coef(boot, parm = "Direct", type = "boot"),
                    mean(boot$reps$t[, 8]))
  expect_equivalent(coef(boot, parm = "Total", type = "boot"),
                    mean(boot$reps$t[, 11]))
  expect_equivalent(coef(boot, parm = "ab", type = "boot"),
                    boot$ab)

  # effects computed on original sample
  expect_equivalent(coef(boot, parm = "a", type = "data"),
                    boot$fit$a)
  expect_equivalent(coef(boot, parm = "b", type = "data"),
                    boot$fit$b)
  expect_equivalent(coef(boot, parm = "Direct", type = "data"),
                    boot$fit$direct)
  expect_equivalent(coef(boot, parm = "Total", type = "data"),
                    boot$fit$total)
  expect_equivalent(coef(boot, parm = "ab", type = "data"),
                    boot$fit$a * boot$fit$b)

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
  expect_s3_class(summary_boot$summary$fit_mx, "summary_lm")
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
  expect_named(summary_boot$summary$fit_ymx$F_test,
               c("statistic", "df", "p_value"))
  df_test_boot <- summary_boot$summary$fit_ymx$F_test$df
  expect_identical(df_test_boot[1], 4L)
  expect_identical(df_test_boot[2], summary_boot$summary$fit_ymx$s$df)
  expect_type(summary_data$summary$fit_ymx$F_test, "list")
  expect_named(summary_data$summary$fit_ymx$F_test,
               c("statistic", "df", "p_value"))
  df_test_data <- summary_data$summary$fit_ymx$F_test$df
  expect_identical(df_test_data[1], 4L)
  expect_identical(df_test_data[2], summary_data$summary$fit_ymx$s$df)
  # no plot is created
  expect_null(summary_boot$plot)
  expect_null(summary_data$plot)

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
  expect_identical(summary_boot$summary$m, "M")
  expect_identical(summary_boot$summary$covariates, c("C1", "C2"))
  expect_identical(summary_data$summary$x, "X")
  expect_identical(summary_data$summary$y, "Y")
  expect_identical(summary_data$summary$m, "M")
  expect_identical(summary_data$summary$covariates, c("C1", "C2"))

})

test_that("effect summaries have correct names", {

  # a path
  expect_identical(dim(summary_boot$summary$fit_mx$coefficients),
                   c(4L, 5L))
  expect_identical(rownames(summary_boot$summary$fit_mx$coefficients),
                   mx_names)
  expect_identical(colnames(summary_boot$summary$fit_mx$coefficients)[1:2],
                   c("Data", "Boot"))
  expect_identical(dim(summary_data$summary$fit_mx$coefficients),
                   c(4L, 4L))
  expect_identical(rownames(summary_data$summary$fit_mx$coefficients),
                   mx_names)
  expect_identical(colnames(summary_data$summary$fit_mx$coefficients)[1],
                   "Estimate")
  # b path
  expect_identical(dim(summary_boot$summary$fit_ymx$coefficients),
                   c(5L, 5L))
  expect_identical(rownames(summary_boot$summary$fit_ymx$coefficient),
                   ymx_names)
  expect_identical(colnames(summary_boot$summary$fit_ymx$coefficient)[1:2],
                   c("Data", "Boot"))
  expect_identical(dim(summary_data$summary$fit_ymx$coefficient),
                   c(5L, 4L))
  expect_identical(rownames(summary_data$summary$fit_ymx$coefficient),
                   ymx_names)
  expect_identical(colnames(summary_data$summary$fit_ymx$coefficient)[1],
                   "Estimate")
  # direct effect
  expect_identical(dim(summary_boot$summary$direct),
                   c(1L, 5L))
  expect_identical(rownames(summary_boot$summary$direct),
                   "X")
  expect_identical(colnames(summary_boot$summary$direct)[1:2],
                   c("Data", "Boot"))
  expect_identical(dim(summary_data$summary$direct),
                   c(1L, 4L))
  expect_identical(rownames(summary_data$summary$direct),
                   "X")
  expect_identical(colnames(summary_data$summary$direct)[1],
                   "Estimate")
  # total effect
  expect_identical(dim(summary_boot$summary$total),
                   c(1L, 5L))
  expect_identical(rownames(summary_boot$summary$total),
                   "X")
  expect_identical(colnames(summary_boot$summary$total)[1:2],
                   c("Data", "Boot"))
  expect_identical(dim(summary_data$summary$total),
                   c(1L, 4L))
  expect_identical(rownames(summary_data$summary$total),
                   "X")
  expect_identical(colnames(summary_data$summary$total)[1],
                   "Estimate")

})

test_that("effect summaries contain correct coefficient values", {

  # effects computed on original sample
  expect_equivalent(summary_boot$summary$fit_mx$coefficients[2, "Data"],
                    boot$fit$a)
  expect_identical(summary_boot$summary$fit_ymx$coefficients[2, "Data"],
                   boot$fit$b)
  expect_identical(summary_boot$summary$direct["X", "Data"],
                   boot$fit$direct)
  expect_identical(summary_boot$summary$total["X", "Data"],
                   boot$fit$total)
  expect_equivalent(summary_data$summary$fit_mx$coefficients[2, "Estimate"],
                    boot$fit$a)
  expect_identical(summary_data$summary$fit_ymx$coefficients[2, "Estimate"],
                   boot$fit$b)
  expect_identical(summary_data$summary$direct["X", "Estimate"],
                   boot$fit$direct)
  expect_identical(summary_data$summary$total["X", "Estimate"],
                   boot$fit$total)

  # bootstrapped effects
  expect_equivalent(summary_boot$summary$fit_mx$coefficients[2, "Boot"],
                    mean(boot$reps$t[, 3]))
  expect_equivalent(summary_boot$summary$fit_ymx$coefficients[2, "Boot"],
                    mean(boot$reps$t[, 7]))
  expect_equal(summary_boot$summary$direct["X", "Boot"],
               mean(boot$reps$t[, 8]))
  expect_equal(summary_boot$summary$total["X", "Boot"],
               mean(boot$reps$t[, 11]))

})

test_that("output of retest() has correct structure", {

  # bootstrap test
  expect_identical(class(boot_less), class(boot))
  expect_identical(class(boot_greater), class(boot))
  expect_identical(class(boot_perc), class(boot))
  # regression fit
  expect_identical(boot_less$fit, boot$fit)
  expect_identical(boot_greater$fit, boot$fit)
  expect_identical(boot_perc$fit, boot$fit)
  # bootstrap replicates
  expect_identical(boot_less$reps, boot$reps)
  expect_identical(boot_greater$reps, boot$reps)
  expect_identical(boot_perc$reps, boot$reps)

})

test_that("arguments of retest() are correctly passed", {

  # alternative hypothesis
  expect_identical(boot_less$alternative, "less")
  expect_identical(boot_greater$alternative, "greater")
  expect_identical(boot_perc$alternative, "twosided")
  # confidence level
  expect_identical(boot_less$level, level[2])
  expect_identical(boot_greater$level, level[2])
  expect_identical(boot_perc$level, level[2])
  # type of confidence intervals
  expect_identical(boot_less$type, "bca")
  expect_identical(boot_greater$type, "bca")
  expect_identical(boot_perc$type, "perc")
  # indirect effect
  expect_identical(boot_less$ab, boot$ab)
  expect_identical(boot_greater$ab, boot$ab)
  expect_identical(boot_perc$ab, boot$ab)
  # confidence interval for alternative = "less"
  expect_identical(length(boot_less$ci), length(boot$ci))
  expect_equivalent(boot_less$ci[1], -Inf)
  expect_equal(boot_less$ci[2], boot$ci[2])
  expect_identical(colnames(confint(boot_less)), c("Lower", "Upper"))
  # confidence interval for alternative = "greater"
  expect_identical(length(boot_greater$ci), length(boot$ci))
  expect_equal(boot_greater$ci[1], boot$ci[1])
  expect_equivalent(boot_greater$ci[2], Inf)
  expect_identical(colnames(confint(boot_greater)), c("Lower", "Upper"))
  # confidence interval for type = "perc"
  expect_identical(length(boot_perc$ci), length(boot$ci))

})

test_that("output of p_value() method has correct attributes", {

  digits <- 3
  p_boot <- p_value(boot_perc, type = "boot", digits = digits)
  p_data <- p_value(boot_perc, type = "data", digits = digits)
  # bootstrapped effects
  expect_length(p_boot, 5L)
  expect_named(p_boot, coef_names)
  expect_equal(p_boot["ab"], round(p_boot["ab"], digits = digits))
  # effects computed on original sample
  expect_length(p_data, 5L)
  expect_named(p_data, coef_names)
  expect_equal(p_data["ab"], round(p_data["ab"], digits = digits))

})

test_that("objects returned by setup_xxx_plot() have correct structure", {

  ## ci plot without p-value
  # check data frame for confidence interval
  expect_s3_class(ci$ci, "data.frame")
  # check dimensions
  expect_identical(dim(ci$ci), c(2L, 4L))
  # check column names
  expect_named(ci$ci, c("Effect", "Estimate", "Lower", "Upper"))
  # check that direct effect and indirect effect are plotted by default
  effect_names <- c("Direct", "ab")
  effect_factor <- factor(effect_names, levels = effect_names)
  expect_identical(ci$ci$Effect, effect_factor)
  # check confidence level
  expect_identical(ci$level, level[1])
  # check logical for multiple methods
  expect_false(ci$have_methods)

  ## ci plot with p-value
  # check data frame for confidence interval and p-value
  expect_s3_class(ci_perc$ci, "data.frame")
  expect_s3_class(ci_perc$p_value, "data.frame")
  # check dimensions
  expect_identical(dim(ci_perc$ci), c(2L, 5L))
  expect_identical(dim(ci_perc$p_value), c(2L, 3L))
  # check column names
  expect_named(ci_perc$ci, c("Label", "Effect", "Estimate", "Lower", "Upper"))
  expect_named(ci_perc$p_value, c("Label", "Effect", "Value"))
  # check that labels are correct
  label_names <- c("Confidence interval", "p-Value")
  expect_identical(ci_perc$ci$Label,
                   factor(rep.int(label_names[1], 2), levels = label_names))
  expect_identical(ci_perc$p_value$Label,
                   factor(rep.int(label_names[2], 2), levels = label_names))
  # check that direct effect and indirect effect are plotted by default
  effect_names <- c("Direct", "ab")
  effect_factor <- factor(effect_names, levels = effect_names)
  expect_identical(ci_perc$ci$Effect, effect_factor)
  expect_identical(ci_perc$p_value$Effect, effect_factor)
  # check confidence level
  expect_identical(ci$level, level[1])
  # check logical for multiple methods
  expect_false(ci$have_methods)

  ## density plot
  # check data frame for confidence interval
  expect_s3_class(density$density, "data.frame")
  # check dimensions
  expect_identical(ncol(density$density), 2L)
  expect_gt(nrow(density$density), 0L)
  # check column names
  expect_named(density$density, c("ab", "Density"))
  # check data frame confidence interval
  expect_s3_class(density$ci, "data.frame")
  # check dimensions
  expect_identical(dim(density$ci), c(1L, 3L))
  # check column names
  expect_named(density$ci, c("Estimate", "Lower", "Upper"))
  # check type of test
  expect_identical(density$test, "boot")
  # check confidence level
  expect_identical(density$level, level[1])
  # check logical for multiple effects
  expect_false(density$have_effect)
  # check logical for multiple methods
  expect_false(density$have_methods)

  ## ellipse plot
  expect_identical(ellipse, setup_ellipse_plot(boot$fit))

  ## weight plot
  expect_error(setup_weight_plot(boot))

})


# run mediation analysis through formula interface with data argument
set.seed(seed)
boot_f1 <- test_mediation(Y ~ m(M) + X + covariates(C1, C2),
                          data = test_data, test = "boot", R = R,
                          level = 0.9, type = "bca", method = "regression",
                          robust = FALSE)
# run mediation analysis through formula interface without data argument
set.seed(seed)
boot_f2 <- test_mediation(Y ~ m(M) + X + covariates(C1, C2),
                          test = "boot", R = R, level = 0.9, type = "bca",
                          method = "regression", robust = FALSE)
# define mediator and covariates outside formula
med <- m(M)
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
