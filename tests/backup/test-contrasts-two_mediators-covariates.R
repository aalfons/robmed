context("contrasts: two mediators, covariates")


## load package
library("robmed", quietly = TRUE)

## control parameters
n <- 250          # number of observations
a <- c <- 0.2     # true effects
b <- 0            # true effect
R <- 1000         # number of bootstrap samples
seed <- 20210407  # seed for the random number generator

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

## control parameters
level <- 0.95
ctrl <- reg_control(max_iterations = 500)

## run bootstrap tests
# without contrasts
set.seed(seed)
boot_without <- test_mediation(test_data, x = "X", y = "Y", m = c("M1", "M2"),
                               covariates = c("C1", "C2"), test = "boot",
                               R = R, level = level, type = "perc",
                               method = "regression", robust = TRUE,
                               contrast = FALSE, control = ctrl)
# re-test with contrasts
reboot_estimates <- retest(boot_without, contrast = "estimates")
reboot_absolute <- retest(boot_without, contrast = "absolute")
# test differences of indirect effect
set.seed(seed)
boot_estimates <- test_mediation(test_data, x = "X", y = "Y", m = c("M1", "M2"),
                                 covariates = c("C1", "C2"), test = "boot",
                                 R = R, level = level, type = "perc",
                                 method = "regression", robust = TRUE,
                                 contrast = "estimates", control = ctrl)
# test differences of absolute values of indirect effect
set.seed(seed)
boot_absolute <- test_mediation(test_data, x = "X", y = "Y", m = c("M1", "M2"),
                                covariates = c("C1", "C2"), test = "boot",
                                R = R, level = level, type = "perc",
                                method = "regression", robust = TRUE,
                                contrast = "absolute", control = ctrl)

## create data for plotting
ci_without <- setup_ci_plot(boot_without)
ci_estimates <- setup_ci_plot(boot_estimates)
ci_absolute <- setup_ci_plot(boot_absolute)
ci_p_without <- setup_ci_plot(boot_without, p_value = TRUE)
ci_p_estimates <- setup_ci_plot(boot_estimates, p_value = TRUE)
ci_p_absolute <- setup_ci_plot(boot_absolute, p_value = TRUE)
density_without <- setup_density_plot(boot_without)
density_estimates <- setup_density_plot(boot_estimates)
density_absolute <- setup_density_plot(boot_absolute)
ellipse_without <- setup_ellipse_plot(boot_without)
ellipse_estimates <- setup_ellipse_plot(boot_estimates)
ellipse_absolute <- setup_ellipse_plot(boot_absolute)
weight_without <- setup_weight_plot(boot_without)
weight_estimates <- setup_weight_plot(boot_estimates)
weight_absolute <- setup_weight_plot(boot_absolute)

## stuff needed to check correctness
indirect_names <- c("Total", "M1", "M2", "Contrast")
ab_names <- paste("ab", indirect_names, sep = "_")


## run tests

test_that("output has correct structure", {

  # bootstrap test
  expect_s3_class(boot_estimates, "boot_test_mediation")
  expect_s3_class(boot_estimates, "test_mediation")
  expect_s3_class(boot_absolute, "boot_test_mediation")
  expect_s3_class(boot_absolute, "test_mediation")

  # regression fit
  expect_s3_class(boot_estimates$fit, "reg_fit_mediation")
  expect_s3_class(boot_estimates$fit, "fit_mediation")
  expect_s3_class(boot_absolute$fit, "reg_fit_mediation")
  expect_s3_class(boot_absolute$fit, "fit_mediation")

})

test_that("model fit works correctly", {

  # extract model fits
  fit_without <- boot_without$fit
  fit_estimates <- boot_estimates$fit
  fit_absolute <- boot_absolute$fit

  # relevant components are the same as when no contrasts are computed
  keep <- setdiff(names(fit_without), c("ab", "contrast"))
  expect_equal(fit_estimates[keep], fit_without[keep])
  expect_equal(fit_absolute[keep], fit_without[keep])

  # multiple indirect effects have correct dimensions and names
  expect_length(fit_estimates$ab, 4L)
  expect_identical(names(fit_estimates$ab), indirect_names)
  expect_length(fit_absolute$ab, 4L)
  expect_identical(names(fit_absolute$ab), indirect_names)

  # contrasts computed correctly
  expect_equivalent(fit_estimates$ab["Contrast"],
                    fit_estimates$ab["M1"] - fit_estimates$ab["M2"])
  expect_equivalent(fit_absolute$ab["Contrast"],
                    abs(fit_absolute$ab["M1"]) - abs(fit_absolute$ab["M2"]))

  # type of contrasts stored correctly
  expect_identical(fit_estimates$contrast, "estimates")
  expect_identical(fit_absolute$contrast, "absolute")

})

test_that("bootstrap test works correctly", {

  # relevant components are the same as when no contrasts are computed
  # ("boot" object should also be the same but needs to be checked separately
  # because of the function that is stored)
  keep <- setdiff(names(boot_without), c("ab", "ci", "fit", "reps"))
  expect_equal(boot_estimates[keep], boot_without[keep])
  expect_equal(boot_absolute[keep], boot_without[keep])

  # check that "boot" object is the same as when no contrasts are computed
  keep <- setdiff(names(boot_without$reps), c("statistic"))
  expect_equal(boot_estimates$reps[keep], boot_without$reps[keep])
  expect_equal(boot_absolute$reps[keep], boot_without$reps[keep])

  # multiple indirect effects have correct dimensions and names
  expect_length(boot_estimates$ab, 4L)
  expect_named(boot_estimates$ab, indirect_names)
  expect_length(boot_absolute$ab, 4L)
  expect_named(boot_absolute$ab, indirect_names)

  # confidence intervals have correct dimensions and names
  expect_identical(dim(boot_estimates$ci), c(4L, 2L))
  expect_identical(rownames(boot_estimates$ci), indirect_names)
  expect_identical(dim(boot_absolute$ci), c(4L, 2L))
  expect_identical(rownames(boot_absolute$ci), indirect_names)

})

test_that("retest() works correctly for two mediators", {

  # object should be the same as when argument is supplied from the start
  # ("boot" object should also be the same but needs to be checked separately
  # because of the function that is stored)
  keep <- setdiff(names(boot_estimates), "reps")
  expect_equal(reboot_estimates[keep], boot_estimates[keep])
  keep <- setdiff(names(boot_absolute), "reps")
  expect_equal(reboot_absolute[keep], boot_absolute[keep])

  # check that "boot" object is the same
  keep <- setdiff(names(boot_estimates$reps), "statistic")
  expect_equal(reboot_estimates$reps[keep], boot_estimates$reps[keep])
  keep <- setdiff(names(boot_absolute$reps), "statistic")
  expect_equal(reboot_absolute$reps[keep], boot_absolute$reps[keep])

})


test_that("objects returned by setup_xxx_plot() have correct structure", {

  ## ci plot without p-value
  # relevant components are the same as when no contrasts are computed
  keep <- setdiff(names(ci_without), "ci")
  expect_equal(ci_estimates[keep], ci_without[keep])
  expect_equal(ci_absolute[keep], ci_without[keep])
  # check data frame for confidence interval
  expect_s3_class(ci_estimates$ci, "data.frame")
  expect_s3_class(ci_absolute$ci, "data.frame")
  # check dimensions
  expect_identical(dim(ci_estimates$ci), c(5L, 4L))
  expect_identical(dim(ci_absolute$ci), c(5L, 4L))
  # check that direct effect and indirect effect are plotted by default
  effect_names <- c("Direct", ab_names)
  effect_factor <- factor(effect_names, levels = effect_names)
  expect_identical(ci_estimates$ci$Effect, effect_factor)
  expect_identical(ci_absolute$ci$Effect, effect_factor)
  # check that values for direct and indirect effects are correct
  expect_equivalent(ci_estimates$ci[-5L, ], ci_without$ci)
  expect_equivalent(ci_absolute$ci[-5L, ], ci_without$ci)


  ## ci plot with p-value
  # relevant components are the same as when no contrasts are computed
  keep <- setdiff(names(ci_without), c("ci", "p_value"))
  expect_equal(ci_estimates[keep], ci_without[keep])
  expect_equal(ci_absolute[keep], ci_without[keep])
  # check data frame for confidence interval and p-value
  expect_s3_class(ci_p_estimates$ci, "data.frame")
  expect_s3_class(ci_p_estimates$p_value, "data.frame")
  expect_s3_class(ci_p_absolute$ci, "data.frame")
  expect_s3_class(ci_p_absolute$p_value, "data.frame")
  # check dimensions
  expect_identical(dim(ci_p_estimates$ci), c(5L, 5L))
  expect_identical(dim(ci_p_estimates$p_value), c(5L, 3L))
  expect_identical(dim(ci_p_absolute$ci), c(5L, 5L))
  expect_identical(dim(ci_p_absolute$p_value), c(5L, 3L))
  # check that labels are correct
  label_names <- c("Confidence interval", "p-Value")
  expect_identical(ci_p_estimates$ci$Label,
                   factor(rep.int(label_names[1], 5), levels = label_names))
  expect_identical(ci_p_estimates$p_value$Label,
                   factor(rep.int(label_names[2], 5), levels = label_names))
  expect_identical(ci_p_absolute$ci$Label,
                   factor(rep.int(label_names[1], 5), levels = label_names))
  expect_identical(ci_p_absolute$p_value$Label,
                   factor(rep.int(label_names[2], 5), levels = label_names))
  # check that direct effect and indirect effect are plotted by default
  effect_names <- c("Direct", ab_names)
  effect_factor <- factor(effect_names, levels = effect_names)
  expect_identical(ci_p_estimates$ci$Effect, effect_factor)
  expect_identical(ci_p_estimates$p_value$Effect, effect_factor)
  expect_identical(ci_p_absolute$ci$Effect, effect_factor)
  expect_identical(ci_p_absolute$p_value$Effect, effect_factor)
  # check that values for direct and indirect effects are correct
  expect_equivalent(ci_p_estimates$ci[-5L, ], ci_p_without$ci)
  expect_equivalent(ci_p_estimates$p_value[-5L, ], ci_p_without$p_value)
  expect_equivalent(ci_p_absolute$ci[-5L, ], ci_p_without$ci)
  expect_equivalent(ci_p_absolute$p_value[-5L, ], ci_p_without$p_value)

  ## density plot
  # relevant components are the same as when no contrasts are computed
  keep <- setdiff(names(density_without), c("density", "ci"))
  expect_equal(density_estimates[keep], density_without[keep])
  expect_equal(density_absolute[keep], density_without[keep])
  # check data frame for confidence interval
  expect_s3_class(density_estimates$density, "data.frame")
  expect_s3_class(density_absolute$density, "data.frame")
  # check dimensions
  expect_identical(ncol(density_estimates$density), 3L)
  expect_gt(nrow(density_estimates$density), 0L)
  expect_identical(ncol(density_absolute$density), 3L)
  expect_gt(nrow(density_absolute$density), 0L)
  # check column names
  expect_named(density_estimates$density, c("Effect", "ab", "Density"))
  expect_named(density_absolute$density, c("Effect", "ab", "Density"))
  # check that contrast is included by default
  expect_true("Contrast" %in% density_estimates$density$Effect)
  expect_true("Contrast" %in% density_absolute$density$Effect)
  # check that values for density estimates are correct
  remove_estimates <- which(density_estimates$density$Effect == "Contrast")
  expect_equivalent(density_estimates$density[-remove_estimates, ],
                    density_without$density)
  remove_absolute <- which(density_absolute$density$Effect == "Contrast")
  expect_equivalent(density_absolute$density[-remove_absolute, ],
                    density_without$density)
  # check data frame confidence interval
  expect_s3_class(density_estimates$ci, "data.frame")
  expect_s3_class(density_absolute$ci, "data.frame")
  # check dimensions
  expect_identical(dim(density_estimates$ci), c(4L, 4L))
  expect_identical(dim(density_absolute$ci), c(4L, 4L))
  # check column names
  expect_named(density_estimates$ci, c("Effect", "Estimate", "Lower", "Upper"))
  expect_named(density_absolute$ci, c("Effect", "Estimate", "Lower", "Upper"))
  # check that values for confidence intervals are correct
  expect_equivalent(density_estimates$ci[-4L, ], density_without$ci)
  expect_equivalent(density_absolute$ci[-4L, ], density_without$ci)

  ## ellipse plot
  expect_identical(ellipse_estimates, ellipse_without)
  expect_identical(ellipse_absolute, ellipse_without)

  ## weight plot
  expect_identical(weight_estimates, weight_without)
  expect_identical(weight_absolute, weight_without)

})
