context("contrasts in mediation models with multiple mediators")


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
M3 <- rnorm(n)
Y <- b * M1 + c * X + rnorm(n)
C1 <- rnorm(n)
C2 <- rnorm(n)
test_data <- data.frame(X, Y, M1, M2, M3, C1, C2)

## control parameters
level <- 0.95
ctrl <- reg_control(max_iterations = 500)

## bootstrap test with a single mediator
# without contrasts
set.seed(seed)
boot1_without <- test_mediation(test_data, x = "X", y = "Y", m = "M1",
                                covariates = c("C1", "C2"), test = "boot",
                                R = R, level = level, type = "perc",
                                method = "regression", robust = TRUE,
                                contrast = FALSE, control = ctrl)
# test differences of indirect effect
set.seed(seed)
boot1_estimates <- test_mediation(test_data, x = "X", y = "Y", m = "M1",
                                  covariates = c("C1", "C2"), test = "boot",
                                  R = R, level = level, type = "perc",
                                  method = "regression", robust = TRUE,
                                  contrast = "estimates", control = ctrl)
# test differences of absolute values of indirect effect
set.seed(seed)
boot1_absolute <- test_mediation(test_data, x = "X", y = "Y", m = "M1",
                                 covariates = c("C1", "C2"), test = "boot",
                                 R = R, level = level, type = "perc",
                                 method = "regression", robust = TRUE,
                                 contrast = "absolute", control = ctrl)

## bootstrap test with two mediators
# without contrasts
set.seed(seed)
boot2_without <- test_mediation(test_data, x = "X", y = "Y", m = c("M1", "M2"),
                                covariates = c("C1", "C2"), test = "boot",
                                R = R, level = level, type = "perc",
                                method = "regression", robust = TRUE,
                                contrast = FALSE, control = ctrl)
# re-test with contrasts
reboot2_estimates <- retest(boot2_without, contrast = "estimates")
reboot2_absolute <- retest(boot2_without, contrast = "absolute")
# test differences of indirect effect
set.seed(seed)
boot2_estimates <- test_mediation(test_data, x = "X", y = "Y", m = c("M1", "M2"),
                                  covariates = c("C1", "C2"), test = "boot",
                                  R = R, level = level, type = "perc",
                                  method = "regression", robust = TRUE,
                                  contrast = "estimates", control = ctrl)
# test differences of absolute values of indirect effect
set.seed(seed)
boot2_absolute <- test_mediation(test_data, x = "X", y = "Y", m = c("M1", "M2"),
                                 covariates = c("C1", "C2"), test = "boot",
                                 R = R, level = level, type = "perc",
                                 method = "regression", robust = TRUE,
                                 contrast = "absolute", control = ctrl)

## bootstrap test with three mediators
# without contrasts
set.seed(seed)
boot3_without <- test_mediation(test_data, x = "X", y = "Y",
                                m = c("M1", "M2", "M3"),
                                covariates = c("C1", "C2"),
                                test = "boot", R = R, level = level,
                                type = "perc", method = "regression",
                                robust = TRUE, contrast = FALSE,
                                control = ctrl)
# re-test with contrasts
reboot3_estimates <- retest(boot3_without, contrast = "estimates")
reboot3_absolute <- retest(boot3_without, contrast = "absolute")
# test differences of indirect effect
set.seed(seed)
boot3_estimates <- test_mediation(test_data, x = "X", y = "Y",
                                  m = c("M1", "M2", "M3"),
                                  covariates = c("C1", "C2"),
                                  test = "boot", R = R, level = level,
                                  type = "perc", method = "regression",
                                  robust = TRUE, contrast = "estimates",
                                  control = ctrl)
# test differences of absolute values of indirect effect
set.seed(seed)
boot3_absolute <- test_mediation(test_data, x = "X", y = "Y",
                                 m = c("M1", "M2", "M3"),
                                 covariates = c("C1", "C2"),
                                 test = "boot", R = R, level = level,
                                 type = "perc", method = "regression",
                                 robust = TRUE, contrast = "absolute",
                                 control = ctrl)


## unit tests

test_that("single mediator yields no contrasts", {

  # there should be no contrasts, so objects should be the same
  expect_equal(boot1_estimates, boot1_without)
  expect_equal(boot1_absolute, boot1_without)

  # retest() should produce warning
  expect_warning(retest(boot1_without, contrast = "estimates"))
  expect_warning(retest(boot1_without, contrast = "absolute"))
})


test_that("contrasts for two mediators have correct structure", {

  # "boot" object is unchanged
  expect_equal(boot2_estimates$reps, boot2_without$reps)
  expect_equal(boot2_absolute$reps, boot2_without$reps)

})

test_that("retest() works correctly for two mediators", {
  # object should be the same as when argument is supplied from the start
  expect_equal(reboot2_estimates, boot2_estimates)
  expect_equal(reboot2_absolute, boot2_absolute)
})



test_that("retest() works correctly for three mediators", {
  # object should be the same as when argument is supplied from the start
  expect_equal(reboot3_estimates, boot3_estimates)
  expect_equal(reboot3_absolute, boot3_absolute)
})
