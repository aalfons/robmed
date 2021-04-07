context("contrasts: single mediator, no covariates")


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
M <- a * X + rnorm(n)
Y <- b * M + c * X + rnorm(n)
test_data <- data.frame(X, Y, M)

## control parameters
level <- 0.95
ctrl <- reg_control(max_iterations = 500)

## run bootstrap tests
# without contrasts
set.seed(seed)
boot_without <- test_mediation(test_data, x = "X", y = "Y", m = "M",
                               test = "boot", R = R, level = level,
                               type = "perc", method = "regression",
                               robust = TRUE, contrast = FALSE,
                               control = ctrl)
# test differences of indirect effect
set.seed(seed)
boot_estimates <- test_mediation(test_data, x = "X", y = "Y", m = "M",
                                 test = "boot", R = R, level = level,
                                 type = "perc", method = "regression",
                                 robust = TRUE, contrast = "estimates",
                                 control = ctrl)
# test differences of absolute values of indirect effect
set.seed(seed)
boot_absolute <- test_mediation(test_data, x = "X", y = "Y", m = "M",
                                test = "boot", R = R, level = level,
                                type = "perc", method = "regression",
                                robust = TRUE, contrast = "absolute",
                                control = ctrl)


## run tests

test_that("single mediator yields no contrasts", {

  # there should be no contrasts, so objects should be the same
  expect_equal(boot_estimates, boot_without)
  expect_equal(boot_absolute, boot_without)

  # retest() should produce warning
  expect_warning(retest(boot_without, contrast = "estimates"))
  expect_warning(retest(boot_without, contrast = "absolute"))
})
