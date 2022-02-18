context("contrasts of indirect effects: no covariates")


## load package
library("robmed", quietly = TRUE)

## control parameters for data generation
n <- 250          # number of observations
seed <- 20210407  # seed for the random number generator

## generate data
set.seed(seed)
X1 <- rnorm(n)
X2 <- rnorm(n)
C1 <- rnorm(n)
C2 <- rnorm(n)
M1 <- rnorm(n)
M2 <- rnorm(n)
Y <- rnorm(n)
test_data <- data.frame(X1, X2, C1, C2, M1, M2, Y)

## control parameters for bootstrap tests
R <- 100                                   # number of bootstrap samples
level <- 0.95                              # confidence level
type <- "perc"                             # type of confidence interval
ctrl <- reg_control(max_iterations = 500)  # for MM-regression estimator


# ----------------
# simple mediation
# ----------------

## arguments for variables
x <- "X1"
y <- "Y"
m <- "M1"
covariates <- c("C1", "C2")

## perform bootstrap test without formula interface
set.seed(seed)
boot_without <- test_mediation(test_data, x = x, y = y, m = m,
                               covariates = covariates, test = "boot",
                               R = R, level = level, type = type,
                               robust = TRUE, control = ctrl)

## perform bootstrap tests with formula interface
# formula interface with data argument
set.seed(seed)
boot_f1 <- test_mediation(Y ~ m(M1) + X1 + covariates(C1, C2),
                          data = test_data, test = "boot", R = R,
                          level = level, type = type, robust = TRUE,
                          control = ctrl)
# formula interface without data argument
set.seed(seed)
boot_f2 <- test_mediation(Y ~ m(M1) + X1 + covariates(C1, C2), test = "boot",
                          R = R, level = level, type = type, robust = TRUE,
                          control = ctrl)
# define mediator and control variables outside formula
med <- m(M1)
cov <- covariates(C1, C2)
set.seed(seed)
boot_f3 <- test_mediation(Y ~ med + X1 + cov, data = test_data, test = "boot",
                          R = R, level = level, type = type, robust = TRUE,
                          control = ctrl)


## run tests

test_that("formula interface works correctly", {

  ## check that results are the same as with default method
  # ("boot" object should also be the same but needs to be checked separately
  # because of the function that is stored)
  keep <- setdiff(names(boot_without), "reps")
  expect_equal(boot_f1[keep], boot_without[keep])
  expect_equal(boot_f2[keep], boot_without[keep])
  expect_equal(boot_f3[keep], boot_without[keep])
  # check that "boot" object is the same
  keep <- setdiff(names(boot_without$reps), "statistic")
  expect_equal(boot_f1$reps[keep], boot_without$reps[keep])
  expect_equal(boot_f2$reps[keep], boot_without$reps[keep])
  expect_equal(boot_f3$reps[keep], boot_without$reps[keep])

})


# -------------------------
# two independent variables
# -------------------------

## arguments for variables
x <- c("X1", "X2")
y <- "Y"
m <- "M1"
covariates <- c("C1", "C2")

## perform bootstrap test without formula interface
set.seed(seed)
boot_without <- test_mediation(test_data, x = x, y = y, m = m,
                               covariates = covariates, test = "boot",
                               R = R, level = level, type = type,
                               robust = TRUE, control = ctrl)

## perform bootstrap tests with formula interface
# formula interface with data argument
set.seed(seed)
boot_f1 <- test_mediation(Y ~ m(M1) + X1 + X2 + covariates(C1, C2),
                          data = test_data, test = "boot", R = R,
                          level = level, type = type, robust = TRUE,
                          control = ctrl)
# formula interface without data argument
set.seed(seed)
boot_f2 <- test_mediation(Y ~ m(M1) + X1 + X2 + covariates(C1, C2),
                          test = "boot", R = R, level = level, type = type,
                          robust = TRUE, control = ctrl)
# define mediator and control variables outside formula
med <- m(M1)
cov <- covariates(C1, C2)
set.seed(seed)
boot_f3 <- test_mediation(Y ~ med + X1 + X2 + cov, data = test_data,
                          test = "boot", R = R, level = level, type = type,
                          robust = TRUE, control = ctrl)


## run tests

test_that("formula interface works correctly", {

  ## check that results are the same as with default method
  # ("boot" object should also be the same but needs to be checked separately
  # because of the function that is stored)
  keep <- setdiff(names(boot_without), "reps")
  expect_equal(boot_f1[keep], boot_without[keep])
  expect_equal(boot_f2[keep], boot_without[keep])
  expect_equal(boot_f3[keep], boot_without[keep])
  # check that "boot" object is the same
  keep <- setdiff(names(boot_without$reps), "statistic")
  expect_equal(boot_f1$reps[keep], boot_without$reps[keep])
  expect_equal(boot_f2$reps[keep], boot_without$reps[keep])
  expect_equal(boot_f3$reps[keep], boot_without$reps[keep])

})


# ----------------------
# two parallel mediators
# ----------------------

## arguments for variables
x <- "X1"
y <- "Y"
m <- c("M1", "M2")
covariates <- c("C1", "C2")

## perform bootstrap test without formula interface
set.seed(seed)
boot_without <- test_mediation(test_data, x = x, y = y, m = m,
                               covariates = covariates, model = "parallel",
                               test = "boot", R = R, level = level,
                               type = type, robust = TRUE, control = ctrl)

## perform bootstrap tests with formula interface
# formula interface with data argument
set.seed(seed)
boot_f1 <- test_mediation(Y ~ parallel_m(M1, M2) + X1 + covariates(C1, C2),
                          data = test_data, test = "boot", R = R,
                          level = level, type = type, robust = TRUE,
                          control = ctrl)
# formula interface without data argument
set.seed(seed)
boot_f2 <- test_mediation(Y ~ parallel_m(M1, M2) + X1 + covariates(C1, C2),
                          test = "boot", R = R, level = level, type = type,
                          robust = TRUE, control = ctrl)
# define mediators and control variables outside formula
med <- parallel_m(M1, M2)
cov <- covariates(C1, C2)
set.seed(seed)
boot_f3 <- test_mediation(Y ~ med + X1 + cov, data = test_data, test = "boot",
                          R = R, level = level, type = type, robust = TRUE,
                          control = ctrl)


## run tests

test_that("formula interface works correctly", {

  ## check that results are the same as with default method
  # ("boot" object should also be the same but needs to be checked separately
  # because of the function that is stored)
  keep <- setdiff(names(boot_without), "reps")
  expect_equal(boot_f1[keep], boot_without[keep])
  expect_equal(boot_f2[keep], boot_without[keep])
  expect_equal(boot_f3[keep], boot_without[keep])
  # check that "boot" object is the same
  keep <- setdiff(names(boot_without$reps), "statistic")
  expect_equal(boot_f1$reps[keep], boot_without$reps[keep])
  expect_equal(boot_f2$reps[keep], boot_without$reps[keep])
  expect_equal(boot_f3$reps[keep], boot_without$reps[keep])

})


# --------------------
# two serial mediators
# --------------------

## arguments for variables
x <- "X1"
y <- "Y"
m <- c("M1", "M2")
covariates <- c("C1", "C2")

## perform bootstrap test without formula interface
set.seed(seed)
boot_without <- test_mediation(test_data, x = x, y = y, m = m,
                               covariates = covariates, model = "serial",
                               test = "boot", R = R, level = level,
                               type = type, robust = TRUE, control = ctrl)

## perform bootstrap tests with formula interface
# formula interface with data argument
set.seed(seed)
boot_f1 <- test_mediation(Y ~ serial_m(M1, M2) + X1 + covariates(C1, C2),
                          data = test_data, test = "boot", R = R,
                          level = level, type = type, robust = TRUE,
                          control = ctrl)
# formula interface without data argument
set.seed(seed)
boot_f2 <- test_mediation(Y ~ serial_m(M1, M2) + X1 + covariates(C1, C2),
                          test = "boot", R = R, level = level, type = type,
                          robust = TRUE, control = ctrl)
# define mediators and control variables outside formula
med <- serial_m(M1, M2)
cov <- covariates(C1, C2)
set.seed(seed)
boot_f3 <- test_mediation(Y ~ med + X1 + cov, data = test_data, test = "boot",
                          R = R, level = level, type = type, robust = TRUE,
                          control = ctrl)


## run tests

test_that("formula interface works correctly", {

  ## check that results are the same as with default method
  # ("boot" object should also be the same but needs to be checked separately
  # because of the function that is stored)
  keep <- setdiff(names(boot_without), "reps")
  expect_equal(boot_f1[keep], boot_without[keep])
  expect_equal(boot_f2[keep], boot_without[keep])
  expect_equal(boot_f3[keep], boot_without[keep])
  # check that "boot" object is the same
  keep <- setdiff(names(boot_without$reps), "statistic")
  expect_equal(boot_f1$reps[keep], boot_without$reps[keep])
  expect_equal(boot_f2$reps[keep], boot_without$reps[keep])
  expect_equal(boot_f3$reps[keep], boot_without$reps[keep])

})
