context("contrasts of indirect effects: covariates")


## load package
library("robmed", quietly = TRUE)

## control parameters for data generation
n <- 250               # number of observations
a <- c <- c(0.2, 0.2)  # true effects
d <- 0.2               # true effects
b <- c(0, 0.2)         # true effect
seed <- 20210407       # seed for the random number generator

## control parameters for bootstrap tests
R <- 100                                   # number of bootstrap samples
level <- 0.95                              # confidence level
ctrl <- reg_control(max_iterations = 500)  # for MM-regression estimator


# ----------------
# simple mediation
# ----------------

## set seed for reproducibility
set.seed(seed)

## generate data
X <- rnorm(n)
M <- a[1] * X + rnorm(n)
Y <- b[1] * M + c[1] * X + rnorm(n)
C1 <- rnorm(n)
C2 <- rnorm(n)
test_data <- data.frame(X, Y, M, C1, C2)

## variables for bootstrap tests
x <- "X"                     # independent variable
y <- "Y"                     # dependent variable
m <- "M"                     # mediator variable
covariates <- c("C1", "C2")  # control variables

## perform bootstrap tests
# without contrasts
set.seed(seed)
boot_without <- test_mediation(test_data, x = x, y = y, m = m,
                               covariates = covariates, test = "boot",
                               R = R, level = level, type = "perc",
                               method = "regression", robust = TRUE,
                               contrast = FALSE, control = ctrl)
# test differences of indirect effect
set.seed(seed)
boot_estimates <- test_mediation(test_data, x = x, y = y, m = m,
                                 covariates = covariates, test = "boot",
                                 R = R, level = level, type = "perc",
                                 method = "regression", robust = TRUE,
                                 contrast = "estimates", control = ctrl)
# test differences of absolute values of indirect effect
set.seed(seed)
boot_absolute <- test_mediation(test_data, x = x, y = y, m = m,
                                covariates = covariates, test = "boot",
                                R = R, level = level, type = "perc",
                                method = "regression", robust = TRUE,
                                contrast = "absolute", control = ctrl)


## run tests

test_that("simple mediation yields no contrasts", {

  ## there should be no contrasts, so objects should be the same
  # ("boot" object should also be the same but needs to be checked separately
  # because of the function that is stored)
  keep <- setdiff(names(boot_without), "reps")
  expect_equal(boot_estimates[keep], boot_without[keep])
  expect_equal(boot_absolute[keep], boot_without[keep])
  # check that "boot" object is the same
  keep <- setdiff(names(boot_without$reps), "statistic")
  expect_equal(boot_estimates$reps[keep], boot_without$reps[keep])
  expect_equal(boot_absolute$reps[keep], boot_without$reps[keep])

  ## retest() should produce warning
  expect_warning(retest(boot_without, contrast = "estimates"))
  expect_warning(retest(boot_without, contrast = "absolute"))
})


# -------------------------
# two independent variables
# -------------------------

## set seed for reproducibility
set.seed(seed)

## generate data
X1 <- rnorm(n)
X2 <- rnorm(n)
M <- a[1] * X1 + a[2] * X2 + rnorm(n)
Y <- b[1] * M + c[1] * X1 + c[2] * X2 + rnorm(n)
C1 <- rnorm(n)
C2 <- rnorm(n)
test_data <- data.frame(X1, X2, Y, M, C1, C2)

## variables for bootstrap tests
x <- c("X1", "X2")           # independent variables
y <- "Y"                     # dependent variable
m <- "M"                     # mediator variable
covariates <- c("C1", "C2")  # control variables

## perform bootstrap tests
# without contrasts
set.seed(seed)
boot_without <- test_mediation(test_data, x = x, y = y, m = m,
                               covariates = covariates, test = "boot",
                               R = R, level = level, type = "perc",
                               method = "regression", robust = TRUE,
                               contrast = FALSE, control = ctrl)
# with contrasts
boot_list <- list(
  estimates = {
    # test differences of indirect effect
    set.seed(seed)
    test_mediation(test_data, x = x, y = y, m = m,
                   covariates = covariates, test = "boot",
                   R = R, level = level, type = "perc",
                   method = "regression", robust = TRUE,
                   contrast = "estimates", control = ctrl)
  },
  absolute = {
    # test differences of absolute values of indirect effect
    set.seed(seed)
    test_mediation(test_data, x = x, y = y, m = m,
                   covariates = covariates, test = "boot",
                   R = R, level = level, type = "perc",
                   method = "regression", robust = TRUE,
                   contrast = "absolute", control = ctrl)
  }
)

## stuff needed to check correctness
indirect_names <- c(x, "Contrast")


## run tests

# loop over types of contrasts
contrasts <- names(boot_list)
for (contrast in contrasts) {

  # extract bootstrap test with current contrasts
  boot <- boot_list[[contrast]]


  # run tests

  test_that("output has correct structure for two independent variables", {

    # bootstrap test
    expect_s3_class(boot, "boot_test_mediation")
    expect_s3_class(boot, "test_mediation")

    # regression fit
    expect_s3_class(boot$fit, "reg_fit_mediation")
    expect_s3_class(boot$fit, "fit_mediation")

  })

  test_that("model fit works correctly for two independent variables", {

    # extract model fits
    fit_without <- boot_without$fit
    fit <- boot$fit

    # relevant components are the same as when no contrasts are computed
    keep <- setdiff(names(fit_without), c("indirect", "ab", "contrast"))
    expect_equal(fit[keep], fit_without[keep])

    # multiple indirect effects have correct dimensions and names
    expect_length(fit$indirect, 3L)
    expect_identical(names(fit$indirect), indirect_names)

    # contrasts computed correctly
    if (contrast == "absolute") {
      expect_equivalent(fit$indirect["Contrast"],
                        abs(fit$indirect[x[1]]) - abs(fit$indirect[x[2]]))
    } else {
      expect_equivalent(fit$indirect["Contrast"],
                        fit$indirect[x[1]] - fit$indirect[x[2]])
    }

    # type of contrasts stored correctly
    expect_identical(fit$contrast, contrast)

  })

  test_that("bootstrap test works correctly for two independent variables", {

    # relevant components are the same as when no contrasts are computed
    # ("boot" object should also be the same but needs to be checked separately
    # because of the function that is stored)
    keep <- setdiff(names(boot), c("indirect", "ab", "ci", "fit", "reps"))
    expect_equal(boot[keep], boot_without[keep])

    # check that "boot" object is the same as when no contrasts are computed
    keep <- setdiff(names(boot$reps), c("statistic"))
    expect_equal(boot$reps[keep], boot_without$reps[keep])

    # multiple indirect effects have correct dimensions and names
    expect_length(boot$indirect, 3L)
    expect_named(boot$indirect, indirect_names)

    # confidence intervals have correct dimensions and names
    expect_identical(dim(boot$ci), c(3L, 2L))
    expect_identical(rownames(boot$ci), indirect_names)
    expect_identical(colnames(boot$ci), colnames(boot_without$ci))

  })

  test_that("retest() works correctly for two independent variables", {

    # compute contrasts with retest()
    reboot <- retest(boot_without, contrast = contrast)
    # object should be the same as when argument is supplied from the start
    # ("boot" object should also be the same but needs to be checked separately
    # because of the function that is stored)
    keep <- setdiff(names(boot), "reps")
    expect_equal(reboot[keep], boot[keep])

    # check that "boot" object is the same
    keep <- setdiff(names(boot$reps), "statistic")
    expect_equal(reboot$reps[keep], boot$reps[keep])

  })

}


# ----------------------
# two parallel mediators
# ----------------------

## set seed for reproducibility
set.seed(seed)

## generate data
X <- rnorm(n)
M1 <- a[1] * X + rnorm(n)
M2 <- a[2] * X + rnorm(n)
Y <- b[1] * M1 + b[2] * M1 + c[1] * X + rnorm(n)
C1 <- rnorm(n)
C2 <- rnorm(n)
test_data <- data.frame(X, Y, M1, M2, C1, C2)

## variables for bootstrap tests
x <- "X"                     # independent variables
y <- "Y"                     # dependent variable
m <- c("M1", "M2")           # mediator variable
covariates <- c("C1", "C2")  # control variables

## perform bootstrap tests
# without contrasts
set.seed(seed)
boot_without <- test_mediation(test_data, x = x, y = y, m = m,
                               covariates = covariates, model = "parallel",
                               test = "boot", R = R, level = level,
                               type = "perc", method = "regression",
                               robust = TRUE, contrast = FALSE, control = ctrl)
# with contrasts
boot_list <- list(
  estimates = {
    # test differences of indirect effect
    set.seed(seed)
    test_mediation(test_data, x = x, y = y, m = m, covariates = covariates,
                   model = "parallel", test = "boot", R = R, level = level,
                   type = "perc", method = "regression", robust = TRUE,
                   contrast = "estimates", control = ctrl)
  },
  absolute = {
    # test differences of absolute values of indirect effect
    set.seed(seed)
    test_mediation(test_data, x = x, y = y, m = m, covariates = covariates,
                   model = "parallel", test = "boot", R = R, level = level,
                   type = "perc", method = "regression", robust = TRUE,
                   contrast = "absolute", control = ctrl)
  }
)

## stuff needed to check correctness
indirect_names <- c("Total", m, "Contrast")


## run tests

# loop over types of contrasts
contrasts <- names(boot_list)
for (contrast in contrasts) {

  # extract bootstrap test with current contrasts
  boot <- boot_list[[contrast]]


  # run tests

  test_that("output has correct structure for two parallel mediators", {

    # bootstrap test
    expect_s3_class(boot, "boot_test_mediation")
    expect_s3_class(boot, "test_mediation")

    # regression fit
    expect_s3_class(boot$fit, "reg_fit_mediation")
    expect_s3_class(boot$fit, "fit_mediation")

  })

  test_that("model fit works correctly for two parallel mediators", {

    # extract model fits
    fit_without <- boot_without$fit
    fit <- boot$fit

    # relevant components are the same as when no contrasts are computed
    keep <- setdiff(names(fit_without), c("indirect", "ab", "contrast"))
    expect_equal(fit[keep], fit_without[keep])

    # multiple indirect effects have correct dimensions and names
    expect_length(fit$indirect, 4L)
    expect_identical(names(fit$indirect), indirect_names)

    # contrasts computed correctly
    if (contrast == "absolute") {
      expect_equivalent(fit$indirect["Contrast"],
                        abs(fit$indirect[m[1]]) - abs(fit$indirect[m[2]]))
    } else {
      expect_equivalent(fit$indirect["Contrast"],
                        fit$indirect[m[1]] - fit$indirect[m[2]])
    }

    # type of contrasts stored correctly
    expect_identical(fit$contrast, contrast)

  })

  test_that("bootstrap test works correctly for two parallel mediators", {

    # relevant components are the same as when no contrasts are computed
    # ("boot" object should also be the same but needs to be checked separately
    # because of the function that is stored)
    keep <- setdiff(names(boot), c("indirect", "ab", "ci", "fit", "reps"))
    expect_equal(boot[keep], boot_without[keep])

    # check that "boot" object is the same as when no contrasts are computed
    keep <- setdiff(names(boot$reps), c("statistic"))
    expect_equal(boot$reps[keep], boot_without$reps[keep])

    # multiple indirect effects have correct dimensions and names
    expect_length(boot$indirect, 4L)
    expect_named(boot$indirect, indirect_names)

    # confidence intervals have correct dimensions and names
    expect_identical(dim(boot$ci), c(4L, 2L))
    expect_identical(rownames(boot$ci), indirect_names)
    expect_identical(colnames(boot$ci), colnames(boot_without$ci))

  })

  test_that("retest() works correctly for two parallel mediators", {

    # compute contrasts with retest()
    reboot <- retest(boot_without, contrast = contrast)
    # object should be the same as when argument is supplied from the start
    # ("boot" object should also be the same but needs to be checked separately
    # because of the function that is stored)
    keep <- setdiff(names(boot), "reps")
    expect_equal(reboot[keep], boot[keep])

    # check that "boot" object is the same
    keep <- setdiff(names(boot$reps), "statistic")
    expect_equal(reboot$reps[keep], boot$reps[keep])

  })

}


# --------------------
# two serial mediators
# --------------------

## set seed for reproducibility
set.seed(seed)

## generate data
X <- rnorm(n)
M1 <- a[1] * X + rnorm(n)
M2 <- d * M1 + a[2] * X + rnorm(n)
Y <- b[1] * M1 + b[2] * M2 + c[1] * X + rnorm(n)
C1 <- rnorm(n)
C2 <- rnorm(n)
test_data <- data.frame(X, Y, M1, M2, C1, C2)

## variables for bootstrap tests
x <- "X"                     # independent variables
y <- "Y"                     # dependent variable
m <- c("M1", "M2")           # mediator variable
covariates <- c("C1", "C2")  # control variables

## perform bootstrap tests
# without contrasts
set.seed(seed)
boot_without <- test_mediation(test_data, x = x, y = y, m = m,
                               covariates = covariates, model = "serial",
                               test = "boot", R = R, level = level,
                               type = "perc", method = "regression",
                               robust = TRUE, contrast = FALSE, control = ctrl)
# with contrasts
boot_list <- list(
  estimates = {
    # test differences of indirect effect
    set.seed(seed)
    test_mediation(test_data, x = x, y = y, m = m, covariates = covariates,
                   model = "serial", test = "boot", R = R, level = level,
                   type = "perc", method = "regression", robust = TRUE,
                   contrast = "estimates", control = ctrl)
  },
  absolute = {
    # test differences of absolute values of indirect effect
    set.seed(seed)
    test_mediation(test_data, x = x, y = y, m = m, covariates = covariates,
                   model = "serial", test = "boot", R = R, level = level,
                   type = "perc", method = "regression", robust = TRUE,
                   contrast = "absolute", control = ctrl)
  }
)

## stuff needed to check correctness
collapsed_m <- paste(m, collapse = "->")
indirect_names <- c("Total", m, collapsed_m, paste0("Contrast", 1:3))


## run tests

# loop over types of contrasts
contrasts <- names(boot_list)
for (contrast in contrasts) {

  # extract bootstrap test with current contrasts
  boot <- boot_list[[contrast]]


  # run tests

  test_that("output has correct structure for two parallel mediators", {

    # bootstrap test
    expect_s3_class(boot, "boot_test_mediation")
    expect_s3_class(boot, "test_mediation")

    # regression fit
    expect_s3_class(boot$fit, "reg_fit_mediation")
    expect_s3_class(boot$fit, "fit_mediation")

  })

  test_that("model fit works correctly for two parallel mediators", {

    # extract model fits
    fit_without <- boot_without$fit
    fit <- boot$fit

    # relevant components are the same as when no contrasts are computed
    keep <- setdiff(names(fit_without), c("indirect", "ab", "contrast"))
    expect_equal(fit[keep], fit_without[keep])

    # multiple indirect effects have correct dimensions and names
    expect_length(fit$indirect, 7L)
    expect_identical(names(fit$indirect), indirect_names)

    # contrasts computed correctly
    if (contrast == "absolute") {
      expect_equivalent(fit$indirect["Contrast1"],
                        abs(fit$indirect[m[1]]) - abs(fit$indirect[m[2]]))
      expect_equivalent(fit$indirect["Contrast2"],
                        abs(fit$indirect[m[1]]) - abs(fit$indirect[collapsed_m]))
      expect_equivalent(fit$indirect["Contrast3"],
                        abs(fit$indirect[m[2]]) - abs(fit$indirect[collapsed_m]))
    } else {
      expect_equivalent(fit$indirect["Contrast1"],
                        fit$indirect[m[1]] - fit$indirect[m[2]])
      expect_equivalent(fit$indirect["Contrast2"],
                        fit$indirect[m[1]] - fit$indirect[collapsed_m])
      expect_equivalent(fit$indirect["Contrast3"],
                        fit$indirect[m[2]] - fit$indirect[collapsed_m])
    }

    # type of contrasts stored correctly
    expect_identical(fit$contrast, contrast)

  })

  test_that("bootstrap test works correctly for two parallel mediators", {

    # relevant components are the same as when no contrasts are computed
    # ("boot" object should also be the same but needs to be checked separately
    # because of the function that is stored)
    keep <- setdiff(names(boot), c("indirect", "ab", "ci", "fit", "reps"))
    expect_equal(boot[keep], boot_without[keep])

    # check that "boot" object is the same as when no contrasts are computed
    keep <- setdiff(names(boot$reps), c("statistic"))
    expect_equal(boot$reps[keep], boot_without$reps[keep])

    # multiple indirect effects have correct dimensions and names
    expect_length(boot$indirect, 7L)
    expect_named(boot$indirect, indirect_names)

    # confidence intervals have correct dimensions and names
    expect_identical(dim(boot$ci), c(7L, 2L))
    expect_identical(rownames(boot$ci), indirect_names)
    expect_identical(colnames(boot$ci), colnames(boot_without$ci))

  })

  test_that("retest() works correctly for two parallel mediators", {

    # compute contrasts with retest()
    reboot <- retest(boot_without, contrast = contrast)
    # object should be the same as when argument is supplied from the start
    # ("boot" object should also be the same but needs to be checked separately
    # because of the function that is stored)
    keep <- setdiff(names(boot), "reps")
    expect_equal(reboot[keep], boot[keep])

    # check that "boot" object is the same
    keep <- setdiff(names(boot$reps), "statistic")
    expect_equal(reboot$reps[keep], boot$reps[keep])

  })

}
