context("covariance fit")


## load package
library("robmed", quietly = TRUE)

## control parameters
n <- 250            # number of observations
a <- c <- 0.2       # true effects
b <- 0              # true effect
seed <- 20150601    # seed for the random number generator

## set seed for reproducibility
set.seed(seed)

## generate data
X1 <- rnorm(n)
X2 <- rnorm(n)
M1 <- a * X1 + rnorm(n)
M2 <- rnorm(n)
Y <- b * M1 + c * X1 + rnorm(n)
C1 <- rnorm(n)
C2 <- rnorm(n)
test_data <- data.frame(X1, X2, Y, M1, M2, C1, C2)

## arguments for bootstrap test
R <- 100                         # number of bootstrap samples
level <- 0.9                     # confidence level
type <- "perc"                   # type of confidence intervals
ctrl <- cov_control(prob = 0.9)  # control parameters for MM-regression

## fit mediation models
boot_list <- list(
  winsorized = {
    set.seed(seed)
    test_mediation(test_data, x = "X1", y = "Y", m = "M1",
                   test = "boot", R = R, level = level, type = type,
                   method = "covariance", robust = TRUE, control = ctrl)
  },
  standard = {
    set.seed(seed)
    test_mediation(test_data, x = "X1", y = "Y", m = "M1",
                   test = "boot", R = R, level = level, type = type,
                   method = "covariance", robust = FALSE)
  }
)

## compute summaries
summary_boot_list <- lapply(boot_list, summary, type = "boot")
summary_data_list <- lapply(boot_list, summary, type = "data")

## stuff needed to check correctness
coef_names <- c("a", "b", "Total", "Direct", "Indirect")


## loop over methods
methods <- names(boot_list)
for (method in methods) {

  ## extract information for current method
  boot <- boot_list[[method]]
  summary_boot <- summary_boot_list[[method]]
  summary_data <- summary_data_list[[method]]


  ## run tests

  test_that("output has correct structure", {

    # bootstrap test
    expect_s3_class(boot, "boot_test_mediation")
    expect_s3_class(boot, "test_mediation")
    # regression fit
    expect_s3_class(boot$fit, "cov_fit_mediation")
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
    expect_identical(boot$type, type)
    # variable names
    expect_identical(boot$fit$x, "X1")
    expect_identical(boot$fit$y, "Y")
    expect_identical(boot$fit$m, "M1")
    expect_identical(boot$fit$covariates, character())
    # robust or nonrobust fit and test
    if (method == "winsorized") {
      expect_true(boot$fit$robust)
      expect_equal(boot$fit$control, ctrl)
    } else {
      expect_false(boot$fit$robust)
      expect_null(boot$fit$control)
    }

  })

  test_that("dimensions are correct", {

    # effects are scalars
    expect_length(boot$indirect, 1L)
    # only one indirect effect, so only one confidence interval
    expect_length(boot$ci, 2L)
    # dimensions of bootstrap replicates
    expect_identical(dim(boot$reps$t), c(as.integer(R), 7L))

  })

  test_that("values of coefficients are correct", {

    expect_equivalent(boot$indirect, mean(boot$reps$t[, 1]))

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
                      mean(boot$reps$t[, 5]))
    expect_equivalent(coef(boot, parm = "Total", type = "boot"),
                      mean(boot$reps$t[, 7]))
    expect_equivalent(coef(boot, parm = "Direct", type = "boot"),
                      mean(boot$reps$t[, 6]))
    expect_equivalent(coef(boot, parm = "Indirect", type = "boot"),
                      boot$indirect)

    # effects computed on original sample
    expect_equivalent(coef(boot, parm = "a", type = "data"),
                      boot$fit$a)
    expect_equivalent(coef(boot, parm = "b", type = "data"),
                      boot$fit$b)
    expect_equivalent(coef(boot, parm = "Total", type = "data"),
                      boot$fit$total)
    expect_equivalent(coef(boot, parm = "Direct", type = "data"),
                      boot$fit$direct)
    expect_equivalent(coef(boot, parm = "Indirect", type = "data"),
                      boot$fit$indirect)

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
    expect_equivalent(confint(boot, parm = "Indirect", type = "boot"), boot$ci)
    # confidence intervals based on theory (except for indirect effect)
    expect_equivalent(confint(boot, parm = "Indirect", type = "data"), boot$ci)

  })

  test_that("summary has correct structure", {

    # summary
    expect_s3_class(summary_boot, "summary_test_mediation")
    expect_s3_class(summary_data, "summary_test_mediation")
    # original output of test for indirect effect
    expect_identical(summary_boot$object, boot)
    expect_identical(summary_data$object, boot)
    # summary of the model fit
    expect_s3_class(summary_boot$summary, "summary_cov_fit_mediation")
    expect_s3_class(summary_boot$summary, "summary_fit_mediation")
    expect_s3_class(summary_data$summary, "summary_cov_fit_mediation")
    expect_s3_class(summary_data$summary, "summary_fit_mediation")
    # regression standard error for model y ~ m + x
    expect_null(summary_boot$summary$s)
    expect_null(summary_data$summary$s)
    # R-squared for model y ~ m + x
    expect_null(summary_boot$summary$R2)
    expect_null(summary_data$summary$R2)
    # F-test for model y ~ m + x
    expect_null(summary_boot$summary$F_test)
    expect_null(summary_data$summary$F_test)
    # no plot is created
    expect_null(summary_boot$plot)
    expect_null(summary_data$plot)

  })

  test_that("attributes are correctly passed through summary", {

    # robustness
    if (method == "winsorized") {
      expect_true(summary_boot$summary$robust)
      expect_true(summary_data$summary$robust)
    } else {
      expect_false(summary_boot$summary$robust)
      expect_false(summary_data$summary$robust)
    }
    # number of observations
    expect_identical(summary_boot$summary$n, as.integer(n))
    expect_identical(summary_data$summary$n, as.integer(n))
    # variable names
    expect_identical(summary_boot$summary$x, "X1")
    expect_identical(summary_boot$summary$y, "Y")
    expect_identical(summary_boot$summary$m, "M1")
    expect_null(summary_boot$summary$covariates)
    expect_identical(summary_data$summary$x, "X1")
    expect_identical(summary_data$summary$y, "Y")
    expect_identical(summary_data$summary$m, "M1")
    expect_null(summary_data$summary$covariates)

  })

  test_that("effect summaries have correct names", {

    # a path
    expect_identical(dim(summary_boot$summary$a), c(1L, 5L))
    expect_identical(rownames(summary_boot$summary$a), "X1")
    expect_identical(colnames(summary_boot$summary$a)[1:2], c("Data", "Boot"))
    expect_identical(dim(summary_data$summary$a), c(1L, 4L))
    expect_identical(rownames(summary_data$summary$a), "X1")
    expect_identical(colnames(summary_data$summary$a)[1], "Estimate")
    # b path
    expect_identical(dim(summary_boot$summary$b), c(1L, 5L))
    expect_identical(rownames(summary_boot$summary$b), "M1")
    expect_identical(colnames(summary_boot$summary$b)[1:2], c("Data", "Boot"))
    expect_identical(dim(summary_data$summary$b), c(1L, 4L))
    expect_identical(rownames(summary_data$summary$b), "M1")
    expect_identical(colnames(summary_data$summary$b)[1], "Estimate")
    # total effect
    expect_identical(dim(summary_boot$summary$total),
                     c(1L, 5L))
    expect_identical(rownames(summary_boot$summary$total),
                     "X1")
    expect_identical(colnames(summary_boot$summary$total)[1:2],
                     c("Data", "Boot"))
    expect_identical(dim(summary_data$summary$total),
                     c(1L, 4L))
    expect_identical(rownames(summary_data$summary$total),
                     "X1")
    expect_identical(colnames(summary_data$summary$total)[1],
                     "Estimate")
    # direct effect
    expect_identical(dim(summary_boot$summary$direct),
                     c(1L, 5L))
    expect_identical(rownames(summary_boot$summary$direct),
                     "X1")
    expect_identical(colnames(summary_boot$summary$direct)[1:2],
                     c("Data", "Boot"))
    expect_identical(dim(summary_data$summary$direct),
                     c(1L, 4L))
    expect_identical(rownames(summary_data$summary$direct),
                     "X1")
    expect_identical(colnames(summary_data$summary$direct)[1],
                     "Estimate")
    # no model fits
    expect_null(summary_boot$summary$fit_mx)
    expect_null(summary_boot$summary$fit_ymx)
    expect_null(summary_data$summary$fit_mx)
    expect_null(summary_data$summary$fit_ymx)

  })

  test_that("effect summaries contain correct coefficient values", {

    # effects computed on original sample
    expect_identical(summary_boot$summary$a["X1", "Data"],
                     boot$fit$a)
    expect_identical(summary_boot$summary$b["M1", "Data"],
                     boot$fit$b)
    expect_identical(summary_boot$summary$total["X1", "Data"],
                     boot$fit$total)
    expect_identical(summary_boot$summary$direct["X1", "Data"],
                     boot$fit$direct)
    expect_identical(summary_data$summary$a["X1", "Estimate"],
                     boot$fit$a)
    expect_identical(summary_data$summary$b["M1", "Estimate"],
                     boot$fit$b)
    expect_identical(summary_data$summary$total["X1", "Estimate"],
                     boot$fit$total)
    expect_identical(summary_data$summary$direct["X1", "Estimate"],
                     boot$fit$direct)

    # bootstrapped effects
    expect_equal(summary_boot$summary$a["X1", "Boot"],
                 mean(boot$reps$t[, 3]))
    expect_equal(summary_boot$summary$b["M1", "Boot"],
                 mean(boot$reps$t[, 5]))
    expect_equal(summary_boot$summary$total["X1", "Boot"],
                 mean(boot$reps$t[, 7]))
    expect_equal(summary_boot$summary$direct["X1", "Boot"],
                 mean(boot$reps$t[, 6]))

  })

}
