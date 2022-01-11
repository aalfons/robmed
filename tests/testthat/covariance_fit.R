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

## fit mediation models
fit_list <- list(
  winsorized = {
    set.seed(seed)
    fit_mediation(test_data, x = "X1", y = "Y", m = "M1",
                  method = "covariance", robust = TRUE,
                  prob = 0.9)
  },
  standard = {
    set.seed(seed)
    fit_mediation(test_data, x = "X1", y = "Y", m = "M1",
                  method = "covariance", robust = FALSE)
  }
)

## compute summaries
summary_list <- lapply(fit_list, summary)

## other relevant information
methods <- names(fit_list)
classes <- c(winsorized = "cov_Huber", standard = "cov_ML")

## stuff needed to check correctness
coef_names <- c("a", "b", "Total", "Direct", "Indirect")


## loop over methods
for (method in methods) {

  ## extract information for current method
  fit <- fit_list[[method]]
  summary <- summary_list[[method]]


  ## run tests

  test_that("output has correct structure", {

    # covariance fit
    expect_s3_class(fit, "cov_fit_mediation")
    expect_s3_class(fit, "fit_mediation")
    # covariance matrix
    expect_s3_class(fit$cov, classes[method])

  })

  test_that("arguments are correctly passed", {

    # variable names
    expect_identical(fit$x, "X1")
    expect_identical(fit$y, "Y")
    expect_identical(fit$m, "M1")
    expect_identical(fit$covariates, character())
    # robust or nonrobust fit
    if (method == "winsorized") {
      expect_true(fit$robust)
      expect_equal(fit$control, cov_control(prob = 0.9))
    } else {
      expect_false(fit$robust)
      expect_null(fit$control)
    }

  })

  test_that("dimensions are correct", {

    # effects are scalars
    expect_length(fit$a, 1L)
    expect_length(fit$b, 1L)
    expect_length(fit$direct, 1L)
    expect_length(fit$total, 1L)
    expect_length(fit$indirect, 1L)
    # center and covariance matrix
    expect_length(fit$cov$center, 3L)
    expect_identical(dim(fit$cov$cov), c(3L, 3L))
    # dimensions of data
    expect_identical(dim(fit$data), c(as.integer(n), 3L))

  })

  test_that("values of coefficients are correct", {

    expect_equivalent(fit$indirect, fit$a * fit$b)
    expect_equivalent(fit$total, fit$indirect + fit$direct)

  })

  test_that("output of coef() method has correct attributes", {

    coefficients <- coef(fit)
    expect_length(coefficients, 5L)
    expect_named(coefficients, coef_names)

  })

  test_that("coef() method returns correct values of coefficients", {

    expect_equivalent(coef(fit, parm = "a"), fit$a)
    expect_equivalent(coef(fit, parm = "b"), fit$b)
    expect_equivalent(coef(fit, parm = "Total"), fit$total)
    expect_equivalent(coef(fit, parm = "Direct"), fit$direct)
    expect_equivalent(coef(fit, parm = "Indirect"), fit$indirect)

  })

  test_that("summary returns original object", {
    expect_identical(fit, summary)
  })

}


## only implemented for simple mediation without covariates

# loop of methods
for (robust in c(TRUE, FALSE)) {

  test_that("multiple independent variables not implemented", {

    # run regression fit
    set.seed(seed)
    suppressWarnings(
      reg_fit <- fit_mediation(test_data, x = c("X1", "X2"), y = "Y", m = "M1",
                               method = "regression", robust = robust)
    )

    # try to run with multiple independent variables (should give warning)
    set.seed(seed)
    expect_warning(
      cov_fit <- fit_mediation(test_data, x = c("X1", "X2"), y = "Y", m = "M1",
                               method = "covariance", robust = robust)
    )

    # these should be the same
    expect_equal(cov_fit, reg_fit)

  })

  test_that("parallel multiple mediators not implemented", {

    # run regression fit
    set.seed(seed)
    suppressWarnings(
      reg_fit <- fit_mediation(test_data, x = "X1", y = "Y", m = c("M1", "M2"),
                               method = "regression", robust = robust,
                               model = "parallel")
    )

    # try to run with multiple mediators (should give warning)
    set.seed(seed)
    expect_warning(
      cov_fit <- fit_mediation(test_data, x = "X1", y = "Y", m = c("M1", "M2"),
                               method = "covariance", robust = robust,
                               model = "parallel")
    )

    # these should be the same
    expect_equal(cov_fit, reg_fit)

  })

  test_that("serial multiple mediators not implemented", {

    # run regression fit
    set.seed(seed)
    suppressWarnings(
      reg_fit <- fit_mediation(test_data, x = "X1", y = "Y", m = c("M1", "M2"),
                               method = "regression", robust = robust,
                               model = "serial")
    )

    # try to run with multiple mediators (should give warning)
    set.seed(seed)
    expect_warning(
      cov_fit <- fit_mediation(test_data, x = "X1", y = "Y", m = c("M1", "M2"),
                               method = "covariance", robust = robust,
                               model = "serial")
    )

    # these should be the same
    expect_equal(cov_fit, reg_fit)

  })

  test_that("covariates not implemented", {

    # run regression fit
    set.seed(seed)
    suppressWarnings(
      reg_fit <- fit_mediation(test_data, x = "X1", y = "Y", m = "M1",
                               covariates = c("C1", "C2"),
                               method = "regression",
                               robust = robust)
    )

    # try to run with covariates (should give warning)
    set.seed(seed)
    expect_warning(
      cov_fit <- fit_mediation(test_data, x = "X1", y = "Y", m = "M1",
                               covariates = c("C1", "C2"),
                               method = "covariance",
                               robust = robust)
    )

    # these should be the same
    expect_equal(cov_fit, reg_fit)

  })

}
