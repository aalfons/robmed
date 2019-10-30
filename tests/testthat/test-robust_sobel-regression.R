context("robust Sobel test: regression")


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
ctrl <- reg_control(efficiency = 0.95)
sobel <- test_mediation(test_data, x = "X", y = "Y", m = "M1", test = "sobel",
                        method = "regression", robust = TRUE, median = FALSE,
                        control = ctrl)

## compute summary
summary_sobel <- summary(sobel)

## create data for plotting
dot <- fortify(sobel, method = "dot")
density <- fortify(sobel, method = "density")

## stuff needed to check correctness
coef_names <- c("a", "b", "Direct", "Total", "ab")
mx_names <- c("(Intercept)", "X")
ymx_names <- c("(Intercept)", "M1", "X")


## run tests

test_that("output has correct structure", {

  # bootstrap test
  expect_s3_class(sobel, "sobel_test_mediation")
  expect_s3_class(sobel, "test_mediation")
  # regression fit
  expect_s3_class(sobel$fit, "reg_fit_mediation")
  expect_s3_class(sobel$fit, "fit_mediation")
  # indirect effect
  expect_is(sobel$ab, "numeric")
  # standard error
  expect_is(sobel$se, "numeric")
  # test statistic
  expect_is(sobel$statistic, "numeric")
  # p-value
  expect_is(sobel$p_value, "numeric")

})

test_that("arguments are correctly passed", {

  # alternative hypothesis
  expect_identical(sobel$alternative, "twosided")
  # variable names
  expect_identical(sobel$fit$x, "X")
  expect_identical(sobel$fit$y, "Y")
  expect_identical(sobel$fit$m, "M1")
  expect_identical(sobel$fit$covariates, character())
  # robust fit and test
  expect_true(sobel$fit$robust)
  expect_false(sobel$fit$median)
  expect_equal(sobel$fit$control, ctrl)

})

test_that("dimensions are correct", {

  # indirect effect
  expect_length(sobel$ab, 1L)
  # standard error
  expect_length(sobel$se, 1L)
  # test statistic
  expect_length(sobel$statistic, 1L)
  # p-value
  expect_length(sobel$p_value, 1L)

})

test_that("values of coefficients are correct", {

  expect_equivalent(sobel$ab, sobel$fit$a * sobel$fit$b)

})

test_that("output of coef() method has correct attributes", {

  coef_sobel <- coef(sobel)
  expect_length(coef_sobel, 5L)
  expect_named(coef_sobel, coef_names)

})

test_that("coef() method returns correct values of coefficients", {

  expect_equivalent(coef(sobel, parm = "ab"), sobel$ab)

})

test_that("output of confint() method has correct attributes", {

  ci_sobel <- confint(sobel, level = 0.9)
  expect_equal(dim(ci_sobel), c(5L, 2L))
  expect_equal(rownames(ci_sobel), coef_names)
  expect_equal(colnames(ci_sobel), c("5 %", "95 %"))

})

test_that("confint() method returns correct values of confidence intervals", {

  ci_default <- confint(sobel, parm = "ab")  # should be 95%
  ci_90 <- confint(sobel, parm = "ab", level = 0.9)
  # default CI should be wider
  expect_lt(ci_default["ab", 1], ci_90["ab", 1])
  expect_gt(ci_default["ab", 2], ci_90["ab", 2])

})

test_that("summary has correct structure", {

  # summary
  expect_s3_class(summary_sobel, "summary_test_mediation")
  # original output of test for indirect effect
  expect_identical(summary_sobel$object, sobel)
  # summary of the model fit
  expect_s3_class(summary_sobel$summary, "summary_reg_fit_mediation")
  expect_s3_class(summary_sobel$summary, "summary_fit_mediation")
  # summary for model m ~ x
  expect_s3_class(summary_sobel$summary$fit_mx, "summary_lmrob")
  # summary for model y ~ m + x
  expect_s3_class(summary_sobel$summary$fit_ymx, "summary_lmrob")
  # regression standard error for model y ~ m + x
  expect_type(summary_sobel$summary$fit_ymx$s, "list")
  expect_named(summary_sobel$summary$fit_ymx$s, c("value", "df"))
  # R-squared for model y ~ m + x
  expect_type(summary_sobel$summary$fit_ymx$R2, "list")
  expect_named(summary_sobel$summary$fit_ymx$R2, c("R2", "adj_R2"))
  # F-test for model y ~ m + x
  expect_type(summary_sobel$summary$fit_ymx$F_test, "list")
  expect_named(summary_sobel$summary$fit_ymx$F_test, c("statistic", "df", "p_value"))
  df_test <- summary_sobel$summary$fit_ymx$F_test$df
  expect_identical(df_test[1], 2)
  expect_identical(df_test[2], Inf)

})

test_that("attributes are correctly passed through summary", {

  # robustness
  expect_true(summary_sobel$summary$robust)
  expect_false(summary_sobel$summary$median)
  # number of observations
  expect_identical(summary_sobel$summary$n, as.integer(n))
  # variable names
  expect_identical(summary_sobel$summary$x, "X")
  expect_identical(summary_sobel$summary$y, "Y")
  expect_identical(summary_sobel$summary$m, "M1")
  expect_identical(summary_sobel$summary$covariates, character())

})

test_that("effect summaries have correct names", {

  # a path
  expect_identical(dim(summary_sobel$summary$fit_mx$coefficients), c(2L, 4L))
  expect_identical(rownames(summary_sobel$summary$fit_mx$coefficients), mx_names)
  expect_identical(colnames(summary_sobel$summary$fit_mx$coefficients)[1], "Estimate")
  # b path
  expect_identical(dim(summary_sobel$summary$fit_ymx$coefficients), c(3L, 4L))
  expect_identical(rownames(summary_sobel$summary$fit_ymx$coefficient), ymx_names)
  expect_identical(colnames(summary_sobel$summary$fit_ymx$coefficient)[1], "Estimate")
  # direct effect
  expect_identical(dim(summary_sobel$summary$direct), c(1L, 4L))
  expect_identical(rownames(summary_sobel$summary$direct), "X")
  expect_identical(colnames(summary_sobel$summary$direct)[1], "Estimate")
  # total effect
  expect_identical(dim(summary_sobel$summary$total), c(1L, 4L))
  expect_identical(rownames(summary_sobel$summary$total), "X")
  expect_identical(colnames(summary_sobel$summary$total)[1], "Estimate")

})

test_that("effect summaries contain correct coefficient values", {

  expect_identical(summary_sobel$summary$fit_mx$coefficients["X", "Estimate"], sobel$fit$a)
  expect_identical(summary_sobel$summary$fit_ymx$coefficients["M1", "Estimate"], sobel$fit$b)
  expect_identical(summary_sobel$summary$direct["X", "Estimate"], sobel$fit$direct)
  expect_identical(summary_sobel$summary$total["X", "Estimate"], sobel$fit$total)

})

test_that("data returned by fortify() has correct structure", {

  ## dot plot
  # check dimensions
  expect_s3_class(dot, "data.frame")
  expect_identical(dim(dot), c(2L, 4L))
  # check column names
  column_names <- c("Effect", "Point", "Lower", "Upper")
  expect_named(dot, column_names)
  # check that direct effect and indirect effect are plotted by default
  effect_names <- c("Direct", "ab")
  expect_identical(dot$Effect, factor(effect_names, levels = effect_names))

  ## density plot
  # check dimensions
  expect_s3_class(density, "data.frame")
  expect_identical(ncol(density), 2L)
  expect_gt(nrow(density), 0L)
  # check column names
  column_names <- c("ab", "Density")
  expect_named(density, column_names)

})

test_that("data returned by fortify() has correct attributes", {

  ## dot plot
  # check aesthetic mapping
  mapping <- aes_string(x = "Effect", y = "Point",
                        ymin = "Lower", ymax = "Upper")
  expect_equal(attr(dot, "mapping"), mapping)
  # check default geom()
  expect_identical(attr(dot, "geom"), geom_pointrange)
  # check facets
  expect_null(attr(dot, "facets"))
  # check that method is stored correctly
  expect_identical(attr(dot, "method"), "dot")

  ## density plot
  # check aesthetic mapping
  mapping <- aes_string(x = "ab", y = "Density")
  expect_equal(attr(density, "mapping"), mapping)
  # check default geom()
  expect_equal(attr(density, "geom"), function(..., stat) {
    geom_density(..., stat = "identity")
  })
  # check facets
  expect_null(attr(density, "facets"))
  # check title
  expect_identical(attr(density, "main"), "Assumed normal distribution")
  # check confidence interval
  ci <- attr(density, "ci")
  expect_s3_class(ci, "data.frame")
  expect_identical(dim(ci), c(1L, 4L))
  expect_named(ci, c("ab", "Density", "Lower", "Upper"))
  # check that method is stored correctly
  expect_identical(attr(density, "method"), "density")

})
