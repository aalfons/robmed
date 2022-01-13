context("bootstrap test: simple mediation, no covariates")


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
X <- rnorm(n)
M <- a * X + rnorm(n)
Y <- b * M + c * X + rnorm(n)
test_data <- data.frame(X, Y, M)

## arguments for bootstrap tests
R <- 100                                # number of bootstrap samples
level <- 0.9                            # confidence level
type <- "perc"                          # type of confidence intervals
reg_ctrl <- reg_control(efficiency = 0.95)  # for MM-regression estimator
cov_ctrl <- cov_control(prob = 0.9)         # for winsorized covariance matrix

## perform bootstrap tests
boot_list <- list(
  robust = {
    set.seed(seed)
    test_mediation(test_data, x = "X", y = "Y", m = "M",
                   test = "boot", R = R, level = level, type = type,
                   method = "regression", robust = TRUE, control = reg_ctrl)
  },
  median = {
    set.seed(seed)
    test_mediation(test_data, x = "X", y = "Y", m = "M",
                   test = "boot", R = R, level = level, type = type,
                   method = "regression", robust = "median")
  },
  OLS = {
    set.seed(seed)
    test_mediation(test_data, x = "X", y = "Y", m = "M",
                   test = "boot", R = R, level = level, type = type,
                   method = "regression", robust = FALSE, family = "gaussian")
  },
  # student = {
  #   set.seed(seed)
  #   suppressWarnings(
  #     test_mediation(test_data, x = "X", y = "Y", m = "M",
  #                    test = "boot", R = R, level = level, type = type,
  #                    method = "regression", robust = FALSE, family = "student")
  #   )
  # },
  # select = {
  #   set.seed(seed)
  #   suppressWarnings(
  #     test_mediation(test_data, x = "X", y = "Y", m = "M",
  #                    test = "boot", R = R, level = level, type = type,
  #                    method = "regression", robust = FALSE, family = "select")
  #   )
  # },
  winsorized = {
    set.seed(seed)
    test_mediation(test_data, x = "X", y = "Y", m = "M",
                   test = "boot", R = R, level = level, type = type,
                   method = "covariance", robust = TRUE, control = cov_ctrl)
  },
  ML = {
    set.seed(seed)
    test_mediation(test_data, x = "X", y = "Y", m = "M",
                   test = "boot", R = R, level = level, type = type,
                   method = "covariance", robust = FALSE)
  }
)

## compute summaries
summary_boot_list <- lapply(boot_list, summary, type = "boot")
summary_data_list <- lapply(boot_list, summary, type = "data")

## correct values
coef_names <- c("a", "b", "Total", "Direct", "Indirect")


## common tests for all model fits

# loop over methods
methods <- names(boot_list)
for (method in methods) {

  # extract information for current method
  boot <- boot_list[[method]]
  summary_boot <- summary_boot_list[[method]]
  summary_data <- summary_data_list[[method]]


  ## run tests

  test_that("output has correct structure", {

    # bootstrap test
    expect_s3_class(boot, "boot_test_mediation")
    expect_s3_class(boot, "test_mediation")
    # model fit
    expect_s3_class(boot$fit, "fit_mediation")
    # bootstrap replicates
    expect_s3_class(boot$reps, "boot")

  })

  test_that("arguments are correctly passed", {

    # alternative hypothesis
    expect_identical(boot$alternative, "twosided")
    # number of bootstrap replicates
    # (doesn't hold for some methods if there are computational issues)
    expect_identical(boot$R, as.integer(R))
    # confidence level
    expect_identical(boot$level, level)
    # type of confidence intervals
    expect_identical(boot$type, type)
    # variable names
    expect_identical(boot$fit$x, "X")
    expect_identical(boot$fit$y, "Y")
    expect_identical(boot$fit$m, "M")
    expect_identical(boot$fit$covariates, character())
    # robust or nonrobust fit and test
    if (method == "robust") {
      expect_identical(boot$fit$robust, "MM")
      expect_equal(boot$fit$control, reg_ctrl)
    } else if (method == "median") {
      expect_identical(boot$fit$robust, "median")
      expect_null(boot$fit$control)
    } else if (method == "winsorized") {
      expect_true(boot$fit$robust)
      expect_equal(boot$fit$control, cov_ctrl)
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
    expect_identical(names(boot$ci), c("Lower", "Upper"))
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
    expect_equivalent(coef(boot, parm = "total", type = "boot"),
                      mean(boot$reps$t[, 7]))
    expect_equivalent(coef(boot, parm = "direct", type = "boot"),
                      mean(boot$reps$t[, 6]))
    expect_equivalent(coef(boot, parm = "indirect", type = "boot"),
                      boot$indirect)

    # effects computed on original sample
    expect_equivalent(coef(boot, parm = "a", type = "data"),
                      boot$fit$a)
    expect_equivalent(coef(boot, parm = "b", type = "data"),
                      boot$fit$b)
    expect_equivalent(coef(boot, parm = "total", type = "data"),
                      boot$fit$total)
    expect_equivalent(coef(boot, parm = "direct", type = "data"),
                      boot$fit$direct)
    expect_equivalent(coef(boot, parm = "indirect", type = "data"),
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
    expect_equivalent(confint(boot, parm = "indirect", type = "boot"), boot$ci)
    # confidence intervals based on theory (except for indirect effect)
    expect_equivalent(confint(boot, parm = "indirect", type = "data"), boot$ci)

  })

  test_that("summary has correct structure", {

    # summary
    expect_s3_class(summary_boot, "summary_test_mediation")
    expect_s3_class(summary_data, "summary_test_mediation")
    # original output of test for indirect effect
    expect_identical(summary_boot$object, boot)
    expect_identical(summary_data$object, boot)
    # summary of the model fit
    expect_s3_class(summary_boot$summary, "summary_fit_mediation")
    expect_s3_class(summary_data$summary, "summary_fit_mediation")
    # diagnostic plot only for ROBMED
    if (method == "robust") {
      expect_s3_class(summary_boot$plot, "gg_weight_plot")
      expect_s3_class(summary_data$plot, "gg_weight_plot")
    } else {
      expect_null(summary_boot$plot)
      expect_null(summary_data$plot)
    }

  })

  test_that("attributes are correctly passed through summary", {

    # robustness
    if (method == "robust") {
      expect_identical(summary_boot$summary$robust, "MM")
      expect_identical(summary_data$summary$robust, "MM")
    } else if (method == "median") {
      expect_identical(summary_boot$summary$robust, "median")
      expect_identical(summary_data$summary$robust, "median")
    } else if (method == "winsorized") {
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
    expect_identical(summary_boot$summary$x, "X")
    expect_identical(summary_boot$summary$y, "Y")
    expect_identical(summary_boot$summary$m, "M")
    expect_identical(summary_boot$summary$covariates, character())
    expect_identical(summary_data$summary$x, "X")
    expect_identical(summary_data$summary$y, "Y")
    expect_identical(summary_data$summary$m, "M")
    expect_identical(summary_data$summary$covariates, character())

  })

  test_that("effect summaries have correct names", {

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

  })

  test_that("effect summaries contain correct coefficient values", {

    # effects computed on original sample
    expect_identical(summary_boot$summary$total["X", "Data"],
                     boot$fit$total)
    expect_identical(summary_boot$summary$direct["X", "Data"],
                     boot$fit$direct)
    expect_identical(summary_data$summary$total["X", "Estimate"],
                     boot$fit$total)
    expect_identical(summary_data$summary$direct["X", "Estimate"],
                     boot$fit$direct)

    # bootstrapped effects
    expect_equal(summary_boot$summary$total["X", "Boot"],
                 mean(boot$reps$t[, 7]))
    expect_equal(summary_boot$summary$direct["X", "Boot"],
                 mean(boot$reps$t[, 6]))

  })

  test_that("output of p_value() method has correct attributes", {

    digits <- 3
    p_boot <- p_value(boot, type = "boot", digits = digits)
    p_data <- p_value(boot, type = "data", digits = digits)
    # bootstrapped effects
    expect_length(p_boot, 5L)
    expect_named(p_boot, coef_names)
    expect_equal(p_boot["Indirect"], round(p_boot["Indirect"], digits = digits))
    # effects computed on original sample
    expect_length(p_data, 5L)
    expect_named(p_data, coef_names)
    expect_equal(p_data["Indirect"], round(p_data["Indirect"], digits = digits))

  })

}



## additional tests for regression fits

# relevant information
reg_methods <- c("robust", "median", "OLS", "student", "select")
classes <- c(robust = "summary_lmrob", median = "summary_rq",
             OLS = "summary_lm", student = "summary_lmse",
             select = "summary_lm")

# loop over methods
for (method in intersect(methods, reg_methods)) {

  # extract information for current method
  boot <- boot_list[[method]]
  summary_boot <- summary_boot_list[[method]]
  summary_data <- summary_data_list[[method]]

  ## correct values
  class <- classes[method]
  family <- if (method %in% c("student", "select")) method else "gaussian"
  if (method == "student") {
    mx_names <- c("(Intercept.DP)", "X")
    ymx_names <- c("(Intercept.DP)", "M", "X")
  } else {
    mx_names <- c("(Intercept)", "X")
    ymx_names <- c("(Intercept)", "M", "X")
  }


  ## run tests

  test_that("output has correct structure", {

    # regression fit
    expect_s3_class(boot$fit, "reg_fit_mediation")

  })

  test_that("arguments are correctly passed", {

    # assumed error distribution
    expect_identical(boot$fit$family, family)
    # mediation model
    expect_identical(boot$fit$model, "simple")
    # no contrasts
    expect_false(boot$fit$contrast)

  })

  test_that("summary has correct structure", {

    # summary of the model fit
    expect_s3_class(summary_boot$summary, "summary_reg_fit_mediation")
    expect_s3_class(summary_data$summary, "summary_reg_fit_mediation")

  })

  # different combinations of model summaries
  combinations <- expand.grid(summary = c("boot", "data"),
                              fit = c("fit_mx", "fit_ymx"),
                              stringsAsFactors = FALSE)
  # correct degrees of freedom in the numerator for F-test
  df1 <- c(1, 1, 2, 2)

  # loop over model summaries
  for (i in 1:nrow(combinations)) {

    # extract model summary
    fit_name <- combinations[i, "fit"]
    model_summary <- switch(combinations[i, "summary"],
                            boot = summary_boot$summary[[fit_name]],
                            data = summary_data$summary[[fit_name]])

    # perform tests
    test_that("model summary has correct structure", {

      # model summary
      expect_s3_class(model_summary, class)
      # regression standard error, R-squared, and F-test
      if (method %in% c("robust", "OLS", "select")) {
        # regression standard error
        expect_type(model_summary$s, "list")
        expect_named(model_summary$s, c("value", "df"))
        # R-squared
        expect_type(model_summary$R2, "list")
        expect_named(model_summary$R2, c("R2", "adj_R2"))
        # F-test
        expect_type(model_summary$F_test, "list")
        expect_named(model_summary$F_test, c("statistic", "df", "p_value"))
        # degrees of freedom of F-test
        df_test <- model_summary$F_test$df
        if (method == "robust") {
          expect_identical(df_test[1], df1[i])
          expect_identical(df_test[2], Inf)
        } else {
          expect_identical(df_test[1], as.integer(df1[i]))
          expect_identical(df_test[2], model_summary$s$df)
        }
      } else {
        # regression standard error
        expect_null(model_summary$s)
        # R-squared
        expect_null(model_summary$R2)
        # F-test
        expect_null(model_summary$F_test)
      }
      # additional information for robust bootstrap test
      if (method == "robust") {
        fit <- boot$fit[[fit_name]]
        # information on convergence
        expect_type(model_summary$algorithm, "list")
        expect_named(model_summary$algorithm, c("converged", "method"))
        expect_identical(model_summary$algorithm$converged,
                         fit$converged)
        expect_identical(model_summary$algorithm$method,
                         fit$control$method)
        # information on outliers
        expect_type(model_summary$outliers, "list")
        expect_named(model_summary$outliers,
                     c("indices", "weights", "threshold"))
        expect_type(model_summary$outliers$indices, "integer")
        expect_identical(model_summary$outliers$weights,
                         weights(fit, type = "robustness"))
        expect_identical(model_summary$outliers$threshold,
                         summary(fit)$control$eps.outlier)
      }

    })

  }

  test_that("effect summaries have correct names", {

    # regression m ~ x
    expect_identical(dim(summary_boot$summary$fit_mx$coefficients),
                     c(2L, 5L))
    expect_identical(rownames(summary_boot$summary$fit_mx$coefficients),
                     mx_names)
    expect_identical(colnames(summary_boot$summary$fit_mx$coefficients)[1:2],
                     c("Data", "Boot"))
    expect_identical(dim(summary_data$summary$fit_mx$coefficients),
                     c(2L, 4L))
    expect_identical(rownames(summary_data$summary$fit_mx$coefficients),
                     mx_names)
    expect_identical(colnames(summary_data$summary$fit_mx$coefficients)[1],
                     "Estimate")
    # regression y ~ m + x
    expect_identical(dim(summary_boot$summary$fit_ymx$coefficients),
                     c(3L, 5L))
    expect_identical(rownames(summary_boot$summary$fit_ymx$coefficient),
                     ymx_names)
    expect_identical(colnames(summary_boot$summary$fit_ymx$coefficient)[1:2],
                     c("Data", "Boot"))
    expect_identical(dim(summary_data$summary$fit_ymx$coefficient),
                     c(3L, 4L))
    expect_identical(rownames(summary_data$summary$fit_ymx$coefficient),
                     ymx_names)
    expect_identical(colnames(summary_data$summary$fit_ymx$coefficient)[1],
                     "Estimate")

  })

  test_that("effect summaries contain correct coefficient values", {

    # effects computed on original sample
    expect_equivalent(summary_boot$summary$fit_mx$coefficients[2, "Data"],
                      boot$fit$a)
    expect_identical(summary_boot$summary$fit_ymx$coefficients[2, "Data"],
                     boot$fit$b)
    expect_equivalent(summary_data$summary$fit_mx$coefficients[2, "Estimate"],
                      boot$fit$a)
    expect_identical(summary_data$summary$fit_ymx$coefficients[2, "Estimate"],
                     boot$fit$b)

    # bootstrapped effects
    expect_equivalent(summary_boot$summary$fit_mx$coefficients[2, "Boot"],
                      mean(boot$reps$t[, 3]))
    expect_equivalent(summary_boot$summary$fit_ymx$coefficients[2, "Boot"],
                      mean(boot$reps$t[, 5]))

  })

}


## additional tests for covariance fits

# relevant information
cov_methods <- c("winsorized", "ML")

# loop over methods
for (method in intersect(methods, cov_methods)) {

  # extract information for current method
  boot <- boot_list[[method]]
  summary_boot <- summary_boot_list[[method]]
  summary_data <- summary_data_list[[method]]


  ## run tests


  test_that("output has correct structure", {

    # covariance fit
    expect_s3_class(boot$fit, "cov_fit_mediation")

  })

  test_that("summary has correct structure", {

    # summary of the model fit
    expect_s3_class(summary_boot$summary, "summary_cov_fit_mediation")
    expect_s3_class(summary_data$summary, "summary_cov_fit_mediation")
    # no regression model summaries
    expect_null(summary_boot$summary$fit_mx)
    expect_null(summary_boot$summary$fit_ymx)
    expect_null(summary_data$summary$fit_mx)
    expect_null(summary_data$summary$fit_ymx)

  })

  test_that("effect summaries have correct names", {

    # a path
    expect_identical(dim(summary_boot$summary$a), c(1L, 5L))
    expect_identical(rownames(summary_boot$summary$a), "X")
    expect_identical(colnames(summary_boot$summary$a)[1:2], c("Data", "Boot"))
    expect_identical(dim(summary_data$summary$a), c(1L, 4L))
    expect_identical(rownames(summary_data$summary$a), "X")
    expect_identical(colnames(summary_data$summary$a)[1], "Estimate")
    # b path
    expect_identical(dim(summary_boot$summary$b), c(1L, 5L))
    expect_identical(rownames(summary_boot$summary$b), "M")
    expect_identical(colnames(summary_boot$summary$b)[1:2], c("Data", "Boot"))
    expect_identical(dim(summary_data$summary$b), c(1L, 4L))
    expect_identical(rownames(summary_data$summary$b), "M")
    expect_identical(colnames(summary_data$summary$b)[1], "Estimate")

  })

  test_that("effect summaries contain correct coefficient values", {

    # effects computed on original sample
    expect_identical(summary_boot$summary$a["X", "Data"], boot$fit$a)
    expect_identical(summary_boot$summary$b["M", "Data"], boot$fit$b)
    expect_identical(summary_data$summary$a["X", "Estimate"], boot$fit$a)
    expect_identical(summary_data$summary$b["M", "Estimate"], boot$fit$b)

    # bootstrapped effects
    expect_equal(summary_boot$summary$a["X", "Boot"],
                 mean(boot$reps$t[, 3]))
    expect_equal(summary_boot$summary$b["M", "Boot"],
                 mean(boot$reps$t[, 5]))

  })

}
