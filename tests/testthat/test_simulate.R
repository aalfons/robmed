context("simulating data based on estimated mediation model")


## load package and data
library("robmed", quietly = TRUE)
data("BSG2014")

## seed to be used for the random number generator
seed <- 20211117

## arguments for bootstrap tests
R <- 100
type <- "perc"

## arguments for data generation
args <- expand.grid(n = c(NA_integer_, 100L),
                    explanatory = c("sim", "boot"),
                    errors = c("sim", "boot"),
                    stringsAsFactors = FALSE)
seq_args <- seq_len(nrow(args))


# ----------------
# simple mediation
# ----------------

## define variables
x <- "ValueDiversity"
y <- "TeamCommitment"
m <- "TaskConflict"

## set seed of the random number generator
set.seed(seed)
## perform mediation analysis
fit <- fit_mediation(BSG2014, x = x, y = y, m = m)
boot <- test_mediation(fit, R = R, type = type)

## loop over arguments for data generation
for (i in seq_args) {

  # number of observations
  n_sim <- args[i, "n"]
  if (is.na(n_sim)) {
    n_sim <- NULL
    n_true <- nrow(BSG2014)
  } else n_true <- n_sim

  # simulate data based on model fit
  set.seed(seed)
  data_fit <- sim_mediation(fit, n = n_sim,
                            explanatory = args[i, "explanatory"],
                            errors = args[i, "errors"])
  # simulate data based on bootstrap test
  set.seed(seed)
  data_boot <- sim_mediation(boot, n = n_sim,
                             explanatory = args[i, "explanatory"],
                             errors = args[i, "errors"])

  # run tests

  test_that("simulated data has correct attributes", {

    expect_identical(dim(data_fit), c(n_true, 3L))
    expect_named(data_fit, c(x, y, m))

  })

  test_that("method for tests inherits from method for model fits", {

    expect_identical(data_boot, data_fit)

  })

}



# -------------------------
# serial multiple mediators
# -------------------------

## define variables
x <- "ValueDiversity"
y <- "TeamScore"
m <- c("TaskConflict", "TeamCommitment")

## set seed of the random number generator
set.seed(seed)
## perform mediation analysis
fit <- fit_mediation(BSG2014, x = x, y = y, m = m, model = "serial")
boot <- test_mediation(fit, R = R, type = type)

## loop over arguments for data generation
for (i in seq_args) {

  # number of observations
  n_sim <- args[i, "n"]
  if (is.na(n_sim)) {
    n_sim <- NULL
    n_true <- nrow(BSG2014)
  } else n_true <- n_sim

  # simulate data based on model fit
  set.seed(seed)
  data_fit <- sim_mediation(fit, n = n_sim,
                            explanatory = args[i, "explanatory"],
                            errors = args[i, "errors"])
  # simulate data based on bootstrap test
  set.seed(seed)
  data_boot <- sim_mediation(boot, n = n_sim,
                             explanatory = args[i, "explanatory"],
                             errors = args[i, "errors"])

  # run tests

  test_that("simulated data has correct attributes", {

    expect_identical(dim(data_fit), c(n_true, 4L))
    expect_named(data_fit, c(x, y, m))

  })

  test_that("method for tests inherits from method for model fits", {

    expect_identical(data_boot, data_fit)

  })

}


# -------------------------------------------------
# parallel multiple mediators and control variables
# -------------------------------------------------

## define variables
x <- "SharedLeadership"
y <- "TeamPerformance"
m <- c("ProceduralJustice", "InteractionalJustice")
covariates <- c("AgeDiversity", "GenderDiversity")

## set seed of the random number generator
set.seed(seed)
## perform mediation analysis
fit <- fit_mediation(BSG2014, x = x, y = y, m = m, covariates = covariates,
                     model = "parallel")
boot <- test_mediation(fit, R = R, type = type)

## loop over arguments for data generation
for (i in seq_args) {

  # number of observations
  n_sim <- args[i, "n"]
  if (is.na(n_sim)) {
    n_sim <- NULL
    n_true <- nrow(BSG2014)
  } else n_true <- n_sim

  # simulate data based on model fit
  set.seed(seed)
  data_fit <- sim_mediation(fit, n = n_sim,
                            explanatory = args[i, "explanatory"],
                            errors = args[i, "errors"])
  # simulate data based on bootstrap test
  set.seed(seed)
  data_boot <- sim_mediation(boot, n = n_sim,
                             explanatory = args[i, "explanatory"],
                             errors = args[i, "errors"])

  # run tests

  test_that("simulated data has correct attributes", {

    expect_identical(dim(data_fit), c(n_true, 6L))
    expect_named(data_fit, c(x, y, m, covariates))

  })

  test_that("method for tests inherits from method for model fits", {

    expect_identical(data_boot, data_fit)

  })

}
