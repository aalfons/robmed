# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

## compute a robust F-test (which actually uses the chi-squared distribution)
#' @useDynLib robmed

# generic function
robFTest <- function(object, ...) UseMethod("robFTest")

# method for mediation model fits via regression
robFTest.regFitMediation <- function(object, ...) {
  # check if Tukey bisquare loss function is used, otherwise return NULL
  control <- object$control
  if (control$psi != "bisquare") return()
  # extract the full model and fit the null model
  fitFull <- object$fitYMX
  fY <- as.formula(paste(object$y, "~ 1"))
  fitNull <- lmrob(fY, data = object$fit$data, control = control)
  # perform robust F-test
  s <- fitFull$scale
  df <- fitFull$rank - fitNull$rank
  FTau <- (2 / df) * sum(chi.weight(residuals(fitNull)/s, control$tuning.psi) -
                           chi.weight(residuals(fitFull)/s, control$tuning.psi))
  tmp <- lmRob.eff0(control$tuning.psi)
  const <- tmp$alfa / tmp$beta
  # if X ~ F(df1, df2), then Y = df1*X ~ chi-squared(df1)
  statistic <- FTau / (const * df)
  pValue <- pchisq(FTau / const, df = df, lower.tail = FALSE)
  # return results
  list(statistic = statistic, df = c(df, Inf), pValue = pValue)
}


## utility functions

# These are modified from package 'robust'.  Unfortunately, package 'robust'
# does not export those functions.  The original Fortran code is not modified.

# weight function
chi.weight <- function(x, tuning.psi) {
  n <- length(x)
  result <- .Fortran("rlchiam2",  # ../src/lmrobmm.f
           as.integer(n),
           as.double(x),
           fvals = double(n),
           2L,  # this (integer) value indicates to use Tukey bisquare loss
           as.double(tuning.psi),
           PACKAGE = "robmed")
  result$fvals
}

# asymptotic relative efficiency (ARE) of a GM-estimator (Tukey bisquare loss):
lmRob.eff0 <- function(tuning.psi, upper = 10, til = 1e-4, maxit = 150,
                       tol = 1e-5, eps = .Machine$double.eps,
                       uflow = .Machine$double.xmin) {
  index <- c(1, 1, 1,
             2, # this value indicates to use Tukey bisquare loss
             1, 1, 0)
  tuning <- c(0, 0, tuning.psi, 0, tuning.psi, eps, uflow, upper, til)
  result <- .Fortran("rlref0bi",  # ../src/lmrobbi.f
           as.integer(index),
           as.double(tuning),
           xlcnst = -1.0,
           1L,
           1.0,
           as.integer(maxit),
           as.double(tol),
           nit = 0L,
           alfa = 0.0,
           beta = 0.0,
           reff = 0.0,
           PACKAGE = "robmed")
  result[c("nit", "alfa", "beta", "reff")]
}
