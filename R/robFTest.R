# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

## compute a robust F-test (which actually uses the chi-squared distribution)

# generic function
robFTest <- function(object, ...) UseMethod("robFTest")

# method for mediation model fits via regression
robFTest.regFitMediation <- function(object, ...) {
  # check if Tukey bisquare loss function is used, otherwise return NULL
  control <- object$control
  if (control$psi != "bisquare" || is.null(control$efficiency)) return()
  # extract the full model and fit the null model
  fitFull <- object$fitYMX
  fY <- as.formula(paste(object$y, "~ 1"))
  fitNull <- lmrob(fY, data = object$fit$data, control = control)
  # perform robust F-test
  s <- fitFull$scale
  df <- fitFull$rank - fitNull$rank
  chiFull <- Mchi(residuals(fitFull) / s, control$tuning.psi, control$psi)
  chiNull <- Mchi(residuals(fitNull) / s, control$tuning.psi, control$psi)
  tau <- (2 / df) * sum(chiNull - chiFull)
  if(control$efficiency == 0.8)       const <- 0.4140647
  else if(control$efficiency == 0.85) const <- 0.3572601  # default
  else if(control$efficiency == 0.9)  const <- 0.2953806
  else if(control$efficiency == 0.95) const <- 0.2180425
  # Y = df1*X ~ chi-squared(df1) if X ~ F(df1, df2) and df2 goes to infinity
  chi2 <- tau / const
  pValue <- pchisq(chi2, df = df, lower.tail = FALSE)
  # return results
  list(statistic = chi2 / df, df = c(df, Inf), pValue = pValue)
}
