# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

## compute a robust F-test (which actually uses the chi-squared distribution)

# generic function
robFTest <- function(object, ...) UseMethod("robFTest")

# method for mediation model fits via regression
robFTest.regFitMediation <- function(object, ...) {
  # extract the full model and fit the null model
  fitFull <- object$fitYMX
  fY <- as.formula(paste(object$y, "~ 1"))
  control <- object$control
  fitNull <- lmrob(fY, data = object$fit$data, control = control)
  # perform robust F-test
  s <- fitFull$scale
  df <- fitFull$rank - fitNull$rank
  FTau <- (2 / df) * sum(robust:::chi.weight(residuals(fitNull) / s, ips = 2,
                                             xk = control$tuning.psi) -
                           robust:::chi.weight(residuals(fitFull) / s, ips = 2,
                                               xk = control$tuning.psi))
  tmp <- robust:::lmRob.eff0(itype = 1, ta = control$tuning.psi,
                             tc = control$tuning.psi, ipsi = 2)
  const <- tmp$alfa / tmp$beta
  # if X ~ F(df1, df2), then Y = df1*X ~ chi-squared(df1)
  statistic <- FTau / (const * df)
  pValue <- pchisq(FTau / const, df = df, lower.tail = FALSE)
  # return results
  list(statistic = statistic, df = c(df, Inf), pValue = pValue)
}
