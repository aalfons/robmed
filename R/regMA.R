# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

## fit a mediation model based on regression
regMA <- function(x, y, m, data, robust = TRUE, control=lmrob.control(...), 
                  ...) {
  # initializations
  robust <- isTRUE(robust)
  fYX <- as.formula(paste(y, "~", x))
  fMX <- as.formula(paste(m, "~", x))
  fYMX <- as.formula(paste(y, "~", m, "+", x))
  # compute regression models
  if(robust) {
    fitMX <- lmrob(fMX, data=data, control=control, model=FALSE, x=FALSE)
    fitYX <- lmrob(fYX, data=data, control=control, model=FALSE, x=FALSE)
    fitYMX <- lmrob(fYMX, data=data, control=control, model=FALSE, x=FALSE)
  } else {
    fitYX <- lm(fYX, data=data, model=FALSE)
    fitMX <- lm(fMX, data=data, model=FALSE)
    fitYMX <- lm(fYMX, data=data, model=FALSE)
  }
  # extract coefficients of mediation model
  a <- unname(coef(fitMX)[2])
  b <- unname(coef(fitYMX)[2])
  c <- unname(coef(fitYMX)[3])
  cPrime <- unname(coef(fitYX)[2])
  # return results
  result <- list(a=a, b=b, c=c, cPrime=cPrime, robust=robust, 
                 fitYX=fitYX, fitMX=fitMX, fitYMX=fitYMX)
  class(result) <- c("regMA", "fitMA")
  result
}
  