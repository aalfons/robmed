# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

#' @export
print.bootTestMediation <- function(x, digits = max(3, getOption("digits")-3),
                                    ...) {
  cat("Bootstrap results for indirect effect\n")
  # print indirect effect
  cat("\nIndirect effect (ab path):\n")
  a <- x$fit$a
  b <- x$fit$b
  ab <- cbind(a*b, x$ab)
  m <- names(x$fit$data)[3]
  dimnames(ab) <- list(m, c("Data", "Boot"))
  print(ab, digits=digits, ...)
  # print confidence interval
  cat("\n")
  cat(format(100 * x$level), "percent confidence interval:\n")
  ci <- t(x$ci)
  dimnames(ci) <- list(m, c("Lower", "Upper"))
  print(ci, digits=digits, ...)
  # print additional information
  cat(sprintf("\nNumber of bootstrap replicates: %d\n", x$R))
  # return object invisibly
  invisible(x)
}

#' @export
print.covHuber <- function(x, ...) {
  # print estimates
  cat("Huber M-estimator\n")
  cat("\nLocation vector estimate:\n")
  print(x$center, ...)
  cat("\nScatter matrix estimate:\n")
  print(x$cov, ...)
  # return object invisibly
  invisible(x)
}

#' @export
print.covML <- function(x, ...) {
  # print estimates
  cat("Maximum likelihood estimator\n")
  cat("\nMean vector estimate:\n")
  print(x$center, ...)
  cat("\nCovariance matrix estimate:\n")
  print(x$cov, ...)
  # return object invisibly
  invisible(x)
}

#' @export
print.fitMediation <- function(x, ...) {
  # print estimated effects
  prefix <- if(x$robust) "Robust mediation" else "Mediation"
  cat(prefix, "model fit\n")
  cat("\nEffects:\n")
  print(coefficients(x), ...)
  # return object invisibly
  invisible(x)
}

#' @export
print.sobelTestMediation <- function(x, digits = max(3, getOption("digits")-3),
                                     ...) {
  # print indirect effect
  cat("Normal theory test for indirect effect\n")
  cat("\nIndirect effect (ab path):\n")
  ab <- cbind(x$ab, x$se, x$statistic, x$pValue)
  m <- x$fit$m
  cn <- switch(x$alternative, twosided="Pr(>|z|)",
               less="Pr(<z)", greater="Pr(>z)")
  dimnames(ab) <- list(m, c("Estimate", "Std. Error", "z value", cn))
  printCoefmat(ab, digits=digits, ...)
  # return object invisibly
  invisible(x)
}

#' @export
print.summaryFitMediation <- function(x, digits = max(3, getOption("digits")-3),
                                      signif.stars = getOption("show.signif.stars"),
                                      signif.legend = signif.stars, ...) {
  # initializations
  p <- length(x$variables)
  haveCovariates <- p > 3
  # print information on data
  cat("\nIndependent, dependent and proposed mediator variables:\n")
  cat(sprintf("x = %s\n", x$variables[1]))
  cat(sprintf("y = %s\n", x$variables[2]))
  cat(sprintf("m = %s\n", x$variables[3]))
  if(haveCovariates) {
    cat("\nControl variables:\n")
    print(x$variables[4:p], quote=FALSE)
  }
  # print sample size
  cat(sprintf("\nSample size: %d\n", x$n))
  # print effects
  if(haveCovariates) cat("---\nPartial effect of x on m (a path):\n")
  else cat("---\nEffect of x on m (a path):\n")
  printCoefmat(x$a, digits=digits, signif.stars=signif.stars,
               signif.legend=FALSE, ...)
  cat("\nDirect effect of m on y (b path):\n")
  printCoefmat(x$b, digits=digits, signif.stars=signif.stars,
               signif.legend=FALSE, ...)
  cat("\nDirect effect of x on y (c path):\n")
  printCoefmat(x$c, digits=digits, signif.stars=signif.stars,
               signif.legend=FALSE, ...)
  cat("\nTotal effect of x on y (c' path):\n")
  printCoefmat(x$cPrime, digits=digits, signif.stars=signif.stars,
               signif.legend=FALSE, ...)
  if(haveCovariates) {
    cat("\nPartial effects of control variables on y:\n")
    printCoefmat(x$covariates, digits=digits, signif.stars=signif.stars,
                 signif.legend=FALSE, ...)
  }
  # print model summary for y ~ m + x + covariates
  postfix <- if(haveCovariates) " + control variables" else ""
  cat(sprintf("---\nModel summary for y ~ m + x%s\n", postfix))
  if(x$robust) {
    cat("\nRobust residual standard error: ", format(signif(x$s$value, digits)),
        "\n", sep="")
    if(!is.null(x$FTest)) {
      cat("Robust R-squared:  ", formatC(x$FTest$R2, digits=digits),
          ",\tAdjusted robust R-squared:  ", formatC(x$FTest$adjR2, digits=digits),
          "\n", sep="")
    }
  } else {
    postfix <- sprintf(" on %d degrees of freedom", x$s$df)
    cat("\nResidual standard error: ", format(signif(x$s$value, digits)),
        postfix, "\n", sep="")
    if(!is.null(x$FTest)) {
      cat("Multiple R-squared:  ", formatC(x$FTest$R2, digits=digits),
          ",\tAdjusted R-squared:  ", formatC(x$FTest$adjR2, digits=digits),
          "\nF-statistic: ", formatC(x$FTest$statistic, digits=digits), " on ",
          x$FTest$df[1], " and ", x$FTest$df[2], " DF,  p-value: ",
          format.pval(x$FTest$pValue, digits=digits), "\n", sep="")
    }
  }
  # print legend for significance stars
  if(isTRUE(signif.stars) && isTRUE(signif.legend)) printLegend()
  # return object invisibly
  invisible(x)
}

#' @export
print.summaryTestMediation <- function(x, digits = max(3, getOption("digits")-3),
                                       signif.stars = getOption("show.signif.stars"),
                                       signif.legend = signif.stars, ...) {
  # print summary of mediation model fit
  print(x$summary, digits=digits, signif.stars=signif.stars,
        signif.legend=FALSE, ...)
  # print indirect effect
  cat("---\n")
  print(x$object, digits=digits, signif.stars=signif.stars,
        signif.legend=FALSE, ...)
  # print legend for significance stars
  if(isTRUE(signif.stars) && isTRUE(signif.legend)) printLegend()
  # return object invisibly
  invisible(x)
}

## internal function to print legend for significance stars
printLegend <- function() {
  cat("---\nSignif. codes:  0", sQuote("***"), "0.001", sQuote("**"),
      "0.01", sQuote("*"), "0.05", sQuote("."), "0.1", sQuote(" "), "1\n")
}
