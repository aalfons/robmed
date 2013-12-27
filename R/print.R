# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

#' @S3method print bootMA
print.bootMA <- function(x, digits = max(3, getOption("digits")-3), ...) {
  cat("Bootstrap results for indirect effect\n")
  # print indirect effect
  cat("\nIndirect effect (ab path):\n")
  a <- x$fit$a
  b <- x$fit$b
  ab <- cbind(a*b, x$ab)
  m <- names(x$data)[3]
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

#' @S3method print covHuber
print.covHuber <- function(x, ...) {
  cat("Huber-type M-estimator\n")
  cat("\nLocation estimate:\n")
  print(x$center)
  cat("\nScatter matrix estimate:\n")
  print(x$cov)
}

#' @S3method print sobelMA
print.sobelMA <- function(x, digits = max(3, getOption("digits")-3), ...) {
  cat("Normal theory test for indirect effect\n")
  # print indirect effect
  cat("\nIndirect effect (ab path):\n")
  ab <- cbind(x$ab, x$se, x$statistic, x$pValue)
  m <- names(x$data)[3]
  cn <- switch(x$alternative, twosided="Pr(>|z|)", 
               less="Pr(<z)", greater="Pr(>z)")
  dimnames(ab) <- list(m, c("Estimate", "Std. Error", "z value", cn))
  printCoefmat(ab, digits=digits, ...)
  # return object invisibly
  invisible(x)
}


#' @S3method print summaryTestMA
print.summaryTestMA <- function(x, digits = max(3, getOption("digits")-3), 
                            signif.stars = getOption("show.signif.stars"), 
                            signif.legend = signif.stars, ...) {
  # initializations
  s <- x$summary
  x <- x$object
  # print information on data
  cn <- names(x$data)
  cat("\nIndependent, dependent and proposed mediator variables:\n")
  cat(sprintf("x = %s\n", cn[1]))
  cat(sprintf("y = %s\n", cn[2]))
  cat(sprintf("m = %s\n", cn[3]))
  # print sample size
  cat(sprintf("\nSample size: %d\n", s$n))
  # print effects
  cat("---\nEffect of x on m (a path):\n")
  printCoefmat(s$a, digits=digits, signif.stars=signif.stars, 
               signif.legend=FALSE, ...)
  cat("\nDirect effect of m on y (b path):\n")
  printCoefmat(s$b, digits=digits, signif.stars=signif.stars, 
               signif.legend=FALSE, ...)
  cat("\nDirect effect of x on y (c path):\n")
  printCoefmat(s$c, digits=digits, signif.stars=signif.stars, 
               signif.legend=FALSE, ...)
  cat("\nTotal effect of x on y (c' path):\n")
  printCoefmat(s$cPrime, digits=digits, signif.stars=signif.stars, 
               signif.legend=FALSE, ...)
  # print model summary for y ~ m + x
  cat("---\nModel summary for y ~ m + x\n")
  prefix <- if(s$robust) "Robust residual" else "Residual"
  postfix <- sprintf(" on %d degrees of freedom", s$s$df)
  cat("\n", prefix, " standard error: ", format(signif(s$s$value, digits)), 
      postfix, "\n", sep="")
  if(!is.null(s$FTest)) {
    cat("Multiple R-squared:  ", formatC(s$FTest$R2, digits=digits), 
        ",\tAdjusted R-squared:  ", formatC(s$FTest$adjR2, digits=digits), 
        "\nF-statistic: ", formatC(s$FTest$statistic, digits=digits), " on ", 
        s$FTest$df[1], " and ", s$FTest$df[2], " DF,  p-value: ", 
        format.pval(s$FTest$pValue, digits=digits), "\n", sep="")
  }
  # print indirect effect
  cat("---\n")
  print(x, digits=digits, signif.stars=signif.stars, signif.legend=FALSE, ...)
  # print legend for significance stars
  if(isTRUE(signif.stars) && isTRUE(signif.legend)) {
    cat("---\nSignif. codes:  0", sQuote("***"), "0.001", sQuote("**"), 
        "0.01", sQuote("*"), "0.05", sQuote("."), "0.1", sQuote(" "), "1\n")
  }
  # return object invisibly
  invisible(x)
}
