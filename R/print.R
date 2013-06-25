# ------------------------------------
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ------------------------------------

#' @S3method print bootMA
print.bootMA <- function(x, digits = max(3, getOption("digits")-3), ...) {
  cat("Bootstrap results for indirect effect\n")
  # print indirect effect
  cat("\nIndirect effect (ab path)\n")
  a <- unname(coef(x$fitMX)[2])
  b <- unname(coef(x$fitYMX)[2])
  ab <- t(c(a*b, x$ab))
  m <- names(x$data)[3]
  dimnames(ab) <- list(m, c("Data", "Boot"))
  print(ab, digits=digits, ...)
  # print confidence interval
  cat("\n")
  cat(format(100 * (1-x$alpha)), "percent confidence interval:\n")
  ci <- t(x$ci)
  dimnames(ci) <- list(m, c("Lower", "Upper"))
  print(ci, digits=digits, ...)
  # print additional information
  cat(sprintf("\nNumber of bootstrap replicates: %d\n", x$R))
  # return object invisibly
  invisible(x)
}

#' @S3method print sobelMA
print.sobelMA <- function(x, digits = max(3, getOption("digits")-3), ...) {
  cat("Normal theory test for indirect effect\n")
  # print indirect effect
  cat("\nIndirect effect (ab path)\n")
  ab <- t(c(x$ab, x$se, x$statistic, x$pValue))
  m <- names(x$data)[3]
  pn <- switch(x$alternative, twosided="Pr(>|z|)", 
               less="Pr(<z)", greater="Pr(>z)")
  dimnames(ab) <- list(m, c("Estimate", "Std. Error", "z value", pn))
  printCoefmat(ab, digits=digits, ...)
  # return object invisibly
  invisible(x)
}


#' @S3method print summaryMA
print.summaryMA <- function(x, digits = max(3, getOption("digits")-3), 
                            signif.stars = getOption("show.signif.stars"), 
                            signif.legend = signif.stars, ...) {
  # print information on data
  nam <- names(x$object$data)
  cat("\nIndependent, dependent and proposed mediator variables:\n")
  cat(sprintf("x = %s\n", nam[1]))
  cat(sprintf("y = %s\n", nam[2]))
  cat(sprintf("m = %s\n", nam[3]))
  # print sample size
  cat(sprintf("\nSample size: %d\n", nobs(x$object$fitYMX)))
  # print effects
  a <- x$summaryMX$coefficients[-1, , drop=FALSE]
  bc <- x$summaryYMX$coefficients[-1, , drop=FALSE]
  cPrime <- x$summaryYX$coefficients[-1, , drop=FALSE]
  cat("---\nEffect of x on m (a path)\n")
  printCoefmat(a, digits=digits, signif.stars=signif.stars, 
               signif.legend=FALSE, ...)
  cat("\nDirect effect of m on y (b path)\n")
  printCoefmat(bc[1, , drop=FALSE], digits=digits, signif.stars=signif.stars, 
               signif.legend=FALSE, ...)
  cat("\nTotal effect of x on y (c' path)\n")
  printCoefmat(cPrime, digits=digits, signif.stars=signif.stars, 
               signif.legend=FALSE, ...)
  cat("\nDirect effect of x on y (c path)\n")
  printCoefmat(bc[2, , drop=FALSE], digits=digits, signif.stars=signif.stars, 
               signif.legend=FALSE, ...)
  # print model summary for y ~ m + x
  s <- x$summaryYMX
  cat("---\nModel summary for y ~ m + x\n")
  if(x$object$robust && !x$object$transform) {
    cat("\nRobust residual standard error: ", format(signif(s$scale, digits)), 
        "\n", sep="")
  } else {
    cat("\nResidual standard error:", format(signif(s$sigma, digits)), 
        "on", s$df[2], "degrees of freedom\n")
    if(nzchar(mess <- naprint(s$na.action))) cat("  (", mess, ")\n", sep = "")
    if(!is.null(s$fstatistic)) {
      cat("Multiple R-squared: ", formatC(s$r.squared, digits=digits))
      cat(",\tAdjusted R-squared: ", formatC(s$adj.r.squared, digits=digits), 
          "\nF-statistic:", formatC(s$fstatistic[1], digits=digits), "on", 
          s$fstatistic[2], "and", s$fstatistic[3], "DF,  p-value:", 
          format.pval(pf(s$fstatistic[1], s$fstatistic[2], s$fstatistic[3], 
                         lower.tail=FALSE), digits=digits))
      cat("\n")
    }
  }
  # print indirect effect
  cat("---\n")
  print(x$object, digits=digits, signif.stars=signif.stars, 
        signif.legend=FALSE, ...)
  # print legend for significance stars
  if(isTRUE(signif.stars) && isTRUE(signif.legend)) {
    cat("---\nSignif. codes:  0", sQuote("***"), "0.001", sQuote("**"), 
        "0.01", sQuote("*"), "0.05", sQuote("."), "0.1", sQuote(" "), "1\n")
  }
  # return object invisibly
  invisible(x)
}
