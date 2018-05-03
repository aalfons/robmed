# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

#' @export
print.boot_test_mediation <- function(x, digits = max(3, getOption("digits")-3),
                                      ...) {
  cat("Bootstrap results for indirect effect\n")
  # print indirect effect
  cat("\nIndirect effect (ab path):\n")
  m <- x$fit$m
  p_m <- length(m)
  a <- x$fit$a
  b <- x$fit$b
  ab <- a*b
  if(p_m > 1L) ab <- c(Total = sum(ab), ab)
  ab <- cbind(Data = ab, Boot = x$ab)
  if(p_m == 1) rownames(ab) <- m
  print(ab, digits=digits, ...)
  # print confidence interval
  cat("\n")
  cat(format(100 * x$level), "percent confidence interval:\n")
  if(p_m == 1) {
    ci <- t(x$ci)
    dimnames(ci) <- list(m, c("Lower", "Upper"))
  } else {
    ci <- x$ci
    colnames(ci) <- c("Lower", "Upper")
  }
  print(ci, digits=digits, ...)
  # print additional information
  cat(sprintf("\nNumber of bootstrap replicates: %d\n", x$R))
  # return object invisibly
  invisible(x)
}

#' @export
print.cov_Huber <- function(x, ...) {
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
print.cov_ML <- function(x, ...) {
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
print.fit_mediation <- function(x, ...) {
  # print estimated effects
  prefix <- if(x$robust) "Robust mediation" else "Mediation"
  cat(prefix, "model fit\n")
  cat("\nEffects:\n")
  print(coefficients(x), ...)
  # return object invisibly
  invisible(x)
}

#' @export
print.sobel_test_mediation <- function(x, digits = max(3, getOption("digits")-3),
                                       ...) {
  # print indirect effect
  cat("Normal theory test for indirect effect\n")
  cat("\nIndirect effect (ab path):\n")
  ab <- cbind(x$ab, x$se, x$statistic, x$p_value)
  m <- x$fit$m
  cn <- switch(x$alternative, twosided="Pr(>|z|)",
               less="Pr(<z)", greater="Pr(>z)")
  dimnames(ab) <- list(m, c("Estimate", "Std. Error", "z value", cn))
  printCoefmat(ab, digits=digits, ...)
  # return object invisibly
  invisible(x)
}

#' @export
print.summary_fit_mediation <- function(x, digits = max(3, getOption("digits")-3),
                                        signif.stars = getOption("show.signif.stars"),
                                        signif.legend = signif.stars, ...) {
  # initializations
  p_m <- length(x$m)
  have_covariates <- length(x$covariates) > 0L
  # print information on data
  cat("\nIndependent, dependent and proposed mediator variables:\n")
  cat(sprintf("x = %s\n", x$x))
  cat(sprintf("y = %s\n", x$y))
  cat(sprintf("m = %s\n", x$m[1]))
  if(p_m > 1) {
    for(j in seq(2, p_m)) cat(sprintf("    %s\n", x$m[j]))
  }
  if(have_covariates) {
    cat("\nControl variables:\n")
    print(x$covariates, quote=FALSE)
  }
  # print sample size
  cat(sprintf("\nSample size: %d\n", x$n))
  # print effects
  if(have_covariates) cat("---\nPartial effect of x on m (a path):\n")
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
  printCoefmat(x$c_prime, digits=digits, signif.stars=signif.stars,
               signif.legend=FALSE, ...)
  if(have_covariates) {
    cat("\nPartial effects of control variables on y:\n")
    printCoefmat(x$covariate_effects, digits=digits, signif.stars=signif.stars,
                 signif.legend=FALSE, ...)
  }
  # print model summary for y ~ m + x + covariates
  postfix <- if(have_covariates) " + control variables" else ""
  cat(sprintf("---\nModel summary for y ~ m + x%s\n", postfix))
  postfix <- sprintf(" on %d degrees of freedom", x$s$df)
  if(x$robust) {
    cat("\nRobust residual standard error: ", format(signif(x$s$value, digits)),
        postfix, "\n", sep="")
    if(!is.null(x$F_test)) {
      cat("Robust R-squared:  ",
          formatC(x$R2$R2, digits=digits),
          ",\tAdjusted robust R-squared:  ",
          formatC(x$R2$adj_R2, digits=digits),
          "\nRobust F-statistic: ", formatC(x$F_test$statistic, digits=digits),
          " on ", x$F_test$df[1], " and ", x$F_test$df[2], " DF,  p-value: ",
          format.pval(x$F_test$p_value, digits=digits), "\n", sep="")
    }
  } else {
    cat("\nResidual standard error: ", format(signif(x$s$value, digits)),
        postfix, "\n", sep="")
    if(!is.null(x$F_test)) {
      cat("Multiple R-squared:  ", formatC(x$R2$R2, digits=digits),
          ",\tAdjusted R-squared:  ", formatC(x$R2$adj_R2, digits=digits),
          "\nF-statistic: ", formatC(x$F_test$statistic, digits=digits),
          " on ", x$F_test$df[1], " and ", x$F_test$df[2], " DF,  p-value: ",
          format.pval(x$F_test$p_value, digits=digits), "\n", sep="")
    }
  }
  # print legend for significance stars
  if(isTRUE(signif.stars) && isTRUE(signif.legend)) print_legend()
  # return object invisibly
  invisible(x)
}

#' @export
print.summary_test_mediation <- function(x, digits = max(3, getOption("digits")-3),
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
  if(isTRUE(signif.stars) && isTRUE(signif.legend)) print_legend()
  # return object invisibly
  invisible(x)
}

## internal function to print legend for significance stars
print_legend <- function() {
  cat("---\nSignif. codes:  0", sQuote("***"), "0.001", sQuote("**"),
      "0.01", sQuote("*"), "0.05", sQuote("."), "0.1", sQuote(" "), "1\n")
}
