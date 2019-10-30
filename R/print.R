# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

# information to display for different objects in first line of output
print_info <- function(x, ...) UseMethod("print_info")

# information to display for regression fits
print_info.reg_fit_mediation <- function(x, ...) {
  prefix <- if (x$robust && !x$median) "Robust mediation" else "Mediation"
  postfix <- if (x$median) "via median regression" else "via regression"
  cat(sprintf("%s model fit %s\n", prefix, postfix))
}

# information to display for covariance matrix fits
print_info.cov_fit_mediation <- function(x, ...) {
  prefix <- if (x$robust) "Robust mediation" else "Mediation"
  cat(sprintf("%s model fit via covariance matrix\n", prefix))
}

# information to display for bootstrap tests
print_info.boot_test_mediation <- function(x, ...) {
  # type of test and model fit
  fit <- x$fit
  if (inherits(x$fit, "reg_fit_mediation")) {
    prefix <- if (fit$robust && !fit$median) "Robust bootstrap" else "Bootstrap"
    postfix <- if (fit$median) " via median regression" else " via regression"
  } else if (inherits(fit, "cov_fit_mediation")) {
    prefix <- if (fit$robust) "Robust bootstrap" else "Bootstrap"
    postfix <- " via covariance matrix"
  } else {
    prefix <- "Bootstrap"
    postfix <- ""
  }
  # use plural for multiple mediators
  p_m <- length(x$fit$m)
  plural <- if (length(p_m) == 1L) "" else "s"
  # return message
  cat(sprintf("%s test%s for indirect effect%s%s\n",
              prefix, plural, plural, postfix))
}

# information to display for sobel tests
print_info.sobel_test_mediation <- function(x, ...) {
  # type of test and model fit
  fit <- x$fit
  if (inherits(x$fit, "reg_fit_mediation")) {
    prefix <- if (fit$robust && !fit$median) "Robust normal" else "Normal"
    postfix <- if (fit$median) " via median regression" else " via regression"
  } else if (inherits(fit, "cov_fit_mediation")) {
    prefix <- if (fit$robust) "Robust normal" else "Normal"
    postfix <- " via covariance matrix"
  } else {
    prefix <- "Normal"
    postfix <- ""
  }
  # return message
  cat(sprintf("%s test for indirect effect%s\n", prefix, postfix))
}

# information on variables in summary of mediation analysis
print_info.summary_fit_mediation <- function(x, ...) {
  # initializations
  p_m <- length(x$m)
  have_covariates <- length(x$covariates) > 0L
  # print information on variables
  if (p_m == 1L) {
    cat(sprintf("x = %s\n", x$x))
    cat(sprintf("y = %s\n", x$y))
    cat(sprintf("m = %s\n", x$m))
  } else {
    width <- nchar(p_m) + 1L
    cat(sprintf(paste0("%-", width, "s = %s\n"), "x", x$x))
    cat(sprintf(paste0("%-", width, "s = %s\n"), "y", x$y))
    m_labels <- paste0("m", seq_len(p_m))
    cat(sprintf(paste0("%-", width, "s = %s\n"), m_labels, x$m), sep = "")
  }
  if (have_covariates) {
    cat("\nCovariates:\n")
    print(x$covariates, quote = FALSE)
  }
  # print sample size
  cat(sprintf("\nSample size: %d\n", x$n))
}


#' @export
print.fit_mediation <- function(x, info = TRUE, ...) {
  # print information on type of model fit
  if (isTRUE(info)) print_info(x, ...)
  # print estimated effects
  cat("\nEffects:\n")
  print(coef(x), ...)
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
print.boot_test_mediation <- function(x, digits = max(3, getOption("digits")-3),
                                      info = TRUE, ...) {
  # print information on type of model fit
  if (isTRUE(info)) print_info(x, ...)
  # initializations
  m <- x$fit$m
  p_m <- length(m)
  plural <- if (p_m == 1L) "" else "s"
  cat(sprintf("\nIndirect effect%s of x on y:\n", plural))
  # extract indirect effect
  a <- x$fit$a
  b <- x$fit$b
  ab <- a*b
  if (p_m > 1L) ab <- c(Total = sum(ab), ab)
  ab <- cbind(Data = ab, Boot = x$ab)
  if (p_m == 1L) rownames(ab) <- m
  # extract confidence interval
  ci <- if (p_m == 1L) t(x$ci) else x$ci
  colnames(ci) <- c("Lower", "Upper")
  # combine and print
  print(cbind(ab, ci), digits = digits, ...)
  # print additional information
  cat("---\nLevel of confidence: ", format(100 * x$level), " %\n", sep = "")
  cat(sprintf("\nNumber of bootstrap replicates: %d\n", x$R))
  ## return object invisibly
  invisible(x)
}

#' @export
print.sobel_test_mediation <- function(x, digits = max(3, getOption("digits")-3),
                                       info = TRUE, ...) {
  # print information on type of model fit
  if (isTRUE(info)) print_info(x, ...)
  # print indirect effect
  cat("\nIndirect effect of x on y:\n")
  ab <- cbind(x$ab, x$se, x$statistic, x$p_value)
  m <- x$fit$m
  cn <- switch(x$alternative, twosided = "Pr(>|z|)",
               less = "Pr(<z)", greater = "Pr(>z)")
  dimnames(ab) <- list(m, c("Estimate", "Std. Error", "z value", cn))
  printCoefmat(ab, digits = digits, ...)
  # return object invisibly
  invisible(x)
}


#' @export
print.summary_lm <- function(x, digits = max(3, getOption("digits")-3),
                             signif.stars = getOption("show.signif.stars"),
                             signif.legend = signif.stars, ...) {
  # print coefficient matrix
  cat("Coefficients:\n")
  printCoefmat(x$coefficients, digits = digits, signif.stars = signif.stars,
               signif.legend = signif.legend, ...)
  # print model summary
  cat("\nResidual standard error:", format(signif(x$s$value, digits)),
      "on", x$s$df, "degrees of freedom\n")
  cat("Multiple R-squared:  ", formatC(x$R2$R2, digits = digits),
      ",\tAdjusted R-squared:  ", formatC(x$R2$adj_R2, digits = digits),
      "\n", sep = "")
  cat("F-statistic: ", formatC(x$F_test$statistic, digits = digits),
      " on ", x$F_test$df[1], " and ", x$F_test$df[2], " DF,  p-value: ",
      format.pval(x$F_test$p_value, digits = digits), "\n", sep = "")
  # return object invisibly
  invisible(x)
}

#' @export
print.summary_lmrob <- function(x, digits = max(3, getOption("digits")-3),
                                signif.stars = getOption("show.signif.stars"),
                                signif.legend = signif.stars, ...) {
  # print coefficient matrix
  cat("Coefficients:\n")
  printCoefmat(x$coefficients, digits = digits, signif.stars = signif.stars,
               signif.legend = signif.legend, ...)
  # print model summary
  cat("\nRobust residual standard error:", format(signif(x$s$value, digits)),
      "on", x$s$df, "degrees of freedom\n")
  cat("Robust R-squared:  ", formatC(x$R2$R2, digits = digits),
      ",\tAdjusted robust R-squared:  ", formatC(x$R2$adj_R2, digits = digits),
      "\n", sep = "")
  cat("Robust F-statistic: ", formatC(x$F_test$statistic, digits = digits),
      " on ", x$F_test$df[1], " and ", x$F_test$df[2], " DF,  p-value: ",
      format.pval(x$F_test$p_value, digits = digits), "\n", sep = "")
  # return object invisibly
  invisible(x)
}

#' @export
print.summary_rq <- function(x, digits = max(3, getOption("digits")-3),
                             signif.stars = getOption("show.signif.stars"),
                             signif.legend = signif.stars, ...) {
  # print coefficient matrix
  cat("Coefficients:\n")
  printCoefmat(x$coefficients, digits = digits, signif.stars = signif.stars,
               signif.legend = signif.legend, ...)
  # return object invisibly
  invisible(x)
}


#' @export
print.summary_reg_fit_mediation <- function(x, digits = max(3, getOption("digits")-3),
                                            signif.stars = getOption("show.signif.stars"),
                                            signif.legend = signif.stars, ...) {
  # initializations
  p_m <- length(x$m)
  # information on variables
  print_info(x, ...)
  # print summary of regression m ~ x + covariates
  if (p_m == 1L) {
    cat(sprintf("---\nOutcome variable: %s\n\n", x$m))
    print(x$fit_mx, digits = digits, signif.stars = signif.stars,
          signif.legend = FALSE, ...)
  } else {
    for (m in x$m) {
      cat(sprintf("---\nOutcome variable: %s\n\n", m))
      print(x$fit_mx[[m]], digits = digits, signif.stars = signif.stars,
            signif.legend = FALSE, ...)
    }
  }
  # print summary of regression y ~ m + x + covariates
  cat(sprintf("---\nOutcome variable: %s\n\n", x$y))
  print(x$fit_ymx, digits = digits, signif.stars = signif.stars,
        signif.legend = FALSE, ...)
  # print summary of total and direct effects of x on y
  cat("---\nTotal effect of x on y:\n")
  printCoefmat(x$total, digits = digits, signif.stars = signif.stars,
               signif.legend = FALSE, ...)
  cat("\nDirect effect of x on y:\n")
  printCoefmat(x$direct, digits = digits, signif.stars = signif.stars,
               signif.legend = FALSE, ...)
  # return object invisibly
  invisible(x)
}

#' @export
print.summary_cov_fit_mediation <- function(x, digits = max(3, getOption("digits")-3),
                                            signif.stars = getOption("show.signif.stars"),
                                            signif.legend = signif.stars, ...) {
  # information on variables
  print_info(x, ...)
  # print summary of effect of x on m
  cat("---\nEffect of x on m:\n")
  printCoefmat(x$a, digits = digits, signif.stars = signif.stars,
               signif.legend = FALSE, ...)
  # print summary of direct effect of m on y
  cat("\nDirect effect of m on y:\n")
  printCoefmat(x$b, digits = digits, signif.stars = signif.stars,
               signif.legend = FALSE, ...)
  # print summary of total and direct effects of x on y
  cat("\nTotal effect of x on y:\n")
  printCoefmat(x$total, digits = digits, signif.stars = signif.stars,
               signif.legend = FALSE, ...)
  cat("\nDirect effect of x on y:\n")
  printCoefmat(x$direct, digits = digits, signif.stars = signif.stars,
               signif.legend = FALSE, ...)
  # return object invisibly
  invisible(x)
}


#' @export
print.summary_test_mediation <- function(x, digits = max(3, getOption("digits")-3),
                                         signif.stars = getOption("show.signif.stars"),
                                         signif.legend = signif.stars, ...) {
  # print information on type of model fit
  print_info(x$object, ...)
  cat("\n")
  # print summary of mediation model fit
  print(x$summary, digits = digits, signif.stars = signif.stars,
        signif.legend = FALSE, ...)
  # print indirect effect
  print(x$object, info = FALSE, digits = digits, signif.stars = signif.stars,
        signif.legend = FALSE, ...)
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
