# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

# information to display for different objects in first line of output
#' @noRd
print_info <- function(object, ...) UseMethod("print_info")

# information to display for regression fits
#' @noRd
print_info.reg_fit_mediation <- function(object, ...) {
  prefix <- if (is_robust(object)) "Robust mediation" else "Mediation"
  if (object$robust == "median") postfix <- "via median regression"
  else {
    errors <- switch(object$family, student = " with t errors",
                     skewnormal = " with skew-normal errors",
                     skewt = " with skew-t errors",
                     select = " with selection of error distribution")
    postfix <- paste0("via regression", errors)
  }
  cat(sprintf("%s model fit %s\n", prefix, postfix))
}

# information to display for covariance matrix fits
#' @noRd
print_info.cov_fit_mediation <- function(object, ...) {
  prefix <- if (is_robust(object)) "Robust mediation" else "Mediation"
  cat(sprintf("%s model fit via covariance matrix\n", prefix))
}

# information to display for bootstrap tests
#' @noRd
print_info.boot_test_mediation <- function(object, ...) {
  # type of test and model fit
  fit <- object$fit
  if (inherits(object$fit, "reg_fit_mediation")) {
    prefix <- if (fit$robust == "MM") "Robust bootstrap" else "Bootstrap"
    if (fit$robust == "median") postfix <- " via median regression"
    else {
      errors <- switch(fit$family, student = " with t errors",
                       skewnormal = " with skew-normal errors",
                       skewt = " with skew-t errors",
                       select = "\nError distribution is selected via BIC")
      postfix <- paste0(" via regression", errors)
    }
  } else if (inherits(fit, "cov_fit_mediation")) {
    prefix <- if (fit$robust) "Robust bootstrap" else "Bootstrap"
    postfix <- " via covariance matrix"
  } else {
    prefix <- "Bootstrap"
    postfix <- ""
  }
  # use plural for multiple mediators
  plural <- if (is.null(fit$model) || fit$model == "simple") "" else "s"
  # return message
  cat(sprintf("%s test%s for indirect effect%s%s\n",
              prefix, plural, plural, postfix))
}

# information to display for sobel tests
#' @noRd
print_info.sobel_test_mediation <- function(object, ...) {
  # type of test and model fit
  fit <- object$fit
  if (inherits(object$fit, "reg_fit_mediation")) {
    prefix <- if (fit$robust == "MM") "Robust normal" else "Normal"
    if (fit$robust == "median") postfix <- " via median regression"
    else {
      errors <- switch(fit$family, student = " with t errors",
                       skewnormal = " with skew-normal errors",
                       skewt = " with skew-t errors")
      postfix <- paste0(" via regression", errors)
    }
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
#' @noRd
print_info.summary_fit_mediation <- function(object, ...) {
  # initializations
  p_x <- length(object$x)
  p_m <- length(object$m)
  have_simple <- is.null(object$model) || object$model == "simple"
  have_covariates <- length(object$covariates) > 0L
  # in case of multiple mediators, print information on type of mediation model
  if (p_m > 1L) {
    prefix <- switch(object$model, parallel = "Parallel", serial = "Serial")
    cat(prefix, "multiple mediator model\n\n")
  }
  # print information on variables
  if (have_simple) {
    cat(sprintf("x = %s\n", object$x))
    cat(sprintf("y = %s\n", object$y))
    cat(sprintf("m = %s\n", object$m))
  } else {
    width <- max(nchar(p_x), nchar(p_m)) + 1L
    if (p_x == 1L) cat(sprintf(paste0("%-", width, "s = %s\n"), "x", object$x))
    else {
      x_labels <- paste0("x", seq_len(p_x))
      cat(sprintf(paste0("%-", width, "s = %s\n"), x_labels, object$x), sep = "")
    }
    cat(sprintf(paste0("%-", width, "s = %s\n"), "y", object$y))
    if (p_m == 1L) cat(sprintf(paste0("%-", width, "s = %s\n"), "m", object$m))
    else {
      m_labels <- paste0("m", seq_len(p_m))
      cat(sprintf(paste0("%-", width, "s = %s\n"), m_labels, object$m), sep = "")
    }
  }
  if (have_covariates) {
    cat("\nCovariates:\n")
    print(object$covariates, quote = FALSE)
  }
  # print sample size
  cat(sprintf("\nSample size: %d\n", object$n))
}


#' @export
print.fit_mediation <- function(x, info = TRUE, ...) {
  # print information on type of model fit
  if (isTRUE(info)) print_info(x, ...)
  # print estimated effects
  if (is.null(x$model) || x$model == "simple") {
    cat("\nEffects:\n")
    print(coef(x), ...)
  } else {
    # some initializations
    p_x <- length(x$x)
    p_m <- length(x$m)
    # print effects
    cat("\na paths:\n")
    print(x[["a"]], ...)
    cat(sprintf("\nb path%s:\n", if (p_m == 1L) "" else "s"))
    print(x[["b"]], ...)
    if (!is.null(x[["d"]])) {
      cat(sprintf("\nd path%s:\n", if (p_m > 2L) "s" else ""))
      print(x[["d"]], ...)
    }
    cat(sprintf("\nTotal effect%s:\n", if (p_x == 1L) "" else "s"))
    print(x[["total"]], ...)
    cat(sprintf("\nDirect effect%s:\n", if (p_x == 1L) "" else "s"))
    print(x[["direct"]], ...)
    cat("\nIndirect effects:\n")
    print(x[["indirect"]], ...)
  }
  # print information on how contrasts are computed
  # (if there are no contrasts, nothing is printed)
  print_contrast_info(x)
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
  # initializations
  names_x <- x$fit$x
  p_x <- length(names_x)
  names_m <- x$fit$m
  p_m <- length(names_m)
  model <- x$fit$model
  have_simple <- is.null(model) || model == "simple"
  # print information on type of model fit
  if (isTRUE(info)) print_info(x, ...)
  # print indirect effects
  plural <- if (have_simple) "" else "s"
  cat(sprintf("\nIndirect effect%s of x on y:\n", plural))
  # extract indirect effect and confidence interval
  indirect <- cbind(Data = x$fit$indirect, Boot = x$indirect,
                    if (have_simple) t(x$ci) else x$ci)
  if (have_simple) rownames(indirect) <- names_m
  # replace complicated row names with shorter labels
  if (p_m > 1L && (p_x > 1L || model == "serial")) {
    labels <- get_indirect_labels(p_m, model = model)
    if (p_x == 1L) rownames(indirect)[1L + seq_along(labels)] <- labels
    else {
      replace <- grep("->", rownames(indirect))
      replacement_names <- unlist(lapply(names_x, paste, labels, sep = "_"),
                                  use.names = FALSE)
      rownames(indirect)[replace] <- replacement_names
    }
  } else labels <- NULL
  # combine and print
  print(indirect, digits = digits, ...)
  # print information on indirect effect paths for models with multiple
  # independent variables and multiple mediators, or serial multiple mediators
  # (if we have a different model, nothing is printed)
  print_indirect_info(x$fit, labels = labels)
  # print information on how contrasts are computed
  # (if there are no contrasts, nothing is printed)
  print_contrast_info(x$fit, labels = labels)
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
  ab <- cbind(x$fit$indirect, x$se, x$statistic, x$p_value)
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
  # check if algorithm converged and print information accordingly
  if (x$algorithm$converged) cat("Coefficients:\n")
  else if (x$s$value == 0) cat("Exact fit detected\n\nCoefficients:\n")
  else {
    cat("Algorithm did not converge\n\n")
    if (x$algorithm$method == "S") {
      cat("Coefficients of the *initial* S-estimator:\n")
    } else {
      cat(sprintf("Coefficients of the %s-estimator:\n", x$algorithm$method))
    }
  }
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
  # print information on robustness weights
  cat("\nRobustness weights:\n")
  # -----
  # indices <- x$outliers$indices
  # n_outliers <- length(indices)
  # if (n_outliers == 0) {
  #   # print that there is no clear outlier and minimum robustness weight
  #   cat("No potential outliers with weight < ",
  #       formatC(x$outliers$threshold, digits = max(2, digits-3), width = 1),
  #       " detected. The minimum weight is ",
  #       formatC(min(x$outliers$weights), digits = max(2, digits-3), width = 1),
  #       ".\n", sep = "")
  # } else if (n_outliers == 1) {
  #   # robustness weight of potential outlier
  #   outlier_weight <- x$outliers$weights[indices]
  #   # print information on potential outliers
  #   cat("Observation ", indices, " is a potential outlier with weight ",
  #       formatC(outlier_weight, digits = max(2, digits-3), width = 1),
  #       # " < ",
  #       # formatC(x$outliers$threshold, digits = max(2, digits-3), width = 1),
  #       "\n", sep = "")
  # } else {
  #   # largest weight still below the outlier threshold
  #   max_outlier_weight <- max(x$outliers$weights[indices])
  #   # print information on potential outliers
  #   cat(n_outliers, " observations are potential outliers with weight <= ",
  #       formatC(max_outlier_weight, digits = max(2, digits-3), width = 1),
  #       ":\n", sep = "")
  #   print(indices)
  # }
  # -----
  # this is a bit complicated, but allows to print the same information
  # differently in SPSS extension bundle
  outlier_info <- get_outlier_info(x$outliers, digits = digits)
  cat(outlier_info$msg)
  if (!is.null(outlier_info$indices)) print(outlier_info$indices)
  # -----
  # return object invisibly
  invisible(x)
}

# print information on clear outliers and robustness weights
get_outlier_info <- function(outliers, digits = max(3, getOption("digits")-3)) {
  # extract indices and number of clear outliers
  indices <- outliers$indices
  n_outliers <- length(indices)
  # prepare information to return
  if (n_outliers == 0) {
    # information that there is no clear outlier and minimum robustness weight
    msg <- paste("No potential outliers with weight < ",
          formatC(outliers$threshold, digits = max(2, digits-3), width = 1),
          " detected.\nThe minimum weight is ",
          formatC(min(outliers$weights), digits = max(2, digits-3), width = 1),
          ".\n", sep = "")
    indices_to_print <- NULL
  } else if (n_outliers == 1) {
    # robustness weight of potential outlier
    outlier_weight <- outliers$weights[indices]
    # information on potential outlier
    msg <- paste("Observation ", indices,
                 " is a potential outlier with weight ",
                 formatC(outlier_weight, digits = max(2, digits-3), width = 1),
                 # " < ",
                 # formatC(x$outliers$threshold, digits = max(2, digits-3),
                 #         width = 1),
                 "\n", sep = "")
    indices_to_print <- NULL
  } else {
    # largest weight still below the outlier threshold
    max_outlier_weight <- max(outliers$weights[indices])
    # information on potential outliers
    msg <- paste(n_outliers,
                 " observations are potential outliers with weight <= ",
                 formatC(max_outlier_weight, digits = max(2, digits-3),
                         width = 1),
                 ":\n", sep = "")
    indices_to_print <- indices
  }
  # return message and indices to print
  list(msg = msg, indices = indices_to_print)
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
  p_x <- length(x$x)
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
  plural <- if (p_x == 1L) "" else "s"
  cat(sprintf("---\nTotal effect%s of x on y:\n", plural))
  printCoefmat(x$total, digits = digits, signif.stars = signif.stars,
               signif.legend = FALSE, ...)
  cat(sprintf("\nDirect effect%s of x on y:\n", plural))
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
  cat("\nPartial effect of m on y:\n")
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
  # create plot if requested
  p <- x$plot
  if (!is.null(p)) print(p)
  # return object invisibly
  invisible(x)
}


## internal function to print information on labels for indirect effects for
## serial multiple mediator models
# This now expects an object of class "fit_mediation".  In future versions,
# this could also be turned into a generic function if necessary.
print_indirect_info <- function(object, labels = NULL, ...) {
  # initializations
  if (is.null(labels)) return()  # nothing to print if there are no labels
  x <- object$x
  p_x <- length(x)
  m <- object$m
  p_m <- length(m)
  y <- object$y
  model <- object$model  # only implemented for regression fit
  # label for independent variable
  x_label <- if (p_x == 1L) x else "x"
  # get information on indirect effect paths
  if (model == "serial") {
    # currently only implemented for two or three hypothesized mediators
    if (p_m == 2L) {
      # format strings of relevant variable names to be of equal length
      format_m <- format(m)
      format_ym <- format(c(y, m[2L]))
      # generate strings describing the paths
      paths <- c(paste(x_label, format_m, format_ym[1L], sep = " -> "),
                 paste(x_label, format_m[1L], format_ym[2L], y, sep = " -> "))
    } else {
      # format strings of relevant variable names to be of equal length
      format_m <- format(m)
      format_ym23 <- format(c(y, m[-1L]))
      format_ym3 <- format(c(y, m[3L]))
      # generate strings describing the paths
      paths <- c(paste(x_label, format_m, format_ym23[1L],
                       sep = " -> "),
                 paste(x_label, format_m[c(1L, 1L, 2L)],
                       format_ym23[c(2L, 3L, 3L)],
                       format_ym3[1L], sep = " -> "),
                 paste(x_label, format_m[1L], format_ym23[2L],
                       format_ym3[2L], y, sep = " -> "))
    }
  } else {
    # format strings of relevant variable names to be of equal length
    format_m <- format(m)
    # generate strings describing the paths
    paths <- paste(x_label, format_m, y, sep = " -> ")
  }
  # create data frame with information on indirect effect paths
  indirect_info <- data.frame(Label = labels, Path = paths)
  # print information on indirect effect paths
  cat("\nIndirect effect paths:\n")
  print(indirect_info, right = FALSE, row.names = FALSE)
}

## internal function to print information on contrast definitions
# This now expects an object of class "fit_mediation".  In future versions,
# this could also be turned into a generic function if necessary.
print_contrast_info <- function(object, labels = NULL, ...) {
  # initializations
  contrast <- object$contrast              # only implemented for regression fit
  have_contrast <- is.character(contrast)  # but this always works
  # if applicable, print indirect effect contrast definitions
  if (have_contrast) {
    # if no labels are supplied, get names of indirect effects
    if (is.null(labels)) {
      labels <- get_indirect_names(object$x, object$m,
                                   model = object$model)
    }
    nr_indirect <- length(labels)
    # print information on contrast definitions
    plural <- if (nr_indirect > 2L) "s" else ""
    cat(sprintf("\nIndirect effect contrast definition%s:\n", plural))
    contrast_info <- get_contrast_info(labels, type = contrast)
    print(contrast_info, right = FALSE, row.names = FALSE)
  }
}

## internal function to print legend for significance stars
print_legend <- function() {
  cat("---\nSignif. codes:  0", sQuote("***"), "0.001", sQuote("**"),
      "0.01", sQuote("*"), "0.05", sQuote("."), "0.1", sQuote(" "), "1\n")
}
