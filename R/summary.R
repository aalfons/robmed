# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

## summary of a mediation model fit based on a scatter matrix
## (currently doesn't do anything)
#' @export
summary.cov_fit_mediation <- function(object, ...) object

## summary of a mediation model fit based on regression
## (currently doesn't do anything)
#' @export
summary.reg_fit_mediation <- function(object, ...) object


## summary of mediation analysis objects

#' Summary of results from (robust) mediation analysis
#'
#' Summarize results from (robust) mediation analysis for proper interpretation.
#'
#' @name summary.test_mediation
#'
#' @param object  an object inheriting from class
#' \code{"\link{test_mediation}"} containing results from (robust) mediation
#' analysis.
#' @param other  a character string specifying how to summarize the effects
#' other than the indirect effect(s).  Possible values are \code{"boot"} (the
#' default) to compute significance tests using the normal approximation of the
#' bootstrap distribution (i.e., to assume a normal distribution of the
#' corresponding effect with the standard deviation computed from the bootstrap
#' replicates), or \code{"theory"} to compute significance tests via
#' statistical theory (e.g., t-tests if the coefficients are estimated via
#' regression).  Note that this is only relevant for mediation analysis via a
#' bootstrap test, where significance of the indirect effect is always assessed
#' via a percentile-based confidence interval due to the asymmetry of its
#' distribution.
#' @param \dots  additional arguments are currently ignored.
#'
#' @return An object of class \code{"summary_test_mediation"} with the
#' following components:
#' \item{object}{the \code{object} passed to the \code{summary} method, which
#' contains the results from testing the indirect effect.}
#' \item{summary}{an object containing all necessary information to summarize
#' the effects other than the indirect effect.}
#'
#' @author Andreas Alfons
#'
#' @seealso \code{\link{test_mediation}}
#'
#' @examples
#' data("BSG2014")
#' test <- test_mediation(BSG2014,
#'                        x = "ValueDiversity",
#'                        y = "TeamCommitment",
#'                        m = "TaskConflict")
#' summary(test)
#'
#' @keywords utilities

NULL


#' @rdname summary.test_mediation
#' @method summary boot_test_mediation
#' @export

summary.boot_test_mediation <- function(object, other = c("boot", "theory"),
                                        ...) {
  # get significance of effects and summary of model fit
  # component 'boot' only exists for bootstrap test, otherwise NULL
  other <- match.arg(other)
  if(other == "boot") summary <- get_summary(object$fit, boot = object$reps)
  else summary <- get_summary(object$fit)
  # construct return object
  result <- list(object = object, summary = summary)
  class(result) <- "summary_test_mediation"
  result
}


#' @rdname summary.test_mediation
#' @method summary sobel_test_mediation
#' @export

summary.sobel_test_mediation <- function(object, ...) {
  # get significance of effects and summary of model fit
  # component 'se' only exists for Sobel test, otherwise NULL
  summary <- get_summary(object$fit)
  # construct return object
  result <- list(object=object, summary=summary)
  class(result) <- "summary_test_mediation"
  result
}


## internal function to compute significance of effects

get_summary <- function(object, ...) UseMethod("get_summary")

get_summary.cov_fit_mediation <- function(object, boot = NULL, ...) {
  # extract variable names
  x <- object$x
  y <- object$y
  m <- object$m
  # extract coefficients
  a <- object$a
  b <- object$b
  c <- object$c
  c_prime <- object$c_prime
  coefficients <- c(a, b, c, c_prime)
  # extract covariance matrix
  S <- object$cov$cov[c(x, m, y), c(x, m, y)]
  # extract number of observations
  n <- nobs(object$cov)
  # compute standard errors and z-statistics
  if(is.null(boot)) {
    # compute the inverse of the Fisher information matrix for the unique
    # elements of the covariance matrix
    D <- duplication_matrix(3)
    inv_S <- solve(S)
    W <- t(D) %*% kronecker(inv_S, inv_S) %*% D / 2
    Omega_Sigma <- solve(W)
    # parameters in mediation model
    # apply the delta method to obtain the inverse of the Fisher information
    # matrix for the coefficients in the mediation model, see Zu & Yuan (2010)
    s_epsilon_mx <- S[m,m] - a^2 * S[x,x]
    h_dot <- matrix(c(-a/S[x,x], a*c/S[x,x], -a^2*c/s_epsilon_mx-c/S[x,x], 1, a^2, c^2,
                     1/S[x,x], (a*b-c)/s_epsilon_mx, -b/S[x,x]-a*(a*b-c)/s_epsilon_mx, 0, -2*a, 2*b*c,
                     0, -a/s_epsilon_mx, a^2/s_epsilon_mx+1/S[x,x], 0, 0, -2*c,
                     0, -b/s_epsilon_mx, a*b/s_epsilon_mx, 0, 1, b^2,
                     0, 1/s_epsilon_mx, -a/s_epsilon_mx, 0, 0, -2*b,
                     0, 0, 0, 0, 0, 1),
                   nrow=6, ncol=6)
    Omega_Theta <- h_dot %*% Omega_Sigma %*% t(h_dot)
    # total effect
    Omega_Sigma_yx <- Omega_Sigma[c(1, 3, 6), c(1, 3, 6)]
    h_dot_yx <- matrix(c(-c_prime/S[x,x], 1, 0, 1/S[x,x], 0, -1, 0, 0, 1),
                     nrow=3, ncol=3)
    Omega_Theta_yx <- h_dot_yx %*% Omega_Sigma_yx %*% t(h_dot_yx)
    # compute standard errors and z-statistics
    means <- NULL
    se <- sqrt(c(diag(Omega_Theta)[1:3], Omega_Theta_yx[1,1]) / n)
    z <- coefficients / se
    tn <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  } else {
    # compute means, standard errors and z-statistics from bootstrap replicates
    means <- colMeans(boot$t[, -1], na.rm=TRUE)
    se <- apply(boot$t[, -1], 2, sd, na.rm=TRUE)
    z <- means / se
    tn <- c("Data", "Boot", "Std. Error", "z value", "Pr(>|z|)")
  }
  # perform z-tests and combine results
  p_value <- p_value_z(z)
  coefficients <- cbind(coefficients, means, se, z, p_value)
  dimnames(coefficients) <- list(c(x, m, x, x), tn)
  # residual standard error as list (for compatibility with regression method)
  s_epsilon_ymx <- S[y,y] - b^2*S[m,m] - c^2*S[x,x] - 2*b*c*S[m,x]
  s <- list(value=s_epsilon_ymx)
  # return results
  result <- list(a=coefficients[1, , drop=FALSE],
                 b=coefficients[2, , drop=FALSE],
                 c=coefficients[3, , drop=FALSE],
                 c_prime=coefficients[4, , drop=FALSE],
                 robust=object$robust, s=s, n=n,
                 variables=c(x, y, m))
  class(result) <- "summary_fit_mediation"
  result
}

get_summary.reg_fit_mediation <- function(object, boot = NULL, ...) {
  # initializations
  x <- object$x
  y <- object$y
  m <- object$m
  p_m <- length(m)
  covariates <- object$covariates
  robust <- object$robust
  median <- object$median
  have_boot <- !is.null(boot)
  # compute summary of y ~ m + x + covariates
  if (median) summary_ymx <- summary(object$fit_ymx, se = "iid")
  else summary_ymx <- summary(object$fit_ymx)
  # extract number of observations
  n <- nobs(object$fit_ymx)
  # perform tests for significance of effects
  if (have_boot) {
    # extract coefficients and add coefficients of covariates
    coefficients <- c(coefficients(object),
                      coef(object$fit_ymx)[-seq_len(p_m + 2L)])
    # compute standard errors and z-statistics from bootstrap replicates
    remove <- if(p_m == 1L) 1L else seq_len(1L + p_m)
    means <- colMeans(boot$t[, -remove], na.rm = TRUE)
    se <- apply(boot$t[, -remove], 2, sd, na.rm = TRUE)
    z <- means / se
    # perform z-tests and combine results
    p_value <- p_value_z(z)
    coefficients <- cbind(coefficients, means, se, z, p_value)
    rn <- c(if(p_m == 1L) x else paste(m, x, sep = "~"), m, x, x, covariates)
    tn <- c("Data", "Boot", "Std. Error", "z value", "Pr(>|z|)")
    dimnames(coefficients) <- list(rn, tn)
    # split up effect summaries
    a <- coefficients[seq_len(p_m), , drop = FALSE]
    b <- coefficients[p_m + seq_len(p_m), , drop = FALSE]
    c <- coefficients[2L * p_m + 1L, , drop = FALSE]
    c_prime <- coefficients[2L * p_m + 2L, , drop = FALSE]
  } else {
    # compute summaries of regression models and extract t-tests for coefficients
    if (p_m == 1L) {
      if (median) summary_mx <- summary(object$fit_mx, se = "iid")
      else summary_mx <- summary(object$fit_mx)
      a <- summary_mx$coefficients[2L, , drop = FALSE]
      b <- summary_ymx$coefficients[2L, , drop = FALSE]
      c <- summary_ymx$coefficients[3L, , drop = FALSE]
    } else {
      if (median) summary_mx <- lapply(object$fit_mx, summary, se = "iid")
      else summary_mx <- lapply(object$fit_mx, summary)
      a <- lapply(summary_mx, function(s) s$coefficients[2L, , drop = FALSE])
      a <- do.call(rbind, a)
      rownames(a) <- paste(m, x, sep = "~")
      b <- summary_ymx$coefficients[1L + seq_len(p_m), , drop = FALSE]
      c <- summary_ymx$coefficients[p_m + 2L, , drop = FALSE]
    }
    if (robust) {
      # standard errors and t-test not available
      c_prime <- matrix(c(object$c_prime, rep.int(NA_real_, 3L)), nrow = 1L)
      dimnames(c_prime) <- dimnames(c)
    } else {
      summary_yx <- summary(object$fit_yx)
      c_prime <- summary_yx$coefficients[2L, , drop = FALSE]
    }
  }
  # initialize return object
  result <- list(a = a, b = b, c = c, c_prime = c_prime)
  # add partial effects of control variables if they exist
  if (length(covariates) > 0L) {
    if (have_boot) {
      remove <- seq_len(2L * p_m + 2L)
      result$covariate_effects <- coefficients[-remove, , drop = FALSE]
    } else {
      remove <- seq_len(p_m + 2L)
      result$covariate_effects <- summary_ymx$coefficients[-remove, , drop=FALSE]
    }
  }
  # add robustness information
  result$robust <- robust
  result$median <- median
  if (!median) {
    # add residual standard error
    result$s <- list(value = summary_ymx$sigma, df = summary_ymx$df[2L])
    # add (robust) R-squared
    result$R2 <- list(R2 = summary_ymx$r.squared,
                      adj_R2 = summary_ymx$adj.r.squared)
    # add (robust) F-test
    if (robust) {
      # compute robust F-test for robust fit
      result$F_test <- rob_F_test(object)
    } else {
      # add F-test for nonrobust fit
      statistic <- unname(summary_ymx$fstatistic[1L])
      df <- as.integer(summary_ymx$fstatistic[-1L])
      p_value <- pf(statistic, df[1L], df[2L], lower.tail = FALSE)
      result$F_test <- list(statistic = statistic, df = df,
                            p_value = p_value)
    }
  }
  # add number of observations and variable names
  result <- c(result, list(n = n, x = x, y = y, m = m, covariates = covariates))
  ## add class and return results
  class(result) <- "summary_fit_mediation"
  result
}


## compute duplication matrix according to Magnus & Neudecker (1999, p.49)
## (required for computing the Fischer information matrix of a mediation model
## fit based on a scatter matrix)
duplication_matrix <- function(p){
  D <- diag(p)
  index <- seq(p*(p+1)/2)
  D[lower.tri(D, diag=TRUE)] <- index
  D[upper.tri(D)] <- D[lower.tri(D)]
  outer(c(D), index, function(i, j) ifelse(i == j, 1, 0 ))
}
