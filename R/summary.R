# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

## summary of a mediation model fit based on a scatter matrix
## (currently doesn't do anything)
#' @method summary cov_fit_mediation
#' @export
summary.cov_fit_mediation <- function(object, ...) object

## summary of a mediation model fit based on regression
## (currently doesn't do anything)
#' @method summary reg_fit_mediation
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
#' @param type  a character string specifying how to summarize the effects
#' other than the indirect effect(s).  Possible values are \code{"boot"} (the
#' default) to compute significance tests using the normal approximation of the
#' bootstrap distribution (i.e., to assume a normal distribution of the
#' corresponding effect with the standard deviation computed from the bootstrap
#' replicates), or \code{"data"} to compute significance tests via
#' statistical theory based on the original data (e.g., t-tests if the
#' coefficients are estimated via regression).  Note that this is only relevant
#' for mediation analysis via a bootstrap test, where significance of the
#' indirect effect is always assessed via a percentile-based confidence
#' interval due to the asymmetry of its distribution.
#' @param plot  a logical indicating whether to produce a diagnostic plot of
#' robust regression weights (see \code{\link{weight_plot}()}).  This is only
#' used for mediation analysis objects fitted with the robust MM-estimator (see
#' \code{\link{test_mediation}()}).
#' @param \dots  additional arguments are currently ignored.
#'
#' @return An object of class \code{"summary_test_mediation"} with the
#' following components:
#' \item{object}{the \code{object} passed to the \code{summary} method, which
#' contains the results from testing the indirect effect(s).}
#' \item{summary}{an object containing all necessary information to summarize
#' the effects other than the indirect effect(s).}
#'
#' @author Andreas Alfons
#'
#' @references
#' Alfons, A., Ates, N.Y. and Groenen, P.J.F. (2021) A robust bootstrap test
#' for mediation analysis.  \emph{Organizational Research Methods},
#' doi: 10.1177/1094428121999096.
#'
#' @seealso \code{\link{test_mediation}}
#'
#' @examples
#' data("BSG2014")
#'
#' # set seed of the random number generator
#' set.seed(20211117)
#'
#' ## The results in Alfons et al. (2021) were obtained with an
#' ## older version of the random number generator.  To reproduce
#' ## those results, uncomment the two lines below.
#' # RNGversion("3.5.3")
#' # set.seed(20150601)
#'
#' # perform mediation analysis
#' test <- test_mediation(TeamCommitment ~ m(TaskConflict) + ValueDiversity,
#'                        data = BSG2014)
#' summary(test)
#'
#' @keywords utilities

NULL


#' @rdname summary.test_mediation
#' @method summary boot_test_mediation
#' @export

summary.boot_test_mediation <- function(object, type = c("boot", "data"),
                                        plot = TRUE, ...) {
  # get significance of effects and summary of model fit
  # component 'boot' only exists for bootstrap test, otherwise NULL
  type <- match.arg(type)
  fit <- object$fit
  if(type == "boot") summary <- get_summary(fit, boot = object$reps)
  else summary <- get_summary(fit)
  # construct return object
  result <- list(object = object, summary = summary)
  # if requested, add diagnostic plot (only ROBMED)
  have_robmed <- inherits(fit, "reg_fit_mediation") && fit$robust == "MM"
  if (have_robmed && isTRUE(plot)) {
    p <- weight_plot(fit) +
      scale_color_manual("", values = c("black", "#00BFC4")) +
      theme(legend.position = "top")
    result$plot <- p
  }
  # add class and return object
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
  result <- list(object = object, summary = summary)
  class(result) <- "summary_test_mediation"
  result
}


## internal function to compute significance of effects

get_summary <- function(object, ...) UseMethod("get_summary")

# ensure that the summary of NULL is NULL
get_summary.NULL <- function(object, ...) NULL

# for a list: get summary for each list element
# -----
# FIXME: I think the code below only works when the function is exported
# -----
# get_summary.list <- function(object, ...) lapply(object, get_summary, ...)
# -----
# dirty hack:
get_summary.list <- function(object, ...) {
  lapply(object, function(x, ...) {
    # FIXME: this throws error when MM-estimator doesn't converge
    #        The reason is that the (unconverged) object has class "lmrob.S"
    #        instead of class "lmrob".
    if (inherits(x, "lmrob")) get_summary.lmrob(x, ...)
    else if (inherits(x, "lm")) get_summary.lm(x, ...)
    else if (inherits(x, "rq")) get_summary.rq(x, ...)
    else if (inherits(x, "lmse")) get_summary.lmse(x, ...)
    else stop("not implemented yet")
  })
}
# -----

# for a linear model: return coefficient matrix, regression standard error,
# R-squared and F-test
get_summary.lm <- function(object, ...) {
  # compute the usual summary and extract coefficient matrix
  summary <- summary(object)
  coefficients <- coef(summary)
  # reformat residual standard error
  s <- list(value = summary$sigma, df = summary$df[2L])
  # reformat R-squared and adjusted R-squared
  R2 <- list(R2 = summary$r.squared, adj_R2 = summary$adj.r.squared)
  # reformat F-test
  statistic <- unname(summary$fstatistic[1L])
  df <- as.integer(summary$fstatistic[-1L])
  p_value <- pf(statistic, df[1L], df[2L], lower.tail = FALSE)
  F_test <- list(statistic = statistic, df = df, p_value = p_value)
  # return results
  result <- list(coefficients = coefficients, s = s, R2 = R2, F_test = F_test)
  class(result) <- "summary_lm"
  result
}

# for a robust linear model: return coefficient matrix, robust regression
# standard error, robust R^2 and robust F-test
get_summary.lmrob <- function(object, ...) {
  # compute the usual summary
  summary <- summary(object)
  # extract information on algorithm
  algorithm <- list(converged = summary$converged,
                    method = summary$control$method)
  # extract coefficient matrix
  coefficients <- coef(summary)
  # reformat residual standard error
  s <- list(value = summary$sigma, df = summary$df[2L])
  # reformat R-squared and adjusted R-squared
  R2 <- list(R2 = summary$r.squared, adj_R2 = summary$adj.r.squared)
  # compute robust F-test for robust fit
  F_test <- rob_F_test(object)
  # detected outliers
  robustness_weights <- summary$rweights
  threshold <- summary$control$eps.outlier
  outliers <- list(indices = which(unname(robustness_weights) < threshold),
                   weights = robustness_weights, threshold = threshold)
  # return results
  result <- list(algorithm = algorithm, coefficients = coefficients, s = s,
                 R2 = R2, F_test = F_test, outliers = outliers)
  class(result) <- "summary_lmrob"
  result
}

# for median regression: return list that contains only coefficient matrix
get_summary.rq <- function(object, ...) {
  # compute the usual summary and extract coefficient matrix
  summary <- summary(object, se = "iid")
  coefficients <- coef(summary)
  colnames(coefficients)[1] <- "Estimate"  # be consistent with other models
  # return results
  result <- list(coefficients = coefficients)
  class(result) <- "summary_rq"
  result
}

get_summary.reg_fit_mediation <- function(object, boot = NULL, ...) {
  ## initializations
  x <- object$x
  y <- object$y
  m <- object$m
  covariates <- object$covariates
  p_x <- length(x)
  p_m <- length(m)
  p_covariates <- length(covariates)
  have_robust <- is_robust(object)
  robust <- object$robust
  family <- object$family
  model <- object$model
  have_simple <- is.null(model) || model == "simple"
  have_yx <- !is.null(object$fit_yx)
  have_boot <- !is.null(boot)
  # extract number of observations
  n <- nobs(object$fit_ymx)
  ## compute bootstrap inference if available
  if (have_boot) {
    # extract coefficients
    if (p_m == 1L) coef_mx <- coef(object$fit_mx)
    else coef_mx <- unlist(unname(lapply(object$fit_mx, coef)))
    total <- object$total
    if (p_x == 1L) names(total) <- x
    coefficients <- c(coef_mx, coef(object$fit_ymx), total)
    # get list of index vectors that indicate which columns of the bootstrap
    # replicates correspond to the respective models
    index_list <- get_index_list(p_x, p_m, p_covariates, model = model,
                                 fit_yx = family != "gaussian" && have_yx)
    # extract bootstrap replicates of regression coefficients
    if (length(index_list$fit_yx) == 0) boot_coefficients <- boot$t
    else boot_coefficients <- boot$t[, -index_list$fit_yx, drop = FALSE]
    # if applicable, add bootstrap replicates for total effect
    indices_total <- ncol(boot_coefficients) + seq_len(p_x)
    boot_list <- extract_boot(object, boot = boot, index_list = index_list)
    boot_coefficients <- cbind(boot_coefficients, boot_list$total)
    # compute standard errors and z-statistics from bootstrap replicates
    means <- colMeans(boot_coefficients, na.rm = TRUE)
    se <- apply(boot_coefficients, 2, sd, na.rm = TRUE)
    z <- means / se
    # perform z-tests and combine results
    p_value <- p_value_z(z)
    coefficients <- cbind(Data = coefficients, Boot = means,
                          "Std. Error" = se, "z value" = z,
                          "Pr(>|z|)" = p_value)
  }
  ## compute summary of m ~ x + covariates
  # robust F test requires that response variable is stored in "lmrob" object
  if (robust == "MM") {
    if (p_m == 1L) {
      object$fit_mx$response <- object$data[, object$m, drop = TRUE]
    } else {
      object$fit_mx <- mapply(function(fit, m) {
        fit$response <- object$data[, m, drop = TRUE]
        fit
      }, fit = object$fit_mx, m = object$m, SIMPLIFY = FALSE)
    }
  }
  # compute summary of model
  summary_mx <- get_summary(object$fit_mx)
  # if bootstrap inference is requested, replace the usual coefficient matrix
  # with one that has bootstrap tests
  if (have_boot) {
    if (p_m == 1L) {
      keep <- index_list$fit_mx
      summary_mx$coefficients <- coefficients[keep, , drop = FALSE]
    } else {
      summary_mx <- mapply(function(summary, keep) {
        summary$coefficients <- coefficients[keep, , drop = FALSE]
        summary
      }, summary = summary_mx, keep = index_list$fit_mx, SIMPLIFY = FALSE)
    }
  }
  ## compute summary of y ~ m + x + covariates
  # robust F test requires that response variable is stored in "lmrob" object
  if (robust == "MM") object$fit_ymx$response <- object$data[, object$y]
  # compute summary of model
  summary_ymx <- get_summary(object$fit_ymx)
  # if bootstrap inference is requested, replace the usual coefficient matrix
  # with one that has bootstrap tests
  if (have_boot) {
    keep <- index_list$fit_ymx
    summary_ymx$coefficients <- coefficients[keep, , drop = FALSE]
  }
  ## extract direct effect of x on y
  direct <- summary_ymx$coefficients[x, , drop = FALSE]
  # summary of total effect of x on y
  if (have_boot) {
    total <- coefficients[indices_total, , drop = FALSE]
  } else if (have_yx) {
    # compute summary of y ~ x + covariates and extract summary of total effect
    summary_yx <- get_summary(object$fit_yx)
    total <- coef(summary_yx)[x, , drop = FALSE]
  } else {
    # standard errors and t-test not available
    total <- cbind(object$total, matrix(NA_real_, nrow = p_x, ncol = 3L))
    dimnames(total) <- dimnames(direct)
  }
  # return results
  result <- list(fit_mx = summary_mx, fit_ymx = summary_ymx, total = total,
                 direct = direct, x = x, y = y, m = m, covariates = covariates,
                 n = n, robust = robust, model = model)
  class(result) <- c("summary_reg_fit_mediation", "summary_fit_mediation")
  result
}

get_summary.cov_fit_mediation <- function(object, boot = NULL, ...) {
  # extract variable names
  x <- object$x
  y <- object$y
  m <- object$m
  # extract coefficients
  a <- object$a
  b <- object$b
  total <- object$total
  direct <- object$direct
  # the order of how coefficients are displayed changed in version 0.10.0, but
  # the computation of the standard error via the delta method (based on the
  # original data) requires the old order of coefficients
  coefficients <- c(a, b, direct, total)
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
    h_dot <- matrix(c(-a/S[x,x], a*direct/S[x,x], -a^2*direct/s_epsilon_mx-direct/S[x,x], 1, a^2, direct^2,
                      1/S[x,x], (a*b-direct)/s_epsilon_mx, -b/S[x,x]-a*(a*b-direct)/s_epsilon_mx, 0, -2*a, 2*b*direct,
                      0, -a/s_epsilon_mx, a^2/s_epsilon_mx+1/S[x,x], 0, 0, -2*direct,
                      0, -b/s_epsilon_mx, a*b/s_epsilon_mx, 0, 1, b^2,
                      0, 1/s_epsilon_mx, -a/s_epsilon_mx, 0, 0, -2*b,
                      0, 0, 0, 0, 0, 1),
                    nrow = 6, ncol = 6)
    Omega_Theta <- h_dot %*% Omega_Sigma %*% t(h_dot)
    # total effect
    Omega_Sigma_yx <- Omega_Sigma[c(1, 3, 6), c(1, 3, 6)]
    h_dot_yx <- matrix(c(-total/S[x,x], 1, 0, 1/S[x,x], 0, -1, 0, 0, 1),
                       nrow = 3, ncol = 3)
    Omega_Theta_yx <- h_dot_yx %*% Omega_Sigma_yx %*% t(h_dot_yx)
    # compute standard errors and z-statistics
    means <- NULL
    se <- sqrt(c(diag(Omega_Theta)[1:3], Omega_Theta_yx[1,1]) / n)
    z <- coefficients / se
    tn <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  } else {
    # compute means, standard errors and z-statistics from bootstrap replicates
    keep <- c(1, 2, 4, 3)
    means <- colMeans(boot$t[, keep], na.rm = TRUE)
    se <- apply(boot$t[, keep], 2, sd, na.rm = TRUE)
    z <- means / se
    tn <- c("Data", "Boot", "Std. Error", "z value", "Pr(>|z|)")
  }
  # perform z-tests and combine results
  p_value <- p_value_z(z)
  coefficients <- cbind(coefficients, means, se, z, p_value)
  dimnames(coefficients) <- list(c(x, m, x, x), tn)
  # # residual standard error as list (for compatibility with regression method)
  # s_epsilon_ymx <- S[y,y] - b^2*S[m,m] - direct^2*S[x,x] - 2*b*direct*S[m,x]
  # s <- list(value = s_epsilon_ymx)
  # return results
  result <- list(a = coefficients[1, , drop = FALSE],
                 b = coefficients[2, , drop = FALSE],
                 total = coefficients[4, , drop = FALSE],
                 direct = coefficients[3, , drop = FALSE],
                 x = x, y = y, m = m, covariates = object$covariates,
                 n = n, robust = object$robust)
  class(result) <- c("summary_cov_fit_mediation", "summary_fit_mediation")
  result
}


## compute duplication matrix according to Magnus & Neudecker (1999, p.49)
## (required for computing the Fisher information matrix of a mediation model
## fit based on a scatter matrix)
duplication_matrix <- function(p){
  D <- diag(p)
  index <- seq_len(p*(p+1)/2)
  D[lower.tri(D, diag=TRUE)] <- index
  D[upper.tri(D)] <- D[lower.tri(D)]
  outer(c(D), index, function(i, j) ifelse(i == j, 1, 0 ))
}
