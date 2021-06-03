# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

#' (Robust) mediation analysis
#'
#' Perform (robust) mediation analysis via a (fast and robust) bootstrap test
#' or Sobel's test.
#'
#' With \code{method = "regression"}, and \code{robust = TRUE} or
#' \code{robust = "MM"}, the tests are based on robust regressions with the
#' MM-estimator from \code{\link[robustbase]{lmrob}()}.  The bootstrap test is
#' thereby performed via the fast and robust bootstrap.  This is the default
#' behavior.
#'
#' Note that the MM-estimator of regression implemented in
#' \code{\link[robustbase]{lmrob}()} can be seen as weighted least squares
#' estimator, where the weights are dependent on how much an observation is
#' deviating from the rest.  The trick for the fast and robust bootstrap is
#' that on each bootstrap sample, first a weighted least squares estimator
#' is computed (using those robustness weights from the original sample)
#' followed by a linear correction of the coefficients.  The purpose of this
#' correction is to account for the additional uncertainty of obtaining the
#' robustness weights.
#'
#' With \code{method = "regression"} and \code{robust = "median"}, the tests
#' are based on median regressions with \code{\link[quantreg]{rq}()}.  Note
#' that the bootstrap test is performed via the standard bootstrap, as the
#' fast and robust bootstrap is not applicable.  Unlike the robust regressions
#' described above, median regressions are not robust against outliers in
#' the explanatory variables, and the standard bootstrap can suffer from
#' oversampling of outliers in the bootstrap samples.
#'
#' With \code{method = "covariance"} and \code{robust = TRUE}, the
#' tests are based on a Huber M-estimator of location and scatter.  For the
#' bootstrap test, the M-estimates are used to first clean the data via a
#' transformation.  Then the standard bootstrap is performed with the cleaned
#' data.  Note that this covariance-based approach is less robust than the
#' approach based on robust regressions described above.  Furthermore, the
#' bootstrap does not account for the variability from cleaning the data.
#'
#' \code{robmed()} is a wrapper function for performing robust mediation
#' analysis via regressions and the fast and robust bootstrap.
#'
#' @aliases print.boot_test_mediation print.sobel_test_mediation
#'
#' @param object  the first argument will determine the method of the generic
#' function to be dispatched.  For the default method, this should be a data
#' frame containing the variables.  There is also a method for a mediation
#' model fit as returned by \code{\link{fit_mediation}()}.
#' @param formula  	an object of class "formula" (or one that can be coerced to
#' that class): a symbolic description of the model to be fitted.  Hypothesized
#' mediator variables should be wrapped in a call to \code{\link{m}()} (see
#' examples), and any optional control variables should be wrapped in a call to
#' \code{\link{covariates}()}.
#' @param data  for the \code{formula} method, a data frame containing the
#' variables.
#' @param x  a character, integer or logical vector specifying the columns of
#' \code{object} containing the independent variables.
#' @param y  a character string, an integer or a logical vector specifying the
#' column of \code{object} containing the dependent variable.
#' @param m  a character, integer or logical vector specifying the columns of
#' \code{object} containing the hypothesized mediator variables.
#' @param covariates  optional; a character, integer or logical vector
#' specifying the columns of \code{object} containing additional covariates to
#' be used as control variables.
#' @param test  a character string specifying the test to be performed for
#' the indirect effects.  Possible values are \code{"boot"} (the default) for
#' the bootstrap, or \code{"sobel"} for Sobel's test.  Currently, Sobel's test
#' is not implemented for models with multiple indirect effects.
#' @param alternative  a character string specifying the alternative hypothesis
#' in the test for the indirect effects.  Possible values are \code{"twosided"}
#' (the default), \code{"less"} or \code{"greater"}.
#' @param R  an integer giving the number of bootstrap replicates.  The default
#' is to use 5000 bootstrap replicates.
#' @param level  numeric; the confidence level of the confidence interval in
#' the bootstrap test.  The default is to compute a 95\% confidence interval.
#' @param type  a character string specifying the type of confidence interval
#' to be computed in the bootstrap test.  Possible values are \code{"bca"} (the
#' default) for the bias-corrected and accelerated bootstrap, or \code{"perc"}
#' for the percentile bootstrap.
#' @param order  a character string specifying the order of approximation of
#' the standard error in Sobel's test.  Possible values are \code{"first"}
#' (the default) for a first-order approximation, and \code{"second"} for a
#' second-order approximation.
#' @param method  a character string specifying the method of estimation for
#' the mediation model.  Possible values are \code{"regression"} (the default)
#' to estimate the effects via regressions, or \code{"covariance"} to estimate
#' the effects via the covariance matrix.  Note that the effects are
#' always estimated via regressions if more than one independent variable or
#' hypothesized mediator is specified, or if control variables are supplied.
#' @param robust  a logical indicating whether to perform a robust test
#' (defaults to \code{TRUE}).  For estimation via regressions
#' (\code{method = "regression"}), this can also be a character string, with
#' \code{"MM"} specifying the MM-estimator of regression, and \code{"median"}
#' specifying median regression.
#' @param family  a character string specifying the error distribution to be
#' used in maximum likelihood estimation of regression models.  Possible values
#' are \code{"gaussian"} for a normal distribution (the default),
#' \code{skewnormal} for a skew-normal distribution, \code{"student"} for
#' Student's t distribution, \code{"skewt"} for a skew-t distribution, or
#' \code{"select"} to select among these four distributions via BIC (see
#' \code{\link{fit_mediation}()} for details).  This is only relevant if
#' \code{method = "regression"} and \code{robust = FALSE}.
#' @param contrast  a logical indicating whether to compute pairwise contrasts
#' of the indirect effects (defaults to \code{FALSE}).  This can also be a
#' character string, with \code{"estimates"} for computing the pairwise
#' differences of the indirect effects (such that it is tested whether two
#' indirect effects are equal), and \code{"absolute"} for computing the
#' pairwise differences of the absolute values of the indirect effects
#' (such that it is tested whether two indirect effects are equal in
#' magnitude).  This is only relevant for models with multiple indirect
#' effects, which are currently only implemented for estimation via
#' regressions (\code{method = "regression"}).
#' @param fit_yx  a logical indicating whether to fit the regression model
#' \code{y ~ x + covariates} to estimate the total effect (the default is
#' \code{TRUE}).  This is only relevant if \code{method = "regression"} and
#' \code{robust = FALSE}.
#' @param control  a list of tuning parameters for the corresponding robust
#' method.  For robust regression (\code{method = "regression"}, and
#' \code{robust = TRUE} or \code{robust = "MM"}), a list of tuning
#' parameters for \code{\link[robustbase]{lmrob}()} as generated by
#' \code{\link{reg_control}()}.  For winsorized covariance matrix estimation
#' (\code{method = "covariance"} and \code{robust = TRUE}), a list of tuning
#' parameters for \code{\link{cov_Huber}()} as generated by
#' \code{\link{cov_control}()}.  No tuning parameters are necessary for median
#' regression (\code{method = "regression"} and \code{robust = "median"}).
#' @param \dots  additional arguments to be passed down.  For the bootstrap
#' tests, those can be used to specify arguments of \code{\link[boot]{boot}()},
#' for example for parallel computing.
#'
#' @return An object inheriting from class \code{"test_mediation"} (class
#' \code{"boot_test_mediation"} if \code{test = "boot"} or
#' \code{"sobel_test_mediation"} if \code{test = "sobel"}) with the
#' following components:
#' \item{a}{a numeric vector containing the bootstrap point estimates of the
#' effects of the independent variables on the proposed mediator variables
#' (only \code{"boot_test_mediation"}).}
#' \item{b}{a numeric vector containing the bootstrap point estimates of the
#' direct effects of the proposed mediator variables on the dependent variable
#' (only \code{"boot_test_mediation"}).}
#' \item{direct}{a numeric vector containing the bootstrap point estimates of
#' the direct effects of the independent variables on the dependent variable
#' (only \code{"boot_test_mediation"}).}
#' \item{total}{a numeric vector containing the bootstrap point estimates of
#' the total effects of the independent variables on the dependent variable
#' (only \code{"boot_test_mediation"}).}
#' \item{ab}{a numeric vector containing the bootstrap point estimates of the
#' indirect effects (only \code{"boot_test_mediation"}).}
#' \item{ci}{a numeric vector of length two or a matrix of two columns
#' containing the bootstrap confidence intervals for the indirect effects
#' (only \code{"boot_test_mediation"}).}
#' \item{reps}{an object of class \code{"\link[boot]{boot}"} containing
#' the bootstrap replicates of the effects (only \code{"boot_test_mediation"}).}
#' \item{se}{numeric; the standard error of the indirect effect according
#' to Sobel's formula (only \code{"sobel_test_mediation"}).}
#' \item{statistic}{numeric; the test statistic for Sobel's test (only
#' \code{"sobel_test_mediation"}).}
#' \item{p_value}{numeric; the p-value from Sobel's test (only
#' \code{"sobel_test_mediation"}).}
#' \item{alternative}{a character string specifying the alternative
#' hypothesis in the test for the indirect effects.}
#' \item{R}{an integer giving the number of bootstrap replicates (only
#' \code{"boot_test_mediation"}).}
#' \item{level}{numeric; the confidence level of the bootstrap confidence
#' interval (only \code{"boot_test_mediation"}).}
#' \item{type}{a character string specifying the type of bootstrap
#' confidence interval (only \code{"boot_test_mediation"}).}
#' \item{fit}{an object inheriting from class
#' \code{"\link{fit_mediation}"} containing the estimation results of the
#' mediation model on the original data.}
#'
#' @note For the fast and robust bootstrap, the simpler correction of
#' Salibian-Barrera & Van Aelst (2008) is used rather than the originally
#' proposed correction of Salibian-Barrera & Zamar (2002).
#'
#' The formula interface is still experimental and may change in future
#' versions.
#'
#' @author Andreas Alfons
#'
#' @references
#' Alfons, A., Ates, N.Y. and Groenen, P.J.F. (2021) A robust bootstrap test
#' for mediation analysis.  \emph{Organizational Research Methods},
#' \doi{10.1177/1094428121999096}.
#'
#' Azzalini, A. and Arellano-Valle, R. B. (2013) Maximum penalized likelihood
#' estimation for skew-normal and skew-t distributions.  \emph{Journal of
#' Statistical Planning and Inference}, \bold{143}(2), 419--433.
#'
#' Preacher, K.J. and Hayes, A.F. (2004) SPSS and SAS procedures for estimating
#' indirect effects in simple mediation models. \emph{Behavior Research Methods,
#' Instruments, & Computers}, \bold{36}(4), 717--731.
#'
#' Preacher, K.J. and Hayes, A.F. (2008) Asymptotic and resampling strategies
#' for assessing and comparing indirect effects in multiple mediator models.
#' \emph{Behavior Research Methods}, \bold{40}(3), 879--891.
#'
#' Salibian-Barrera, M. and Van Aelst, S. (2008) Robust model selection using
#' fast and robust bootstrap. \emph{Computational Statistics & Data Analysis},
#' \bold{52}(12), 5121--5135
#'
#' Salibian-Barrera, M. and Zamar, R. (2002) Bootstrapping robust estimates of
#' regression. \emph{The Annals of Statistics}, \bold{30}(2), 556--582.
#'
#' Sobel, M.E. (1982) Asymptotic confidence intervals for indirect effects in
#' structural equation models. \emph{Sociological Methodology}, \bold{13},
#' 290--312.
#'
#' Yuan, Y. and MacKinnon, D.P. (2014) Robust mediation analysis based on
#' median regression. \emph{Psychological Methods}, \bold{19}(1),
#' 1--20.
#'
#' Zu, J. and Yuan, K.-H. (2010) Local influence and robust procedures for
#' mediation analysis. \emph{Multivariate Behavioral Research}, \bold{45}(1),
#' 1--44.
#'
#' @seealso \code{\link{fit_mediation}()}
#'
#' \code{\link[=coef.test_mediation]{coef}()},
#' \code{\link[=confint.test_mediation]{confint}()} and
#' \code{\link[=plot-methods]{plot}()} methods, \code{\link{p_value}()}
#'
#' \code{\link[boot]{boot}()}, \code{\link[robustbase]{lmrob}()},
#' \code{\link[stats]{lm}()}, \code{\link{cov_Huber}()}, \code{\link{cov_ML}()}
#'
#' @examples
#' data("BSG2014")
#'
#' ## The results in Alfons et al. (2021) were obtained with an
#' ## older version of the random number generator.  To reproduce
#' ## those results, uncomment the call to RNGversion() below.
#'
#' # RNGversion("3.5.3")
#' seed <- 20150601
#'
#' # formula interface
#' set.seed(seed)
#' test1 <- test_mediation(TeamCommitment ~ m(TaskConflict) + ValueDiversity,
#'                         data = BSG2014)
#' summary(test1)
#'
#' # default method
#' set.seed(seed)
#' test2 <- test_mediation(BSG2014,
#'                         x = "ValueDiversity",
#'                         y = "TeamCommitment",
#'                         m = "TaskConflict")
#' summary(test2)
#'
#' @keywords multivariate
#'
#' @import boot
#' @import robustbase
#' @importFrom quantreg rq.fit
#' @export

test_mediation <- function(object, ...) UseMethod("test_mediation")


#' @rdname test_mediation
#' @method test_mediation formula
#' @export

test_mediation.formula <- function(formula, data, test = c("boot", "sobel"),
                                   alternative = c("twosided", "less", "greater"),
                                   R = 5000, level = 0.95,
                                   type = c("bca", "perc"),
                                   order = c("first", "second"),
                                   method = c("regression", "covariance"),
                                   robust = TRUE, family = "gaussian",
                                   contrast = FALSE, fit_yx = TRUE,
                                   control = NULL, ...) {
  ## fit mediation model
  if (missing(data)) {
    fit <- fit_mediation(formula, method = method, robust = robust,
                         family = family, contrast = contrast,
                         fit_yx = fit_yx, control = control)
  } else {
    fit <- fit_mediation(formula, data = data, method = method,
                         robust = robust, family = family,
                         contrast = contrast, fit_yx = fit_yx,
                         control = control)
  }
  ## call method for fitted model
  test_mediation(fit, test = test, alternative = alternative, R = R,
                 level = level, type = type, order = order, ...)
}


#' @rdname test_mediation
#' @method test_mediation default
#' @export

test_mediation.default <- function(object, x, y, m, covariates = NULL,
                                   test = c("boot", "sobel"),
                                   alternative = c("twosided", "less", "greater"),
                                   R = 5000, level = 0.95,
                                   type = c("bca", "perc"),
                                   order = c("first", "second"),
                                   method = c("regression", "covariance"),
                                   robust = TRUE, family = "gaussian",
                                   contrast = FALSE, fit_yx = TRUE,
                                   control = NULL, ...) {
  ## fit mediation model
  fit <- fit_mediation(object, x = x, y = y, m = m, covariates = covariates,
                       method = method, robust = robust, family = family,
                       contrast = contrast, fit_yx = fit_yx, control = control)
  ## call method for fitted model
  test_mediation(fit, test = test, alternative = alternative, R = R,
                 level = level, type = type, order = order, ...)
}


#' @rdname test_mediation
#' @method test_mediation fit_mediation
#' @export

test_mediation.fit_mediation <- function(object, test = c("boot", "sobel"),
                                         alternative = c("twosided", "less", "greater"),
                                         R = 5000, level = 0.95,
                                         type = c("bca", "perc"),
                                         order = c("first", "second"),
                                         ...) {
  ## initializations
  test <- match.arg(test)
  alternative <- match.arg(alternative)
  p_x <- length(object$x)
  p_m <- length(object$m)
  if ((p_x > 1L || p_m > 1L) && test == "sobel") {
    test <- "boot"
    warning("Sobel test not available with multiple independent variables or ",
            "multiple mediators; using bootstrap test")
  }
  ## perform mediation analysis
  if (test == "boot") {
    # further inizializations
    level <- rep(as.numeric(level), length.out = 1)
    if (is.na(level) || level < 0 || level > 1) level <- formals()$level
    type <- match.arg(type)
    # perform bootstrap test
    boot_test_mediation(object, alternative = alternative, R = R,
                        level = level, type = type, ...)
  } else if (test == "sobel") {
    # further inizializations
    order <- match.arg(order)
    # perform Sobel test
    sobel_test_mediation(object, alternative = alternative, order = order)
  } else stop("test not implemented")
}


#' @rdname test_mediation
#' @export

robmed <- function(..., test = "boot", method = "regression", robust = TRUE) {
  test_mediation(..., test = "boot", method = "regression", robust = TRUE)
}


## internal function for bootstrap test
boot_test_mediation <- function(fit,
                                alternative = c("twosided", "less", "greater"),
                                R = 5000, level = 0.95, type = c("bca", "perc"),
                                ...) {

  # initializations
  p_x <- length(fit$x)                     # number of independent variables
  p_m <- length(fit$m)                     # number of mediators
  p_covariates <- length(fit$covariates)   # number of covariates
  nr_indirect <- p_x * p_m                 # number of indirect effects
  contrast <- fit$contrast                 # only implemented for regression fit
  have_contrast <- is.character(contrast)  # but this always works
  # useful index sequences
  seq_x <- seq_len(p_x)
  seq_m <- seq_len(p_m)

  # check whether we have a regression fit or a covariance fit
  if (inherits(fit, "reg_fit_mediation")) {

    # indices of independent variables in data matrix to be used in bootstrap
    j_x <- match(fit$x, names(fit$data)) + 1L
    # indices of the dependent variable in data matrix to be used in bootstrap
    j_y <- match(fit$y, names(fit$data)) + 1L
    # indices of mediators in data matrix to be used in bootstrap
    j_m <- match(fit$m, names(fit$data)) + 1L
    # indices of covariates in data matrix to be used in bootstrap
    j_covariates <- match(fit$covariates, names(fit$data)) + 1L
    # combine data
    n <- nrow(fit$data)
    z <- cbind(rep.int(1, n), as.matrix(fit$data))

    # perform (fast and robust) bootstrap
    robust <- fit$robust
    family <- fit$family
    if (robust == "MM") {

      # This implementation uses the simpler approximation of
      # Salibian-Barrera & Van Aelst (2008) rather than that of
      # Salibian-Barrera & Zamar (2002).  The difference is that
      # the latter also requires a correction of the residual scale.

      # extract regression models
      fit_mx <- fit$fit_mx
      fit_ymx <- fit$fit_ymx
      # extract control object from robust regressions
      # (necessary to compute correction matrices)
      psi_control <- get_psi_control(fit_ymx)  # the same for all model fits

      # extract (square root of) robustness weights and combine data
      if (p_m == 1) w_m <- sqrt(weights(fit_mx, type = "robustness"))
      else w_m <- sqrt(sapply(fit_mx, weights, type = "robustness"))
      w_y <- sqrt(weights(fit_ymx, type = "robustness"))
      # compute matrices for linear corrections and
      # extract coefficients from full sample
      if (p_m == 1) {
        corr_m <- correction_matrix(z[, c(1L, j_x, j_covariates)],
                                    weights = w_m,
                                    residuals = residuals(fit_mx),
                                    scale = fit_mx$scale,
                                    control = psi_control)
        coef_m <- coef(fit_mx)
      } else {
        corr_m <- lapply(fit$m, function(m, z) {
          correction_matrix(z, weights = w_m[, m],
                            residuals = residuals(fit_mx[[m]]),
                            scale = fit_mx[[m]]$scale,
                            control = psi_control)
        }, z = z[, c(1L, j_x, j_covariates)])
        coef_m <- lapply(fit_mx, coef)
      }
      corr_y <- correction_matrix(z[, c(1L, j_m, j_x, j_covariates)],
                                  weights = w_y,
                                  residuals = residuals(fit_ymx),
                                  scale = fit_ymx$scale,
                                  control = psi_control)
      coef_y <- coef(fit_ymx)
      # number of variables in predictor matrices of regression models
      d_m <- 1L + p_x + p_covariates
      d_y <- 1L + p_m + p_x + p_covariates

      # define function for fast and robust bootstrap
      if (p_m == 1) {
        # one mediator
        robust_bootstrap <- function(z, i, w_m, corr_m, coef_m,
                                     w_y, corr_y, coef_y) {
          # extract bootstrap sample from the data
          z_i <- z[i, , drop = FALSE]
          w_m_i <- w_m[i]
          w_y_i <- w_y[i]
          # check whether there are enough observations with nonzero weights
          if(sum(w_m_i > 0) <= d_m || sum(w_y_i > 0) <= d_y) return(NA_real_)
          # compute coefficients from weighted regression m ~ x + covariates
          weighted_x_i <- w_m_i * z_i[, c(1L, j_x, j_covariates)]
          weighted_m_i <- w_m_i * z_i[, j_m]
          coef_m_i <- solve(crossprod(weighted_x_i)) %*%
            crossprod(weighted_x_i, weighted_m_i)
          # compute coefficients from weighted regression y ~ m + x + covariates
          weighted_mx_i <- w_y_i * z_i[, c(1L, j_m, j_x, j_covariates)]
          weighted_y_i <- w_y_i * z_i[, j_y]
          coef_y_i <- solve(crossprod(weighted_mx_i)) %*%
            crossprod(weighted_mx_i, weighted_y_i)
          # compute corrected coefficients
          coef_m_i <- unname(drop(coef_m + corr_m %*% (coef_m_i - coef_m)))
          coef_y_i <- unname(drop(coef_y + corr_y %*% (coef_y_i - coef_y)))
          # compute effects
          a <- coef_m_i[1L + seq_x]
          b <- coef_y_i[2L]
          direct <- coef_y_i[2L + seq_x]
          ab <- a * b
          sum_ab <- if (p_x > 1L) sum(ab)
          total <- ab + direct
          # return effects
          c(sum_ab, ab, coef_m_i, coef_y_i, total)
        }
      } else {
        # multiple mediators
        robust_bootstrap <- function(z, i, w_m, corr_m, coef_m,
                                     w_y, corr_y, coef_y) {
          # extract bootstrap sample from the data
          z_i <- z[i, , drop = FALSE]
          w_m_i <- w_m[i, , drop = FALSE]
          w_y_i <- w_y[i]
          # check whether there are enough observations with nonzero weights
          if(any(colSums(w_m_i > 0) <= d_m) || sum(w_y_i > 0) <= d_y) {
            return(NA_real_)
          }
          # compute coefficients from weighted regression m ~ x + covariates
          coef_m_i <- lapply(fit$m, function(m, x_i) {
            w_i <- w_m_i[, m]
            weighted_x_i <- w_i * x_i
            weighted_m_i <- w_i * z_i[, m]
            solve(crossprod(weighted_x_i)) %*%
              crossprod(weighted_x_i, weighted_m_i)
          }, x_i = z_i[, c(1L, j_x, j_covariates)])
          # compute coefficients from weighted regression y ~ m + x + covariates
          weighted_mx_i <- w_y_i * z_i[, c(1L, j_m, j_x, j_covariates)]
          weighted_y_i <- w_y_i * z_i[, j_y]
          coef_y_i <- solve(crossprod(weighted_mx_i)) %*%
            crossprod(weighted_mx_i, weighted_y_i)
          # compute corrected coefficients
          coef_m_i <- unname(mapply(function(coef_m, coef_m_i, corr_m) {
            drop(coef_m + corr_m %*% (coef_m_i - coef_m))
          }, coef_m = coef_m, coef_m_i = coef_m_i, corr_m = corr_m))
          coef_y_i <- unname(drop(coef_y + corr_y %*% (coef_y_i - coef_y)))
          # compute effects
          a <- coef_m_i[1L + seq_x, ]
          b <- coef_y_i[1L + seq_m]
          direct <- coef_y_i[1L + p_m + seq_x]
          ab <- if (p_x == 1L) a * b else sweep(a, 2, b, FUN = "*")
          sum_ab <- sum(ab)
          total <- if (p_x == 1L) sum_ab + direct else rowSums(ab) + direct
          # return effects
          c(sum_ab, ab, coef_m_i, coef_y_i, total)
        }
      }

      # perform fast and robust bootstrap
      bootstrap <- local_boot(z, robust_bootstrap, R = R, w_m = w_m,
                              corr_m = corr_m, coef_m = coef_m, w_y = w_y,
                              corr_y = corr_y, coef_y = coef_y, ...)
      R <- colSums(!is.na(bootstrap$t))  # adjust number of replicates for NAs

    } else if (robust == "median") {

      # define function for standard bootstrap with median regression
      # (the fast and robust bootstrap does not work for median regression)
      median_bootstrap <- function(z, i) {
        # extract bootstrap sample from the data
        z_i <- z[i, , drop = FALSE]
        # compute coefficients from regressions m ~ x + covariates
        x_i <- z_i[, c(1L, j_x, j_covariates)]
        m_i <- z_i[, j_m]
        if (p_m == 1L) {
          coef_m_i <- unname(rq.fit(x_i, m_i, tau = 0.5)$coefficients)
        } else {
          coef_m_i <- unname(apply(m_i, 2, function(current_m_i) {
            rq.fit(x_i, current_m_i, tau = 0.5)$coefficients
          }))
        }
        # compute coefficients from regression y ~ m + x + covariates
        mx_i <- z_i[, c(1L, j_m, j_x, j_covariates)]
        y_i <- z_i[, j_y]
        coef_y_i <- unname(drop(solve(crossprod(mx_i)) %*% crossprod(mx_i, y_i)))
        # compute effects
        if (p_m == 1L) a <- coef_m_i[1L + seq_x]
        else a <- coef_m_i[1L + seq_x, ]
        b <- coef_y_i[1L + seq_m]
        direct <- coef_y_i[1L + p_m + seq_x]
        if (p_x > 1L && p_m > 1L) ab <- sweep(a, 2, b, FUN = "*")
        else ab <- a * b
        sum_ab <- if (p_x > 1L || p_m > 1L) sum(ab)
        if (p_m == 1L) total <- ab + direct
        else if (p_x == 1L) total <- sum_ab + direct
        else total <- rowSums(ab) + direct
        # return effects
        c(sum_ab, ab, coef_m_i, coef_y_i, total)
      }

      # perform standard bootstrap with median regression
      bootstrap <- local_boot(z, median_bootstrap, R = R, ...)
      R <- nrow(bootstrap$t)  # make sure that number of replicates is correct

    } else if (family == "gaussian") {

      # define function for standard bootstrap mediation test
      standard_bootstrap <- function(z, i) {
        # extract bootstrap sample from the data
        z_i <- z[i, , drop = FALSE]
        # compute coefficients from regressions m ~ x + covariates
        x_i <- z_i[, c(1L, j_x, j_covariates)]
        m_i <- z_i[, j_m]
        coef_m_i <- unname(drop(solve(crossprod(x_i)) %*% crossprod(x_i, m_i)))
        # compute coefficients from regression y ~ m + x + covariates
        mx_i <- z_i[, c(1L, j_m, j_x, j_covariates)]
        y_i <- z_i[, j_y]
        coef_y_i <- unname(drop(solve(crossprod(mx_i)) %*% crossprod(mx_i, y_i)))
        # compute effects
        if (p_m == 1L) a <- coef_m_i[1L + seq_x]
        else a <- coef_m_i[1L + seq_x, ]
        b <- coef_y_i[1L + seq_m]
        direct <- coef_y_i[1L + p_m + seq_x]
        if (p_x > 1L && p_m > 1L) ab <- sweep(a, 2, b, FUN = "*")
        else ab <- a * b
        sum_ab <- if (p_x > 1L || p_m > 1L) sum(ab)
        if (p_m == 1L) total <- ab + direct
        else if (p_x == 1L) total <- sum_ab + direct
        else total <- rowSums(ab) + direct
        # return effects
        c(sum_ab, ab, coef_m_i, coef_y_i, total)
      }

      # perform standard bootstrap
      bootstrap <- local_boot(z, standard_bootstrap, R = R, ...)
      R <- nrow(bootstrap$t)  # make sure that number of replicates is correct

    } else if (family == "select") {

      # check whether to perform the regression for the total effect
      estimate_yx <- !is.null(fit$fit_yx)
      # obtain parameters as required for package 'sn'
      control <- list(method = "MLE")

      # use values from full sample as starting values for optimization
      if (p_m == 1L) {
        start <- list(mx = fit$fit_mx$start, ymx = fit$fit_ymx$start,
                      yx = fit$fit_yx$start)
      } else {
        start <- list(mx = lapply(fit$fit_mx, "[[", "start"),
                      ymx = fit$fit_ymx$start,
                      yx = fit$fit_yx$start)
      }

      # define function for bootstrap with selection of error distribution
      select_bootstrap <- function(z, i, start, control, estimate_yx) {
        # extract bootstrap sample from the data
        z_i <- z[i, , drop = FALSE]
        # skew-t distribution can be unstable on bootstrap samples
        tryCatch({
          # compute coefficients from regression m ~ x + covariates
          x_i <- z_i[, c(1L, j_x, j_covariates)]
          if (p_m == 1L) {
            m_i <- z_i[, j_m]
            coef_mx_i <- lmselect_boot(x_i, m_i, start = start$mx,
                                       control = control)
          } else {
            coef_mx_i <- mapply(function(j, start) {
              m_i <- z_i[, j]
              coef_mx_i <- lmselect_boot(x_i, m_i, start = start,
                                         control = control)
            }, j = j_m, start = start$mx)

          }
          # compute coefficients from regression y ~ m + x + covariates
          mx_i <- z_i[, c(1L, j_m, j_x, j_covariates)]
          y_i <- z_i[, j_y]
          coef_ymx_i <- lmselect_boot(mx_i, y_i, start = start$ymx,
                                      control = control)
          # compute coefficients from regression y ~ x + covariates
          if (estimate_yx) {
            coef_yx_i <- lmselect_boot(x_i, y_i, start = start$yx,
                                       control = control)
            total <- coef_yx_i[1L + seq_x]
          } else total <- rep.int(NA_real_, p_x)
          # compute indirect effects
          if (p_m == 1L) a <- coef_mx_i[1L + seq_x]
          else a <- coef_mx_i[1L + seq_x, ]
          b <- coef_ymx_i[1L + seq_m]
          if (p_x > 1L && p_m > 1L) ab <- sweep(a, 2, b, FUN = "*")
          else ab <- a * b
          sum_ab <- if (p_x > 1L || p_m > 1L) sum(ab)
          # return effects
          c(sum_ab, ab, coef_mx_i, coef_ymx_i, total)
        }, error = function(condition) NA_real_)
      }

      # perform bootstrap with selection of error distribution
      bootstrap <- local_boot(z, select_bootstrap, R = R, start = start,
                              control = control, estimate_yx = estimate_yx,
                              ...)
      R <- colSums(!is.na(bootstrap$t))  # adjust number of replicates for NAs

    } else {

      # check whether to perform the regression for the total effect
      estimate_yx <- !is.null(fit$fit_yx)
      # obtain parameters as required for package 'sn'
      selm_args <- get_selm_args(family)
      control <- list(method = "MLE")

      # use values from full sample as starting values for optimization
      if (family == "skewnormal") {
        # starting values in centered parametrization
        if (p_m == 1L) {
          start <- list(mx = get_cp(fit$fit_mx), ymx = get_cp(fit$fit_ymx),
                        yx = get_cp(fit$fit_yx))
        } else {
          start <- list(mx = lapply(fit$fit_mx, get_cp),
                        ymx = get_cp(fit$fit_ymx),
                        yx = get_cp(fit$fit_yx))
        }
      } else {
        # starting values in direct parametrization
        if (p_m == 1L) {
          start <- list(mx = get_dp(fit$fit_mx), ymx = get_dp(fit$fit_ymx),
                        yx = get_dp(fit$fit_yx))

        } else {
          start <- list(mx = lapply(fit$fit_mx, get_dp),
                        ymx = get_dp(fit$fit_ymx),
                        yx = get_dp(fit$fit_yx))
        }
      }

      # define function for bootstrap with skew-elliptical errors
      selm_bootstrap <- function(z, i, family, start, fixed.param,
                                 control, estimate_yx) {
        # extract bootstrap sample from the data
        z_i <- z[i, , drop = FALSE]
        # skew-t distribution can be unstable on bootstrap samples
        tryCatch({
          # compute coefficients from regression m ~ x + covariates
          x_i <- z_i[, c(1L, j_x, j_covariates)]
          if (p_m == 1L) {
            m_i <- z_i[, j_m]
            fit_mx_i <- selm.fit(x_i, m_i, family = family,
                                 start = start$mx,
                                 fixed.param = fixed.param,
                                 selm.control = control)
            coef_mx_i <- unname(get_coef(fit_mx_i$param, family))
          } else {
            coef_mx_i <- mapply(function(j, start) {
              m_i <- z_i[, j]
              fit_mx_i <- selm.fit(x_i, m_i, family = family,
                                   start = start,
                                   fixed.param = fixed.param,
                                   selm.control = control)
              unname(get_coef(fit_mx_i$param, family))
            }, j = j_m, start = start$mx)
          }
          # compute coefficients from regression y ~ m + x + covariates
          mx_i <- z_i[, c(1L, j_m, j_x, j_covariates)]
          y_i <- z_i[, j_y]
          fit_ymx_i <- selm.fit(mx_i, y_i, family = family,
                                start = start$ymx,
                                fixed.param = fixed.param,
                                selm.control = control)
          coef_ymx_i <- unname(get_coef(fit_ymx_i$param, family))
          # compute coefficients from regression y ~ x + covariates
          if (estimate_yx) {
            fit_yx_i <- selm.fit(x_i, y_i, family = family,
                                 start = start$yx,
                                 fixed.param = fixed.param,
                                 selm.control = control)
            coef_yx_i <- unname(get_coef(fit_yx_i$param, family))
            total <- coef_yx_i[1L + seq_x]
          } else total <- rep.int(NA_real_, p_x)
          # compute indirect effects
          if (p_m == 1L) a <- coef_mx_i[1L + seq_x]
          else a <- coef_mx_i[1L + seq_x, ]
          b <- coef_ymx_i[1L + seq_m]
          if (p_x > 1L && p_m > 1L) ab <- sweep(a, 2, b, FUN = "*")
          else ab <- a * b
          sum_ab <- if (p_x > 1L || p_m > 1L) sum(ab)
          # return effects
          c(sum_ab, ab, coef_mx_i, coef_ymx_i, total)
        }, error = function(condition) NA_real_)
      }

      # perform bootstrap with skew-elliptical errors
      bootstrap <- local_boot(z, selm_bootstrap, R = R,
                              family = selm_args$family, start = start,
                              fixed.param = selm_args$fixed.param,
                              control = control, estimate_yx = estimate_yx,
                              ...)
      R <- colSums(!is.na(bootstrap$t))  # adjust number of replicates for NAs

    }

  } else if(inherits(fit, "cov_fit_mediation")) {

    # extract data and variable names
    x <- fit$x
    y <- fit$y
    m <- fit$m
    data <- fit$data
    # check if the robust transformation of Zu & Yuan (2010) should be applied
    if(fit$robust) {
      cov <- fit$cov
      data[] <- mapply("-", data, cov$center, SIMPLIFY=FALSE, USE.NAMES=FALSE)
      data <- weights(cov, type="consistent") * data
    }

    # perform bootstrap
    bootstrap <- local_boot(data, function(z, i) {
      # extract bootstrap sample from the data
      z_i <- z[i, , drop=FALSE]
      # compute MLE of covariance matrix on bootstrap sample
      S <- cov_ML(z_i)$cov
      # compute effects
      a <- S[m, x] / S[x, x]
      det <- S[x, x] * S[m, m] - S[m, x]^2
      b <- (-S[m, x] * S[y, x] + S[x, x] * S[y, m]) / det
      direct <- (S[m, m] * S[y, x] - S[m, x] * S[y, m]) / det
      total <- S[y, x] / S[x, x]
      c(a*b, NA_real_, a, NA_real_, b, direct, total)
    }, R=R, ...)
    R <- nrow(bootstrap$t)  # make sure that number of replicates is correct

  } else stop("method not implemented")

  # get indices of columns of bootstrap replicates that that correspond to
  # the respective models
  index_list <- get_index_list(p_x, p_m, p_covariates)
  # extract bootstrap estimates of effects other than the indirect effects
  if (p_m == 1L) indices_a <- index_list$fit_mx[1L + seq_x]
  else indices_a <- sapply(index_list$fit_mx, "[", 1L + seq_x)
  a <- colMeans(bootstrap$t[, indices_a, drop = FALSE], na.rm = TRUE)
  indices_b <- index_list$fit_ymx[1L + seq_m]
  b <- colMeans(bootstrap$t[, indices_b, drop = FALSE], na.rm = TRUE)
  indices_direct <- index_list$fit_ymx[1L + p_m + seq_x]
  direct <- colMeans(bootstrap$t[, indices_direct, drop = FALSE], na.rm = TRUE)
  # -----
  # colMeans() may return NaN instead of NA when all values are NA.  Here this
  # is the case for regression with skew-elliptical errors when the regression
  # for the total effect is not performed.
  # -----
  # total <- colMeans(bootstrap$t[, index_list$total, drop = FALSE], na.rm = TRUE)
  # -----
  bootstrap_total <- bootstrap$t[, index_list$total, drop = FALSE]
  if (all(is.na(bootstrap_total))) total <- rep.int(NA_real_, p_x)
  else total <- colMeans(bootstrap_total, na.rm = TRUE)
  # -----
  # extract bootstrap estimates and confidence intervals of indirect effects
  indices_ab <- index_list$ab
  ab <- colMeans(bootstrap$t[, indices_ab, drop = FALSE], na.rm = TRUE)
  if (nr_indirect == 1L) {
    ci <- confint(bootstrap, parm = indices_ab, level = level,
                  alternative = alternative, type = type)
  } else {
    ci <- lapply(indices_ab, function(j) {
      confint(bootstrap, parm = j, level = level,
              alternative = alternative, type = type)
    })
    ci <- do.call(rbind, ci)
    # if requested, compute contrasts of indirect effects
    if (have_contrast) {
      # list of all combinations of indices of the relevant indirect effects
      combinations <- combn(indices_ab[-1L], 2, simplify = FALSE)
      # prepare "boot" object for the calculation of the confidence intervals
      contrast_bootstrap <- bootstrap
      contrast_bootstrap$t0 <- get_contrasts(bootstrap$t0, combinations,
                                             type = contrast)
      contrast_bootstrap$t <- get_contrasts(bootstrap$t, combinations,
                                            type = contrast)
      # compute bootstrap estimates of the contrasts
      contrasts <- colMeans(contrast_bootstrap$t, na.rm = TRUE)
      # compute confidence intervals of contrasts
      indices_contrasts <- seq_along(combinations)
      contrast_ci <- lapply(indices_contrasts, function(j) {
        confint(contrast_bootstrap, parm = j, level = level,
                alternative = alternative, type = type)
      })
      contrast_ci <- do.call(rbind, contrast_ci)
      # add contrasts to results for the indirect effects
      ab <- c(ab, contrasts)
      ci <- rbind(ci, contrast_ci)
    }
    # add names to effects other than the indirect effect
    names(a) <- names(fit$a)
    if (p_m > 1L) names(b) <- names(fit$b)
    if (p_x > 1L) {
      names(direct) <- names(fit$direct)
      names(total) <- names(fit$total)
    }
    # add names for indirect effects and confidence intervals
    names(ab) <- rownames(ci) <- names(fit$ab)
  }
  # construct return object
  result <- list(a = a, b = b, direct = direct, total = total, ab = ab,
                 ci = ci, reps = bootstrap, alternative = alternative,
                 R = as.integer(R[1L]), level = level, type = type, fit = fit)
  class(result) <- c("boot_test_mediation", "test_mediation")
  result

}


## internal function for sobel test
sobel_test_mediation <- function(fit,
                                 alternative = c("twosided", "less", "greater"),
                                 order = c("first", "second"), ...) {
  # extract coefficients
  a <- fit$a
  b <- fit$b
  ab <- fit$ab
  # compute standard errors
  summary <- get_summary(fit)
  if (inherits(fit, "reg_fit_mediation")) {
    sa <- summary$fit_mx$coefficients[2L, 2L]
    sb <- summary$fit_ymx$coefficients[2L, 2L]
  } else if (inherits(fit, "cov_fit_mediation")) {
    sa <- summary$a[, 2L]
    sb <- summary$b[, 2L]
  } else stop("not implemented for this type of model fit")
  # compute test statistic and p-Value
  if (order == "first") se <- sqrt(b^2 * sa^2 + a^2 * sb^2)
  else se <- sqrt(b^2 * sa^2 + a^2 * sb^2 + sa^2 * sb^2)
  z <- ab / se
  p_value <- p_value_z(z, alternative = alternative)
  # construct return item
  result <- list(se = se, statistic = z, p_value = p_value,
                 alternative = alternative, order = order,
                 fit = fit)
  class(result) <- c("sobel_test_mediation", "test_mediation")
  result
}


# ## wrapper function for boot() that ignores unused arguments, but allows
# ## arguments for parallel computing to be passed down
# local_boot <- function(..., sim, stype, L, m, ran.gen, mle) boot(...)

## get control arguments for psi function as used in a given model fit
get_psi_control <- function(object) object$control[c("tuning.psi", "psi")]

## compute matrix for linear correction
# (see Salibian-Barrera & Van Aelst, 2008)
# The definition of the weigths in Salibian-Barrera & Van Aelst (2008) does not
# include the residual scale, whereas the robustness weights in lmrob() do.
# Hence the residual scale shows up in Equation (16) of Salibian-Barrera & Van
# Aelst (2008), but here the residual scale is already included in the weights.
correction_matrix <- function(X, weights, residuals, scale, control) {
  tmp <- Mpsi(residuals/scale, cc=control$tuning.psi, psi=control$psi, deriv=1)
  solve(crossprod(X, tmp * X)) %*% crossprod(weights * X)
}

## internal function to compute p-value based on normal distribution
p_value_z <- function(z, alternative = c("twosided", "less", "greater")) {
  # initializations
  alternative <- match.arg(alternative)
  # compute p-value
  switch(alternative, twosided = 2 * pnorm(abs(z), lower.tail = FALSE),
         less = pnorm(z), greater = pnorm(z, lower.tail = FALSE))
}


# The function for bootstrap replicates is required to return a vector.  This
# means that the columns of the bootstrap replicates contain coefficient
# estimates from different models.  This utility function returns the indices
# that correspond to the respective models, which makes it easier to extract
# the desired coefficients.
get_index_list <- function(p_x, p_m, p_covariates, indirect = TRUE) {
  # initializations
  nr_indirect <- p_x * p_m
  # numbers of coefficients in different models
  if (indirect) p_ab <- if (nr_indirect == 1L) 1L else 1L + nr_indirect
  else p_ab <- 0L
  p_mx <- rep.int(1L + p_x + p_covariates, p_m)
  p_ymx <- 1L + p_m + p_x + p_covariates
  p_total <- p_x
  p_all <- sum(p_ab, p_mx, p_ymx, p_total)
  # indices of vector for each bootstrap replicate
  indices <- seq_len(p_all)
  # the first columns correspond to indirect effect(s) of x on y
  first <- 1L
  indices_ab <- if (indirect) seq.int(first, length.out = p_ab) else integer()
  # the next columns correspond to models m ~ x + covariates
  first <- first + p_ab
  if (p_m == 1L) {
    indices_mx <- seq.int(from = first, length.out = p_mx)
  } else {
    first_mx <- first + c(0L, cumsum(p_mx[-1L]))
    indices_mx <- mapply(seq.int, from = first_mx, length.out = p_mx,
                         SIMPLIFY = FALSE, USE.NAMES = FALSE)
  }
  # the next columns correspond to model y ~ m + x + covariates
  first <- first + sum(p_mx)
  indices_ymx <- seq.int(from = first, length.out = p_ymx)
  # the last column corresponds to the total effect of x on y
  first <- first + p_ymx
  index_total <- seq.int(from = first, length.out = p_total)
  # return list of indices
  list(ab = indices_ab, fit_mx = indices_mx, fit_ymx = indices_ymx,
       total = index_total)
}


## internal functions to compute contrasts

# compute contrasts
get_contrasts <- function(x, combinations = NULL, type = "estimates") {
  # compute combinations if not supplied
  if (is.null(combinations)) {
    indices <- get_contrast_indices(x)
    combinations <- combn(indices, 2, simplify = FALSE)
  }
  # define function to compute contrasts
  if (type == "estimates") fun <- get_original_contrast
  else if (type == "absolute") fun <- get_absolute_contrast
  else stop(sprintf("%s contrasts not implemented", type))
  # compute contrasts
  sapply(combinations, fun, x = x)
}

# obtain indices to be used for computing contrasts
get_contrast_indices <- function(x) {
  if (is.matrix(x) || is.data.frame(x)) seq_len(ncol(x))
  else seq_along(x)
}

# obtain names for contrasts
get_contrast_names <- function(n) {
  if (n > 1) paste0("Contrast", seq_len(n))
  else "Contrast"
}

# compute contrasts as differences of values
get_original_contrast <- function(j, x) {
  if (is.matrix(x) || is.data.frame(x)) x[, j[1]] - x[, j[2]]
  else x[j[1]] - x[j[2]]
}

# compute contrasts as differences of absolute values
get_absolute_contrast <- function(j, x) {
  if (is.matrix(x) || is.data.frame(x)) abs(x[, j[1]]) - abs(x[, j[2]])
  else abs(x[j[1]]) - abs(x[j[2]])
}

# obtain information on how contrasts are computed
get_contrast_info <- function(x, m, type = "estimates", prefix = FALSE) {
  # initializations
  p_x <- length(x)
  p_m <- length(m)
  # names used for indirect effects
  if (p_x > 1 && p_m > 1) names <- sapply(m, paste, x, sep = ".")
  else if (p_x > 1) names <- x
  else if (p_m > 1) names <- m
  else {
    # should not happen
    stop("contrasts are only applicable in case of multiple indirect effects")
  }
  # compute combinations of names
  combinations <- combn(names, 2, simplify = FALSE)
  n_contrasts <- length(combinations)
  # obtain labels for contrasts
  labels <- get_contrast_names(n_contrasts)
  if (prefix) labels <- paste("ab", labels, sep = "_")
  # obtain information on contrasts
  if (type == "estimates") {
    fun <- function(names) paste(paste("ab", names, sep = "_"), collapse = " - ")
  } else if (type == "absolute") {
    fun <- function(names) paste(paste0("|ab_", names, "|"), collapse = " - ")
  } else stop(sprintf("%s contrasts not implemented", type))
  # return information on contrasts
  data.frame(Label = labels, Definition = sapply(combinations, fun))
}
