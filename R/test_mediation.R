# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

#' (Robust) mediation analysis
#'
#' Perform (robust) mediation analysis via a (fast-and-robust) bootstrap test
#' or Sobel's test.
#'
#' With \code{method = "regression"}, and \code{robust = TRUE} or
#' \code{robust = "MM"}, the tests are based on robust regressions with the
#' MM-estimator from \code{\link[robustbase]{lmrob}()}.  The bootstrap test is
#' thereby performed via the fast-and-robust bootstrap.  This is the default
#' behavior.
#'
#' Note that the MM-estimator of regression implemented in
#' \code{\link[robustbase]{lmrob}()} can be seen as weighted least squares
#' estimator, where the weights are dependent on how much an observation is
#' deviating from the rest.  The trick for the fast-and-robust bootstrap is
#' that on each bootstrap sample, first a weighted least squares estimator
#' is computed (using those robustness weights from the original sample)
#' followed by a linear correction of the coefficients.  The purpose of this
#' correction is to account for the additional uncertainty of obtaining the
#' robustness weights.
#'
#' With \code{method = "regression"} and \code{robust = "median"}, the tests
#' are based on median regressions with \code{\link[quantreg]{rq}()}.  Note
#' that the bootstrap test is performed via the standard bootstrap, as the
#' fast-and-robust bootstrap is not applicable.  Unlike the robust regressions
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
#' analysis via regressions and the fast-and-robust bootstrap.
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
#' \code{object} containing the independent variables of interest.
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
#' to be computed in the bootstrap test.  Possible values are \code{"perc"}
#' (the default) for the percentile bootstrap, or \code{"bca"} for the
#' bias-corrected and accelerated (BCa) bootstrap.
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
#' @param model  a character string specifying the type of model in case of
#' multiple mediators.  Possible values are \code{"parallel"} (the default) for
#' the parallel multiple mediator model, or \code{"serial"} for the serial
#' multiple mediator model.  This is only relevant for models with multiple
#' hypothesized mediators, which are currently only implemented for estimation
#' via regressions (\code{method = "regression"}).
#' @param contrast  a logical indicating whether to compute pairwise contrasts
#' of the indirect effects (defaults to \code{FALSE}).  This can also be a
#' character string, with \code{"estimates"} for computing the pairwise
#' differences of the indirect effects (such that it is tested whether two
#' indirect effects are equal), and \code{"absolute"} for computing the
#' pairwise differences of the absolute values of the indirect effects
#' (such that it is tested whether two indirect effects are equal in
#' magnitude).  This is only relevant for models with multiple indirect
#' effects, which are currently only implemented for estimation via
#' regressions (\code{method = "regression"}).  For models with multiple
#' independent variables of interest and multiple hypothesized mediators,
#' contrasts are only computed between indirect effects corresponding to
#' the same independent variable.
#' @param fit_yx  a logical indicating whether to fit the regression model
#' \code{y ~ x + covariates} to estimate the total effect (the default is
#' \code{TRUE}).  This is only relevant if \code{method = "regression"} and
#' \code{robust = FALSE}.
#' @param control  a list of tuning parameters for the corresponding robust
#' method.  For the robust MM-estimator of regression
#' (\code{method = "regression"}, and \code{robust = TRUE} or
#' \code{robust = "MM"}), a list of tuning parameters for
#' \code{\link[robustbase]{lmrob}()} as generated by
#' \code{\link{MM_reg_control}()}.  For median regression
#' (\code{method = "regression"}, and \code{robust = "median"}), a list of
#' tuning parameters for \code{\link[quantreg]{rq}()} as generated by
#' \code{\link{median_reg_control}()}.  For winsorized covariance matrix
#' estimation (\code{method = "covariance"} and \code{robust = TRUE}), a list
#' of tuning parameters for \code{\link{cov_Huber}()} as generated by
#' \code{\link{cov_control}()}.
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
#' \item{d}{in case of a serial multiple mediator model, a numeric vector
#' containing the bootstrap point estimates of the effects of proposed mediator
#' variables on other mediator variables occurring later in the sequence (only
#' \code{"boot_test_mediation"} if applicable.}
#' \item{total}{a numeric vector containing the bootstrap point estimates of
#' the total effects of the independent variables on the dependent variable
#' (only \code{"boot_test_mediation"}).}
#' \item{direct}{a numeric vector containing the bootstrap point estimates of
#' the direct effects of the independent variables on the dependent variable
#' (only \code{"boot_test_mediation"}).}
#' \item{indirect}{a numeric vector containing the bootstrap point estimates of
#' the indirect effects (only \code{"boot_test_mediation"}).}
#' \item{ab}{for back-compatibility with versions <0.10.0, the bootstrap
#' point estimates of the indirect effects are also included here (only
#' \code{"boot_test_mediation"}).  \bold{This component is deprecated and
#' may be removed as soon as the next version.}}
#' \item{ci}{a numeric vector of length two or a matrix of two columns
#' containing the bootstrap confidence intervals for the indirect effects
#' (only \code{"boot_test_mediation"}).}
#' \item{reps}{an object of class \code{"\link[boot]{boot}"} containing
#' the bootstrap replicates (only \code{"boot_test_mediation"}).  For
#' regression model fits, bootstrap replicates of the coefficients in the
#' individual regression models are stored.}
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
#' @inheritSection fit_mediation Mediation models
#'
#' @note For the fast-and-robust bootstrap, the simpler correction of
#' Salibian-Barrera & Van Aelst (2008) is used rather than the originally
#' proposed correction of Salibian-Barrera & Zamar (2002).
#'
#' The default method takes a data frame its first argument so that it can
#' easily be used with the pipe operator (\R's built-in \code{|>} or
#' \pkg{magrittr}'s \code{\%>\%}).
#'
#' Since version 1.1.0, bias-corrected and accelerated (BCa) bootstrap
#' confidence intervals are no longer recommended, and the option
#' \code{type = "bca"} may disappear in the future.
#'
#' @author Andreas Alfons
#'
#' @references
#' Alfons, A., Ates, N.Y. and Groenen, P.J.F. (2022a) A Robust Bootstrap Test
#' for Mediation Analysis.  \emph{Organizational Research Methods},
#' \bold{25}(3), 591--617.  doi:10.1177/1094428121999096.
#'
#' Alfons, A., Ates, N.Y. and Groenen, P.J.F. (2022b) Robust Mediation Analysis:
#' The \R Package \pkg{robmed}.  \emph{Journal of Statistical Software},
#' \bold{103}(13), 1--45.  doi:10.18637/jss.v103.i13.
#'
#' Alfons, A. and Schley, D.R. (2025) \emph{Robust Mediation Analysis: What We
#' Talk About When We Talk About Robustness}.  PsyArXiV,
#' doi:10.31234/osf.io/2hqdy.
#'
#' Azzalini, A. and Arellano-Valle, R. B. (2013) Maximum Penalized Likelihood
#' Estimation for Skew-Normal and Skew-t Distributions.  \emph{Journal of
#' Statistical Planning and Inference}, \bold{143}(2), 419--433.
#' doi:10.1016/j.jspi.2012.06.022.
#'
#' Preacher, K.J. and Hayes, A.F. (2004) SPSS and SAS Procedures for Estimating
#' Indirect Effects in Simple Mediation Models.  \emph{Behavior Research
#' Methods, Instruments, & Computers}, \bold{36}(4), 717--731.
#' doi:10.3758/bf03206553.
#'
#' Preacher, K.J. and Hayes, A.F. (2008) Asymptotic and Resampling Strategies
#' for Assessing and Comparing Indirect Effects in Multiple Mediator Models.
#' \emph{Behavior Research Methods}, \bold{40}(3), 879--891.
#' doi:10.3758/brm.40.3.879.
#'
#' Salibian-Barrera, M. and Van Aelst, S. (2008) Robust Model Selection Using
#' Fast and Robust Bootstrap.  \emph{Computational Statistics & Data Analysis},
#' \bold{52}(12), 5121--5135.  doi:10.1016/j.csda.2008.05.007.
#'
#' Salibian-Barrera, M. and Zamar, R. (2002) Bootstrapping Robust Estimates of
#' Regression.  \emph{The Annals of Statistics}, \bold{30}(2), 556--582.
#' doi:10.1214/aos/1021379865.
#'
#' Sobel, M.E. (1982) Asymptotic Confidence Intervals for Indirect Effects in
#' Structural Equation Models.  \emph{Sociological Methodology}, \bold{13},
#' 290--312.  doi:10.2307/270723.
#'
#' Yuan, Y. and MacKinnon, D.P. (2014) Robust Mediation Analysis Based on
#' Median Regression.  \emph{Psychological Methods}, \bold{19}(1), 1--20.
#' doi:10.1037/a0033820.
#'
#' Zu, J. and Yuan, K.-H. (2010) Local Influence and Robust Procedures for
#' Mediation Analysis.  \emph{Multivariate Behavioral Research}, \bold{45}(1),
#' 1--44.  doi:10.1080/00273170903504695.
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
#' ## seed to be used for the random number generator
#' seed <- 20150601
#'
#' ## simple mediation
#' set.seed(seed)
#' boot_simple <- test_mediation(TeamCommitment ~
#'                                 m(TaskConflict) +
#'                                   ValueDiversity,
#'                               data = BSG2014)
#' summary(boot_simple)
#'
#' # depending on the seed of the random number generator, one
#' # may get a p value slightly below or above the arbitrary
#' # 5% threshold
#' p_value(boot_simple, parm = "indirect")
#'
#' # The results in Alfons et al. (2022a) were obtained with an
#' # older version of the random number generator and with BCa
#' # bootstrap intervals (which are no longer recommended).
#' # To reproduce those results, uncomment the lines below.
#' # RNGversion("3.5.3")
#' # set.seed(seed)
#' # boot_simple <- test_mediation(TeamCommitment ~
#' #                                 m(TaskConflict) +
#' #                                   ValueDiversity,
#' #                               data = BSG2014,
#' #                               type = "bca")
#' # summary(boot_simple)
#'
#' \donttest{
#' ## serial multiple mediators
#' set.seed(seed)
#' boot_serial <- test_mediation(TeamScore ~
#'                                 serial_m(TaskConflict,
#'                                          TeamCommitment) +
#'                                 ValueDiversity,
#'                               data = BSG2014,
#'                               level = 0.9)
#' summary(boot_serial)
#'
#' ## parallel multiple mediators and control variables
#' set.seed(seed)
#' boot_parallel <- test_mediation(TeamPerformance ~
#'                                   parallel_m(ProceduralJustice,
#'                                              InteractionalJustice) +
#'                                   SharedLeadership +
#'                                   covariates(AgeDiversity,
#'                                              GenderDiversity),
#'                                 data = BSG2014)
#' summary(boot_parallel)
#' }
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
                                   type = c("perc", "bca"),
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
                                   type = c("perc", "bca"),
                                   order = c("first", "second"),
                                   method = c("regression", "covariance"),
                                   robust = TRUE, family = "gaussian",
                                   model = c("parallel", "serial"),
                                   contrast = FALSE, fit_yx = TRUE,
                                   control = NULL, ...) {
  ## fit mediation model
  fit <- fit_mediation(object, x = x, y = y, m = m, covariates = covariates,
                       method = method, robust = robust, family = family,
                       model = model, contrast = contrast, fit_yx = fit_yx,
                       control = control)
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
                                         type = c("perc", "bca"),
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
    # check confidence level and number of bootstrap samples
    level <- rep(as.numeric(level), length.out = 1L)
    if (is.na(level) || level < 0 || level > 1) {
      level <- formals()$level
      warning("confidence level must be between 0 and 1; using ",
              format(100 * level), " %")
    }
    R <- rep(as.integer(R), length.out = 1L)
    if (is.na(R) || R <= 0L) {
      R <- formals()$R
      warning("number of bootstrap samples must be a positive integer; ",
              sprintf("using %d samples", R))
    }
    # check type of confidence intervals
    type <- match.arg(type)
    # check type if BCa confidence intervals can be computed
    if (type == "bca") {
      n <- nrow(object$data)
      if (R < n) {
        type <- "perc"
        warning(sprintf("number of bootstrap samples must be at least %d ", n),
                "to compute BCa confidence intervals; using percentile ",
                "confidence intervals")
      }
    }
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

# generic function
#' @noRd
boot_test_mediation <- function(fit, ...) UseMethod("boot_test_mediation")

# method for regression fits
#' @noRd
boot_test_mediation.reg_fit_mediation <- function(fit,
                                                  alternative = "twosided",
                                                  R = 5000, level = 0.95,
                                                  type = "perc", ...) {

  # initializations
  p_x <- length(fit$x)                     # number of independent variables
  p_m <- length(fit$m)                     # number of mediators
  p_covariates <- length(fit$covariates)   # number of covariates
  model <- fit$model                       # only implemented for regression fit
  have_serial <- !is.null(model) && model == "serial"  # but this always works

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

  # in case of a serial multiple mediator model, construct list containing
  # indices of predictor variables
  if (have_serial) {
    j_mx <- vector("list", length = p_m)
    j_mx[[1L]] <- c(1L, j_x, j_covariates)
    for (j in seq_len(p_m-1L)) {
      j_mx[[j+1L]] <- c(1L, j_m[seq_len(j)], j_x, j_covariates)
    }
  }

  # perform (fast-and-robust) bootstrap
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
    psi_control <- get_psi_control(fit)

    # extract (square root of) robustness weights and combine data
    if (p_m == 1L) w_m <- sqrt(weights(fit_mx, type = "robustness"))
    else w_m <- sqrt(sapply(fit_mx, weights, type = "robustness"))
    w_y <- sqrt(weights(fit_ymx, type = "robustness"))
    # for regression m ~ x + covariates: compute matrices for linear
    # corrections and extract coefficients from full sample
    if (p_m == 1L) {
      ## single mediator
      corr_m <- correction_matrix(z[, c(1L, j_x, j_covariates)],
                                  weights = w_m,
                                  residuals = residuals(fit_mx),
                                  scale = fit_mx$scale,
                                  control = psi_control)
      coef_m <- coef(fit_mx)
    } else {
      ## multiple mediators
      # compute matrices for linear corrections
      if (have_serial) {
        corr_m <- mapply(function(m, current_j_mx) {
          correction_matrix(z[, current_j_mx], weights = w_m[, m],
                            residuals = residuals(fit_mx[[m]]),
                            scale = fit_mx[[m]]$scale,
                            control = psi_control)
        }, m = fit$m, current_j_mx = j_mx,
        SIMPLIFY = FALSE, USE.NAMES = FALSE)
      } else {
        corr_m <- lapply(fit$m, function(m, z) {
          correction_matrix(z, weights = w_m[, m],
                            residuals = residuals(fit_mx[[m]]),
                            scale = fit_mx[[m]]$scale,
                            control = psi_control)
        }, z = z[, c(1L, j_x, j_covariates)])
      }
      # extract coefficients
      coef_m <- lapply(fit_mx, coef)
    }
    # for regression y ~ m + x + covariates: compute matrices for linear
    # corrections and extract coefficients from full sample
    corr_y <- correction_matrix(z[, c(1L, j_m, j_x, j_covariates)],
                                weights = w_y,
                                residuals = residuals(fit_ymx),
                                scale = fit_ymx$scale,
                                control = psi_control)
    coef_y <- coef(fit_ymx)
    # number of variables in predictor matrices of regression models
    d_m <- if (have_serial) sapply(j_mx, length) else 1L + p_x + p_covariates
    d_y <- 1L + p_m + p_x + p_covariates

    # define function for fast-and-robust bootstrap
    if (p_m == 1L) {

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
        # return coefficients
        c(coef_m_i, coef_y_i)
      }

    } else if (have_serial) {

      # multiple serial mediators
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
        coef_m_i <- mapply(function(m, current_j_mx) {
          w_i <- w_m_i[, m]
          weighted_mx_i <- w_i * z_i[, current_j_mx]
          weighted_m_i <- w_i * z_i[, m]
          solve(crossprod(weighted_mx_i)) %*%
            crossprod(weighted_mx_i, weighted_m_i)
        }, m = fit$m, current_j_mx = j_mx,
        SIMPLIFY = FALSE, USE.NAMES = FALSE)
        # compute coefficients from weighted regression y ~ m + x + covariates
        weighted_mx_i <- w_y_i * z_i[, c(1L, j_m, j_x, j_covariates)]
        weighted_y_i <- w_y_i * z_i[, j_y]
        coef_y_i <- solve(crossprod(weighted_mx_i)) %*%
          crossprod(weighted_mx_i, weighted_y_i)
        # compute corrected coefficients
        coef_m_i <- mapply(function(coef_m, coef_m_i, corr_m) {
          unname(drop(coef_m + corr_m %*% (coef_m_i - coef_m)))
        }, coef_m = coef_m, coef_m_i = coef_m_i, corr_m = corr_m)
        coef_y_i <- unname(drop(coef_y + corr_y %*% (coef_y_i - coef_y)))
        # make sure coefficients are vectors
        coef_m_i <- unlist(coef_m_i, use.names = FALSE)
        # return coefficients
        c(coef_m_i, coef_y_i)
      }

    } else {

      # multiple parallel mediators
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
        # return coefficients
        c(coef_m_i, coef_y_i)
      }

    }

    # perform fast-and-robust bootstrap
    bootstrap <- local_boot(z, robust_bootstrap, R = R, w_m = w_m,
                            corr_m = corr_m, coef_m = coef_m, w_y = w_y,
                            corr_y = corr_y, coef_y = coef_y, ...)
    R <- colSums(!is.na(bootstrap$t))  # adjust number of replicates for NAs

  } else if (robust == "median") {

    # extract algorithm object from model fit
    control <- fit$control
    # define function for standard bootstrap with median regression
    # (the fast-and-robust bootstrap does not work for median regression)
    median_bootstrap <- function(z, i, control) {
      # extract bootstrap sample from the data
      z_i <- z[i, , drop = FALSE]
      # compute coefficients from regressions m ~ x + covariates
      if (have_serial) {
        # compute coefficients
        coef_m_i <- mapply(function(current_j_m, current_j_mx) {
          mx_i <- z_i[, current_j_mx]
          m_i <- z_i[, current_j_m]
          unname(rq.fit(mx_i, m_i, tau = 0.5,
                        method = control$method)$coefficients)
        }, current_j_m = j_m, current_j_mx = j_mx,
        SIMPLIFY = FALSE, USE.NAMES = FALSE)
        # make sure coefficients are vectors
        coef_m_i <- unlist(coef_m_i, use.names = FALSE)
      } else {
        # extract predictor matrix and mediators
        x_i <- z_i[, c(1L, j_x, j_covariates)]
        m_i <- z_i[, j_m]
        # compute coefficients
        if (p_m == 1L) {
          coef_m_i <- unname(rq.fit(x_i, m_i, tau = 0.5,
                                    method = control$method)$coefficients)
        } else {
          coef_m_i <- unname(apply(m_i, 2, function(current_m_i) {
            rq.fit(x_i, current_m_i, tau = 0.5,
                   method = control$method)$coefficients
          }))
        }
      }
      # compute coefficients from regression y ~ m + x + covariates
      mx_i <- z_i[, c(1L, j_m, j_x, j_covariates)]
      y_i <- z_i[, j_y]
      coef_y_i <- unname(rq.fit(mx_i, y_i, tau = 0.5,
                                method = control$method)$coefficients)
      # return coefficients
      c(coef_m_i, coef_y_i)
    }

    # perform standard bootstrap with median regression
    bootstrap <- local_boot(z, median_bootstrap, R = R, control = control, ...)
    R <- nrow(bootstrap$t)  # make sure that number of replicates is correct

  } else if (family == "gaussian") {

    # define function for standard bootstrap mediation test
    standard_bootstrap <- function(z, i) {
      # extract bootstrap sample from the data
      z_i <- z[i, , drop = FALSE]
      # compute coefficients from regressions m ~ x + covariates
      if (have_serial) {
        # compute coefficients
        coef_m_i <- mapply(function(current_j_m, current_j_mx) {
          mx_i <- z_i[, current_j_mx]
          m_i <- z_i[, current_j_m]
          coef_m_i <- unname(drop(solve(crossprod(mx_i)) %*% crossprod(mx_i, m_i)))
        }, current_j_m = j_m, current_j_mx = j_mx,
        SIMPLIFY = FALSE, USE.NAMES = FALSE)
        # make sure coefficients are vectors
        coef_m_i <- unlist(coef_m_i, use.names = FALSE)
      } else {
        # compute coefficients
        x_i <- z_i[, c(1L, j_x, j_covariates)]
        m_i <- z_i[, j_m]
        coef_m_i <- unname(drop(solve(crossprod(x_i)) %*% crossprod(x_i, m_i)))
      }
      # compute coefficients from regression y ~ m + x + covariates
      mx_i <- z_i[, c(1L, j_m, j_x, j_covariates)]
      y_i <- z_i[, j_y]
      coef_y_i <- unname(drop(solve(crossprod(mx_i)) %*% crossprod(mx_i, y_i)))
      # return effects
      c(coef_m_i, coef_y_i)
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
        # extract predictor matrix
        if (model != "serial" || estimate_yx) {
          x_i <- z_i[, c(1L, j_x, j_covariates)]
        }
        # compute coefficients from regression m ~ x + covariates
        if (have_serial) {
          # compute coefficients
          coef_mx_i <- mapply(function(current_j_m, current_j_mx, start) {
            mx_i <- z_i[, current_j_mx]
            m_i <- z_i[, current_j_m]
            lmselect_boot(mx_i, m_i, start = start, control = control)
          }, current_j_m = j_m, current_j_mx = j_mx, start = start$mx,
          SIMPLIFY = FALSE, USE.NAMES = FALSE)
          # make sure coefficients are vectors
          coef_mx_i <- unlist(coef_mx_i, use.names = FALSE)
        } else {
          # compute coefficients
          if (p_m == 1L) {
            m_i <- z_i[, j_m]
            coef_mx_i <- lmselect_boot(x_i, m_i, start = start$mx,
                                       control = control)
          } else {
            coef_mx_i <- mapply(function(j, start) {
              m_i <- z_i[, j]
              lmselect_boot(x_i, m_i, start = start, control = control)
            }, j = j_m, start = start$mx)
          }
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
        } else coef_yx_i <- NULL
        # return coefficients
        c(coef_mx_i, coef_ymx_i, coef_yx_i)
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
        # extract predictor matrix
        if (model != "serial" || estimate_yx) {
          x_i <- z_i[, c(1L, j_x, j_covariates)]
        }
        # compute coefficients from regression m ~ x + covariates
        if (have_serial) {
          # compute coefficients
          coef_mx_i <- mapply(function(current_j_m, current_j_mx, start) {
            mx_i <- z_i[, current_j_mx]
            m_i <- z_i[, current_j_m]
            fit_mx_i <- selm.fit(mx_i, m_i, family = family,
                                 start = start,
                                 fixed.param = fixed.param,
                                 selm.control = control)
            unname(get_coef(fit_mx_i$param, family))
          }, current_j_m = j_m, current_j_mx = j_mx, start = start$mx,
          SIMPLIFY = FALSE, USE.NAMES = FALSE)
          # make sure coefficients are vectors
          coef_mx_i <- unlist(coef_mx_i, use.names = FALSE)
        } else {
          # compute coefficients
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
        } else coef_yx_i <- NULL
        # return effects
        c(coef_mx_i, coef_ymx_i, coef_yx_i)
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

  # extract bootstrap replicates for the different effects
  boot_list <- extract_boot(fit, boot = bootstrap)
  # compute bootstrap estimates of the different effects
  estimate_list <- lapply(boot_list, colMeans, na.rm = TRUE)
  # colMeans() may return NaN instead of NA when all values are NA.  Here
  # this is the case for regression with skew-elliptical errors when the
  # regression for the total effect is not performed.
  fix_total <- is.nan(estimate_list$total)
  if (any(fix_total)) estimate_list$total[fix_total] <- NA_real_
  # compute confidence intervals of indirect effects
  ci <- boot_ci(fit$indirect, boot_list$indirect, object = bootstrap,
                alternative = alternative, level = level, type = type)

  # return results
  result <- list(ab = estimate_list$indirect,  # for back-compatibility, will be removed
                 ci = ci, reps = bootstrap, alternative = alternative,
                 R = as.integer(R[1L]), level = level, type = type,
                 fit = fit)
  result <- c(estimate_list, result)
  class(result) <- c("boot_test_mediation", "test_mediation")
  result

}

# method for regression fits
#' @noRd
boot_test_mediation.cov_fit_mediation <- function(fit,
                                                  alternative = "twosided",
                                                  R = 5000, level = 0.95,
                                                  type = "perc", ...) {

  # extract variable names and data
  x <- fit$x
  y <- fit$y
  m <- fit$m
  data <- fit$data
  # check if the robust transformation of Zu & Yuan (2010) should be applied
  if (fit$robust) {
    cov <- fit$cov
    data[] <- mapply("-", data, cov$center, SIMPLIFY = FALSE, USE.NAMES = FALSE)
    data <- weights(cov, type = "consistent") * data
  }

  # perform bootstrap
  bootstrap <- local_boot(data, function(z, i) {
    # extract bootstrap sample from the data
    z_i <- z[i, , drop = FALSE]
    # compute MLE of covariance matrix on bootstrap sample
    S <- cov_ML(z_i)$cov
    # compute effects
    det <- S[x, x] * S[m, m] - S[m, x]^2
    a <- S[m, x] / S[x, x]
    b <- (-S[m, x] * S[y, x] + S[x, x] * S[y, m]) / det
    total <- S[y, x] / S[x, x]
    direct <- (S[m, m] * S[y, x] - S[m, x] * S[y, m]) / det
    c(a, b, total, direct, a * b)
  }, R = R, ...)
  R <- nrow(bootstrap$t)  # make sure that number of replicates is correct

  # compute bootstrap estimates of effects
  a <- mean(bootstrap$t[, 1L], na.rm = TRUE)
  b <- mean(bootstrap$t[, 2L], na.rm = TRUE)
  total <- mean(bootstrap$t[, 3L], na.rm = TRUE)
  direct <- mean(bootstrap$t[, 4L], na.rm = TRUE)
  indirect <- mean(bootstrap$t[, 5L], na.rm = TRUE)
  # compute confidence interval for indirect effect
  ci <- extract_ci(parm = 5L, object = bootstrap, alternative = alternative,
                   level = level, type = type)

  # return results
  result <- list(a = a, b = b, total = total, direct = direct,
                 indirect = indirect,
                 ab = indirect,  # for back-compatibility, will be removed
                 ci = ci, reps = bootstrap, alternative = alternative,
                 R = R, level = level, type = type, fit = fit)
  class(result) <- c("boot_test_mediation", "test_mediation")
  result

}


## internal function for sobel test
sobel_test_mediation <- function(fit, alternative = "twosided",
                                 order = "first", ...) {
  # extract coefficients
  a <- fit$a
  b <- fit$b
  indirect <- fit$indirect
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
  z <- indirect / se
  p_value <- p_value_z(z, alternative = alternative)
  # return results
  result <- list(se = se, statistic = z, p_value = p_value,
                 alternative = alternative, order = order,
                 fit = fit)
  class(result) <- c("sobel_test_mediation", "test_mediation")
  result
}
