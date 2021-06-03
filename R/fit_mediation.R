# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

#' (Robustly) fit a mediation model
#'
#' (Robustly) estimate the effects in a mediation model.
#'
#' With \code{method = "regression"}, and \code{robust = TRUE} or
#' \code{robust = "MM"}, the effects are computed via the robust MM-estimator
#' of regression from \code{\link[robustbase]{lmrob}()}.  This is the default
#' behavior.
#'
#' With \code{method = "regression"} and \code{robust = "median"}, the effects
#' are estimated via median regressions with \code{\link[quantreg]{rq}()}.
#' Unlike the robust MM-regressions above, median regressions are not robust
#' against outliers in the explanatory variables.
#'
#' With \code{method = "regression"}, \code{robust = FALSE} and
#' \code{family = "select"}, the error distribution to be used in maximum
#' likelihood estimation of the regression models is selected via BIC.  The
#' following error distributions are included in the selection procedure: a
#' normal distribution, a skew-normal distribution, Student's t distribution,
#' and a skew-t distribution.  Note that the parameters of those distributions
#' are estimated as well.  The skew-normal and skew-t distributions thereby
#' use a centered parametrization such that the residuals are (approximately)
#' centered around 0.  Moreover, the skew-t distribution is only evaluated in
#' the selection procedure if both the skew-normal and Student's t distribution
#' yield an improvement in BIC over the normal distribution.  Otherwise the
#' estimation with a skew-t error distribution can be unstable.  Furthermore,
#' this saves a considerable amount of computation time in a bootstrap test,
#' as estimation with those error distributions is orders of magnitude slower
#' than any other estimation procedure in package \pkg{robmed}.
#'
#' With \code{method = "covariance"} and \code{robust = TRUE}, the effects are
#' estimated based on a Huber M-estimator of location and scatter.  Note that
#' this covariance-based approach is less robust than the approach based on
#' robust MM-regressions described above.
#'
#' @aliases print.fit_mediation summary.reg_fit_mediation
#' summary.cov_fit_mediation
#'
#' @param object  the first argument will determine the method of the generic
#' function to be dispatched.  For the default method, this should be a data
#' frame containing the variables.
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
#' @param method  a character string specifying the method of
#' estimation.  Possible values are \code{"regression"} (the default)
#' to estimate the effects via regressions, or \code{"covariance"} to
#' estimate the effects via the covariance matrix.  Note that the effects are
#' always estimated via regressions if more than one independent variable or
#' hypothesized mediator is specified, or if control variables are supplied.
#' @param robust  a logical indicating whether to robustly estimate the effects
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
#' \sQuote{Details}).  This is only relevant if \code{method = "regression"}
#' and \code{robust = FALSE}.
#' @param contrast  a logical indicating whether to compute pairwise contrasts
#' of the indirect effects (defaults to \code{FALSE}).  This can also be a
#' character string, with \code{"estimates"} for computing the pairwise
#' differences of the indirect effects, and \code{"absolute"} for computing the
#' pairwise differences of the absolute values of the indirect effects.  This
#' is only relevant for models with multiple indirect effects, which are
#' currently only implemented for estimation via regressions
#' (\code{method = "regression"}).
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
#' @param \dots  additional arguments to be passed down.  For the default
#' method, this can be used to specify tuning parameters directly instead
#' of via \code{control}.
#'
#' @return An object inheriting from class \code{"fit_mediation"} (class
#' \code{"reg_fit_mediation"} if \code{method = "regression"} or
#' \code{"cov_fit_mediation"} if \code{method = "covariance"}) with
#' the following components:
#' \item{a}{a numeric vector containing the point estimates of the effects of
#' the independent variables on the proposed mediator variables.}
#' \item{b}{a numeric vector containing the point estimates of the direct
#' effects of the proposed mediator variables on the dependent variable.}
#' \item{direct}{a numeric vector containing the point estimates of the direct
#' effects of the independent variables on the dependent variable.}
#' \item{total}{a numeric vector containing the point estimates of the total
#' effects of the independent variables on the dependent variable.}
#' \item{ab}{a numeric vector containing the point estimates of the indirect
#' effects.}
#' \item{fit_mx}{an object of class \code{"\link[robustbase]{lmrob}"},
#' \code{"\link[quantreg]{rq}"}, \code{"\link[stats]{lm}"} or \code{"lmse"}
#' containing the estimation results from the regression of the proposed
#' mediator variable on the independent variables, or a list of such objects
#' in case of more than one hypothesized mediator (only
#' \code{"reg_fit_mediation"}).}
#' \item{fit_ymx}{an object of class \code{"\link[robustbase]{lmrob}"},
#' \code{"\link[quantreg]{rq}"}, \code{"\link[stats]{lm}"} or \code{"lmse"}
#' containing the estimation results from the regression of the dependent
#' variable on the proposed mediator and independent variables (only
#' \code{"reg_fit_mediation"}).}
#' \item{fit_yx}{an object of class \code{"\link[stats]{lm}"} or \code{"lmse"}
#' containing the estimation results from the regression of the dependent
#' variable on the independent variables (only \code{"reg_fit_mediation"}
#' if arguments \code{robust = FALSE} and \code{fit_yx = TRUE} were used).}
#' \item{cov}{an object of class \code{"\link{cov_Huber}"} or
#' \code{"\link{cov_ML}"} containing the covariance matrix estimates
#' (only \code{"cov_fit_mediation"}).}
#' \item{x, y, m, covariates}{character vectors specifying the respective
#' variables used.}
#' \item{data}{a data frame containing the independent, dependent and
#' proposed mediator variables, as well as covariates.}
#' \item{robust}{either a logical indicating whether the effects were estimated
#' robustly, or one of the character strings \code{"MM"} and \code{"median"}
#' specifying the type of robust regressions.}
#' \item{contrast}{either a logical indicating whether contrasts of the
#' indirect effects were computed, or one of the character strings
#' \code{"estimates"} and \code{"absolute"} specifying the type of contrasts
#' of the indirect effects (only \code{"reg_fit_mediation"}).}
#' \item{control}{a list of tuning parameters used (if applicable).}
#'
#' @note The formula interface is still experimental and may change in future
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
#' Yuan, Y. and MacKinnon, D.P. (2014) Robust mediation analysis based on
#' median regression. \emph{Psychological Methods}, \bold{19}(1),
#' 1--20.
#'
#' Zu, J. and Yuan, K.-H. (2010) Local influence and robust procedures for
#' mediation analysis. \emph{Multivariate Behavioral Research}, \bold{45}(1),
#' 1--44.
#'
#' @seealso \code{\link{test_mediation}()}
#'
#' \code{\link[robustbase]{lmrob}()}, \code{\link[stats]{lm}()},
#' \code{\link{cov_Huber}()}, \code{\link{cov_ML}()}
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
#' fit1 <- fit_mediation(TeamCommitment ~ m(TaskConflict) + ValueDiversity,
#'                       data = BSG2014)
#' test1 <- test_mediation(fit1)
#' summary(test1)
#'
#' # default method
#' set.seed(seed)
#' fit2 <- fit_mediation(BSG2014,
#'                       x = "ValueDiversity",
#'                       y = "TeamCommitment",
#'                       m = "TaskConflict")
#' test2 <- test_mediation(fit2)
#' summary(test2)
#'
#' @keywords multivariate
#'
#' @import boot
#' @import robustbase
#' @importFrom quantreg rq.fit
#' @importFrom utils combn
#' @export

fit_mediation <- function(object, ...) UseMethod("fit_mediation")


#' @rdname fit_mediation
#' @method fit_mediation formula
#' @export

fit_mediation.formula <- function(formula, data, ...) {
  ## prepare model frame
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  ## make sure that there are no interaction terms
  ord <- attr(mt, "order")
  if (length(ord) > 0 && any(ord > 1)) {
    stop("interaction terms are not implemented for mediation models")
  }
  # make sure that dependent variable is specified
  if (attr(mt, "response") == 0) {
    stop("no dependent variable specified in formula")
  }
  d <- dim(mf[[1]])
  if (!is.null(d) && d[2] > 1) {
    stop("only one dependent variable allowed in formula")
  }
  index_y <- 1
  y <- names(mf)[index_y]
  # make sure that mediators are specified
  index_m <- which(sapply(mf, inherits, "parallel_mediators"))
  if (length(index_m) == 0) {
    stop("mediators must be specified using m() in the formula")
  } else if (length(index_m) > 1) {
    stop("use m() only once in the formula to specify all mediators")
  }
  m <- colnames(mf[[index_m]])
  # check if covariates are specified
  index_covariates <- which(sapply(mf, inherits, "covariates"))
  if (length(index_covariates) > 1) {
    stop("use covariates() only once in the formula to specify all covariates")
  }
  have_covariates <- length(index_covariates) > 0
  covariates <- if (have_covariates) colnames(mf[[index_covariates]])
  # make sure that independent variable is specified
  if (length(mf) == (2 + have_covariates)) {
    stop("no independent variable specified in formula")
  }
  index_x <- setdiff(seq(2, length(mf)), c(index_m, index_covariates))
  x <- names(mf)[index_x]
  # prepare to rebuild data frame
  mf <- as.list(mf)
  names(mf)[c(index_m, index_covariates)] <- ""  # ensure the correct names
  mf[[index_m]] <- as.data.frame(mf[[index_m]])
  if (have_covariates) {
    mf[[index_covariates]] <- as.data.frame(mf[[index_covariates]])
  }
  # add additional arguments to be passed to data.frame()
  mf$check.names <- FALSE
  mf$stringsAsFactors <- TRUE
  # rebuild data frame
  data <- do.call(data.frame, mf)
  # call default method
  fit_mediation(data, x = x, y = y, m = m, covariates = covariates, ...)
}


#' @rdname fit_mediation
#' @method fit_mediation default
#' @export

fit_mediation.default <- function(object, x, y, m, covariates = NULL,
                                  method = c("regression", "covariance"),
                                  robust = TRUE, family = "gaussian",
                                  contrast = FALSE, fit_yx = TRUE,
                                  control = NULL, ...) {
  ## initializations
  # prepare data set
  data <- as.data.frame(object)
  # check independent variable
  x <- data[, x, drop = FALSE]
  p_x <- ncol(x)
  if (p_x == 0L) stop("at least one independent variable required")
  convert_x <- !all(sapply(x, is.numeric))
  # check dependent variable
  y <- data[, y, drop = FALSE]
  p_y <- ncol(y)
  if (p_y != 1L) stop("exactly one dependent variable required")
  if (!is.numeric(y[, 1L])) {
    stop("currently only implemented for a numeric dependent variable")
  }
  # check hypothesized mediator variables
  m <- data[, m, drop = FALSE]
  p_m <- ncol(m)
  if (p_m == 0L) stop("at least one hypothesized mediator variable required")
  if (!all(sapply(m, is.numeric))) {
    stop("currently only implemented for numeric hypothesized mediators")
  }
  # extract covariates
  covariates <- data[, covariates, drop = FALSE]
  p_covariates <- ncol(covariates)
  have_covariates <- p_covariates > 0L
  convert_covariates <- have_covariates && !all(sapply(covariates, is.numeric))
  # reorder columns of data frame
  data <- cbind(x, y, m, covariates)
  # extract names
  cn <- names(data)
  x <- cn[seq_len(p_x)]
  y <- cn[p_x + 1L]
  m <- cn[p_x + 1L + seq_len(p_m)]
  covariates <- cn[-(seq_len(p_x + 1L + p_m))]
  # remove incomplete observations
  data <- data[complete.cases(data), ]
  # if necessary, convert non-numeric independent variables
  if (convert_x) {
    # construct variables for design matrix as usual
    x <- data[, x, drop = FALSE]
    x <- model.matrix(~ ., data = x)[, -1L, drop = FALSE]
    p_x <- ncol(x)
    # replace independent variable in data frame with converted one
    data <- cbind(x, data[, c(y, m, covariates), drop = FALSE])
    # update variable name
    x <- colnames(x)
  }
  # if necessary, convert non-numeric covariates
  if (convert_covariates) {
    # construct variables for design matrix as usual
    covariates <- data[, covariates, drop = FALSE]
    covariates <- model.matrix(~ ., data = covariates)[, -1L, drop = FALSE]
    # replace covariates in data frame with converted ones
    data <- cbind(data[, c(x, y, m), drop = FALSE], covariates)
    # update number of covariates and variable names
    p_covariates <- ncol(covariates)
    covariates <- colnames(covariates)
  }
  # check if there are enough observations
  d <- dim(data)
  if (d[1L] <= d[2L]) stop("not enough observations")
  # check other arguments
  method <- match.arg(method)
  if ((p_x > 1L || p_m > 1L || have_covariates) && method == "covariance") {
    method <- "regression"
    warning("covariance method not available with multiple independent ",
            "variables, multiple mediators, or covariates; using regression ",
            "method")
  }
  # regression method requires some more checks, covariance method is simpler
  if (method == "regression") {
    # check for robust method
    if (is.logical(robust)) {
      robust <- isTRUE(robust)
      if (robust) robust <- "MM"
    } else robust <- match.arg(robust, choices = c("MM", "median"))
    if (robust == "MM" && is.null(control)) control <- reg_control(...)
    # check error distribution
    if (is.character(robust)) {
      family <- "gaussian"
      fit_yx <- FALSE  # this is actually ignored for robust estimation
    } else {
      families <- c("gaussian", "student", "skewnormal", "skewt", "select")
      family <- match.arg(family, choices = families)
      fit_yx <- isTRUE(fit_yx)
    }
    # check whether to compute differences of indirect effect(s)
    if (p_x == 1L && p_m == 1L) contrast <- FALSE
    else if (is.logical(contrast)) {
      contrast <- isTRUE(contrast)
      if (contrast) contrast <- "estimates"
    } else contrast <- match.arg(contrast, choices = c("estimates", "absolute"))
    # estimate effects
    reg_fit_mediation(data, x = x, y = y, m = m, covariates = covariates,
                      robust = robust, family = family, contrast = contrast,
                      fit_yx = fit_yx, control = control)
  } else {
    # check for robust method
    robust <- isTRUE(robust)
    if (robust && is.null(control)) control <- cov_control(...)
    # estimate effects
    cov_fit_mediation(data, x = x, y = y, m = m, robust = robust,
                      control = control)
  }
}


## estimate the effects in a mediation model via regressions
reg_fit_mediation <- function(data, x, y, m, covariates = character(),
                              robust = "MM", family = "gaussian",
                              contrast = FALSE, fit_yx = TRUE,
                              control = reg_control()) {
  # number of indirect effects
  p_x <- length(x)
  p_m <- length(m)
  nr_indirect <- p_x * p_m
  # construct predictor matrices for regression models
  n <- nrow(data)
  predictors_mx <- as.matrix(data[, c(x, covariates), drop = FALSE])
  predictors_ymx <- as.matrix(data[, c(m, x, covariates), drop = FALSE])
  # other initializations
  have_robust <- is.character(robust)
  have_contrast <- is.character(contrast)
  estimate_yx <- fit_yx  # avoid name conflicts
  # compute regression models
  if (have_robust) {
    # for the robust methods, the total effect is estimated as c' = ab + c
    # to satisfy this relationship
    # (For median regression, I believe this only makes sense if the
    # conditional distribution is symmetric.  This is ok, because the robust
    # methods assume a normal distribution, which is reflected by setting
    # family = "gaussian".)
    if (robust == "MM") {
      # MM-estimator for robust regression
      if (p_m == 1L) {
        fit_mx <- lmrob_fit(predictors_mx, data[, m], control = control)
      } else {
        fit_mx <- lapply(m, function(m_j) {
          lmrob_fit(predictors_mx, data[, m_j], control = control)
        })
        names(fit_mx) <- m
      }
      fit_ymx <- lmrob_fit(predictors_ymx, data[, y], control = control)
    } else if (robust == "median") {
      # LAD-estimator for median regression
      if (p_m == 1L) fit_mx <- rq_fit(predictors_mx, data[, m], tau = 0.5)
      else {
        fit_mx <- lapply(m, function(m_j) {
          rq_fit(predictors_mx, data[, m_j], tau = 0.5)
        })
        names(fit_mx) <- m
      }
      fit_ymx <- rq_fit(predictors_ymx, data[, y], tau = 0.5)
    } else {
      # shouldn't happen
      stop(sprintf('robust method "%s" not implemented', robust))
    }
    # neither method fits the direct path
    fit_yx <- NULL
  } else if (family == "gaussian") {
    # For OLS, there is not much additional cost in performing the regression
    # for the total effect.  But it is up to the user whether this is done.
    if (p_m == 1L) fit_mx <- lm_fit(predictors_mx, data[, m])
    else {
      fit_mx <- lapply(m, function(m_j) lm_fit(predictors_mx, data[, m_j]))
      names(fit_mx) <- m
    }
    fit_ymx <- lm_fit(predictors_ymx, data[, y])
    fit_yx <- if (estimate_yx) lm_fit(predictors_mx, data[, y])
  } else if (family == "select") {
    # select among normal, skew-normal, t and skew-t errors
    if (p_m == 1L) fit_mx <- lmselect_fit(predictors_mx, data[, m])
    else {
      fit_mx <- lapply(m, function(m_j) {
        lmselect_fit(predictors_mx, data[, m_j])
      })
      names(fit_mx) <- m
    }
    fit_ymx <- lmselect_fit(predictors_ymx, data[, y])
    fit_yx <- if (estimate_yx) lmselect_fit(predictors_mx, data[, y])
  } else {
    # obtain parameters as required for package 'sn'
    selm_args <- get_selm_args(family)
    # perform regression with skew-elliptical errors
    if (p_m == 1L) {
      fit_mx <- selm_fit(predictors_mx, data[, m], family = selm_args$family,
                         fixed.param = selm_args$fixed.param)
    } else {
      fit_mx <- lapply(m, function(m_j) {
        selm_fit(predictors_mx, data[, m_j], family = selm_args$family,
                 fixed.param = selm_args$fixed.param)
      })
      names(fit_mx) <- m
    }
    # The relationship ab + direct = total doesn't hold, so we need to
    # perform the regression for the total effect.  But it is up to the
    # user whether this is done, which can save computation time if the
    # user is not interested in the total effect.
    fit_ymx <- selm_fit(predictors_ymx, data[, y], family = selm_args$family,
                        fixed.param = selm_args$fixed.param)
    fit_yx <- if (estimate_yx) {
      selm_fit(predictors_mx, data[, y], family = selm_args$family,
               fixed.param = selm_args$fixed.param)
    }
  }
  # extract effects a and b and compute indirect effect(s)
  if (p_m == 1L) {
    # extract effects a and b
    if (p_x == 1L) a <- unname(coef(fit_mx)[2L])
    else a <- coef(fit_mx)[1L + seq_len(p_x)]
    b <- unname(coef(fit_ymx)[2L])
    # compute indirect effect(s)
    ab <- a * b
  } else {
    # extract effects a and b
    if (p_x == 1L) {
      a <- sapply(fit_mx, function(fit) unname(coef(fit)[2L]))
    } else {
      a <- lapply(fit_mx, function(fit) coef(fit)[1L + seq_len(p_x)])
    }
    b <- coef(fit_ymx)[1L + seq_len(p_m)]
    # compute indirect effect(s)
    ab_list <- mapply(function(current_a, current_b) current_a * current_b,
                      current_a = a, current_b = b, SIMPLIFY = FALSE)
    # make sure we have vectors
    if (p_x > 1L) a <- unlist(a, use.names = TRUE)
    ab <- unlist(ab_list, use.names = TRUE)
  }
  # if applicable, compute total indirect effect and contrasts
  if (nr_indirect > 1L) {
    # if requested, compute contrasts
    if (have_contrast) {
      contrasts <- get_contrasts(ab, type = contrast)
      n_contrasts <- length(contrasts)
      names(contrasts) <- get_contrast_names(n_contrasts)
    } else contrasts <- NULL
    # compute total indirect effect
    ab_total <- sum(ab)
  }
  # extract direct effect(s)
  if (p_x == 1L) direct <- unname(coef(fit_ymx)[2L + p_m])
  else direct <- coef(fit_ymx)[1L + p_m + seq_len(p_x)]
  # extract total effect(s)
  if (have_robust || (family == "gaussian" && !estimate_yx)) {
    if (p_x == 1L) {
      total <- if (p_m == 1L) ab + direct else ab_total + direct
    } else {
      if (p_m == 1) total <- ab + direct
      else {
        total <- sapply(x, function(current_x) {
          current_ab <- sapply(ab_list, "[", current_x)
          sum(current_ab) + unname(direct[current_x])
        })
      }
    }
  } else if (estimate_yx) {
    if (p_x == 1) total <- unname(coef(fit_yx)[2L])
    else total <- coef(fit_yx)[1L + seq_len(p_x)]
  } else {
    if (p_x == 1) total <- NA_real_
    else {
      total <- rep.int(NA_real_, p_x)
      names(total) <- x
    }
  }
  # combine total indirect effect, individual indirect effects, and contrasts
  if (nr_indirect > 1L) ab <- c(Total = ab_total, ab, contrasts)
  # return results
  result <- list(a = a, b = b, direct = direct, total = total, ab = ab,
                 fit_mx = fit_mx, fit_ymx = fit_ymx, fit_yx = fit_yx,
                 x = x, y = y, m = m, covariates = covariates, data = data,
                 robust = robust, family = family, contrast = contrast)
  if (robust == "MM") result$control <- control
  class(result) <- c("reg_fit_mediation", "fit_mediation")
  result
}


## estimate the effects in a mediation model via the covariance matrix
cov_fit_mediation <- function(data, x, y, m, robust = TRUE,
                              control = cov_control()) {
  # compute scatter matrix (Huber M-estimator or MLE of covariance matrix)
  cov <- if (robust) cov_Huber(data, control = control) else cov_ML(data)
  S <- cov$cov
  # compute coefficients of mediation model
  a <- S[m, x] / S[x, x]
  det <- S[x, x] * S[m, m] - S[m, x]^2
  b <- (-S[m, x] * S[y, x] + S[x, x] * S[y, m]) / det
  direct <- (S[m, m] * S[y, x] - S[m, x] * S[y, m]) / det
  total <- S[y, x] / S[x, x]
  # return results
  result <- list(a = a, b = b, direct = direct, total = total, ab = a * b,
                 cov = cov, x = x, y = y, m = m, covariates = character(),
                 data = data, robust = robust)
  if(robust) result$control <- control
  class(result) <- c("cov_fit_mediation", "fit_mediation")
  result
}


## model fitting functions that make summary() work
# This will allow fit_mediation() to fit the regression models without a
# formula, which is necessary to make the formula interface work.  Otherwise
# formulas that contain log(.) or similar terms would cause errors when fitting
# the models in fit_mediation().

# OLS regression
lm_fit <- function(x, y, intercept = TRUE) {
  # if requested, add constant for intercept
  if (intercept) {
    n <- nrow(x)
    x <- cbind("(Intercept)" = rep.int(1, n), x)
  }
  # fit the linear model
  fit <- lm.fit(x, y)
  # Add a dummy formula as a 'terms' component to make summary() method work.
  # This 'terms' component needs to have an attribute that specifies whether
  # the model has an intercept, that's all.  It's not used in any other way.
  f <- as.formula(NULL)
  attr(f, "intercept") <- as.integer(intercept)
  fit$terms <- f
  # add the class and return the model fit
  class(fit) <- "lm"
  fit
}

# MM regression
lmrob_fit <- function(x, y, intercept = TRUE, control = reg_control()) {
  # if requested, add constant for intercept
  if (intercept) {
    n <- nrow(x)
    x <- cbind("(Intercept)" = rep.int(1, n), x)
  }
  # fit the linear model
  fit <- lmrob.fit(x, y, control = control)
  # Add a dummy formula as a 'terms' component to make summary() method work.
  # This 'terms' component needs to have an attribute that specifies whether
  # the model has an intercept, that's all.  It's not used in any other way.
  f <- as.formula(NULL)
  attr(f, "intercept") <- as.integer(intercept)
  fit$terms <- f
  # class is already added in lmrob.fit(), so just return the model fit
  fit
}

# quantile (median) regression
rq_fit <- function(x, y, intercept = TRUE, tau = 0.5) {
  # summary() method requires data frame and formula to extract response and
  # predictor matrix
  data <- data.frame(y, x, check.names = FALSE)
  formula <- if (intercept) y ~ . else y ~ 0 + .
  # if requested, add constant for intercept
  if (intercept) {
    n <- nrow(x)
    x <- cbind("(Intercept)" = rep.int(1, n), x)
  }
  # fit the linear model
  fit <- rq.fit(x, y, tau = tau, method = "br")
  fit$method <- "br"
  # construct terms object from formula that specifies whether there is a
  # response and an intercept
  terms <- formula
  attr(terms, "response") <- 1L
  attr(terms, "intercept") <- as.integer(intercept)
  # summary() method requires that the terms object is also an attribute of the
  # data frame
  attr(data, "terms") <- terms
  # add information to model fit
  fit$formula <- formula
  fit$terms <- terms
  fit$model <- data
  fit$tau <- tau
  # add the class and return the model fit
  class(fit) <- "rq"
  fit
}
