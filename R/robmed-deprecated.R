# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

#' Deprecated functions in package \pkg{robmed}
#'
#' These functions are provided for compatibility with older versions only, and
#' may be defunct as soon as the next release.
#'
#' @name robmed-deprecated
#'
#' @param x  either a numeric vector containing the independent variable, or
#' (if \code{data} is supplied) a character string, an integer or a logical
#' vector specifying the corresponding column of \code{data}.
#' @param y  either a numeric vector containing the dependent variable, or
#' (if \code{data} is supplied) a character string, an integer or a logical
#' vector specifying the corresponding column of \code{data}.
#' @param m  either a numeric vector or data frame containing the proposed
#' mediator variables, or (if \code{data} is supplied) a character, integer or
#' logical vector specifying the corresponding columns of \code{data}.
#' @param covariates  optional; either a numeric vector or data frame
#' containing additional covariates to be used as control variables, or (if
#' \code{data} is supplied) a character, integer or logical vector specifying
#' the corresponding columns of \code{data}.
#' @param data  an optional \code{data.frame} containing the variables.
#' @param efficiency  a numeric value giving the desired efficiency (defaults
#' to 0.85 for 85\% efficiency).
#' @param maxIterations  an integer giving the maximum number of iterations in
#' iterative algorithms.
#' @param tol  a small positive numeric value to be used to determine
#' convergence in various parts of the algorithm.
#' @param prob  numeric; probability for the quantile of the
#' \eqn{\chi^{2}}{chi-squared} distribution to be used as cutoff point in the
#' Huber weight function (defaults to 0.95).
#' @param seed  optional initial seed for the random number generator (see
#' \code{\link{.Random.seed}}).
#' @param \dots  arguments to be passed down
#'
#' @author Andreas Alfons
#'
#' @seealso \code{\link[base]{Deprecated}}
#'
#' \code{\link{fit_mediation}}, \code{\link{test_mediation}},
#' \code{\link{plot_mediation}}
#'
#' \code{\link{reg_control}}, \code{\link{cov_control}}
#'
#' \code{\link{cov_Huber}}, \code{\link{cov_ML}}
#'
#' @keywords multivariate

NULL


#' @rdname robmed-deprecated
#' @export

fitMediation <- function(x, y, m, covariates = NULL, data, ...) {
  .Deprecated("fit_mediation")
  # initializations
  if(missing(data)) {
    # prepare data frame containing all variables with original names
    x <- substitute(x)
    y <- substitute(y)
    p_m <- ncol(m)
    if(is.null(p_m)) p_m <- 1
    m <- substitute(m)
    if(is.null(covariates)) {
      data <- eval.parent(call("data.frame", x, y, m))
    } else if(is.data.frame(covariates)) {
      data <- cbind(eval.parent(call("data.frame", x, y, m)), covariates)
    } else {
      covariates <- substitute(covariates)
      data <- eval.parent(call("data.frame", x, y, m, covariates))
    }
  } else {
    # prepare data set
    data <- as.data.frame(data)
    x <- data[, x, drop=FALSE]
    y <- data[, y, drop=FALSE]
    m <- data[, m, drop=FALSE]
    p_m <- ncol(m)
    covariates <- data[, covariates, drop=FALSE]
    data <- cbind(x, y, m, covariates)
  }
  # extract names
  cn <- names(data)
  x <- cn[1]
  y <- cn[2]
  m <- cn[2L + seq_len(p_m)]
  covariates <- cn[-(seq_len(2L + p_m))]
  # call new function
  fit_mediation(data, x = x, y = y, m = m, covariates = covariates, ...)
}


#' @rdname robmed-deprecated
#' @export

regControl <- function(efficiency = 0.85, maxIterations = 200,
                       tol = 1e-07, seed = NULL) {
  .Deprecated("reg_control")
  reg_control(efficiency = efficiency, max_iterations = maxIterations,
              tol = tol, seed = seed)
}


#' @rdname robmed-deprecated
#' @export

covControl <- function(prob = 0.95, maxIterations = 200, tol = 1e-07) {
  .Deprecated("cov_control")
  cov_control(prob = prob, max_iterations = maxIterations, tol = tol)
}


#' @rdname robmed-deprecated
#' @export

covHuber <- function(...) {
  .Deprecated("cov_Huber")
  cov_Huber(...)
}


#' @rdname robmed-deprecated
#' @export

covML <- function(...) {
  .Deprecated("cov_ML")
  cov_ML(...)
}


#' @rdname robmed-deprecated
#' @export

testMediation <- function(x, ...) {
  .Deprecated("test_mediation")
  UseMethod("testMediation")
}


#' @rdname robmed-deprecated
#' @method testMediation default
#' @export

testMediation.default <- function(x, y, m, covariates = NULL, data, ...) {
  # initializations
  if(missing(data)) {
    # prepare data frame containing all variables with original names
    x <- substitute(x)
    y <- substitute(y)
    p_m <- ncol(m)
    if(is.null(p_m)) p_m <- 1
    m <- substitute(m)
    if(is.null(covariates)) {
      data <- eval.parent(call("data.frame", x, y, m))
    } else if(is.data.frame(covariates)) {
      data <- cbind(eval.parent(call("data.frame", x, y, m)), covariates)
    } else {
      covariates <- substitute(covariates)
      data <- eval.parent(call("data.frame", x, y, m, covariates))
    }
  } else {
    # prepare data set
    data <- as.data.frame(data)
    x <- data[, x, drop=FALSE]
    y <- data[, y, drop=FALSE]
    m <- data[, m, drop=FALSE]
    p_m <- ncol(m)
    covariates <- data[, covariates, drop=FALSE]
    data <- cbind(x, y, m, covariates)
  }
  # extract names
  cn <- names(data)
  x <- cn[1]
  y <- cn[2]
  m <- cn[2L + seq_len(p_m)]
  covariates <- cn[-(seq_len(2L + p_m))]
  # call new function
  test_mediation(data, x = x, y = y, m = m, covariates = covariates, ...)
}


#' @rdname robmed-deprecated
#' @method testMediation fit_mediation
#' @export

testMediation.fit_mediation <- function(x, ...) test_mediation(x, ...)


#' @rdname robmed-deprecated
#' @export

plotMediation <- function(...) {
  .Deprecated("plot_mediation")
  plot_mediation(...)
}
