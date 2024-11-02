# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------


#' Set up a diagnostic plot with a tolerance ellipse
#'
#' Extract the relevant information for a diagnostic plot with a tolerance
#' ellipse from results of (robust) mediation analysis.
#'
#' This function is used internally by \code{\link{ellipse_plot}()}.  It may
#' also be useful for users who want to produce a similar plot, but who want
#' more control over what information to display or how to display that
#' information.
#'
#' @param object  an object inheriting from class \code{"\link{fit_mediation}"}
#' or \code{"\link{test_mediation}"} containing results from (robust) mediation
#' analysis, or a list of such objects.
#' @param horizontal  a character string specifying the variable to be
#' plotted on the horizontal axis.  If the dependent variable is chosen for
#' the vertical axis, a hypothsized mediator or an independent variable must
#' be selected for the horizontal axis.  If a hypothesized mediator is chosen
#' for the vertical axis, an independent variable must be selected for the
#' horizontal axis (in case of a serial multiple mediator model, a hypothesized
#' mediator occurring earlier in the sequence is also allowed).  The default is
#' to plot the first independent variable on the horizontal axis.
#' @param vertical  a character string specifying the variable to be
#' plotted on the vertical axis: the dependent variable or a hypothesized
#' mediator.  The default is to plot the first hypothesized mediator on the
#' vertical axis.
#' @param partial  a logical indicating whether to extract the observed values
#' of the selected variable for the vertical axis (\code{FALSE}), or the
#' partial residuals with respect to the variable on the horizontal axis
#' (\code{TRUE}).  The latter allows to display the corresponding regression
#' coefficient by a line.
#' @param level  numeric; the confidence level of the tolerance ellipse.  It
#' gives the percentage of observations that are expected to lie within the
#' ellipse under the assumption of a normal distribution, and therefore it
#' controls the size of the ellipse.  The default is such that the ellipse is
#' expected to contain 97.5\% of the observations.
#' @param npoints  the number of grid points used to evaluate the ellipse.  The
#' default is to use 100 grid points.
#' @param \dots  additional arguments to be passed down.
#'
#' @return An object of class \code{"setup_ellipse_plot"} with the following
#' components:
#' \item{data}{a data frame containing the coordinates of the data points to
#' be plotted on the horizontal axis (column \code{x}) and the coordinates
#' on the vertical axis (column \code{y}).  For robust methods that assign
#' outlyingness weights to each data point, those weights are given in column
#' \code{Weight}.  If a list of objects has been supplied and there are
#' multiple objects from such robust methods, or if partial residuals are to be
#' plotted on the vertical axis, there is also a column \code{Method}, which
#' takes the names or indices of the list elements to indicate the different
#' methods.}
#' \item{ellipse}{a data frame containing the coordinates of the tolerance
#' ellipse on the horizontal axis (column \code{x}) and on the vertical axis
#' (column \code{y}).  If a list of objects has been supplied, there is also a
#' column \code{Method}, which takes the names or indices of the list elements
#' to indicate the different methods.}
#' \item{line}{a data frame with columns \code{intercept} and \code{slope}
#' containing the intercept and slope, respectively, of the regression line to
#' be plotted.  If  a list of objects has been supplied, there is also a column
#' \code{Method}, which takes the names or indices of the list elements to
#' indicate the different methods. This is only returned if
#' \code{partial = TRUE}, or in case of a simple mediation model (without
#' control variables) when the hypothesized mediator is plotted on the vertical
#' axis and the independent variable is plotted on the horizontal axis.}
#' \item{horizontal}{a character string giving the variable to be plotted on
#' the horizontal axis.}
#' \item{vertical}{a character string giving the variable to be plotted on the
#' vertical axis}
#' \item{partial}{a logical indicating whether the values to be plotted on the
#' vertical axis correspond to the observed values of the selected variable
#' (\code{FALSE}), or the partial residuals with respect to the variable on the
#' horizontal axis (\code{TRUE}).}
#' \item{robust}{a logical indicating whether the object contains results from
#' a robust method, or a vector of such logicals if a list of objects has been
#' supplied.}
#' \item{have_methods}{a logical indicating whether a list of objects has been
#' supplied.}
#'
#' @author Andreas Alfons
#'
#' @references
#' Alfons, A., Ates, N.Y. and Groenen, P.J.F. (2022) Robust Mediation Analysis:
#' The \R Package \pkg{robmed}.  \emph{Journal of Statistical Software},
#' \bold{103}(13), 1--45.  doi:10.18637/jss.v103.i13.
#'
#' @seealso
#' \code{\link{fit_mediation}()}, \code{\link{test_mediation}()},
#' \code{\link{ellipse_plot}()}
#'
#' @examples
#' data("BSG2014")
#'
#' # fit mediation model
#' fit <- fit_mediation(BSG2014,
#'                      x = "ValueDiversity",
#'                      y = "TeamCommitment",
#'                      m = "TaskConflict")
#'
#' # set up information for plot
#' setup <- setup_ellipse_plot(fit)
#'
#' # plot only data and tolerance ellipse
#' ggplot() +
#'   geom_path(aes(x = x, y = y), data = setup$ellipse,
#'             color = "#00BFC4") +
#'   geom_point(aes(x = x, y = y, fill = Weight),
#'              data = setup$data, shape = 21) +
#'   scale_fill_gradient(limits = 0:1, low = "white",
#'                       high = "black") +
#'   labs(x = setup$horizontal, y = setup$vertical)
#'
#' @keywords hplot
#'
#' @import ggplot2
#' @export

setup_ellipse_plot <- function(object, ...) UseMethod("setup_ellipse_plot")


#' @rdname setup_ellipse_plot
#' @method setup_ellipse_plot test_mediation
#' @export

setup_ellipse_plot.test_mediation <- function(object, ...) {
  # simply call methods for mediation model fit
  setup_ellipse_plot(object$fit, ...)
}


#' @rdname setup_ellipse_plot
#' @method setup_ellipse_plot reg_fit_mediation
#' @export

setup_ellipse_plot.reg_fit_mediation <- function(object,
                                                 horizontal = NULL,
                                                 vertical = NULL,
                                                 partial = FALSE,
                                                 level = 0.975,
                                                 npoints = 100,
                                                 ...) {
  # initializations
  have_robust <- is_robust(object)
  if (have_robust && object$robust == "median") {
    # weights for weighted covariance matrix can be infinite for median
    # regression (if there is a residual that is exactly 0)
    stop("tolerance ellipse not meaningful for median regression")
  }
  if (!have_robust && object$family != "gaussian") {
    stop("tolerance ellipse only implemented for normal distribution")
  }
  # extract variable names
  x <- object$x
  p_x <- length(x)
  y <- object$y
  m <- object$m
  p_m <- length(m)
  covariates <- object$covariates
  p_covariates <- length(covariates)
  model <- object$model
  # check variable on vertical axis
  if (is.null(vertical)) vertical <- m[1L]
  else {
    if (!(is.character(vertical) && length(vertical) == 1L)) {
      stop("only one variable allowed for the vertical axis")
    }
    if (!(vertical %in% m || vertical == y)) {
      stop("variable on the vertical axis must be ",
           "the dependent variable or a mediator")
    }
  }
  # check variable on horizontal axis
  if (is.null(horizontal)) horizontal <- x[1L]
  else {
    if (!(is.character(horizontal) && length(horizontal) == 1L)) {
      stop("only one variable allowed for the horizontal axis")
    }
    # check horizontal axis if a mediator is on vertical axis
    if (model == "serial") {
      # serial multiple mediators: independent variable and mediators earlier
      # in the sequence can be on horizontal axis
      if (vertical %in% m) {
        pos <- which(vertical == m)
        if (pos == 1L && !horizontal %in% x) {
          stop("variable on the horizontal axis must be an independent variable")
        }
        if (pos > 1L && !(horizontal %in% c(x, m[seq_len(pos-1L)]))) {
          stop("variable on the horizontal axis must be an independent variable ",
               "or a mediator occuring earlier in the serial mediator model")
        }
      }
    } else {
      # single mediator or parallel multiple mediators: only independent
      # variable can be on horizontal axis
      if (vertical %in% m && !(horizontal %in% x)) {
        stop("variable on the horizontal axis must be an independent variable")
      }
    }
    # check horizontal axis if dependent variable is on vertical axis
    if (vertical == y && !(horizontal %in% m || horizontal %in% x)) {
      stop("variable on the horizontal axis must be ",
           "an independent variable or a mediator")
    }
  }
  # other initializations
  partial <- isTRUE(partial)
  vertical_mx <- if (model == "serial") vertical == m[1L] else vertical %in% m
  have_mx <- vertical_mx && p_x == 1L && horizontal == x && p_covariates == 0L
  robust <- object$robust == "MM"
  # extract model fit
  if (partial || have_mx || robust) {
    if (vertical %in% m) {
      fit <- if (p_m > 1L) object$fit_mx[[vertical]] else object$fit_mx
    } else fit <- object$fit_ymx
    coefficients <- coef(fit)
  }
  # if applicable, extract residuals
  if (partial || robust) residuals <- residuals(fit)
  # if applicable, extract intercept and slope
  if (partial || have_mx) {
    if (have_mx && !partial) intercept <- unname(coefficients["(Intercept)"])
    else intercept <- 0
    slope <- unname(coefficients[horizontal])
    line <- data.frame(intercept = intercept, slope = slope)
  } else line <- NULL
  # extract data to plot
  if (partial) {
    x <- object$data[, horizontal]
    y <- residuals + slope * x
    data <- data.frame(x = x, y = y)
  } else {
    data <- data.frame(x = object$data[, horizontal],
                       y = object$data[, vertical])
  }
  # obtain location and shape of ellipse
  if (robust) {
    # The weighted covariance matrix with weights from robust regression is
    # underestimated in the direction of the residuals, but the submatrix that
    # involves only the explanatory variables is correctly estimated under the
    # model.  We therefore start with the corrected variance of the residuals
    # and the weighted covariance matrix of the predictor variables, and then
    # transform to obtain the variance of the response variable and the
    # covariances of the response with the explanatory variables.
    # -----
    # # compute correction factor for variance of the residuals
    # control <- fit$control
    # # the following integrals compute the result at the model
    # integrand1 <- function(y) {
    #   Mwgt(y, cc=control$tuning.psi, psi=control$psi) * y^2 * dnorm(y)
    # }
    # numerator <- integrate(integrand1, lower = -10, upper = 10,
    #                        rel.tol = .Machine$double.eps^0.5)
    # integrand2 <- function(y) {
    #   Mwgt(y, cc=control$tuning.psi, psi=control$psi) * dnorm(y)
    # }
    # denominator <- integrate(integrand2, lower = -10, upper = 10,
    #                          rel.tol = .Machine$double.eps^0.5)
    # # the correction factor is the inverse of the result at the model
    # correction <- denominator$value / numerator$value
    # -----
    # the computation of the correction factor is commented out since we
    # can use the residual scale from the "lmrob" object, which is already
    # corrected for downweighting observations
    # -----
    # extract weights in case of robust regression
    w <- weights(fit, type = "robustness")
    if (partial) {
      # compute center
      center_x <- weighted.mean(data$x, w = w)
      center_y <- slope * center_x
      center <- c(x = center_x, y = center_y)
      # compute covariance matrix
      cov_xx <- weighted.var(data$x, w = w, center = center_x)
      cov_xy <- slope * cov_xx
      cov_yy <- slope^2 * cov_xx + fit$scale^2
      cov <- cbind(rbind(cov_xx, cov_xy), rbind(cov_xy, cov_yy))
      dimnames(cov) <- replicate(2, c("x", "y"), simplify = FALSE)
    } else {
      # compute weighted mean and weighted covariance matrix of all explanatory
      # variables, as this is necessary for proper correction
      if (vertical %in% m) {
        if (model == "serial") {
          pos <- which(vertical == m)
          predictors <- c(m[seq_len(pos-1)], x, covariates)
        } else predictors <- c(x, covariates)
      } else predictors <- c(m, x, covariates)
      predictor_data <- object$data[, predictors, drop = FALSE]
      m_x <- sapply(predictor_data, weighted.mean, w = w)
      S_xx <- weighted.cov(predictor_data, w = w, center = m_x)
      # extract regression coefficients
      alpha_hat <- coefficients[1L]
      beta_hat <- coefficients[-1L]
      # reconstruct center estimate
      center_x <- unname(m_x[horizontal])
      center_y <- alpha_hat + crossprod(m_x, beta_hat)
      center <- c(x = center_x, y = center_y)
      # reconstruct estimate of covariance matrix
      cov_xx <- S_xx[horizontal, horizontal, drop = FALSE]
      cov_xy <- S_xx[horizontal, , drop = FALSE] %*% beta_hat
      cov_yy <- t(beta_hat) %*% S_xx %*% beta_hat + fit$scale^2
      cov <- cbind(rbind(cov_xx, cov_xy), rbind(cov_xy, cov_yy))
      dimnames(cov) <- replicate(2, c("x", "y"), simplify = FALSE)
    }
    # add weights to data frame
    data$Weight <- w
  } else {
    # compute mean and covariance matrix
    center <- colMeans(data)
    cov <- cov(data)
  }
  # compute ellipse
  ellipse <- ellipse(center, cov, level = level, npoints = npoints)
  # return data and ellipse
  out <- list(data = data, ellipse = as.data.frame(ellipse), line = line,
              horizontal = horizontal, vertical = vertical, partial = partial,
              robust = robust, have_methods = FALSE)
  class(out) <- "setup_ellipse_plot"
  out
}


#' @rdname setup_ellipse_plot
#' @method setup_ellipse_plot cov_fit_mediation
#' @export

setup_ellipse_plot.cov_fit_mediation <- function(object,
                                                 horizontal = NULL,
                                                 vertical = NULL,
                                                 partial = FALSE,
                                                 level = 0.975,
                                                 npoints = 100,
                                                 ...) {
  # extract variable names
  x <- object$x
  y <- object$y
  m <- object$m
  # check variable on vertical axis
  if (is.null(vertical)) vertical <- m
  else {
    if (!(is.character(vertical) && length(vertical) == 1L)) {
      stop("only one variable allowed for the vertical axis")
    }
    if (vertical != m && vertical != y) {
      stop("variable on the vertical axis must be ",
           "the dependent variable or the mediator")
    }
  }
  # check variable on horizontal axis
  if (is.null(horizontal)) horizontal <- x
  else {
    if (!(is.character(horizontal) && length(horizontal) == 1L)) {
      stop("only one variable allowed for the horizontal axis")
    }
    if (vertical == m && horizontal != x) {
      stop("variable on the horizontal axis must be the independent variable")
    } else if (vertical == y && (horizontal != m && horizontal != x)) {
      stop("variable on the horizontal axis must be ",
           "the independent variable or a mediator")
    }
  }
  # other initializations
  partial <- isTRUE(partial)
  select <- c(horizontal, vertical)
  # extract covariance fit
  fit <- object$cov
  # extract information to be plotted
  if (vertical == m) {
    # extract location and shape of ellipse
    center <- fit$center[select]
    cov <- fit$cov[select, select]
    # compute intercept and slope
    slope <- unname(cov[vertical, horizontal] / cov[horizontal, horizontal])
    intercept <- unname(center[vertical] - slope * center[horizontal])
    # extract data to plot
    data <- data.frame(x = object$data[, horizontal],
                       y = object$data[, vertical])
    # if requested, make correction for partial residual plot
    if (partial) {
      data$y <- data$y - intercept
      center[vertical] <- center[vertical] - intercept
      intercept <- 0
    }
    # combine information for line representing (partial) effect
    line <- data.frame(intercept = intercept, slope = slope)
  } else {
    # dependent variable on the vertical axis
    if (partial) {
      # extract all means and full covariance matrix
      center <- fit$center
      cov <- fit$cov
      # compute regression coefficient
      det <- cov[x, x] * cov[m, m] - cov[m, x]^2
      b <- (-cov[m, x] * cov[y, x] + cov[x, x] * cov[y, m]) / det
      direct <- (cov[m, m] * cov[y, x] - cov[m, x] * cov[y, m]) / det
      i_y <- unname(center[y] - b * center[m] - direct * center[x])
      # compute data to plot
      if (horizontal == x) {
        # compute partial residuals
        residuals <- object$data[, vertical] - i_y - b * object$data[, m]
        # adjust location and shape of ellipse
        center[vertical] <- center[vertical] - i_y - b * center[m]
        cov[y, y] <- cov[y, y] + b^2 * cov[m, m] + 2 * b * cov[y, m]
        cov[y, x] <- cov[y, x] - b * cov[m, x]
      } else {
        # compute partial residuals
        residuals <- object$data[, vertical] - i_y - direct * object$data[, x]
        # adjust location and shape of ellipse
        center[vertical] <- center[vertical] - i_y - direct * center[x]
        cov[y, y] <- cov[y, y] + direct^2 * cov[x, x] + 2 * direct * cov[y, x]
        cov[y, m] <- cov[y, m] - direct * cov[x, m]
      }
      data <- data.frame(x = object$data[, horizontal], y = residuals)
      # extract location and shape of ellipse
      center <- center[select]
      cov <- cov[select, select]
      # compute intercept and slope
      intercept <- 0
      slope <- if (horizontal == x) unname(direct) else unname(b)
      # combine information for line representing (partial) effect
      line <- data.frame(intercept = intercept, slope = slope)
    } else {
      # extract location and shape of ellipse
      center <- fit$center[select]
      cov <- fit$cov[select, select]
      # extract data to plot
      data <- data.frame(x = object$data[, horizontal],
                         y = object$data[, vertical])
      # do not plot line for effect
      line <- NULL
    }
  }
  # in case of robust covariance matrix, add interpretable robustness weights
  if (object$robust) data$Weight <- weights(fit, type = "relative")
  # compute ellipse
  ellipse <- ellipse(center, cov, level = level, npoints = npoints)
  colnames(ellipse) <- c("x", "y")
  # return data and ellipse
  out <- list(data = data, ellipse = as.data.frame(ellipse), line = line,
              horizontal = horizontal, vertical = vertical, partial = partial,
              robust = object$robust, have_methods = FALSE)
  class(out) <- "setup_ellipse_plot"
  out
}


#' @rdname setup_ellipse_plot
#' @method setup_ellipse_plot list
#' @export

setup_ellipse_plot.list <- function(object, ...) {
  # initializations
  is_test <- sapply(object, inherits, "test_mediation")
  is_fit <- sapply(object, inherits, "fit_mediation")
  object <- object[is_test | is_fit]
  if (length(object) == 0L) {
    stop('no objects inheriting from class "test_mediation" or "fit_mediation"')
  }
  # check that variables are the same
  components <- c("x", "y", "m", "covariates")
  variables <- lapply(object, function(x) x$fit[components])
  all_identical <- all(sapply(variables[-1L], identical, variables[[1L]]))
  if (!isTRUE(all_identical)) {
    stop("all mediation objects must use the same variables")
  }
  # check names of list elements
  methods <- names(object)
  if (is.null(methods)) methods <- seq_along(object)
  else {
    replace <- methods == "" | duplicated(methods)
    methods[replace] <- seq_along(object)[replace]
  }
  # compute tolerance ellipse for each list element
  tol_ellipse_list <- lapply(object, setup_ellipse_plot, ...)
  # check for properties of tolerance ellipses
  partial <- tol_ellipse_list[[1]]$partial
  robust <- sapply(tol_ellipse_list, "[[", "robust")
  sum_robust <- sum(robust)
  # combine information of tolerance ellipses
  # TODO: should there be an argument to define what to return?
  if (partial || sum_robust > 1L) {
    # combine data sets
    # For a partial residuals plot, those partial residuals are different for
    # every method.  Otherwise if there are multiple robust methods, there are
    # different weights.
    if (sum_robust == 0L) {
      # add column identifying the method to data for scatter plots
      data_list <- mapply(function(ellipse, method) {
        data.frame(Method = method, ellipse$data, stringsAsFactors = TRUE)
      }, ellipse = tol_ellipse_list, method = methods,
      SIMPLIFY = FALSE, USE.NAMES = FALSE)
    } else {
      # if there is a robust method, add column with weights for nonrobust
      # methods as well
      data_list <- mapply(function(ellipse, method) {
        if (ellipse$robust) {
          data.frame(Method = method, ellipse$data,
                     stringsAsFactors = TRUE)
        } else {
          data.frame(Method = method, ellipse$data, Weight = 1,
                     stringsAsFactors = TRUE)
        }
      }, ellipse = tol_ellipse_list, method = methods,
      SIMPLIFY = FALSE, USE.NAMES = FALSE)
    }
    # combine data sets for scatter plots
    data <- do.call(rbind, data_list)
  } else if (sum_robust == 1L) data <- tol_ellipse_list[[which(robust)]]$data
  else data <- tol_ellipse_list[[1]]$data
  # combine data for ellipses
  ellipse_list <- mapply(function(ellipse, method) {
    data.frame(Method = method, ellipse$ellipse, stringsAsFactors = TRUE)
  }, ellipse = tol_ellipse_list, method = methods,
  SIMPLIFY = FALSE, USE.NAMES = FALSE)
  ellipses <- do.call(rbind, ellipse_list)
  # combine data for lines representing (partial) effects
  line_list <- mapply(function(ellipse, method) {
    line <- ellipse$line
    if (!is.null(line)) {
      data.frame(Method = method, line, stringsAsFactors = TRUE)
    }
  }, ellipse = tol_ellipse_list, method = methods,
  SIMPLIFY = FALSE, USE.NAMES = FALSE)
  lines <- do.call(rbind, line_list)
  # return results
  out <- list(data = data, ellipse = ellipses, line = lines,
              horizontal = tol_ellipse_list[[1]]$horizontal,
              vertical = tol_ellipse_list[[1]]$vertical,
              partial = partial, robust = robust,
              have_methods = TRUE)
  class(out) <- "setup_ellipse_plot"
  out
}


## workhorse function to compute ellipse based on center and covariance matrix
ellipse <- function(center, cov, level = 0.975, npoints = 100) {
  # extract scales and correlation
  scale <- sqrt(diag(cov))
  r <- cov[1, 2] / prod(scale)
  # compute ellipse
  d <- acos(r)
  a <- seq(0, 2 * pi, length.out = npoints)
  q <- sqrt(qchisq(level, df = 2))
  x <- q * scale[1] * cos(a + d/2) + center[1]
  y <- q * scale[2] * cos(a - d/2) + center[2]
  xy <- cbind(x, y)
  # add names of variables and return ellipse
  colnames(xy) <- names(center)
  xy
}


## utility functions

# weighted variance
weighted.var <- function(x, w, center = NULL, ..., na.rm = TRUE) {
  na.rm <- isTRUE(na.rm)
  if (missing(w)) var(x, na.rm = na.rm)
  else {
    # initial checks
    x <- as.numeric(x)
    w <- as.numeric(w)
    if (length(w) != length(x)) stop("'w' must have the same length as 'x'")
    if (is.null(center)) center <- weighted.mean(x, w = w, na.rm = na.rm)
    else if (length(center) != 1) stop("'center' must have length 1")
    # if requested, remove missing values
    if (na.rm) {
      select <- !is.na(x)
      x <- x[select]
      w <- w[select]
    }
    # denominator is chosen such that it reduces to unbiased estimator if all
    # weights are equal to 1
    if (length(x) <= 1 || sum(w > 0) <= 1) NA_real_
    else sum((x - center)^2 * w) / (sum(w) - 1)
  }
}

# weighted covariance matrix
weighted.cov <- function(x, w, center = NULL, ...) {
  x <- as.matrix(x)
  if (missing(w)) cov(x, use = "pairwise.complete.obs")
  else {
    # initial checks
    if (length(w) != nrow(x)) {
      stop("length of 'w' must equal the number of rows in 'x'")
    }
    if (is.null(center)) center <- apply(x, 2, weighted.mean, w = w)
    else if (length(center) != ncol(x)) {
      stop("length of 'center' must equal the number of columns in 'x'")
    }
    # sweep out center estimates and compute weighted cross product
    x <- sweep(x, 2, center, check.margin = FALSE)
    weighted.crossprod(x, w = w)
  }
}

# weighted cross product
weighted.crossprod <- function(x, w) {
  ci <- colnames(x)
  p <- ncol(x)
  if (p == 0) matrix(numeric(), nrow = 0, ncol = 0)
  else if (p == 1) {
    # remove missing values
    select <- !is.na(x)
    x <- x[select, ]
    w <- w[select]
    # denominator is chosen such that it reduces to unbiased estimator if all
    # weights are equal to 1
    matrix(sum(w * x^2) / (sum(w) - 1), nrow = 1, ncol = 1,
           dimnames = replicate(2, ci, simplify = FALSE))
  } else {
    if (is.null(ci)) ci <- seq_len(p)
    sapply(ci, function(j) sapply(ci, function(i) {
      # use pairwise complete observations
      select <- !is.na(x[, i]) & !is.na(x[, j])
      xi <- x[select, i]
      xj <- x[select, j]
      w <- w[select]
      # denominator is chosen such that it reduces to unbiased estimator if all
      # weights are equal to 1
      sum(xi * xj * w) / (sum(w) - 1)
    }))
  }
}
