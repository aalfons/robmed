# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------


#' Set up a diagnostic plot of robust regression weights
#'
#' Extract the relevant information for a diagnostic plot of the regression
#' weights from robust mediation analysis.  This plot allows to easily detect
#' deviations from normality assumptions such as skewness or heavy tails.
#'
#' This function is used internally by \code{\link{weight_plot}()}.  It may
#' also be useful for users who want to produce a similar plot, but who want
#' more control over what information to display or how to display that
#' information.
#'
#' @param object  an object inheriting from class \code{"\link{fit_mediation}"}
#' or \code{"\link{test_mediation}"} containing results from robust mediation
#' analysis.  Only mediation analysis objects fitted with the robust
#' MM-estimator are supported.
#' @param outcome  a character vector specifying the outcome variables of the
#' regressions to be included in the plot.  This must be a subset of the
#' hypothesized mediators and the dependent variable, or \code{NULL} (the
#' default) to include all regressions of the mediation model.
#' @param npoints  the number of grid points used to evaluate the expected
#' percentages.  The default is to use 1000 grid points.
#' @param \dots  additional arguments to be passed down.
#'
#' @return An object of class \code{"setup_weight_plot"} with the following
#' components:
#' \item{data}{a data frame containing the following information: the outcome
#' variable of the regression (column \code{Outcome}; only if multiple
#' regressions are to be included in the plot), whether the row corresponds to
#' the negative or the positive tail of the residual distribution (column
#' \code{Tail}), whether the row corresponds to the expected (under the normal
#' distribution) or the empirical weights (column \code{Weights}), the weight
#' thresholds (column \code{Threshold}), and the corresponding percentage of
#' observations that have a weight below this threshold (column
#' \code{Percentage}).}
#' \item{outcome}{a character vector containing the outcome variables of the
#' regressions to be included in the plot.}
#'
#' @author Andreas Alfons
#'
#' @seealso
#' \code{\link{fit_mediation}()}, \code{\link{test_mediation}()},
#' \code{\link{weight_plot}()}
#'
#' @examples
#' data("BSG2014")
#'
#' # run fast and robust bootstrap test
#' robust_boot <- test_mediation(BSG2014,
#'                               x = "ValueDiversity",
#'                               y = "TeamCommitment",
#'                               m = "TaskConflict",
#'                               robust = TRUE)
#'
#' # set up information for plot
#' setup <- setup_weight_plot(robust_boot)
#' # create diagnostic plot of robust regression weights
#' weight_plot(setup) +
#'   scale_color_manual("", values = c("black", "#00BFC4")) +
#'   theme(legend.position = "top")
#'
#' @keywords hplot
#'
#' @import ggplot2
#' @export

setup_weight_plot <- function(object, ...) UseMethod("setup_weight_plot")


#' @rdname setup_weight_plot
#' @method setup_weight_plot test_mediation
#' @export

setup_weight_plot.test_mediation <- function(object, ...) {
  # simply call methods for mediation model fit
  setup_weight_plot(object$fit, ...)
}


#' @rdname setup_weight_plot
#' @method setup_weight_plot reg_fit_mediation
#' @export

setup_weight_plot.reg_fit_mediation <- function(object,
                                                outcome = NULL,
                                                npoints = 1000,
                                                ...) {
  # initializations
  have_robust <- is_robust(object)
  if (!(have_robust && object$robust == "MM")) {
    stop("weight plot only meaningful for MM-regression")
  }
  # extract relevant variable names
  y <- object$y
  m <- object$m
  # check outcome variable that defines which regression model to plot
  if (is.null(outcome)) outcome <- c(m, y)
  else {
    if (!(is.character(outcome) && length(outcome) > 0L)) {
      stop("outcome variables must be specified as character strings")
    }
    if (!all(outcome %in% c(m, y))) {
      stop("outcome variables must be the dependent variable or a mediator")
    }
  }
  # call workhose function
  if (length(outcome) == 1L) {
    data <- get_weight_percentages(object, outcome = outcome, npoints = npoints)
  } else {
    tmp <- lapply(outcome, function(current) {
      current_data <- get_weight_percentages(object, outcome = current,
                                             npoints = npoints)
      cbind(Outcome = current, current_data)
    })
    data <- do.call(rbind, tmp)
  }
  # return object
  out <- list(data = data, outcome = outcome)
  class(out) <- "setup_weight_plot"
  out
}

# workhorse function to get the relevant data frame for a single regression fit
get_weight_percentages <- function(object, outcome, npoints = 1000) {
  # extract relevant variable names
  y <- object$y
  m <- object$m
  p_m <- length(m)
  # extract selected regression model
  if (outcome == y) fit <- object$fit_ymx
  else if (p_m == 1L) fit <- object$fit_mx
  else fit <- object$fit_mx[[outcome]]
  # extract residuals and weights
  residuals <- residuals(fit)
  weights <- weights(fit, type = "robustness")
  n <- length(weights)
  # extract control parameters of MM-estimator
  psi <- fit$control$psi
  tuning <- fit$control$tuning.psi
  # define grid of thresholds for weights
  thresholds <- seq(0, 1, length.out = npoints)
  # compute expected percentages of observations with weights below thresholds
  # (for positive errors only due to symmetry)
  expected <- sapply(thresholds, function(threshold) {
    # compute point where weight function is equal to threshold
    x <- tuning * sqrt(1 - sqrt(threshold))
    # compute the probability that we exceed this point in a standard normal
    # distribution
    pnorm(x, lower.tail = FALSE)
  })
  # compute empirical percentages of observations with weights below thresholds
  # (for negative and positive residuals separately)
  # FIXME: this should be done as a proper step function (thresholds given
  # by weights, as those are the only points where something happens)
  tails <- c("negative", "positive")
  out_list <- lapply(tails, function(tail) {
    if (tail == "negative") {
      in_tail <- residuals <= 0
    } else {
      in_tail <- residuals > 0
    }
    # subset of residuals and weights
    r <- residuals[in_tail]
    w <- weights[in_tail]
    # compute percentages of observations with weights below thresholds
    empirical <- sapply(thresholds, function(threshold) {
      sum(w <= threshold) / n
    })
    # return data frame
    if (tail == "negative") {
      rbind(
        data.frame(Tail = "Negative residuals", Weights = "Expected",
                   Threshold = thresholds, Percentage = expected),
        data.frame(Tail = "Negative residuals", Weights = "Empirical",
                   Threshold = thresholds, Percentage = empirical)
      )
    } else {
      rbind(
        data.frame(Tail = "Positive residuals", Weights = "Expected",
                   Threshold = thresholds, Percentage = expected),
        data.frame(Tail = "Positive residuals", Weights = "Empirical",
                   Threshold = thresholds, Percentage = empirical)
      )
    }
  })
  # combine everything into one data frame
  do.call(rbind, out_list)
}

# @rdname setup_weight_plot
# @method setup_weight_plot cov_fit_mediation
# @export

# setup_weight_plot.cov_fit_mediation <- function(object, ...) {
#   stop("not yet implemented for robust covariance fits")
# }


## utility functions
