# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------


#' Set up a diagnostic plot of regression weights
#'
#' Extract the relevant information for a diagnostic plot of the regression
#' weights from (robust) mediation analysis.  This plot allows to easily detect
#' deviations from normality assumptions such as skewness or heavy tails.
#'
#' This function is used internally by \code{\link{weight_plot}()}.  It may
#' also be useful for users who want to produce a similar plot, but who want
#' more control over what information to display or how to display that
#' information.
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
  p_m <- length(m)
  # check outcome variable that defines which regression model to plot
  if (is.null(outcome)) outcome <- y
  else {
    if (!is.character(outcome) && length(outcome) == 1L) {
      stop("only one outcome variable allowed")
    }
    if (!(outcome %in% m || outcome == y)) {
      stop("outcome variable must be the dependent variable or a mediator")
    }
  }
  # extract selected regression model
  if (outcome == y) fit <- object$fit_ymx
  else if (length(m) == 1) fit <- object$fit_mx
  else fit <- object$fit_mx[[m]]
  # extract residuals and weights
  residuals <- residuals(fit)
  weights <- weights(fit, type = "robustness")
  n <- length(weights)
  # extract control parameters of MM-estimator
  psi <- fit$control$psi
  tuning <- fit$control$tuning.psi
  # define grid of thresholds for weights
  thresholds <- seq(0, 1, length.out = npoints)
  # compute expected number of observations with weights below threshold
  # (for positive errors only due to symmetry)
  expected <- sapply(thresholds, function(threshold) {
    # compute point where weight function is equal to threshold
    x <- tuning * sqrt(1 - sqrt(threshold))
    # compute the probability that we exceed this point in a standard normal
    # distribution
    pnorm(x, lower.tail = FALSE)
  })
  # compute empirical percentages of weights smaller than threshold
  # (for negative and positive residuals separately)
  # FIXME: this should be done as a proper step function (thresholds given
  # by weights, as those are the only points where something happens)
  tails <- c("negative", "positive")
  tmp <- lapply(tails, function(tail) {
    if (tail == "negative") {
      in_tail <- residuals <= 0
    } else {
      in_tail <- residuals > 0
    }
    # subset of residuals and weights
    r <- residuals[in_tail]
    w <- weights[in_tail]
    # compute percentages of weights below thresholds
    empirical <- sapply(thresholds, function(threshold) {
      sum(w <= threshold) / n
    })
    # return data frame
    if (tail == "negative") {
      rbind(
        data.frame(Tail = "Negative residuals", Type = "Empirical",
                   Weight = thresholds, Percentage = empirical),
        data.frame(Tail = "Negative residuals", Type = "Expected",
                   Weight = thresholds, Percentage = expected)
      )
    } else {
      rbind(
        data.frame(Tail = "Positive residuals", Type = "Empirical",
                   Weight = thresholds, Percentage = empirical),
        data.frame(Tail = "Positive residuals", Type = "Expected",
                   Weight = thresholds, Percentage = expected)
      )
    }
  })
  data <- do.call(rbind, tmp)
  # return object
  out <- list(data = data, outcome = outcome)
  class(out) <- "setup_weight_plot"
  out
}


#' @rdname setup_weight_plot
#' @method setup_weight_plot cov_fit_mediation
#' @export

setup_weight_plot.cov_fit_mediation <- function(object, ...) {
  stop("not yet implemented for robust covariance fits")
}


## utility functions
