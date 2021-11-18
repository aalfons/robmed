# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

#' Plot (robust) mediation analysis results
#'
#' Visualize results from (robust) mediation analysis.
#'
#' The \code{"\link{fit_mediation}"} method calls \code{\link{ellipse_plot}()}
#' or \code{\link{weight_plot}()}, depending on the argument \code{which}.
#'
#' The \code{"\link{test_mediation}"} method calls \code{\link{ci_plot}()},
#' \code{\link{density_plot}()}, \code{\link{ellipse_plot}()}, or
#' \code{\link{weight_plot}()}, depending on the argument \code{which}.
#'
#' @name plot-methods
#'
#' @param object,x  an object inheriting from class
#' \code{"\link{fit_mediation}"} or \code{"\link{test_mediation}"} containing
#' results from (robust) mediation analysis.
#' @param which  a character string specifying which plot to produce.
#' Possible values are \code{"ci"} for a dot plot of selected effects
#' together with confidence intervals (see \code{\link{ci_plot}()}),
#' \code{"density"} for a density plot of the indirect effect(s) (see
#' \code{\link{density_plot}()}), \code{"ellipse"} for a diagnostic plot
#' of the data together with a tolerance ellipse (see
#' \code{\link{ellipse_plot}()}), and \code{"weight"} for a diagnostic plot
#' of robust regression weights (see \code{\link{weight_plot}()}).
#' @param \dots  additional arguments to be passed down.
#'
#' @return An object of class \code{"\link[ggplot2]{ggplot}"}.
#'
#' @author Andreas Alfons
#'
#' @seealso
#' \code{\link{fit_mediation}()}, \code{\link{test_mediation}()}
#'
#' \code{\link{ci_plot}()}, \code{\link{density_plot}()},
#' \code{\link{ellipse_plot}()}, \code{\link{weight_plot}()}
#'
#' @examples
#' data("BSG2014")
#'
#' # run fast-and-robust bootstrap test
#' test <- test_mediation(BSG2014,
#'                        x = "ValueDiversity",
#'                        y = "TeamCommitment",
#'                        m = "TaskConflict")
#'
#' # create plots for robust bootstrap test
#' plot(test, which = "ci")
#' plot(test, which = "density")
#' plot(test, which = "ellipse")
#' plot(test, which = "weight")
#'
#' @keywords hplot
#'
#' @import ggplot2

NULL


## internal function for plotting mediation fits
plot_internal_fit <- function(object, which = c("ellipse", "weight"), ...) {
  # initializations
  which <- match.arg(which)
  # call selected plot function
  if (which == "ellipse") ellipse_plot(object, ...)
  else if (which == "weight") weight_plot(object, ...)
  else stop("type of plot not implemented")  # shouldn't happen
}

## internal function for plotting mediation tests
plot_internal_test <- function(object,
                               which = c("ci", "density", "ellipse", "weight"),
                               ...) {
  # initializations
  which <- match.arg(which)
  # call selected plot function
  if (which == "ci") ci_plot(object, ...)
  else if (which == "density") density_plot(object, ...)
  else if (which == "ellipse") ellipse_plot(object, ...)
  else if (which == "weight") weight_plot(object, ...)
  else stop("type of plot not implemented")  # shouldn't happen
}


#' @rdname plot-methods
#' @method autoplot fit_mediation
#' @export

autoplot.fit_mediation <- function(object, which = c("ellipse", "weight"),
                                   ...) {
  plot_internal_fit(object, which = which, ...)
}


#' @rdname plot-methods
#' @method autoplot test_mediation
#' @export

autoplot.test_mediation <- function(object,
                                    which = c("ci", "density", "ellipse", "weight"),
                                    ...) {
  plot_internal_test(object, which = which, ...)
}


#' @rdname plot-methods
#' @method plot fit_mediation
#' @export

plot.fit_mediation <- function(x, which = c("ellipse", "weight"), ...) {
  plot_internal_fit(x, which = which, ...)
}


#' @rdname plot-methods
#' @method plot test_mediation
#' @export

plot.test_mediation <- function(x,
                                which = c("ci", "density", "ellipse", "weight"),
                                ...) {
  plot_internal_test(x, which = which, ...)
}
