# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------


#' Diagnostic plot with a tolerance ellipse
#'
#' Produce a scatter plot of two variables used in (robust) mediation analysis
#' together with a tolerance ellipse.  Exploiting the relationship between the
#' regression coefficients and the covariance matrix, that tolerance ellipse
#' illustrates how well the regression results represent the data. In addition,
#' a line that visualizes the estimated regression coefficient is added when
#' relevant.
#'
#' A line to visualize the corresponding regression coefficient is added if
#' \code{partial = TRUE}, or in case of a simple mediation model
#' (without control variables) when the hypothesized mediator is plotted on
#' the vertical axis and the independent variable is plotted on the horizontal
#' axis.
#'
#' For robust estimation methods that return outlyingness weights for each
#' data point, those weights are visualized by coloring the points on a grey
#' scale.  If a list of objects has been supplied and there are multiple
#' objects from such robust methods, each method is placed in a separate panel.
#'
#' Methods first call \code{\link{setup_ellipse_plot}()} to extract all
#' necessary information to produce the plot, then the
#' \code{"setup_ellipse_plot"} method is called to produce the plot.
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
#' @param partial  a logical indicating whether the vertical axis should
#' display the observed values of the selected variable (\code{FALSE}), or
#' the partial residuals with respect to the variable on the horizontal axis
#' (\code{TRUE}).  The latter allows to display the corresponding regression
#' coefficient by a line.
#' @param level  numeric; the confidence level of the tolerance ellipse.  It
#' gives the percentage of observations that are expected to lie within the
#' ellipse under the assumption of a normal distribution, and therefore it
#' controls the size of the ellipse.  The default is such that the ellipse is
#' expected to contain 97.5\% of the observations.
#' @param npoints  the number of grid points used to evaluate and draw the
#' ellipse.  The default is to use 100 grid points.
#' @param \dots  additional arguments to be passed down.
#'
#' @return An object of class \code{"\link[ggplot2]{ggplot}"}.
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
#' \code{\link{setup_ellipse_plot}()}
#'
#' \code{\link{ci_plot}()}, \code{\link{density_plot}()},
#' \code{\link{weight_plot}()}, \code{\link[=plot-methods]{plot}()}
#'
#' @examples
#' data("BSG2014")
#'
#' # obtain robust fit of mediation model
#' robust_fit <- fit_mediation(BSG2014,
#'                             x = "ValueDiversity",
#'                             y = "TeamCommitment",
#'                             m = "TaskConflict",
#'                             robust = TRUE)
#'
#' # create plot for robust model fit
#' ellipse_plot(robust_fit)
#'
#' # original data and partial residuals
#' ellipse_plot(robust_fit, horizontal = "TaskConflict",
#'              vertical = "TeamCommitment")
#' ellipse_plot(robust_fit, horizontal = "TaskConflict",
#'              vertical = "TeamCommitment", partial = TRUE)
#'
#' # obtain OLS fit of mediation model
#' ols_fit <- fit_mediation(BSG2014,
#'                          x = "ValueDiversity",
#'                          y = "TeamCommitment",
#'                          m = "TaskConflict",
#'                          robust = FALSE)
#'
#' # compare robust and OLS model fits
#' fit_list <- list("OLS bootstrap" = ols_fit,
#'                  "ROBMED" = robust_fit)
#' ellipse_plot(fit_list)
#'
#' # the plot can be customized in the usual way
#' ellipse_plot(fit_list) + theme_bw() +
#'   labs(title = "OLS vs robust estimation")
#'
#' @keywords hplot
#'
#' @import ggplot2
#' @export

ellipse_plot <- function(object, ...) UseMethod("ellipse_plot")


#' @rdname ellipse_plot
#' @method ellipse_plot default
#' @export

ellipse_plot.default <- function(object, horizontal = NULL, vertical = NULL,
                                 partial = FALSE, level = 0.975, npoints = 100,
                                 ...) {
  # compute tolerance ellipse
  setup <- setup_ellipse_plot(object, horizontal = horizontal,
                              vertical = vertical, partial = partial,
                              level = level, npoints = npoints)
  # call method for tolerance ellipse objects
  ellipse_plot(setup, ...)
}


#' @rdname ellipse_plot
#' @method ellipse_plot setup_ellipse_plot
#' @export

ellipse_plot.setup_ellipse_plot <- function(object, ...) {
  # initializations
  robust <- any(object$robust)
  have_line <- !is.null(object$line)
  have_methods <- object$have_methods
  tmp <- "Method" %in% names(object$data)
  use_color <- have_methods && !tmp
  use_facets <- have_methods && tmp
  # define aesthetic mapping for plotting ellipses and lines
  if (use_color) {
    aes_ellipse <- aes_string(x = "x", y = "y", color = "Method")
    if (have_line) {
      aes_line <- aes_string(intercept = "intercept", slope = "slope",
                             color = "Method")
    }
  } else {
    aes_ellipse <- aes_string(x = "x", y = "y")
    if (have_line) {
      aes_line <- aes_string(intercept = "intercept", slope = "slope")
    }
  }
  # define aesthetic mapping for plotting points
  if (robust) {
    aes_data <- aes_string(x = "x", y = "y", fill = "Weight")
  } else aes_data <- aes_string(x = "x", y = "y")
  # create plot
  p <- ggplot() +
    geom_path_ellipse(aes_ellipse, data = object$ellipse, ...) +
    geom_point_data(aes_data, data = object$data, ...)
  # add line representing (partial) effect
  if (have_line) {
    p <- p +
      geom_abline_effect(aes_line, data = object$line, ...)
  }
  # add nice labels
  if (object$partial) ylab <- paste("Partial residuals of", object$vertical)
  else ylab <- object$vertical
  p <- p + labs(x = object$horizontal, y = ylab)
  # add color gradient for weights
  if (robust) {
    p <- p + scale_fill_gradient(limits = 0:1, low = "white", high = "black")
  }
  # add facets in case of multiple methods
  if (use_facets) p <- p + facet_wrap(~ Method)
  # return plot
  p
}


## custom geom for data: avoid passing unknown argument 'linetype'
geom_point_data <- function(..., linetype, lty, lwd) {
  # extract argument names
  arguments <- list(...)
  argument_names <- names(arguments)
  # replace argument names with standardized ones
  standardized_names <- standardise_aes_names(argument_names)
  names(arguments) <- standardized_names
  # check plot symbol
  shape <- arguments$shape
  if (is.null(shape)) arguments$shape <- 21
  else if (!isTRUE(shape %in% 21:25)) {
    arguments$shape <- 21
    warning("only plot symbols 21 to 25 allowed, using default symbol 21",
            call. = FALSE)
  }
  # use black fill color as default for nonrobust fits
  # (if there are robust fits, the fill color is specified in the aesthetic
  # mapping, otherwise it would be possible that a user overrides the fill
  # color by passing additional arguments via ...)
  aes <- arguments[[1]]
  if (is.null(aes$fill) && is.null(arguments$fill)) arguments$fill <- "black"
  # call existing geom function
  do.call(geom_point, arguments)
}

## custom geom for ellipse: avoid passing unknown argument 'fill'
geom_path_ellipse <- function(..., fill, bg, shape, pch, cex) geom_path(...)

## custom geom for (partial) effect: avoid passing unknown argument 'fill'
geom_abline_effect <- function(..., fill, bg, shape, pch, cex,
                               show.legend = FALSE) {
  geom_abline(..., show.legend = show.legend)
}
