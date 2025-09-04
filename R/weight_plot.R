# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------


#' Diagnostic plot of robust regression weights
#'
#' Produce a diagnostic plot of the regression weights from robust mediation
#' analysis.  This plot allows to easily detect deviations from normality
#' assumptions such as skewness or heavy tails.
#'
#' The horizontal axis contains different weight thresholds, and the vertical
#' axis displays the percentage of observations that have a weight below this
#' threshold.  For comparison, a reference line is drawn for the expected
#' percentages under normally distributed errors.  Observations with negative
#' and positive residuals are shown separately to make it possible to
#' distinguish between symmetric and asymmetric deviations from normality.
#'
#' If the plot reveals more downweighted observations than expected, but
#' roughly the same amounts in both tails, the residual distribution is
#' symmetric but with heavy tails.  If the plot shows that observations in one
#' tail are downweighted more heavily than those in the other tail, the
#' residual distribution is skewed.
#'
#' @param object  an object inheriting from class \code{"\link{fit_mediation}"}
#' or \code{"\link{test_mediation}"} containing results from robust mediation
#' analysis.  Only mediation analysis objects fitted with the robust
#' MM-estimator are supported.
#' @param outcome  a character vector specifying the outcome variables of the
#' regressions to be included in the plot.  This must be a subset of the
#' hypothesized mediators and the dependent variable, or \code{NULL} (the
#' default) to include all regressions of the mediation model.
#' @param npoints  the number of grid points used to evaluate and draw the
#' expected percentages.  The default is to use 1000 grid points.
#' @param \dots  additional arguments to be passed down.
#'
#' @return An object inheriting from class \code{"\link[ggplot2]{ggplot}"}.
#'
#' @note The current implementation is a slight hack of \pkg{ggplot2} and the
#' \pkg{grid} graphics system in order to revert the horizontal axis only in
#' the panels for observations with postive residuals.  It is therefore not
#' possible to change the horizontal axis with
#' \code{\link[ggplot2]{scale_x_continuous}()}.
#'
#' The implementation may change in the future if the required functionality
#' becomes available in \pkg{ggplot2}.
#'
#' @author Andreas Alfons
#'
#' @references
#' Alfons, A., Ates, N.Y. and Groenen, P.J.F. (2022) Robust Mediation Analysis:
#' The \R Package \pkg{robmed}.  \emph{Journal of Statistical Software},
#' \bold{103}(13), 1--45.  doi:10.18637/jss.v103.i13.
#'
#' Alfons, A. and Schley, D.R. (2025) \emph{Robust Mediation Analysis: What We
#' Talk About When We Talk About Robustness}.  PsyArXiV,
#' doi:10.31234/osf.io/2hqdy.
#'
#' @seealso
#' \code{\link{fit_mediation}()}, \code{\link{test_mediation}()},
#' \code{\link{setup_weight_plot}()}
#'
#' \code{\link{ci_plot}()}, \code{\link{density_plot}()},
#' \code{\link{ellipse_plot}()}, \code{\link[=plot-methods]{plot}()}
#'
#' @examples
#' data("BSG2014")
#'
#' # obtain robust fit of mediation model
#' fit <- fit_mediation(BSG2014,
#'                      x = "ValueDiversity",
#'                      y = "TeamCommitment",
#'                      m = "TaskConflict")
#'
#' # create diagnostic plot of robust regression weights
#' weight_plot(fit) +
#'   scale_color_manual("", values = c("black", "#00BFC4")) +
#'   theme(legend.position = "top")
#'
#' # plot only the regression model for the hypothesized mediator
#' weight_plot(fit, outcome = "TaskConflict") +
#'   scale_color_manual("", values = c("black", "#00BFC4")) +
#'   theme(legend.position = "top")
#'
#' @keywords hplot
#'
#' @import ggplot2
#' @export

weight_plot <- function(object, ...) UseMethod("weight_plot")


#' @rdname weight_plot
#' @method weight_plot default
#' @export

weight_plot.default <- function(object, outcome = NULL, npoints = 1000, ...) {
  # compute percentages of observations with weights below thresholds
  setup <- setup_weight_plot(object, outcome = outcome, npoints = npoints)
  # call method for resulting object
  weight_plot(setup, ...)
}


#' @rdname weight_plot
#' @method weight_plot setup_weight_plot
#' @export

weight_plot.setup_weight_plot <- function(object, ...) {
  # initializations
  outcome <- object$outcome
  p_outcomes <- length(outcome)
  xlab <- "Weight threshold"
  ylab <- "Percentage of observations with weight lower than threshold"
  # create plot
  # For the vertical line to separate the panels, the default backround color
  # of the panel strips in ggplot2 is "gray85".  But that is not dark enough
  # and barely visible.  So "darkgray" is used instead.
  p <- ggplot() +
    geom_vline(xintercept = 1, color = "darkgray") +
    geom_line(aes_string(x = "Threshold", y = "Percentage", color = "Weights"),
              data = object$data)
  if (p_outcomes == 1L) {
    p <- p +
      facet_grid(. ~ Tail) +
      labs(title = paste("Outcome variable:", outcome), x = xlab, y = ylab)
  } else {
    p <- p +
      facet_grid(Outcome ~ Tail) +
      labs(x = xlab, y = ylab)
  }
  # dirty hack: add new class with print() method to reverse x-axis in one panel
  class(p) <- c("gg_weight_plot", class(p))
  p
}


## dirty hack: print() method to reverse x-axis in panel for positive residuals
#' @importFrom grid grid.newpage grid.draw seekViewport pushViewport upViewport
#' @export
print.gg_weight_plot <- function (x, newpage = is.null(vp), vp = NULL, ...) {
  # make sure that plot is registered to be fetched by last_plot()
  set_last_plot(x)
  # enforce that there is no horizontal space between panels
  p <- x + theme(panel.spacing.x = unit(0, "points"))
  # prepare two versions of the plot with different x-axis
  p1 <- p + scale_x_continuous(expand = expansion(mult = c(0.05, 0)),
                               trans = "identity")
  p2 <- p + scale_x_continuous(expand = expansion(mult = c(0, 0.05)),
                               trans = "reverse")
  # set up grid graphics system (this must be done before the grobs are
  # created, otherwise an empty page will be created before the plot)
  if (newpage) grid.newpage()
  # build plot and convert to grob
  g1 <- ggplot_gtable(ggplot_build(p1))
  g2 <- ggplot_gtable(ggplot_build(p2))
  # replace the panels from g1 with the ones from g2
  g1$grobs[find_right_panels(g1$layout)] <-
    g2$grobs[find_right_panels(g2$layout)]
  # also replace the x-axis corresponding to those panels
  g1$grobs[find_right_axis(g1$layout)] <-
    g2$grobs[find_right_axis(g2$layout)]
  # draw plot
  if (is.null(vp)) grid.draw(g1)
  else {
    if (is.character(vp)) seekViewport(vp)
    else pushViewport(vp)
    grid.draw(g1)
    upViewport()
  }
  # return plot invisibly
  invisible(x)
}


## utility functions

#find right hand side panels in plot layout
find_right_panels <- function(layout) {
  panels <- grep("panel", layout$name)
  r <- layout[panels, "r"]
  panels[r == max(r)]
}

#find right hand side bottom axis in plot layout
find_right_axis <- function(layout) {
  axes <- grep("axis-b", layout$name)
  r <- layout[axes, "r"]
  axes[r == max(r)]
}
