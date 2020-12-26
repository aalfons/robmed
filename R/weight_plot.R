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
#' @return An object of class \code{"\link[ggplot2]{ggplot}"}.
#'
#' @author Andreas Alfons
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
#' # run fast and robust bootstrap test
#' robust_boot <- test_mediation(BSG2014,
#'                               x = "ValueDiversity",
#'                               y = "TeamCommitment",
#'                               m = "TaskConflict",
#'                               robust = TRUE)
#'
#' # create plot for robust bootstrap test
#' weight_plot(robust_boot)
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
  # build plot and convert to grob
  g1 <- ggplot_gtable(ggplot_build(p1))
  g2 <- ggplot_gtable(ggplot_build(p2))
  # # replace the panels from g1 with the ones from g2
  # g1$grobs[grep("panel-2-1", g1$layout$name)] <-
  #   g2$grobs[grep("panel-2-1", g2$layout$name)]
  # # also replace the x-axis corresponding to those panels
  # g1$grobs[grep('axis-b-2', g1$layout$name)] <-
  #   g2$grobs[grep('axis-b-2', g2$layout$name)]
  # replace the panels from g1 with the ones from g2
  g1$grobs[find_right_panels(g1$layout)] <-
    g2$grobs[find_right_panels(g2$layout)]
  # also replace the x-axis corresponding to those panels
  g1$grobs[find_right_axis(g1$layout)] <-
    g2$grobs[find_right_axis(g2$layout)]
  # set up grid graphics system
  if (newpage) grid.newpage()
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
