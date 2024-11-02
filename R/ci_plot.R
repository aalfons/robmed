# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------


#' Dot plot with confidence intervals
#'
#' Produce a dot plot with confidence intervals of selected effects from
#' (robust) mediation analysis.  In addition to confidence intervals, p-values
#' of the selected effects can be plotted as well.
#'
#' Methods first call \code{\link{setup_ci_plot}()} to extract all necessary
#' information to produce the plot, then the \code{"setup_ci_plot"}
#' method is called to produce the plot.
#'
#' @param object  an object inheriting from class
#' \code{"\link{test_mediation}"} containing results from
#' (robust) mediation analysis, or a list of such objects.
#' @param parm  an integer, character or logical vector specifying which
#' effects to include in the plot.  In case of a character vector, possible
#' values are \code{"a"}, \code{"b"}, \code{"d"} (only serial multiple mediator
#' models), \code{"total"}, \code{"direct"}, and \code{"indirect"}.  The
#' default is to include the direct and the indirect effect(s).
#' @param type  a character string specifying which point estimates and
#' confidence intervals to plot: those based on the bootstrap distribution
#' (\code{"boot"}; the default), or those based on the original data
#' (\code{"data"}).  If \code{"boot"}, the confidence intervals of effects
#' other than the indirect effect(s) are computed using a normal approximation
#' (i.e., assuming a normal distribution of the corresponding effect with the
#' standard deviation computed from the bootstrap replicates).  If
#' \code{"data"}, the confidence intervals of effects other than the indirect
#' effect(s) are computed via statistical theory based on the original data
#' (e.g., based on a t-distribution if the coefficients are estimated via
#' regression).  Note that this is only relevant for mediation analysis via a
#' bootstrap test, where the confidence interval of the indirect effect is
#' always computed via a percentile-based method due to the asymmetry of its
#' distribution.
#' @param p_value  a logical indicating whether to include dot plots of the
#' p-values in addition to those with confidence intervals.  The default is
#' \code{FALSE}.
#' @param digits  an integer determining how many digits to compute for
#' bootstrap p-values of the indirect effects (see \code{\link{p_value}()}).
#' The default is to compute 4 digits after the comma.  This is only relevant
#' if \code{p_value = TRUE}.
#' @param level  numeric;  the confidence level of the confidence intervals
#' from Sobel's test.  The default is to include 95\% confidence intervals.
#' Note that this is not used for bootstrap tests, as those require to specify
#' the confidence level already in \code{\link{test_mediation}()}.
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
#' \code{\link{test_mediation}()}, \code{\link{setup_ci_plot}()}
#'
#' \code{\link{density_plot}()}, \code{\link{ellipse_plot}()},
#' \code{\link{weight_plot}()}, \code{\link[=plot-methods]{plot}()}
#'
#' @examples
#' data("BSG2014")
#'
#' # run fast-and-robust bootstrap test
#' robust_boot <- test_mediation(BSG2014,
#'                               x = "ValueDiversity",
#'                               y = "TeamCommitment",
#'                               m = "TaskConflict",
#'                               level = 0.9,
#'                               robust = TRUE)
#'
#' # create plot for robust bootstrap test
#' ci_plot(robust_boot)
#' ci_plot(robust_boot, color = "#00BFC4")
#'
#' # run OLS bootstrap test
#' ols_boot <- test_mediation(BSG2014,
#'                            x = "ValueDiversity",
#'                            y = "TeamCommitment",
#'                            m = "TaskConflict",
#'                            level = 0.9,
#'                            robust = FALSE)
#'
#' # compare robust and OLS bootstrap tests
#' boot_list <- list("OLS bootstrap" = ols_boot,
#'                   "ROBMED" = robust_boot)
#' ci_plot(boot_list)
#'
#' # the plot can be customized in the usual way
#' ci_plot(boot_list) +
#'   geom_hline(yintercept = 0, color = "darkgrey") +
#'   coord_flip() + theme_bw() +
#'   labs(title = "OLS bootstrap vs ROBMED")
#'
#' @keywords hplot
#'
#' @import ggplot2
#' @export

ci_plot <- function(object, ...) UseMethod("ci_plot")


#' @rdname ci_plot
#' @method ci_plot default
#' @export

ci_plot.default <- function(object, parm = c("direct", "indirect"), ...) {
  # extract information
  setup <- setup_ci_plot(object, parm = parm, ...)
  # call method for corresponding objects
  ci_plot(setup, ...)
}


#' @rdname ci_plot
#' @method ci_plot boot_test_mediation
#' @export

ci_plot.boot_test_mediation <- function(object, parm = c("direct", "indirect"),
                                        type = c("boot", "data"),
                                        p_value = FALSE, digits = 4L,
                                        ...) {
  # extract information
  setup <- setup_ci_plot(object, parm = parm, type = type, p_value = p_value,
                         digits = digits, ...)
  # call method for corresponding objects
  ci_plot(setup, ...)
}


#' @rdname ci_plot
#' @method ci_plot sobel_test_mediation
#' @export

ci_plot.sobel_test_mediation <- function(object, parm = c("direct", "indirect"),
                                         level = 0.95, p_value = FALSE, ...) {
  # extract information
  setup <- setup_ci_plot(object, parm = parm, level = level,
                         p_value = p_value, ...)
  # call method for corresponding objects
  ci_plot(setup, ...)
}


#' @rdname ci_plot
#' @method ci_plot list
#' @export

ci_plot.list <- function(object, parm = c("direct", "indirect"),
                         type = c("boot", "data"), level = 0.95,
                         p_value = FALSE, digits = 4L, ...) {
  # extract information
  setup <- setup_ci_plot(object, parm = parm, type = type, level = level,
                         p_value = p_value, digits = digits, ...)
  # call method for corresponding objects
  ci_plot(setup, ...)
}


#' @rdname ci_plot
#' @method ci_plot setup_ci_plot
#' @export

ci_plot.setup_ci_plot <- function(object, ...) {
  # initializations
  have_methods <- object$have_methods
  have_p_value <- !is.null(object$p_value)
  # define aesthetic mapping for confidence interval
  if (have_methods) {
    mapping_ci <- aes_string(x = "Method", y = "Estimate",
                          ymin="Lower", ymax="Upper")
  } else {
    mapping_ci <- aes_string(x = "\"\"", y = "Estimate",
                          ymin="Lower", ymax="Upper")
  }
  # if requested, plot p-values in the bottom row
  if (have_p_value) {
    # define mapping
    if (have_methods) mapping_p_value <- aes_string(x = "Method", y = "Value")
    else mapping_p_value <- aes_string(x = "\"\"", y = "Value")
    # generate plot
    p <- ggplot() +
      geom_point_p_value(mapping_p_value, data = object$p_value, ...) +
      geom_pointrange_ci(mapping_ci, data = object$ci, ...)
  } else {
    # generate plot
    p <- ggplot() +
      geom_pointrange_ci(mapping_ci, data = object$ci, ...)
  }
  # remove annotation
  p <- p + labs(title = NULL, x = NULL, y = NULL)
  # add facets
  if (have_p_value) {
    p <- p + facet_wrap(~ Label + Effect, nrow = 2, scales = "free_y")
  } else p <- p + facet_wrap(~ Effect, scales = "free")
  # return plot
  p
}


## custom geom for confidence intervals: ignore certain arguments
geom_pointrange_ci <- function(..., stat, position, orientation) {
  geom_pointrange(...)
}

## custom geom for p_values:
# 1) handle size of points in the same way as geom_pointrange()
# 2) avoid passing unknown arguments
geom_point_p_value <- function(..., fatten = 4, linetype,
                               lty, lwd, orientation) {
  # extract argument names
  arguments <- list(...)
  argument_names <- names(arguments)
  # replace argument names with standardized ones
  standardized_names <- standardise_aes_names(argument_names)
  names(arguments) <- standardized_names
  # check size of points
  size <- arguments$size
  if (is.null(size)) size <- 0.5
  arguments$size <- fatten * size
  # check stroke of points (i.e., border thickness for point shapes 21-25)
  stroke <- arguments$stroke
  if (is.null(stroke)) stroke <- 1
  arguments$stroke <- stroke
  # call geom_point()
  do.call(geom_point, arguments)
}
