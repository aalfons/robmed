# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------


#' Density plot of the indirect effect(s)
#'
#' Produce a density plot of the indirect effect(s) from (robust) mediation
#' analysis.  In addition to the density, a vertical line representing the
#' point estimate and a shaded area representing the confidence interval are
#' drawn.
#'
#' Methods first call \code{\link{setup_density_plot}()} to extract all
#' necessary information to produce the plot, then the
#' \code{"setup_density_plot"} method is called to produce the plot.
#'
#' @param object  an object inheriting from class
#' \code{"\link{test_mediation}"} containing results from
#' (robust) mediation analysis, or a list of such objects.
#' @param grid  an optional numeric vector containing the values at which to
#' evaluate the assumed normal density from Sobel's test.  The default is to
#' take 512 equally spaced points between the estimated indirect effect
#' \eqn{\pm}{+/-} three times the standard error according to Sobel's formula.
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
#' \code{\link{test_mediation}()}, \code{\link{setup_density_plot}()}
#'
#' \code{\link{ci_plot}()}, \code{\link{ellipse_plot}()},
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
#' density_plot(robust_boot)
#' density_plot(robust_boot, color = "#00BFC4", fill = "#00BFC4")
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
#' density_plot(boot_list)
#'
#' # the plot can be customized in the usual way
#' density_plot(boot_list) + theme_bw() +
#'   labs(title = "OLS bootstrap vs ROBMED")
#'
#' @keywords hplot
#'
#' @import ggplot2
#' @export

density_plot <- function(object, ...) UseMethod("density_plot")


#' @rdname density_plot
#' @method density_plot default
#' @export

density_plot.default <- function(object, ...) {
  # extract information
  setup <- setup_density_plot(object, ...)
  # call method for corresponding objects
  density_plot(setup, ...)
}


#' @rdname density_plot
#' @method density_plot sobel_test_mediation
#' @export

density_plot.sobel_test_mediation <- function(object, grid = NULL,
                                              level = 0.95, ...) {
  # extract information
  setup <- setup_density_plot(object, grid = grid, level = level, ...)
  # call method for corresponding objects
  density_plot(setup, ...)
}


#' @rdname density_plot
#' @method density_plot list
#' @export

density_plot.list <- function(object, grid = NULL, level = 0.95, ...) {
  # extract information
  setup <- setup_density_plot(object, grid = grid, level = level, ...)
  # call method for corresponding objects
  density_plot(setup, ...)
}


#' @rdname density_plot
#' @method density_plot setup_density_plot
#' @export

density_plot.setup_density_plot <- function(object, ...) {
  # define aesthetic mappings for density estimate, point estimate and
  # confidence interval
  if(object$have_methods) {
    mapping_density <- aes_string(x = "Indirect", y = "Density",
                                  color = "Method")
    mapping_estimate <- aes_string(xintercept = "Estimate", color = "Method")
    mapping_ci <- aes_string(xmin = "Lower", xmax = "Upper",
                             ymin = -Inf, ymax = Inf,
                             fill = "Method")
  } else {
    mapping_density <- aes_string(x = "Indirect", y = "Density")
    mapping_estimate <- aes_string(xintercept = "Estimate")
    mapping_ci <- aes_string(xmin = "Lower", xmax = "Upper",
                             ymin = -Inf, ymax = Inf)
  }
  # define default title
  if (all(object$test == "boot")) title <- "Bootstrap distribution"
  else if (all(object$test == "sobel")) title <- "Assumed normal distribution"
  else title <- NULL
  # generate plot
  p <- ggplot() +
    geom_density_indirect(mapping_density, data = object$density, ...) +
    geom_vline_indirect(mapping_estimate, data = object$ci, ...) +
    geom_rect_ci(mapping_ci, data = object$ci, ...) +
    labs(title = title, x = "Indirect effect", y = "Density")
  # split plot into different panels in case of multiple indirect effects
  if(object$have_effects) p <- p + facet_wrap(~ Effect, scales = "free")
  # return plot
  p
}


## custom geom for density estimate to be used in density plot
#  1) always use stat = "identity" because the density is already estimated
#  2) do not allow for a fill color because a filled rectangle is used to
#     display the confidence interval
geom_density_indirect <- function(..., stat, fill, bg, alpha) {
  geom_density(..., stat = "identity")
}

## custom geom for vertical line to be used in density plot:
#  1) avoid passing argument 'alpha' to ensure that line is of the same style
#     as lines of density
#  2) avoid passing unknown argument 'fill'
geom_vline_indirect <- function(..., fill, bg, alpha) geom_vline(...)

## custom geom for confidence intervals to be used in density plot:
#  fix transparant rectangle without edges and avoid duplication of arguments
geom_rect_ci <- function(...) {
  # extract argument names
  arguments <- list(...)
  argument_names <- names(arguments)
  # replace argument names with standardized ones
  standardized_names <- standardise_aes_names(argument_names)
  names(arguments) <- standardized_names
  # make sure that there is no border
  arguments$colour <- NA
  # use default transparency if not specified otherwise
  if (is.null(arguments$alpha)) arguments$alpha <- 0.2
  # call existing geom function
  do.call(geom_rect, arguments)
}
