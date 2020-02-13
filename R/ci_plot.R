# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------


#' Dot plot with confidence intervals
#'
#' Produce a dot plot with confidence intervals of selected effects from
#' (robust) mediation analysis.
#'
#' Methods first call \code{\link{setup_ci_plot}()} to extract all necessary
#' information to produce the plot, then the \code{"setup_ci_plot"}
#' method is called to produce the plot.
#'
#' @param object  an object inheriting from class
#' \code{"\link{test_mediation}"} containing results from
#' (robust) mediation analysis, or a list of such objects.
#' @param parm  a character string specifying the effects to be included
#' in the plot.  The default is to include the direct and the indirect
#' effect(s).
#' @param type  a character string specifying which point estiamates and
#' confidence intervals to plot: those based on the bootstrap distribution
#' (\code{"boot"}; the default), or those based on the original data
#' (\code{"data"}).  If \code{"boot"}, the confidence intervals of effects
#' other than the indirect effect(s) are computed using a normal approximation
#' (i.e., assuming a normal distribution of the corresponding effect with the
#' standard deviation computed from the bootstrap replicates).  If
#' \code{"data"}, the confidence intervals of effects other than the indirect
#' effect(s) are computed via statistical theory based on the original data
#' (e.g., based on a t-distribution the coefficients are estimated via
#' regression).  Note that this is only relevant for mediation analysis via a
#' bootstrap test, where the confidence interval of the indirect effect is
#' always computed via a percentile-based method due to the asymmetry of its
#' distribution.
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
#' @seealso
#' \code{\link{test_mediation}()}, \code{\link{setup_ci_plot}()}
#'
#' \code{\link{density_plot}()}, \code{\link{ellipse_plot}()},
#' \code{\link[=plot-methods]{plot}()}
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
#' ci_plot(robust_boot)
#' ci_plot(robust_boot, color = "#00BFC4")
#'
#' # run standard bootstrap test
#' standard_boot <- test_mediation(BSG2014,
#'                                 x = "ValueDiversity",
#'                                 y = "TeamCommitment",
#'                                 m = "TaskConflict",
#'                                 robust = FALSE)
#'
#' # compare robust and standard tests
#' tests <- list(Standard = standard_boot, Robust = robust_boot)
#' ci_plot(tests)
#'
#' # the plot can be customized in the usual way
#' ci_plot(tests) +
#'   geom_hline(yintercept = 0, color = "darkgrey") +
#'   coord_flip() + theme_bw() +
#'   labs(title = "Standard vs robust bootstrap test")
#'
#' @keywords hplot
#'
#' @import ggplot2
#' @export

ci_plot <- function(object, ...) UseMethod("ci_plot")


#' @rdname ci_plot
#' @method ci_plot default
#' @export

ci_plot.default <- function(object, parm = NULL, ...) {
  # extract information
  setup <- setup_ci_plot(object, parm = parm, ...)
  # call method for corresponding objects
  ci_plot(setup, ...)
}


#' @rdname ci_plot
#' @method ci_plot boot_test_mediation
#' @export

ci_plot.boot_test_mediation <- function(object, parm = NULL,
                                        type = c("boot", "data"),
                                        ...) {
  # extract information
  setup <- setup_ci_plot(object, parm = parm, type = type, ...)
  # call method for corresponding objects
  ci_plot(setup, ...)
}


#' @rdname ci_plot
#' @method ci_plot sobel_test_mediation
#' @export

ci_plot.sobel_test_mediation <- function(object, parm = NULL,
                                         level = 0.95, ...) {
  # extract information
  setup <- setup_ci_plot(object, parm = parm, level = level, ...)
  # call method for corresponding objects
  ci_plot(setup, ...)
}


#' @rdname ci_plot
#' @method ci_plot list
#' @export

ci_plot.list <- function(object, parm = NULL, type = c("boot", "data"),
                         level = 0.95, ...) {
  # extract information
  setup <- setup_ci_plot(object, parm = parm, type = type, level = level, ...)
  # call method for corresponding objects
  ci_plot(setup, ...)
}


#' @rdname ci_plot
#' @method ci_plot setup_ci_plot
#' @export

ci_plot.setup_ci_plot <- function(object, ...) {
  # define aesthetic mapping for confidence interval
  if (object$have_methods) {
    mapping <- aes_string(x = "Method", y = "Estimate",
                          ymin="Lower", ymax="Upper")
  } else {
    mapping <- aes_string(x = "\"\"", y = "Estimate",
                          ymin="Lower", ymax="Upper")
  }
  # generate plot
  ggplot(object$ci, mapping) +
    geom_pointrange(...) +
    labs(title = NULL, x = NULL, y = NULL) +
    facet_wrap(~ Effect, scales = "free")
}