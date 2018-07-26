# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

#' Plot (robust) mediation analysis results
#'
#' Produce dot plots of selected coefficients from regression models computed
#' in (robust) mediation analysis, or density plots of the indirect effect.
#'
#' @param object,x  an object inheriting from class
#' \code{"\link{test_mediation}"} containing results from (robust) mediation
#' analysis.  For \code{plot_mediation}, a list of such objects may be supplied
#' as well.
#' @param data  an optional numeric vector containing the \eqn{x}-values at
#' which to evaluate the assumed normal density from Sobel's test (only used in
#' case of a density plot).  The default is to take 100 equally spaced points
#' between the estimated indirect effect \eqn{\pm}{+/-} three times the
#' standard error according to Sobel's formula.
#' @param method  a character string specifying which plot to produce.
#' Possible values are \code{"dot"} for a dot plot of selected coefficients, or
#' \code{"density"} for a density plot of the indirect effect(s).
#' @param parm  a character string specifying the coefficients to be included
#' in a dot plot.  The default is to include the direct and the indirect
#' effect(s).
#' @param level  numeric;  the confidence level of the confidence intervals
#' from Sobel's test to be included in a dot plot.  The default is to include
#' 95\% confidence intervals.
#' @param mapping  an aesthetic mapping to override the default behavior (see
#' \code{\link[ggplot2]{aes}} or \code{\link[ggplot2]{aes_}}).
#' @param facets  a faceting formula to override the default behavior (only
#' used in case of a dot plot).  If supplied, \code{\link[ggplot2]{facet_wrap}}
#' or \code{\link[ggplot2]{facet_grid}} is called depending on whether the
#' formula is one-sided or two-sided.
#' @param \dots  additional arguments to be passed to and from methods.
#'
#' @return An object of class \code{"ggplot"} (see
#' \code{\link[ggplot2]{ggplot}}).
#'
#' @author Andreas Alfons
#'
#' @seealso \code{\link{test_mediation}},
#' \code{\link[=fortify.test_mediation]{fortify}}
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
#' # create plots for robust bootstrap test
#' plot(robust_boot, method = "dot")
#' plot(robust_boot, method = "density")
#'
#' # run standard bootstrap test
#' standard_boot <- test_mediation(BSG2014,
#'                                 x = "ValueDiversity",
#'                                 y = "TeamCommitment",
#'                                 m = "TaskConflict",
#'                                 robust = FALSE)
#'
#' # compare robust and standard tests
#' tests <- list(Robust = robust_boot, Standard = standard_boot)
#' plot_mediation(tests, method = "dot")
#' plot_mediation(tests, method = "density")
#'
#' @keywords hplot
#'
#' @import ggplot2
#' @export

plot_mediation <- function(object, ...) UseMethod("plot_mediation")


#' @rdname plot_mediation
#' @method plot_mediation boot_test_mediation
#' @export

plot_mediation.boot_test_mediation <- function(object,
                                               method = c("dot", "density"),
                                               parm = NULL, ...) {
  data <- fortify(object, method=method, parm=parm)
  plot_mediation(data, ...)
}


#' @rdname plot_mediation
#' @method plot_mediation sobel_test_mediation
#' @export

plot_mediation.sobel_test_mediation <- function(object, data,
                                                method = c("dot", "density"),
                                                parm = c("c", "ab"),
                                                level = 0.95, ...) {
  data <- fortify(object, data=data, method=method, parm=parm, level=level)
  plot_mediation(data, ...)
}


#' @rdname plot_mediation
#' @method plot_mediation list
#' @export

plot_mediation.list <- function(object, data, method = c("dot", "density"),
                                parm = NULL, level = 0.95, ...) {
  data <- fortify(object, data=data, method=method, parm=parm, level=level)
  plot_mediation(data, ...)
}


#' @rdname plot_mediation
#' @method plot_mediation default
#' @export

plot_mediation.default <- function(object, mapping = attr(object, "mapping"),
                                   facets = attr(object, "facets"), ...) {
  # create selected plot
  if(attr(object, "method") == "dot") dot_plot(object, mapping, facets, ...)
  else density_plot(object, mapping, facets, ...)
}


#' @rdname plot_mediation
#' @method autoplot test_mediation
#' @export

autoplot.test_mediation <- function(object, ...) plot_mediation(object, ...)


#' @rdname plot_mediation
#' @method plot test_mediation
#' @export

plot.test_mediation <- function(x, ...) plot_mediation(x, ...)


## internal function for dot plot
dot_plot <- function(data, mapping, facets, main = NULL,
                     xlab = NULL, ylab = NULL, ...) {
  # generate plot
  geom <- attr(data, "geom")
  p <- ggplot(data, mapping) + geom(...) + labs(title=main, x=xlab, y=ylab)
  if(!is.null(facets)) {
    # split plot into different panels
    if(length(facets) == 2) p <- p + facet_wrap(facets)
    else p <- p + facet_grid(facets)
  }
  p
}


## internal function for density plot
density_plot <- function(data, mapping, facets, main = NULL,
                         xlab = NULL, ylab = NULL, ...) {
  # define default title and axis labels
  if(is.null(main)) main <- attr(data, "main")
  if(is.null(xlab)) xlab <- "Indirect effect"
  if(is.null(ylab)) ylab <- "Density"
  # extract point estimate and confidence interval
  ci <- attr(data, "ci")
  if("Method" %in% names(data)) {
    mapping_line <- aes_string(xintercept = "ab", color = "Method")
    mapping_rect <- aes_string(xmin = "Lower", xmax = "Upper",
                               ymin = -Inf, ymax = Inf,
                               fill = "Method")
  } else {
    mapping_line <- aes_string(xintercept = "ab")
    mapping_rect <- aes_string(xmin = "Lower", xmax = "Upper",
                               ymin = -Inf, ymax = Inf)
  }
  # generate plot
  geom <- attr(data, "geom")
  p <- ggplot(data, mapping) + geom(...) +
    geom_vline(mapping_line, data = ci, ...) +
    geom_rect(mapping_rect, data = ci, color = NA, alpha = 0.2, ...) +
    labs(title = main, x = xlab, y = ylab)
  if(!is.null(facets)) {
    # split plot into different panels
    if(length(facets) == 2L) p <- p + facet_wrap(facets, scales = "free")
    else p <- p + facet_grid(facets, scales = "free")
  }
  p
}
