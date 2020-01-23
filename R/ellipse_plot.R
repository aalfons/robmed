# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------


#' @import ggplot2
#' @export
ellipse_plot <- function(object, ...) UseMethod("ellipse_plot")

#' @export
ellipse_plot.default <- function(object, horizontal = NULL, vertical = NULL,
                                 partial = FALSE, level = 0.975, npoints = 100,
                                 ...) {
  # compute tolerance ellipse
  ellipse <- tol_ellipse(object, horizontal = horizontal, vertical = vertical,
                         partial = partial, level = level, npoints = npoints)
  # call method for tolerance ellipse objects
  ellipse_plot(ellipse, ...)
}

#' @export
ellipse_plot.tol_ellipse <- function(object, ...) {
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
    geom_ellipse(aes_ellipse, data = object$ellipse, ...) +
    geom_scatter(aes_data, data = object$data, ...)
  # add line representing (partial) effect
  if (have_line) {
    p <- p +
      geom_effect(aes_line, data = object$line, ...)
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
geom_scatter <- function(..., linetype, lty, lwd) {
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
geom_ellipse <- function(..., fill, bg, shape, pch, cex) geom_path(...)

## custom geom for (partial) effect: avoid passing unknown argument 'fill'
geom_effect <- function(..., fill, bg, shape, pch, cex, show.legend = FALSE) {
  geom_abline(..., show.legend = show.legend)
}
