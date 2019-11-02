# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------


#' @import ggplot2
#' @export
ellipse_plot <- function(object, ...) UseMethod("ellipse_plot")

#' @export
ellipse_plot.test_mediation <- function(object, ...) {
  # # initial checks
  # if (!inherits(object$fit, "reg_fit_mediation")) {
  #   stop("not implemented for this type of mediation model")
  # }
  # call method for mediation model fit
  ellipse_plot(object$fit, ...)
}

#' @export
ellipse_plot.fit_mediation <- function(object, horizontal = NULL,
                                       vertical = NULL, partial = FALSE,
                                       level = 0.975, npoints = 100, ...) {
  # obtain data to be plotted and tolerance ellipse
  df_list <- tol_ellipse(object, horizontal = horizontal, vertical = vertical,
                         partial = partial, level = level, npoints = npoints)
  # define aesthetic mapping for plotting points
  if (df_list$robust) aes_data <- aes_string(x = "x", y = "y", fill = "Weight")
  else aes_data <- aes_string(x = "x", y = "y")
  # create plot
  p <- ggplot() +
    geom_ellipse(aes_string(x = "x", y = "y"), data = df_list$ellipse, ...) +
    geom_scatter(aes_data, data = df_list$data, ...)
  # add line representing (partial) effect
  if (!is.null(df_list$intercept) && !is.null(df_list$slope)) {
    p <- p +
      geom_effect(intercept = df_list$intercept, slope = df_list$slope, ...)
  }
  # add nice labels
  if (df_list$partial) ylab <- paste("Partial residuals of", df_list$vertical)
  else ylab <- df_list$vertical
  p <- p + labs(x = df_list$horizontal, y = ylab)
  # add color gradient for weights
  if (df_list$robust) {
    p <- p +
      scale_fill_gradient(limits = 0:1, low = "transparent", high = "black")
  }
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
  # call existing geom function
  do.call(geom_point, arguments)
}

## custom geom for ellipse: avoid passing unknown argument 'fill'
geom_ellipse <- function(..., fill, bg, shape, pch, cex) geom_path(...)

## custom geom for (partial) effect: avoid passing unknown argument 'fill'
geom_effect <- function(..., fill, bg, shape, pch, cex) geom_abline(...)
