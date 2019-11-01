# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------


#' @import ggplot2
#' @export
ellipse_plot <- function(object, ...) UseMethod("ellipse_plot")

#' @export
ellipse_plot.boot_test_mediation <- function(object, ...) {
  # initial checks
  if (!inherits(object$fit, "reg_fit_mediation")) {
    stop("not implemented for this type of mediation model")
  }
  # call method for mediation model fit
  ellipse_plot(object$fit, ...)
}

#' @export
ellipse_plot.reg_fit_mediation <- function(object, horizontal = NULL,
                                           vertical = NULL, partial = FALSE,
                                           ...) {
  # obtain data to be plotted and tolerance ellipse
  df_list <- tol_ellipse(object, horizontal = horizontal, vertical = vertical,
                         partial = partial, ...)
  # define aesthetic mapping for plotting points
  if (df_list$robust) aes_data <- aes_string(x = "x", y = "y", fill = "Weight")
  else aes_data <- aes_string(x = "x", y = "y")
  # create plot
  p <- ggplot() +
    geom_path(aes_string(x = "x", y = "y"), data = df_list$ellipse) +
    geom_point(aes_data, data = df_list$data, shape = 21)
  # add line representing (partial) effect
  if (!is.null(df_list$intercept) && !is.null(df_list$slope)) {
    p <- p + geom_abline(intercept = df_list$intercept, slope = df_list$slope)
  }
  # add nice labels
  if (df_list$partial) ylab <- paste("Partial residuals of", df_list$vertical)
  else ylab <- df_list$vertical
  p <- p + labs(x = df_list$horizontal, y = ylab)
  # add color gradient for weights
  if (df_list$robust) {
    p <- p + scale_fill_gradient(low = "transparent", high = "black")
  }
  # return plot
  p
}
