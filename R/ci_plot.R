# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------


#' @import ggplot2
#' @export
ci_plot <- function(object, ...) UseMethod("ci_plot")

#' @export
ci_plot.default <- function(object, parm = NULL, ...) {
  # extract information
  setup <- setup_ci_plot(object, parm = parm, ...)
  # call method for corresponding objects
  ci_plot(setup, ...)
}

#' @export
ci_plot.boot_test_mediation <- function(object, parm = NULL,
                                        type = c("boot", "data"),
                                        other = c("boot", "theory"),
                                        ...) {
  # extract information
  setup <- setup_ci_plot(object, parm = parm, type = type, other = other, ...)
  # call method for corresponding objects
  ci_plot(setup, ...)
}

#' @export
ci_plot.sobel_test_mediation <- function(object, parm = NULL,
                                         level = 0.95, ...) {
  # extract information
  setup <- setup_ci_plot(object, parm = parm, level = level, ...)
  # call method for corresponding objects
  ci_plot(setup, ...)
}

#' @export
ci_plot.list <- function(object, parm = NULL, type = c("boot", "data"),
                         other = c("boot", "theory"), level = 0.95, ...) {
  # extract information
  setup <- setup_ci_plot(object, parm = parm, type = type,
                 other = other,  level = level, ...)
  # call method for corresponding objects
  ci_plot(setup, ...)
}

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


## internal function for dot plot (deprecated)
dot_plot_fortified <- function(data, mapping, facets, main = NULL,
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
