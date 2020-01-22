# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------


#' @import ggplot2
#' @export
density_plot <- function(object, ...) UseMethod("density_plot")

#' @export
density_plot.default <- function(object, ...) {
  # extract information
  density <- get_density(object, ...)
  # call method for corresponding objects
  density_plot(density, ...)
}

#' @export
density_plot.sobel_test_mediation <- function(object, grid = NULL,
                                              level = 0.95, ...) {
  # extract information
  density <- get_density(object, grid = grid, level = level, ...)
  # call method for corresponding objects
  density_plot(density, ...)
}

#' @export
density_plot.list <- function(object, grid = NULL, level = 0.95, ...) {
  # extract information
  density <- get_density(object, grid = grid, level = level, ...)
  # call method for corresponding objects
  density_plot(density, ...)
}

#' @export
density_plot.indirect_density <- function(object, ...) {
  # define aesthetic mappings for density estimate, point estimate and
  # confidence interval
  if(is.null(object$methods)) {
    mapping_density <- aes_string(x = "ab", y = "Density")
    mapping_line <- aes_string(xintercept = "ab")
    mapping_rect <- aes_string(xmin = "Lower", xmax = "Upper",
                               ymin = -Inf, ymax = Inf)
  } else {
    mapping_density <- aes_string(x = "ab", y = "Density", color = "Method")
    mapping_line <- aes_string(xintercept = "ab", color = "Method")
    mapping_rect <- aes_string(xmin = "Lower", xmax = "Upper",
                               ymin = -Inf, ymax = Inf,
                               fill = "Method")
  }
  # define default title
  if (all(object$test == "boot")) title <- "Bootstrap distribution"
  else if (all(object$test == "sobel")) title <- "Assumed normal distribution"
  else title <- NULL
  # generate plot
  p <- ggplot() +
    geom_densityline(mapping_density, data = object$density, ...) +
    geom_indirect(mapping_line, data = object$ci, ...) +
    geom_ci(mapping_rect, data = object$ci, ...) +
    labs(title = title, x = "Indirect effect", y = "Density")
  # split plot into different panels in case of multiple indirect effects
  if(!is.null(object$effects)) p <- p + facet_wrap(~ Effect, scales = "free")
  # return plot
  p
}


## internal function for density plot
density_plot_fortified <- function(data, mapping, facets, main = NULL,
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
    # geom_vline(mapping_line, data = ci, ...) +
    # geom_rect(mapping_rect, data = ci, color = NA, alpha = 0.2, ...) +
    geom_indirect(mapping_line, data = ci, ...) +
    geom_ci(mapping_rect, data = ci, ...) +
    labs(title = main, x = xlab, y = ylab)
  if(!is.null(facets)) {
    # split plot into different panels
    if(length(facets) == 2L) p <- p + facet_wrap(facets, scales = "free")
    else p <- p + facet_grid(facets, scales = "free")
  }
  p
}


## custom geom for density estimate to be used in density plot
#  1) always use stat = "identity" because the density is already estimated
#  2) do not allow for a fill color because a filled rectangle is used to
#     display the confidence interval
geom_densityline <- function(..., stat, fill, bg, alpha) {
  geom_density(..., stat = "identity")
}

## custom geom for vertical line to be used in density plot:
#  1) avoid passing argument 'alpha' to ensure that line is of the same style
#     as lines of density
#  2) avoid passing unknown argument 'fill'
geom_indirect <- function(..., fill, bg, alpha) geom_vline(...)

## custom geom for confidence intervals to be used in density plot:
#  fix transparant rectangle without edges and avoid duplication of arguments
geom_ci <- function(...) {
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
