# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------


## maybe it's better to always facet by effect; and if there are multiple
## methods, put those on the x-axis


## internal function for dot plot
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
