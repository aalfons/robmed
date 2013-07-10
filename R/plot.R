# ------------------------------------
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ------------------------------------

#' @method plot bootMA
#' @import ggplot2
#' @export

plot.bootMA <- function(x, method = c("dot", "density"), 
                        parm = c("c", "ab"), ...) {
  mediatePlot(x, method=method, parm=parm, ...)
}


#' @method plot sobelMA
#' @import ggplot2
#' @export

plot.sobelMA <- function(x, data, method = c("dot", "density"), 
                         parm = c("c", "ab"), level = 0.95, ...) {
  mediatePlot(x, data=data, method=method, parm=parm, level=level, ...)
}


#' @method autoplot bootMA
#' @import ggplot2
#' @export

autoplot.bootMA <- function(object, method = c("dot", "density"), 
                            parm = c("c", "ab"), ...) {
  mediatePlot(object, method=method, parm=parm, ...)
}


#' @method autoplot sobelMA
#' @import ggplot2
#' @export

autoplot.sobelMA <- function(object, data, method = c("dot", "density"), 
                             parm = c("c", "ab"), level = 0.95, ...) {
  mediatePlot(object, data=data, method=method, parm=parm, level=level, ...)
}


## internal function to select plot

mediatePlot <- function(object, ...) UseMethod("mediatePlot")

mediatePlot.bootMA <- function(object, method = c("dot", "density"), 
                               parm = c("c", "ab"), ...) {
  data <- fortify(object, method=method, parm=parm)
  if(attr(data, "method") == "dot") dotPlot(data, ...)
  else densityPlot(data, ...)
}

mediatePlot.sobelMA <- function(object, data, method = c("dot", "density"), 
                                parm = c("c", "ab"), level = 0.95, ...) {
  data <- fortify(object, data=data, method=method, parm=parm, level=level)
  if(attr(data, "method") == "dot") dotPlot(data, ...)
  else densityPlot(data, ...)
}


## internal function for dot plot
dotPlot <- function(data, mapping, main = NULL, xlab = NULL, ylab = NULL, ...) {
  # define aesthetic mapping for density plot
  if(missing(mapping)) mapping <- attr(data, "mapping")
  # generate plot
  geom <- attr(data, "geom")
  ggplot(data, mapping) + geom(...) + labs(title=main, x=xlab, y=ylab)
}


## internal function for density plot
densityPlot <- function(data, mapping, main = NULL, xlab = NULL, ylab = NULL, 
                        ...) {
  # define aesthetic mapping for density plot
  if(missing(mapping)) mapping <- attr(data, "mapping")
  # define default title and axis labels
  if(is.null(main)) main <- attr(data, "main")
  if(is.null(xlab)) xlab <- "Indirect effect"
  if(is.null(ylab)) ylab <- "Density"
  # extract point estimate and confidence interval
  ci <- attr(data, "ci")
  # generate plot
  geom <- attr(data, "geom")
  ggplot(data, mapping) + geom(...) + 
    geom_vline(aes_string(xintercept="ab"), data=ci, ...) + 
    geom_rect(aes_string(xmin="lower", xmax="upper", ymin=-Inf, ymax=Inf), 
              data=ci, col=NA, alpha=0.2, ...) + 
    labs(title=main, x=xlab, y=ylab)
}
