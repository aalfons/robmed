# ------------------------------------
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ------------------------------------

#' Plot (robust) mediation analysis results
#' 
#' Produce dot plots of selected coefficients from regression models computed 
#' in (robust) mediation analysis, or density plots of the indirect effect.
#' 
#' @param object,x  an object of class \code{"bootMA"} or \code{"sobelMA"} 
#' containing results from (robust) mediation analysis, as returned by 
#' \code{\link{mediate}}.
#' @param data  an optional numeric vector containing the \eqn{x}-values at 
#' which to evaluate the assumed normal density from Sobel's test (only used in 
#' case of a density plot).  The default is to take 100 equally spaced points 
#' between the estimated indirect effect \eqn{\pm}{+/-} three times the 
#' standard error according to Sobel's formula.
#' @param method  a character string specifying which plot to produce.  
#' Possible values are \code{"dot"} for a dot plot of selected coefficients, or 
#' \code{"density"} for a density plot of the indirect effect.
#' @param parm  a character string specifying the coefficients to be included 
#' in a dot plot.  The default is to include the direct and the indirect effect.
#' @param level  numeric;  the confidence level of the confidence intervals 
#' from Sobel's test to be included in a dot plot.  The default is to include 
#' 95\% confidence intervals.
#' @param mapping  an aesthetic mapping to override the default behavior (see 
#' \code{\link[ggplot2]{aes}} or \code{\link[ggplot2]{aes_string}}).
#' @param \dots  additional arguments to be passed to and from methods.
#' 
#' @return An object of class \code{"ggplot"} (see 
#' \code{\link[ggplot2]{ggplot}}).
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link{mediate}}, \code{\link[=fortify.bootMA]{fortify}}
#' 
#' @keywords hplot
#' 
#' @import ggplot2
#' @export

mediatePlot <- function(object, ...) UseMethod("mediatePlot")


#' @rdname mediatePlot
#' @method mediatePlot bootMA
#' @export

mediatePlot.bootMA <- function(object, method = c("dot", "density"), 
                               parm = c("c", "ab"), ...) {
  data <- fortify(object, method=method, parm=parm)
  mediatePlot(data, ...)
}


#' @rdname mediatePlot
#' @method mediatePlot sobelMA
#' @export

mediatePlot.sobelMA <- function(object, data, method = c("dot", "density"), 
                                parm = c("c", "ab"), level = 0.95, ...) {
  data <- fortify(object, data=data, method=method, parm=parm, level=level)
  mediatePlot(data, ...)
}


#' @rdname mediatePlot
#' @method mediatePlot default
#' @export

mediatePlot.default <- function(object, mapping, ...) {
  # define aesthetic mapping for selected plot
  if(missing(mapping)) mapping <- attr(object, "mapping")
  # create selected plot
  if(attr(object, "method") == "dot") dotPlot(object, mapping, ...)
  else densityPlot(object, mapping, ...)
}


#' @rdname mediatePlot
#' @method autoplot bootMA
#' @export

autoplot.bootMA <- function(object, ...) mediatePlot(object, ...)


#' @rdname mediatePlot
#' @method autoplot sobelMA
#' @export

autoplot.sobelMA <- function(object, ...) mediatePlot(object, ...)


#' @rdname mediatePlot
#' @method plot bootMA
#' @export

plot.bootMA <- function(x, ...) mediatePlot(x, ...)


#' @rdname mediatePlot
#' @method plot sobelMA
#' @export

plot.sobelMA <- function(x, ...) mediatePlot(x, ...)


## internal function for dot plot
dotPlot <- function(data, mapping, main = NULL, xlab = NULL, ylab = NULL, ...) {
  # generate plot
  geom <- attr(data, "geom")
  ggplot(data, mapping) + geom(...) + labs(title=main, x=xlab, y=ylab)
}


## internal function for density plot
densityPlot <- function(data, mapping, main = NULL, xlab = NULL, ylab = NULL, 
                        ...) {
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
