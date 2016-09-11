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
#' \code{"\link{testMediation}"} containing results from (robust) mediation
#' analysis.  For \code{plotMediation}, a list of such objects may be supplied
#' as well.
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
#' @seealso \code{\link{testMediation}},
#' \code{\link[=fortify.testMediation]{fortify}}
#'
#' @keywords hplot
#'
#' @import ggplot2
#' @export

plotMediation <- function(object, ...) UseMethod("plotMediation")


#' @rdname plotMediation
#' @method plotMediation bootTestMediation
#' @export

plotMediation.bootTestMediation <- function(object,
                                            method = c("dot", "density"),
                                            parm = c("c", "ab"),
                                            ...) {
  data <- fortify(object, method=method, parm=parm)
  plotMediation(data, ...)
}


#' @rdname plotMediation
#' @method plotMediation sobelTestMediation
#' @export

plotMediation.sobelTestMediation <- function(object, data,
                                             method = c("dot", "density"),
                                             parm = c("c", "ab"),
                                             level = 0.95, ...) {
  data <- fortify(object, data=data, method=method, parm=parm, level=level)
  plotMediation(data, ...)
}


#' @rdname plotMediation
#' @method plotMediation list
#' @export

plotMediation.list <- function(object, data, method = c("dot", "density"),
                               parm = c("c", "ab"), level = 0.95, ...) {
  data <- fortify(object, data=data, method=method, parm=parm, level=level)
  plotMediation(data, ...)
}


#' @rdname plotMediation
#' @method plotMediation default
#' @export

plotMediation.default <- function(object, mapping = attr(object, "mapping"),
                                  facets = attr(object, "facets"), ...) {
  # create selected plot
  if(attr(object, "method") == "dot") dotPlot(object, mapping, facets, ...)
  else densityPlot(object, mapping, ...)
}


#' @rdname plotMediation
#' @method autoplot testMediation
#' @export

autoplot.testMediation <- function(object, ...) plotMediation(object, ...)


#' @rdname plotMediation
#' @method plot testMediation
#' @export

plot.testMediation <- function(x, ...) plotMediation(x, ...)


## internal function for dot plot
dotPlot <- function(data, mapping, facets, main = NULL,
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
densityPlot <- function(data, mapping, main = NULL, xlab = NULL, ylab = NULL,
                        ...) {
  # define default title and axis labels
  if(is.null(main)) main <- attr(data, "main")
  if(is.null(xlab)) xlab <- "Indirect effect"
  if(is.null(ylab)) ylab <- "Density"
  # extract point estimate and confidence interval
  ci <- attr(data, "ci")
  if("Method" %in% names(data)) {
    mappingLine <- aes_string(xintercept="ab", color="Method")
    mappingRect <- aes_string(xmin="Lower", xmax="Upper", ymin=-Inf, ymax=Inf,
                              fill="Method")
  } else {
    mappingLine <- aes_string(xintercept="ab")
    mappingRect <- aes_string(xmin="Lower", xmax="Upper", ymin=-Inf, ymax=Inf)
  }
  # generate plot
  geom <- attr(data, "geom")
  ggplot(data, mapping) + geom(...) + geom_vline(mappingLine, data=ci, ...) +
    geom_rect(mappingRect, data=ci, color=NA, alpha=0.2, ...) +
    labs(title=main, x=xlab, y=ylab)
}
