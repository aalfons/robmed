# ------------------------------------
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ------------------------------------

#' Convert (robust) mediation analysis results into a data frame for plotting
#' 
#' Supplement the estimated coefficients with other useful information for 
#' informative visualization of the (robust) mediation analysis results.  It is 
#' thereby possible to construct data frames for dot plots of selected 
#' coefficients, as well as density plots of the indirect effect.
#' 
#' @method fortify bootMA
#' 
#' @param model  an object of class \code{"bootMA"} or \code{"sobelMA"} 
#' containing results from (robust) mediation analysis, as returned by 
#' \code{\link{mediate}}.
#' @param data  for the \code{"bootMA"} method, this is currently ignored.  For 
#' the \code{"sobelMA"} method, this is an optional numeric vector containing 
#' the \eqn{x}-values at which to evaluate the assumed normal density from 
#' Sobel's test (only used in case of a density plot).  The default is to take 
#' 100 equally spaced points between the estimated indirect effect 
#' \eqn{\pm}{+/-} three times the standard error according to Sobel's formula.
#' @param method  a character string specifying for which plot to construct the 
#' data frame.  Possible values are \code{"dot"} for a dot plot of selected 
#' coefficients, or \code{"density"} for a density plot of the indirect effect.
#' @param parm  a character string specifying the coefficients to be included 
#' in a dot plot.  The default is to include the direct and the indirect effect.
#' @param level  numeric;  the confidence level of the confidence intervals 
#' from Sobel's test to be included in a dot plot.  The default is to include 
#' 95\% confidence intervals.
#' @param \dots  additional arguments are currently ignored.
#' 
#' @return A data frame containing the necessary data for the selected plot, as 
#' well as additional information stored in attributes.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link{mediate}}, \code{\link[=mediatePlot]{plot}}
#' 
#' @keywords utilities
#' 
#' @import ggplot2
#' @export

fortify.bootMA <- function(model, data, method = c("dot", "density"), 
                           parm = c("c", "ab"), ...) {
  # initialization
  method <- match.arg(method)
  # construct data fram with relevant information
  if(method == "dot") {
    # extract point estimates
    coef <- coefficients(model, parm=parm)
    # extract confidence intervals
    ci <- confint(model, parm=parm)
    effect <- rownames(ci)
    dimnames(ci) <- list(NULL, c("lower", "upper"))
    # construct data frame
    data <- data.frame(effect=factor(effect, levels=effect), 
                       point=coef[effect], ci)
    rownames(data) <- NULL
    # add additional information as attributes
    attr(data, "mapping") <- aes_string(x="effect", y="point", 
                                        ymin="lower", ymax="upper")
    attr(data, "geom") <- geom_pointrange
  } else {
    # construct data frame containing bootstrap replicates
    data <- data.frame(ab=model$reps$t)
    # extract point estimate and confidence interval
    ab <- model$ab
    ci <- model$ci
    # add additional information as attributes
    attr(data, "mapping") <- aes_string(x="ab")
    attr(data, "geom") <- function(..., stat) geom_density(..., stat="density")
    attr(data, "main") <- "Bootstrap distribution"
    attr(data, "ci") <- data.frame(ab, lower=ci[1], upper=ci[2])
  }
  # return data frame
  attr(data, "method") <- method
  data
}


#' @rdname fortify.bootMA
#' @method fortify sobelMA
#' @import ggplot2
#' @export

fortify.sobelMA <- function(model, data, method = c("dot", "density"), 
                            parm = c("c", "ab"), level = 0.95, ...) {
  # initialization
  method <- match.arg(method)
  level <- rep(as.numeric(level), length.out=1)
  if(is.na(level) || level < 0 || level > 1) level <- formals()$level
  # construct data fram with relevant information
  if(method == "dot") {
    # extract point estimates
    coef <- coefficients(model, parm=parm)
    # extract confidence intervals
    ci <- confint(model, parm=parm, level=level)
    effect <- rownames(ci)
    dimnames(ci) <- list(NULL, c("lower", "upper"))
    # construct data frame
    data <- data.frame(effect=factor(effect, levels=effect), 
                       point=coef[effect], ci)
    rownames(data) <- NULL
    # add additional information as attributes
    attr(data, "mapping") <- aes_string(x="effect", y="point", 
                                        ymin="lower", ymax="upper")
    attr(data, "geom") <- geom_pointrange
  } else {
    # extract point estimate and standard error
    ab <- model$ab
    se <- model$se
    # compute confidence interval
    ci <- confintZ(ab, se, level=level, alternative=model$alternative)
    # x- and y-values for the density
    if(missing(data)) x <- seq(ab - 3 * se, ab + 3 * se, length.out=100)
    else x <- as.numeric(data)
    y <- dnorm(x, mean=ab, sd=se)
    data <- data.frame(ab=x, density=y)
    # add additional information as attributes
    attr(data, "mapping") <- aes_string(x="ab", y="density")
    attr(data, "geom") <- function(..., stat) geom_density(..., stat="identity")
    attr(data, "main") <- "Assumed normal distribution"
    attr(data, "ci") <- data.frame(ab, density=dnorm(ab, mean=ab, sd=se), 
                                   lower=ci[1], upper=ci[2])
  }
  # return data frame
  attr(data, "method") <- method
  data
}
