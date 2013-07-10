# ------------------------------------
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ------------------------------------

#' @method fortify bootMA
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
