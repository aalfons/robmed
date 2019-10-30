# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

#' Convert (robust) mediation analysis results into a data frame for plotting
#'
#' Supplement the estimated coefficients with other useful information for
#' informative visualization of the (robust) mediation analysis results.  It is
#' thereby possible to construct data frames for dot plots of selected
#' coefficients, as well as density plots of the indirect effect.
#'
#' @name fortify.test_mediation
#'
#' @param model  an object inheriting from class \code{"\link{test_mediation}"}
#' containing results from (robust) mediation analysis.
#' @param data  for the \code{"boot_test_mediation"} method, this is currently
#' ignored.  For the \code{"sobel_test_mediation"} method, this is an optional
#' numeric vector containing the \eqn{x}-values at which to evaluate the
#' assumed normal density from Sobel's test (only used in case of a density
#' plot).  The default is to take 100 equally spaced points between the
#' estimated indirect effect \eqn{\pm}{+/-} three times the standard error
#' according to Sobel's formula.
#' @param method  a character string specifying for which plot to construct the
#' data frame.  Possible values are \code{"dot"} for a dot plot of selected
#' coefficients, or \code{"density"} for a density plot of the indirect
#' effect(s).
#' @param parm  a character string specifying the coefficients to be included
#' in a dot plot.  The default is to include the direct and the indirect
#' effect(s).
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
#' @seealso \code{\link{test_mediation}}, \code{\link{plot_mediation}}
#'
#' @examples
#' data("BSG2014")
#'
#' # run fast and robust bootstrap test
#' test <- test_mediation(BSG2014,
#'                        x = "ValueDiversity",
#'                        y = "TeamCommitment",
#'                        m = "TaskConflict")
#'
#' # data for dot plot
#' dot <- fortify(test, method = "dot")
#' plot_mediation(dot)
#'
#' # data for density plot
#' density <- fortify(test, method = "density")
#' plot_mediation(density)
#'
#' @keywords utilities

NULL


# For multiple densities (indirect effects from different mediators), add a
# column that identifies the effect and add a formula for facetting

# For dotplots with multiple mediators, simply add more dots/errorbars in the
# same plot

#' @rdname fortify.test_mediation
#' @method fortify boot_test_mediation
#' @import ggplot2
#' @export

fortify.boot_test_mediation <- function(model, data,
                                        method = c("dot", "density"),
                                        parm = NULL, ...) {
  # initialization
  method <- match.arg(method)
  p_m <- length(model$fit$m)
  # construct data fram with relevant information
  if(method == "dot") {
    if(is.null(parm)) {
      if(p_m == 1L) parm <- c("Direct", "ab")
      else parm <- c("Direct", paste("ab", names(model$ab), sep = "_"))
    }
    # extract point estimates
    coef <- coef(model, parm=parm)
    # extract confidence intervals
    ci <- confint(model, parm=parm)
    effect <- rownames(ci)
    dimnames(ci) <- list(NULL, c("Lower", "Upper"))
    # construct data frame
    data <- data.frame(Effect=factor(effect, levels=effect),
                       Point=coef[effect], ci)
    rownames(data) <- NULL
    # add additional information as attributes
    attr(data, "mapping") <- aes_string(x="Effect", y="Point",
                                        ymin="Lower", ymax="Upper")
    attr(data, "geom") <- geom_pointrange
  } else {
    # extract point estimate and confidence interval
    ab <- model$ab
    ci <- model$ci
    # construct data frame containing bootstrap density
    if(p_m == 1L) {
      pdf <- density(model$reps$t[, 1L])
      data <- data.frame(ab = pdf$x, Density = pdf$y)
    } else {
      pdf <- lapply(seq_len(1L + p_m), function(j) density(model$reps$t[, j]))
      data <- mapply(function(density, effect) {
        data.frame(ab = density$x, Density = density$y, Effect = effect)
      }, density = pdf, effect = names(ab), SIMPLIFY = FALSE, USE.NAMES = FALSE)
      data <- do.call(rbind, data)
    }
    # add additional information as attributes
    attr(data, "mapping") <- aes_string(x = "ab", y = "Density")
    attr(data, "geom") <- function(..., stat) {
      geom_density(..., stat = "identity")
    }
    attr(data, "main") <- "Bootstrap distribution"
    if(p_m == 1L) {
      attr(data, "ci") <- data.frame(ab = ab, Density = NA_real_,
                                     Lower = ci[1L], Upper = ci[2L])
    } else {
      attr(data, "ci") <- data.frame(ab = unname(ab),
                                     Density = NA_real_,
                                     Lower = unname(ci[, 1L]),
                                     Upper = unname(ci[, 2L]),
                                     Effect = names(ab))
      attr(data, "facets") <- ~ Effect  # split plot into different panels
    }
  }
  # return data frame
  attr(data, "method") <- method
  data
}


#' @rdname fortify.test_mediation
#' @method fortify sobel_test_mediation
#' @import ggplot2
#' @export

fortify.sobel_test_mediation <- function(model, data,
                                         method = c("dot", "density"),
                                         parm = NULL,
                                         level = 0.95, ...) {
  # initialization
  method <- match.arg(method)
  level <- rep(as.numeric(level), length.out=1)
  if(is.na(level) || level < 0 || level > 1) level <- formals()$level
  # construct data fram with relevant information
  if(method == "dot") {
    if(is.null(parm)) parm <- c("Direct", "ab")
    # extract point estimates
    coef <- coef(model, parm=parm)
    # extract confidence intervals
    ci <- confint(model, parm=parm, level=level)
    effect <- rownames(ci)
    dimnames(ci) <- list(NULL, c("Lower", "Upper"))
    # construct data frame
    data <- data.frame(Effect=factor(effect, levels=effect),
                       Point=coef[effect], ci)
    rownames(data) <- NULL
    # add additional information as attributes
    attr(data, "mapping") <- aes_string(x="Effect", y="Point",
                                        ymin="Lower", ymax="Upper")
    attr(data, "geom") <- geom_pointrange
  } else {
    # extract point estimate and standard error
    ab <- model$ab
    se <- model$se
    # compute confidence interval
    ci <- confint_z(ab, se, level=level, alternative=model$alternative)
    # x- and y-values for the density
    if(missing(data)) x <- seq(ab - 3 * se, ab + 3 * se, length.out=100)
    else x <- as.numeric(data)
    y <- dnorm(x, mean=ab, sd=se)
    data <- data.frame(ab=x, Density=y)
    # add additional information as attributes
    attr(data, "mapping") <- aes_string(x="ab", y="Density")
    attr(data, "geom") <- function(..., stat) {
      geom_density(..., stat = "identity")
    }
    attr(data, "main") <- "Assumed normal distribution"
    attr(data, "ci") <- data.frame(ab, Density=dnorm(ab, mean=ab, sd=se),
                                   Lower=ci[1], Upper=ci[2])
  }
  # return data frame
  attr(data, "method") <- method
  data
}


#' @rdname fortify.test_mediation
#' @method fortify list
#' @import ggplot2
#' @export

fortify.list <- function(model, data, ...) {
  ## initializations
  is_boot <- sapply(model, inherits, "boot_test_mediation")
  is_sobel <- sapply(model, inherits, "sobel_test_mediation")
  model <- model[is_boot | is_sobel]
  if(length(model) == 0) {
    stop('no objects inheriting from class "test_mediation"')
  }
  # check if all objects use the same number of mediators
  p_m <- sapply(model, function(object) length(object$fit$m))
  if(length(unique(p_m)) > 1) {
    stop("all objects must use the same number of hypothesized mediators")
  }
  # check names of list elements
  methods <- names(model)
  if(is.null(methods)) methods <- seq_along(model)
  else {
    replace <- methods == "" | duplicated(methods)
    methods[replace] <- seq_along(model)[replace]
  }
  ## fortify each list element
  if(missing(data)) data <- lapply(model, fortify, ...)
  else data <- lapply(model, fortify, data = data, ...) # is this actually used?
  ## extract additional information
  info <- lapply(data, function(x) attributes(x)[-(1L:3L)])
  if(info[[1]]$method == "density") ci <- lapply(info, "[[", "ci")
  info <- info[[1]]
  ## combine information into one data frame
  data <- mapply(function(x, method) {
    cbind(Method = rep.int(method, nrow(x)), x, stringsAsFactors = FALSE)
  }, data, methods, SIMPLIFY = FALSE, USE.NAMES = FALSE)
  data <- do.call(rbind, data)
  data$Method <- factor(data$Method, levels = methods)
  ## modify additional information
  if(info$method == "dot") {
    # additional information for dot plot
    if(nlevels(data$Effect) == 1) {
      # only one effect to be plotted, use method on x-axis
      info$mapping <- aes_string(x="Method", y="Point",
                                 ymin="Lower", ymax="Upper")
    } else {
      # split plot into different panels
      info$facets <- ~ Method
    }
  } else {
    # additional information for density plot
    info$mapping <- aes_string(x = "ab", y = "Density", color = "Method")
    if(any(is_boot) && any(is_sobel)) {
      info$geom <- function(..., stat) {
        geom_density(..., stat = "identity")
      }
      info$main <- NULL
    }
    # combine confidence intervals for the methods
    ci <- do.call(rbind, ci)
    ci <- cbind(Method = factor(methods, levels = methods), ci)
    rownames(ci) <- NULL
    info$ci <- ci
  }
  # add additional information and return data frame
  attributes(data) <- c(attributes(data), info)
  data
}
