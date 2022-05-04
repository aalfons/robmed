# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------


#' Set up information for a dot plot with confidence intervals
#'
#' Extract the relevant information for a dot plot with confidence intervals of
#' selected effects from (robust) mediation analysis.  Information on p-values
#' of the selected effects can be included in addition to confidence intervals.
#'
#' This function is used internally by \code{\link{ci_plot}()}.  It may also
#' be useful for users who want to produce a similar plot, but who want more
#' control over what information to display or how to display that information.
#'
#' @param object  an object inheriting from class
#' \code{"\link{test_mediation}"} containing results from
#' (robust) mediation analysis, or a list of such objects.
#' @param parm  an integer, character or logical vector specifying which
#' effects to include in the plot.  In case of a character vector, possible
#' values are \code{"a"}, \code{"b"}, \code{"d"} (only serial multiple mediator
#' models), \code{"total"}, \code{"direct"}, and \code{"indirect"}.  The
#' default is to include the direct and the indirect effect(s).
#' @param type  a character string specifying which point estiamates and
#' confidence intervals to plot: those based on the bootstrap distribution
#' (\code{"boot"}; the default), or those based on the original data
#' (\code{"data"}).  If \code{"boot"}, the confidence intervals of effects
#' other than the indirect effect(s) are computed using a normal approximation
#' (i.e., assuming a normal distribution of the corresponding effect with the
#' standard deviation computed from the bootstrap replicates).  If
#' \code{"data"}, the confidence intervals of effects other than the indirect
#' effect(s) are computed via statistical theory based on the original data
#' (e.g., based on a t-distribution if the coefficients are estimated via
#' regression).  Note that this is only relevant for mediation analysis via a
#' bootstrap test, where the confidence interval of the indirect effect is
#' always computed via a percentile-based method due to the asymmetry of its
#' distribution.
#' @param level  numeric;  the confidence level of the confidence intervals
#' from Sobel's test.  The default is to include 95\% confidence intervals.
#' Note that this is not used for bootstrap tests, as those require to specify
#' the confidence level already in \code{\link{test_mediation}()}.
#' @param p_value  a logical indicating whether to include information on the
#' p-values in addition to the confidence intervals.  The default is
#' \code{FALSE}.
#' @param digits  an integer determining how many digits to compute for
#' bootstrap p-values of the indirect effects (see \code{\link{p_value}()}).
#' The default is to compute 4 digits after the comma.  This is only relevant
#' if \code{p_value = TRUE}.
#' @param \dots  additional arguments to be passed down.
#'
#' @return An object of class \code{"setup_ci_plot"} with the following
#' components:
#' \item{ci}{a data frame consisting of column \code{Effect} indicating the
#' different effects, column \code{Estimate} containing the point estimates,
#' column \code{Lower} for the lower confidence limit, and column \code{Upper}
#' for the upper confidence limit.  If argument \code{p_value = TRUE}, there is
#' an additional column \code{Label} which gives the default facet label for
#' the confidence intervals. If a list of \code{"\link{test_mediation}"}
#' objects has been supplied, there is also a column \code{Method}, which takes
#' the names or indices of the list elements to indicate the different methods.}
#' \item{p_value}{a data frame consisting of column \code{Label} which gives
#' the default facet label for the p-values, column \code{Effect} indicating
#' the different effects, and column \code{Value} containing the p-values.  If
#' a list of \code{"\link{test_mediation}"} objects has been supplied, there is
#' also a column \code{Method}, which takes the names or indices of the list
#' elements to indicate the different methods.  This is only returned if
#' argument \code{p_value = TRUE}.}
#' \item{level}{numeric; the confidence level used for the confidence intervals
#' of the indirect effect(s).}
#' \item{have_methods}{a logical indicating whether a list of
#' \code{"\link{test_mediation}"} objects has been supplied.  If \code{TRUE},
#' the data frame in the \code{ci} component contains a column \code{Method}.}
#'
#' @author Andreas Alfons
#'
#' @seealso
#' \code{\link{test_mediation}()}, \code{\link{ci_plot}()}
#'
#' @examples
#' data("BSG2014")
#'
#' # run fast-and-robust bootstrap test
#' boot <- test_mediation(BSG2014,
#'                        x = "ValueDiversity",
#'                        y = "TeamCommitment",
#'                        m = "TaskConflict")
#'
#' # set up information for plot
#' setup <- setup_ci_plot(boot, parm = "Indirect")
#'
#' # plot only density and confidence interval
#' ggplot() +
#'   geom_hline(yintercept = 0, color = "darkgrey") +
#'   geom_pointrange(aes(x = "Robust bootstrap", y = Estimate,
#'                       ymin = Lower, ymax = Upper),
#'                   data = setup$ci) +
#'   labs(x = NULL, y = "Indirect effect")
#'
#' @keywords hplot
#'
#' @export

setup_ci_plot <- function(object, ...) UseMethod("setup_ci_plot")


#' @rdname setup_ci_plot
#' @method setup_ci_plot boot_test_mediation
#' @export

setup_ci_plot.boot_test_mediation <- function(object,
                                              parm = c("direct", "indirect"),
                                              type = c("boot", "data"),
                                              p_value = FALSE, digits = 4L,
                                              ...) {
  # initializations
  include_p_value <- isTRUE(p_value)
  # extract point estimates
  coefficients <- coef(object, parm = parm, type = type)
  # extract confidence intervals
  ci <- confint(object, parm = parm, type = type)
  # construct data frame
  effect_labels <- rownames(ci)
  effects <- factor(effect_labels, levels = effect_labels)
  ci <- data.frame(Effect = effects,
                   Estimate = unname(coefficients[effects]),
                   Lower = unname(ci[, 1L]),
                   Upper = unname(ci[, 2L]))
  # if requested, compute or extract p-values
  if (include_p_value) {
    # labels for facets
    labels <- c("Confidence interval", "p-Value")
    labels <- factor(labels, levels = labels)
    # add column to CI data frame for facetting
    ci <- cbind(Label = labels[1L], ci)
    # construct data frame for p-values
    p_values <- p_value(object, parm = parm, type = type, digits = digits)
    p_value <- data.frame(Label = labels[2L], Effect = effects,
                          Value = unname(p_values[effects]))
    # return object contains confidence intervals and p-values
    out <- list(ci = ci, p_value = p_value, level = object$level,
                have_methods = FALSE)
  } else {
    # return object contains only confidence intervals
    out <- list(ci = ci, level = object$level, have_methods = FALSE)
  }
  # add class and return object
  class(out) <- "setup_ci_plot"
  out
}


#' @rdname setup_ci_plot
#' @method setup_ci_plot sobel_test_mediation
#' @export

setup_ci_plot.sobel_test_mediation <- function(object,
                                               parm = c("direct", "indirect"),
                                               level = 0.95, p_value = FALSE,
                                               ...) {
  # initializations
  level <- rep(as.numeric(level), length.out = 1L)
  if (is.na(level) || level < 0 || level > 1) level <- formals()$level
  include_p_value <- isTRUE(p_value)
  # extract point estimates
  coefficients <- coef(object, parm = parm)
  # extract confidence intervals
  ci <- confint(object, parm = parm, level = level)
  # construct data frame
  effect_labels <- rownames(ci)
  effects <- factor(effect_labels, levels = effect_labels)
  ci <- data.frame(Effect = effects,
                   Estimate = unname(coefficients[effects]),
                   Lower = unname(ci[, 1L]),
                   Upper = unname(ci[, 2L]))
  # if requested, compute or extract p-values
  if (include_p_value) {
    # labels for facets
    labels <- c("Confidence interval", "p-Value")
    labels <- factor(labels, levels = labels)
    # construct data frame for p-values
    p_values <- p_value(object, parm = parm)
    p_value <- data.frame(Label = labels[2], Effect = effects,
                          Value = unname(p_values[effects]))
    # add column to CI data frame for facetting
    ci <- cbind(Label = labels[1], ci)
    # return object contains confidence intervals and p-values
    out <- list(ci = ci, p_value = p_value, level = level, have_methods = FALSE)
  } else {
    # return object contains only confidence intervals
    out <- list(ci = ci, level = level, have_methods = FALSE)
  }
  # add class and return object
  class(out) <- "setup_ci_plot"
  out
}


#' @rdname setup_ci_plot
#' @method setup_ci_plot list
#' @export

setup_ci_plot.list <- function(object, ...) {
  # initializations
  is_boot <- sapply(object, inherits, "boot_test_mediation")
  is_sobel <- sapply(object, inherits, "sobel_test_mediation")
  object <- object[is_boot | is_sobel]
  if (length(object) == 0L) {
    stop('no objects inheriting from class "test_mediation"')
  }
  # check that variables are the same
  components <- c("x", "y", "m", "covariates")
  variables <- lapply(object, function(x) x$fit[components])
  all_identical <- all(sapply(variables[-1L], identical, variables[[1L]]))
  if (!isTRUE(all_identical)) {
    stop("all mediation objects must use the same variables")
  }
  # check names of list elements
  methods <- names(object)
  if (is.null(methods)) methods <- seq_along(object)
  else {
    replace <- methods == "" | duplicated(methods)
    methods[replace] <- seq_along(object)[replace]
  }
  # extract information for each list element
  tmp <- lapply(object, setup_ci_plot, ...)
  # check that confidence levels are the same for all objects
  level <- unique(sapply(tmp, "[[", "level"))
  if (length(level) != 1L) {
    stop("confidence level must be the same for all mediation objects")
  }
  # reorganize information on the confidence interval in the proper structure
  ci_list <- mapply(function(object, method) {
    data.frame(Method = method, object$ci, stringsAsFactors = TRUE)
  }, object = tmp, method = methods, SIMPLIFY = FALSE, USE.NAMES = FALSE)
  ci <- do.call(rbind, ci_list)
  # if available, reorganize information on the p-value in the proper structure
  have_p_value <- !is.null(tmp[[1]]$p_value)
  if (have_p_value) {
    p_value_list <- mapply(function(object, method) {
      data.frame(Method = method, object$p_value, stringsAsFactors = TRUE)
    }, object = tmp, method = methods, SIMPLIFY = FALSE, USE.NAMES = FALSE)
    p_value <- do.call(rbind, p_value_list)
    # return object contains confidence intervals and p-values
    out <- list(ci = ci, p_value = p_value, level = level, have_methods = TRUE)
  } else {
    # return object contains only confidence intervals
    out <- list(ci = ci, level = level, have_methods = TRUE)
  }
  # add class and return object
  class(out) <- "setup_ci_plot"
  out
}
