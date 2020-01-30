# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

## function to prepare information for confidence interval plot

#' @export
setup_ci_plot <- function(object, ...) UseMethod("setup_ci_plot")

#' @export
setup_ci_plot.boot_test_mediation <- function(object, parm = NULL,
                                              type = c("boot", "data"),
                                              other = c("boot", "theory"),
                                              ...) {
  # initializations
  p_m <- length(object$fit$m)
  if (is.null(parm)) {
    if (p_m == 1L) parm <- c("Direct", "ab")
    else parm <- c("Direct", paste("ab", names(object$ab), sep = "_"))
  }
  # extract point estimates
  coefficients <- coef(object, parm = parm, type = type)
  # extract confidence intervals
  ci <- confint(object, parm = parm, other = other)
  # construct data frame
  effects <- rownames(ci)
  ci <- data.frame(Effect = factor(effects, levels = effects),
                   Estimate = unname(coefficients[effects]),
                   Lower = unname(ci[, 1L]),
                   Upper = unname(ci[, 2L]))
  # return point estimates and confidence interval
  out <- list(ci = ci, level = object$level, have_methods = FALSE)
  class(out) <- "setup_ci_plot"
  out
}

#' @export
setup_ci_plot.sobel_test_mediation <- function(object, parm = NULL,
                                               level = 0.95, ...) {
  # initializations
  if (is.null(parm)) parm <- c("Direct", "ab")
  level <- rep(as.numeric(level), length.out = 1)
  if (is.na(level) || level < 0 || level > 1) level <- formals()$level
  # extract point estimates
  coefficients <- coef(object, parm = parm)
  # extract confidence intervals
  ci <- confint(object, parm = parm, level = level)
  # construct data frame
  effects <- rownames(ci)
  ci <- data.frame(Effect = factor(effects, levels = effects),
                   Estimate = unname(coefficients[effects]),
                   Lower = unname(ci[, 1L]),
                   Upper = unname(ci[, 2L]))
  # return point estimates and confidence interval
  out <- list(ci = ci, level = level, have_methods = FALSE)
  class(out) <- "setup_ci_plot"
  out
}

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
    data.frame(Method = method, object$ci)
  }, object = tmp, method = methods, SIMPLIFY = FALSE, USE.NAMES = FALSE)
  ci <- do.call(rbind, ci_list)
  # return point estimates and confidence interval
  out <- list(ci = ci, level = level, have_methods = TRUE)
  class(out) <- "setup_ci_plot"
  out
}
