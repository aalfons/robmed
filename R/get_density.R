# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

## function to prepare information for density plot

#' @export
get_density <- function(object, ...) UseMethod("get_density")

#' @export
get_density.boot_test_mediation <- function(object, ...) {
  # initialization
  p_m <- length(object$fit$m)
  have_effects <- p_m > 1L
  # extract point estimate and confidence interval
  ab <- object$ab
  ci <- object$ci
  # extract information to be plotted
  if (have_effects) {
    # information on indirect effects
    effects <- names(ab)
    # construct data frame containing bootstrap density
    pdf_list <- lapply(seq_len(1L+p_m), function(j) density(object$reps$t[, j]))
    density_list <- mapply(function(pdf, effect) {
      data.frame(Effect = effect, ab = pdf$x, Density = pdf$y)
    }, pdf = pdf_list, effect = effects, SIMPLIFY = FALSE, USE.NAMES = FALSE)
    density <- do.call(rbind, density_list)
    # construct data frame containing confidence interval
    ci <- data.frame(Effect = effects, ab = unname(ab),
                     Lower = unname(ci[, 1L]),
                     Upper = unname(ci[, 2L]))
  } else {
    # construct data frame containing bootstrap density
    pdf <- density(object$reps$t[, 1L])
    density <- data.frame(ab = pdf$x, Density = pdf$y)
    # construct data frame containing confidence interval
    ci <- data.frame(ab = ab, Lower = ci[1L], Upper = ci[2L])
  }
  # return density and confidence interval
  out <- list(density = density, ci = ci, test = "boot", level = object$level,
              have_effects = have_effects, have_methods = FALSE)
  class(out) <- "indirect_density"
  out
}

#' @import stats
#' @export
get_density.sobel_test_mediation <- function(object, grid = NULL, level = 0.95,
                                             ...) {
  # initializations
  level <- rep(as.numeric(level), length.out = 1)
  if (is.na(level) || level < 0 || level > 1) level <- formals()$level
  # extract point estimate and standard error
  ab <- object$ab
  se <- object$se
  # construct data frame containing x- and y-values for the density
  if (is.null(grid)) {
    n_grid <- formals(density.default)$n  # same number of points as density()
    grid <- seq(ab - 3 * se, ab + 3 * se, length.out = n_grid)
  } else grid <- as.numeric(grid)
  y <- dnorm(grid, mean = ab, sd = se)
  density <- data.frame(ab = grid, Density = y)
  # construct data frame containing confidence interval
  ci <- confint_z(ab, se, level = level, alternative = object$alternative)
  ci <- data.frame(ab, Lower = ci[1], Upper = ci[2])
  # return density and confidence interval
  out <- list(density = density, ci = ci, test = "sobel", level = level,
              have_effects = FALSE, have_methods = FALSE)
  class(out) <- "indirect_density"
  out
}

#' @export
get_density.list <- function(object, ...) {
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
  tmp <- lapply(object, get_density, ...)
  # check that confidence levels are the same for all objects
  level <- unique(sapply(tmp, "[[", "level"))
  if (length(level) != 1L) {
    stop("confidence level must be the same for all mediation objects")
  }
  # reorganize information on the density in the proper structure
  density_list <- mapply(function(object, method) {
    data.frame(Method = method, object$density)
  }, object = tmp, method = methods, SIMPLIFY = FALSE, USE.NAMES = FALSE)
  density <- do.call(rbind, density_list)
  # reorganize information on the confidence interval in the proper structure
  ci_list <- mapply(function(object, method) {
    data.frame(Method = method, object$ci)
  }, object = tmp, method = methods, SIMPLIFY = FALSE, USE.NAMES = FALSE)
  ci <- do.call(rbind, ci_list)
  # extract information on type of test
  test <- sapply(tmp, "[[", "test")
  # all objects have the same variables, so we can take the information
  # whether there are multiple indirect effects from the first one
  have_effects <- tmp[[1L]]$have_effects
  # return density and confidence interval
  out <- list(density = density, ci = ci, test = test, level = level,
              have_effects = have_effects, have_methods = TRUE)
  class(out) <- "indirect_density"
  out
}
