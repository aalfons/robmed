# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------


#' Set up information for a density plot of the indirect effect(s)
#'
#' Extract the relevant information for a density plot of the indirect
#' effect(s) from results of (robust) mediation analysis.
#'
#' This function is used internally by \code{\link{density_plot}()}.  It may
#' also be useful for users who want to produce a similar plot, but who want
#' more control over what information to display or how to display that
#' information.
#'
#' @param object  an object inheriting from class
#' \code{"\link{test_mediation}"} containing results from
#' (robust) mediation analysis, or a list of such objects.
#' @param grid  an optional numeric vector containing the values at which to
#' evaluate the assumed normal density from Sobel's test.  The default is to
#' take 512 equally spaced points between the estimated indirect effect
#' \eqn{\pm}{+/-} three times the standard error according to Sobel's formula.
#' @param level  numeric;  the confidence level of the confidence intervals
#' from Sobel's test.  The default is to include 95\% confidence intervals.
#' Note that this is not used for bootstrap tests, as those require to specify
#' the confidence level already in \code{\link{test_mediation}()}.
#' @param \dots  additional arguments to be passed down.
#'
#' @return An object of class \code{"setup_density_plot"} with the following
#' components:
#' \item{density}{a data frame containing the values of the indirect effect
#' where the density is estimated (column \code{Indirect}), and the estimated
#' density values (column \code{Density}).  In case of a model with multiple
#' indirect effects, there is a column \code{Effect} that indicates the
#' different indirect effects.  If a list of \code{"\link{test_mediation}"}
#' objects has been supplied, there is also a column \code{Method}, which takes
#' the names or indices of the list elements to indicate the different methods.}
#' \item{ci}{a data frame consisting of column \code{Estimate} containing the
#' point estimates, column \code{Lower} for the lower confidence limit, and
#' column \code{Upper} for the upper confidence limit.  In case of a model with
#' multiple indirect effects, there is a column \code{Effect} that indicates
#' the different indirect effects.  If a list of \code{"\link{test_mediation}"}
#' objects has been supplied, there is also a column \code{Method}, which takes
#' the names or indices of the list elements to indicate the different methods.}
#' \item{test}{a character string indicating whether the object contains
#' results from a bootstrap test (\code{"boot"}) or a Sobel test
#' (\code{"sobel"}), or a vector of such character strings if a list of
#' \code{"\link{test_mediation}"} objects has been supplied.}
#' \item{level}{numeric; the confidence level used for the confidence intervals
#' of the indirect effect(s).}
#' \item{have_effects}{a logical indicating whether the mediation model
#' contains multiple indirect effects  If \code{TRUE}, the data frames in the
#' \code{density} and \code{ci} components contain a column \code{Effect}.}
#' \item{have_methods}{a logical indicating whether a list of
#' \code{"\link{test_mediation}"} objects has been supplied.  If \code{TRUE},
#' the data frames in the \code{density} and \code{ci} components contain a
#' column \code{Method}.}
#'
#' @author Andreas Alfons
#'
#' @references
#' Alfons, A., Ates, N.Y. and Groenen, P.J.F. (2022) Robust Mediation Analysis:
#' The \R Package \pkg{robmed}.  \emph{Journal of Statistical Software},
#' \bold{103}(13), 1--45.  doi:10.18637/jss.v103.i13.
#'
#' @seealso
#' \code{\link{test_mediation}()}, \code{\link{density_plot}()}
#'
#' @examples
#' data("BSG2014")
#'
#' # run fast-and-robust bootstrap test
#' boot <- test_mediation(BSG2014,
#'                        x = "ValueDiversity",
#'                        y = "TeamCommitment",
#'                        m = "TaskConflict",
#'                        level = 0.9)
#'
#' # set up information for plot
#' setup <- setup_density_plot(boot)
#'
#' # plot only density and confidence interval
#' ggplot() +
#'   geom_density(aes(x = Indirect, y = Density), data = setup$density,
#'                stat = "identity") +
#'   geom_rect(aes(xmin = Lower, xmax = Upper,
#'                 ymin = -Inf, ymax = Inf),
#'             data = setup$ci, color = NA, alpha = 0.2) +
#'   labs(x = "Indirect effect", y = "Bootstrap density")
#'
#' @keywords hplot
#'
#' @import stats
#' @export

setup_density_plot <- function(object, ...) UseMethod("setup_density_plot")


#' @rdname setup_density_plot
#' @method setup_density_plot boot_test_mediation
#' @export

setup_density_plot.boot_test_mediation <- function(object, ...) {
  # initialization
  fit <- object$fit
  have_effects <- !is.null(fit$model) && fit$model != "simple"
  contrast <- fit$contrast                 # only implemented for regression fit
  have_contrast <- is.character(contrast)  # but this always works
  # extract point estimate and confidence interval
  indirect <- object$indirect
  ci <- object$ci
  # extract bootstrap replicates of the indirect effects
  if (inherits(fit, "reg_fit_mediation")) {
    boot_indirect <- extract_boot(fit, boot = object$reps)$indirect
  } else if (inherits(fit, "cov_fit_mediation")) {
    boot_indirect <- object$reps$t[, 5L, drop = FALSE]
  } else stop("not implemented for this type of model fit")
  # extract information to be plotted
  if (have_effects) {
    # information on indirect effects
    effect_labels <- names(indirect)
    effects <- factor(effect_labels, levels = effect_labels)
    # compute bootstrap densities and construct data frame
    density_list <- lapply(effects, function(effect) {
      pdf <- density(boot_indirect[, effect], na.rm = TRUE)
      data.frame(Effect = effect, Indirect = pdf$x, Density = pdf$y)
    })
    density <- do.call(rbind, density_list)
    # construct data frame containing confidence interval
    ci <- data.frame(Effect = effects, Estimate = unname(indirect),
                     Lower = unname(ci[, 1L]), Upper = unname(ci[, 2L]))
  } else {
    # construct data frame containing bootstrap density
    pdf <- density(boot_indirect, na.rm = TRUE)
    density <- data.frame(Indirect = pdf$x, Density = pdf$y)
    # construct data frame containing confidence interval
    ci <- data.frame(Estimate = indirect, Lower = ci[1L], Upper = ci[2L])
  }
  # return density and confidence interval
  out <- list(density = density, ci = ci, test = "boot", level = object$level,
              have_effects = have_effects, have_methods = FALSE)
  class(out) <- "setup_density_plot"
  out
}


#' @rdname setup_density_plot
#' @method setup_density_plot sobel_test_mediation
#' @export

setup_density_plot.sobel_test_mediation <- function(object, grid = NULL,
                                                    level = 0.95, ...) {
  # initializations
  level <- rep(as.numeric(level), length.out = 1L)
  if (is.na(level) || level < 0 || level > 1) level <- formals()$level
  # extract point estimate and standard error
  indirect <- object$fit$indirect
  se <- object$se
  # construct data frame containing x- and y-values for the density
  if (is.null(grid)) {
    n_grid <- formals(density.default)$n  # same number of points as density()
    grid <- seq(indirect - 3 * se, indirect + 3 * se, length.out = n_grid)
  } else grid <- as.numeric(grid)
  y <- dnorm(grid, mean = indirect, sd = se)
  density <- data.frame(Indirect = grid, Density = y)
  # construct data frame containing confidence interval
  ci <- confint_z(indirect, se, level = level, alternative = object$alternative)
  ci <- data.frame(Estimate = indirect, Lower = ci[1L], Upper = ci[2L])
  # return density and confidence interval
  out <- list(density = density, ci = ci, test = "sobel", level = level,
              have_effects = FALSE, have_methods = FALSE)
  class(out) <- "setup_density_plot"
  out
}


#' @rdname setup_density_plot
#' @method setup_density_plot list
#' @export

setup_density_plot.list <- function(object, ...) {
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
  tmp <- lapply(object, setup_density_plot, ...)
  # check that confidence levels are the same for all objects
  level <- unique(sapply(tmp, "[[", "level"))
  if (length(level) != 1L) {
    stop("confidence level must be the same for all mediation objects")
  }
  # reorganize information on the density in the proper structure
  density_list <- mapply(function(object, method) {
    data.frame(Method = method, object$density, stringsAsFactors = TRUE)
  }, object = tmp, method = methods, SIMPLIFY = FALSE, USE.NAMES = FALSE)
  density <- do.call(rbind, density_list)
  # reorganize information on the confidence interval in the proper structure
  ci_list <- mapply(function(object, method) {
    data.frame(Method = method, object$ci, stringsAsFactors = TRUE)
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
  class(out) <- "setup_density_plot"
  out
}
