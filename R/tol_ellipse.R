# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

## function to compute tolerance ellipses

#' @importFrom ellipse ellipse
#' @export
tol_ellipse <- function(object, ...) UseMethod("tol_ellipse")

#' @export
tol_ellipse.test_mediation <- function(object, ...) {
  # simply call methods for mediation model fit
  tol_ellipse(object$fit, ...)
}

#' @export
tol_ellipse.reg_fit_mediation <- function(object,
                                          variables = c("mx", "ym", "yx"),
                                          level = 0.975, npoints = 100,
                                          ...) {
  # initializations
  if (object$robust && object$median) {
    stop("tolerance ellipse not meaningful for median regression")
  }
  if (length(object$m) + length(object$covariates) > 1) {
    stop("currently only implemented for simple mediation models")
  }
  variables <- match.arg(variables)
  # select which variables to plot
  if (variables == "mx") select <- c(object$x, object$m)
  else if (variables == "ym") select <- c(object$m, object$y)
  else select <- c(object$x, object$y)
  # extract data
  data <- object$data[, select]
  # obtain location and shape of ellipse
  if (object$robust) {
    # extract weights in case of robust regression
    fit <- if (variables == "mx") "fit_mx" else "fit_ymx"
    w <- weights(object[[fit]], type = "robustness")
    # compute weighted mean and weighted covariance matrix
    center <- sapply(data, weighted.mean, w = w)
    cov <- weighted.cov(data, w = w)
  } else {
    # compute mean and covariance matrix
    center <- colMeans(data)
    cov <- cov(data)
  }
  # compute ellipse
  ellipse <- ellipse(cov, centre = center, level = level, npoints = npoints)
  as.data.frame(ellipse)
}

#' @export
tol_ellipse.cov_fit_mediation <- function(object,
                                          variables = c("mx", "ym", "yx"),
                                          level = 0.975, npoints = 100,
                                          ...) {
  # initializations
  if (length(object$m) + length(object$covariates) > 1) {
    stop("currently only implemented for simple mediation models")
  }
  variables <- match.arg(variables)
  # select which variables to plot
  if (variables == "mx") select <- c(object$x, object$m)
  else if (variables == "ym") select <- c(object$m, object$y)
  else select <- c(object$x, object$y)
  # extract covariance matrix
  center <- object$cov$center[select]
  cov <- object$cov$cov[select, select]
  # compute ellipse
  ellipse <- ellipse(cov, centre = center, level = level, npoints = npoints)
  as.data.frame(ellipse)
}

# ## @export
# tol_ellipse.list <- function(object, ...) {
#   # If I add a column 'Method' to the data frame, it could happen that one of
#   # the variables is also called 'Method'
#   # Maybe let the user do some work for combining ellipses?
#   stop("not implemented yet")
# }


## utility functions

# weighted covariance matrix
weighted.cov <- function(x, w, center = NULL, ...) {
  x <- as.matrix(x)
  if (missing(w)) cov(x, use = "pairwise.complete.obs")
  else {
    # initial check
    if (length(w) != nrow(x)) {
      stop("length of 'w' must equal the number of rows in 'x'")
    }
    if (is.null(center)) center <- apply(x, 2, weighted.mean, w = w)
    else if (length(center) != ncol(x)) {
      stop("length of 'center' must equal the number of columns in 'x'")
    }
    # sweep out center estimates and compute weighted cross product
    x <- sweep(x, 2, center, check.margin = FALSE)
    weighted.crossprod(x, w = w)
  }
}

# weighted cross product
weighted.crossprod <- function(x, w) {
  ci <- colnames(x)
  if (is.null(ci)) ci <- seq_len(ncol(x))
  sapply(ci, function(j) sapply(ci, function(i) {
    select <- !is.na(x[, i]) & !is.na(x[, j])
    xi <- x[select, i]
    xj <- x[select, j]
    w <- w[select]
    sum(xi * xj * w) / (sum(w) - 1)
  }))
}


# ## @import ggplot2
# ## @importFrom ellipse ellipse
# ## @export
# ellipse_plot <- function(object, ...) UseMethod("ellipse_plot")
#
# ## @export
# ellipse_plot.boot_test_mediation <- function(object, ...) {
#   # initial checks
#   if (!inherits(object$fit, "reg_fit_mediation")) {
#     stop("not implemented for this type of mediation model")
#   }
#   if (length(object$fit$m) + length(object$fit$covariates) > 1) {
#     stop("currently only implemented for simple mediation models")
#   }
#   # call method for mediation model fit
#   ellipse_plot(object$fit)
# }
#
# ## @export
# ellipse_plot.reg_fit_mediation <- function(object, fit = c("mx", "ymx"),
#                                            level = 0.975, npoints = 100,
#                                            ...) {
#   # selector for requested model fit
#   postfix <- match.arg(fit)
#   fit <- paste("fit", postfix, sep = "_")
#   if (postfix == "mx") {
#     y <- object$m
#     x <- object$x
#   } else {
#     stop("not implemented yet")
#   }
#   # create data frame for plotting points
#   df_points <- object$data[, c(x, y)]
#   # add weights in case of robust regression
#   aes_points <- list(x = x, y = y)
#   if (object$robust) {
#     # FIXME: make sure that neighter variable is called 'Weight'
#     df_points$Weight <- weights(object[[fit]], type = "robustness")
#     aes_points$colour <- "Weight"
#     # aes_points$size <- "Weight"
#   }
#   aes_points <- do.call(aes_string, aes_points)
#   # create data frame for plotting ellipse
#   df_ellipse <- tol_ellipse(object, fit = postfix, level = level,
#                                 npoints = npoints)
#   # extract coefficients of regression line
#   coefficients <- coef(object[[fit]])
#   # create plot
#   ggplot() +
#     geom_path(aes_string(x = x, y = y), data = df_ellipse) +
#     geom_point(aes_points, data = df_points) +
#     geom_abline(intercept = coefficients[1], slope = coefficients[2])
# }
