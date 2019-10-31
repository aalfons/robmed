# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

## function to compute tolerance ellipses

#' @export
tol_ellipse <- function(object, ...) UseMethod("tol_ellipse")

#' @export
tol_ellipse.test_mediation <- function(object, ...) {
  # simply call methods for mediation model fit
  tol_ellipse(object$fit, ...)
}

#' @export
tol_ellipse.reg_fit_mediation <- function(object, horizontal = NULL,
                                          vertical = NULL, partial = FALSE,
                                          level = 0.975, npoints = 100, ...) {
  # initializations
  if (object$robust && object$median) {
    stop("tolerance ellipse not meaningful for median regression")
  }
  # extract variable names
  x <- object$x
  y <- object$y
  m <- object$m
  # check variable on vertical axis
  if (is.null(vertical)) vertical <- m
  else {
    if (!is.character(vertical) && length(vertical) == 1) {
      stop("only one variable allowed for the vertical axis")
    }
    if (!(vertical %in% m || vertical == y)) {
      stop("variable on the vertical axis must be ",
           "the dependent variable or a mediator")
    }
  }
  # check variable on horizontal axis
  if (is.null(horizontal)) horizontal <- x
  else {
    if (!is.character(horizontal) && length(horizontal) == 1) {
      stop("only one variable allowed for the horizontal axis")
    }
    if (vertical %in% m && horizontal != x) {
      stop("variable on the horizontal axis must be the independent variable")
    } else if (vertical == y && !(horizontal %in% m || horizontal == x)) {
      stop("variable on the horizontal axis must be ",
           "the dependent variable or a mediator")
    }
  }
  # other initializations
  partial <- isTRUE(partial)
  have_mx <- vertical == m && horizontal == x
  # extract model fit
  if (partial || have_mx || object$robust) {
    fit <- if (have_mx) object$fit_mx else object$fit_ymx
  }
  # if applicable, extract intercept and slope
  if (partial || have_mx) {
    coefficients <- coef(fit)
    if (have_mx && !partial) intercept <- unname(coefficients["(Intercept)"])
    else intercept <- 0
    slope <- unname(coefficients[horizontal])
  } else {
    intercept <- NULL
    slope <- NULL
  }
  # extract data to plot
  if (partial) {
    x <- object$data[, horizontal]
    y <- residuals(fit) + slope * x
    data <- data.frame(x = x, y = y)
  } else {
    data <- data.frame(x = object$data[, horizontal],
                       y = object$data[, vertical])
  }
  # obtain location and shape of ellipse
  if (object$robust) {
    # extract weights in case of robust regression
    w <- weights(fit, type = "robustness")
    # compute weighted mean and weighted covariance matrix
    center <- sapply(data, weighted.mean, w = w)
    cov <- weighted.cov(data, w = w)
    # add weights to data frame
    data$Weight <- w
  } else {
    # compute mean and covariance matrix
    center <- colMeans(data)
    cov <- cov(data)
  }
  # compute ellipse
  ellipse <- ellipse(center, cov, level = level, npoints = npoints)
  # return data and ellipse
  list(data = data, ellipse = as.data.frame(ellipse),
       intercept = intercept, slope = slope,
       horizontal = horizontal, vertical = vertical,
       partial = partial, robust = object$robust)
}

# ## @export
# tol_ellipse.cov_fit_mediation <- function(object,
#                                           variables = c("mx", "ym", "yx"),
#                                           level = 0.975, npoints = 100,
#                                           ...) {
#   # initializations
#   if (length(object$m) + length(object$covariates) > 1) {
#     stop("currently only implemented for simple mediation models")
#   }
#   variables <- match.arg(variables)
#   # select which variables to plot
#   if (variables == "mx") select <- c(object$x, object$m)
#   else stop("not implemented yet")
#   # extract covariance matrix
#   center <- object$cov$center[select]
#   cov <- object$cov$cov[select, select]
#   # compute ellipse
#   ellipse <- ellipse(center, cov, level = level, npoints = npoints)
#   # TODO: also return data
#   as.data.frame(ellipse)
# }

# function to compute an ellipse based on center and covariance matrix
ellipse <- function(center, cov, level = 0.975, npoints = 100) {
  # extract scales and correlation
  scale <- sqrt(diag(cov))
  r <- cov[1, 2] / prod(scale)
  # compute ellipse
  d <- acos(r)
  a <- seq(0, 2 * pi, length.out = npoints)
  q <- sqrt(qchisq(level, df = 2))
  x <- q * scale[1] * cos(a + d/2) + center[1]
  y <- q * scale[2] * cos(a - d/2) + center[2]
  xy <- cbind(x, y)
  # add names of variables and return ellipse
  colnames(xy) <- names(center)
  xy
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


#' @import ggplot2
#' @export
ellipse_plot <- function(object, ...) UseMethod("ellipse_plot")

#' @export
ellipse_plot.boot_test_mediation <- function(object, ...) {
  # initial checks
  if (!inherits(object$fit, "reg_fit_mediation")) {
    stop("not implemented for this type of mediation model")
  }
  # call method for mediation model fit
  ellipse_plot(object$fit, ...)
}

#' @export
ellipse_plot.reg_fit_mediation <- function(object, horizontal = NULL,
                                           vertical = NULL, partial = FALSE,
                                           ...) {
  # obtain data to be plotted and tolerance ellipse
  df_list <- tol_ellipse(object, horizontal = horizontal, vertical = vertical,
                         partial = partial, ...)
  # define aesthetic mapping for plotting points
  if (df_list$robust) aes_data <- aes_string(x = "x", y = "y", fill = "Weight")
  else aes_data <- aes_string(x = "x", y = "y")
  # create plot
  p <- ggplot() +
    geom_path(aes_string(x = "x", y = "y"), data = df_list$ellipse) +
    geom_point(aes_data, data = df_list$data, shape = 21)
  # add line representing (partial) effect
  if (!is.null(df_list$intercept) && !is.null(df_list$slope)) {
    p <- p + geom_abline(intercept = df_list$intercept, slope = df_list$slope)
  }
  # add nice labels
  if (df_list$partial) ylab <- paste("Partial residuals of", df_list$vertical)
  else ylab <- df_list$vertical
  p <- p + labs(x = df_list$horizontal, y = ylab)
  # add color gradient for weights
  if (df_list$robust) {
    p <- p + scale_fill_gradient(low = "transparent", high = "black")
  }
  # return plot
  p
}


# foo <- tol_ellipse(test, horizontal = "ValueDiversity", vertical = "TaskConflict")
# df_points <- cbind(foo$data, Weight = foo$weights)
# coefficients <- coef(test$fit$fit_mx)
# ggplot() +
#   geom_point(aes(x = x, y = y, color = Weight),
#              data = df_points) +
#   geom_path(aes(x = x, y = y), data = foo$ellipse) +
#   geom_abline(intercept = coefficients[1], slope = coefficients["ValueDiversity"]) +
#   labs(x = "ValueDiversity", y = "TaskConflict")
#
# bar <- tol_ellipse(test, horizontal = "ValueDiversity", vertical = "TeamCommitment", partial = TRUE)
# df_points <- cbind(bar$data, Weight = bar$weights)
# coefficients <- coef(test$fit$fit_ymx)
# ggplot() +
#   geom_point(aes(x = x, y = y, color = Weight),
#              data = df_points) +
#   geom_path(aes(x = x, y = y), data = bar$ellipse) +
#   geom_abline(intercept = 0, slope = coefficients[bar$horizontal]) +
#   labs(x = bar$horizontal, y = paste("Partial residual of", bar$vertical))
#
# bar <- tol_ellipse(test, horizontal = "TaskConflict", vertical = "TeamCommitment", partial = FALSE)
# df_points <- cbind(bar$data, Weight = bar$weights)
# coefficients <- coef(test$fit$fit_ymx)
# ggplot() +
#   geom_point(aes(x = x, y = y, color = Weight),
#              data = df_points) +
#   geom_path(aes(x = x, y = y), data = bar$ellipse) +
#   labs(x = bar$horizontal, y = bar$vertical)
