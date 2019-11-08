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
  have_mx <- vertical == m && horizontal == x && length(object$covariates) == 0
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
    line <- data.frame(intercept = intercept, slope = slope)
  } else line <- NULL
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
    cov <- weighted.cov(data, w = w)  # FIXME: multiply with correction factor
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
  out <- list(data = data, ellipse = as.data.frame(ellipse), line = line,
              horizontal = horizontal, vertical = vertical, partial = partial,
              robust = object$robust)
  class(out) <- "tol_ellipse"
  out
}

#' @export
tol_ellipse.cov_fit_mediation <- function(object, horizontal = NULL,
                                          vertical = NULL, partial = FALSE,
                                          level = 0.975, npoints = 100, ...) {
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
    if (vertical != m && vertical != y) {
      stop("variable on the vertical axis must be ",
           "the dependent variable or the mediator")
    }
  }
  # check variable on horizontal axis
  if (is.null(horizontal)) horizontal <- x
  else {
    if (!is.character(horizontal) && length(horizontal) == 1) {
      stop("only one variable allowed for the horizontal axis")
    }
    if (vertical == m && horizontal != x) {
      stop("variable on the horizontal axis must be the independent variable")
    } else if (vertical == y && (horizontal != m && horizontal != x)) {
      stop("variable on the horizontal axis must be ",
           "the dependent variable or a mediator")
    }
  }
  # other initializations
  partial <- isTRUE(partial)
  select <- c(horizontal, vertical)
  # extract covariance fit
  fit <- object$cov
  # extract information to be plotted
  if (vertical == m) {
    # extract location and shape of ellipse
    center <- fit$center[select]
    cov <- fit$cov[select, select]
    # compute intercept and slope
    slope <- unname(cov[vertical, horizontal] / cov[horizontal, horizontal])
    intercept <- unname(center[vertical] - slope * center[horizontal])
    # extract data to plot
    data <- data.frame(x = object$data[, horizontal],
                       y = object$data[, vertical])
    # if requested, make correction for partial residual plot
    if (partial) {
      data$y <- data$y - intercept
      center[vertical] <- center[vertical] - intercept
      intercept <- 0
    }
    # combine information for line representing (partial) effect
    line <- data.frame(intercept = intercept, slope = slope)
  } else {
    # dependent variable on the vertical axis
    if (partial) {
      # extract all means and full covariance matrix
      center <- fit$center
      cov <- fit$cov
      # compute regression coefficient
      det <- cov[x, x] * cov[m, m] - cov[m, x]^2
      b <- (-cov[m, x] * cov[y, x] + cov[x, x] * cov[y, m]) / det
      direct <- (cov[m, m] * cov[y, x] - cov[m, x] * cov[y, m]) / det
      i_y <- unname(center[y] - b * center[m] - direct * center[x])
      # compute data to plot
      if (horizontal == x) {
        # compute partial residuals
        residuals <- object$data[, vertical] - i_y - b * object$data[, m]
        # adjust location and shape of ellipse
        center[vertical] <- center[vertical] - i_y - b * center[m]
        cov[y, y] <- cov[y, y] + b^2 * cov[m, m] + 2 * b * cov[y, m]
        cov[y, x] <- cov[y, x] - b * cov[m, x]
      } else {
        # compute partial residuals
        residuals <- object$data[, vertical] - i_y - direct * object$data[, x]
        # adjust location and shape of ellipse
        center[vertical] <- center[vertical] - i_y - direct * center[x]
        cov[y, y] <- cov[y, y] + direct^2 * cov[x, x] + 2 * direct * cov[y, x]
        cov[y, m] <- cov[y, m] - direct * cov[x, m]
      }
      data <- data.frame(x = object$data[, horizontal], y = residuals)
      # extract location and shape of ellipse
      center <- center[select]
      cov <- cov[select, select]
      # compute intercept and slope
      intercept <- 0
      slope <- if (horizontal == x) unname(direct) else unname(b)
      # combine information for line representing (partial) effect
      line <- data.frame(intercept = intercept, slope = slope)
    } else {
      # extract location and shape of ellipse
      center <- fit$center[select]
      cov <- fit$cov[select, select]
      # extract data to plot
      data <- data.frame(x = object$data[, horizontal],
                         y = object$data[, vertical])
      # do not plot line for effect
      line <- NULL
    }
  }
  # in case of robust covariance matrix, add interpretable robustness weights
  if (object$robust) data$Weight <- weights(fit, type = "relative")
  # compute ellipse
  ellipse <- ellipse(center, cov, level = level, npoints = npoints)
  colnames(ellipse) <- c("x", "y")
  # return data and ellipse
  out <- list(data = data, ellipse = as.data.frame(ellipse), line = line,
              horizontal = horizontal, vertical = vertical, partial = partial,
              robust = object$robust)
  class(out) <- "tol_ellipse"
  out
}

#' @export
tol_ellipse.list <- function(object, ...) {
  ## initializations
  is_test <- sapply(object, inherits, "test_mediation")
  is_fit <- sapply(object, inherits, "fit_mediation")
  object <- object[is_test | is_fit]
  if(length(object) == 0L) {
    stop('no objects inheriting from class "test_mediation" or "fit_mediation"')
  }
  # TODO: check that variables are the same
  # check names of list elements
  methods <- names(object)
  if(is.null(methods)) methods <- seq_along(object)
  else {
    replace <- methods == "" | duplicated(methods)
    methods[replace] <- seq_along(object)[replace]
  }
  # compute tolerance ellipse for each list element
  tol_ellipse_list <- lapply(object, tol_ellipse, ...)
  # check for properties of tolerance ellipses
  partial <- tol_ellipse_list[[1]]$partial
  robust <- sapply(tol_ellipse_list, "[[", "robust")
  sum_robust <- sum(robust)
  # combine information of tolerance ellipses
  # TODO: should there be an argument to define what to return?
  if (partial || sum_robust > 1L) {
    # combine data sets
    # For a partial residuals plot, those partial residuals are different for
    # every method.  Otherwise if there are multiple robust methods, there are
    # different weights.
    if (sum_robust == 0L) {
      # add column identifying the method to data for scatter plots
      data_list <- mapply(function(ellipse, method) {
        cbind(Method = method, ellipse$data)
      }, ellipse = tol_ellipse_list, method = methods,
      SIMPLIFY = FALSE, USE.NAMES = FALSE)
    } else {
      # if there is a robust method, add column with weights for nonrobust
      # methods as well
      data_list <- mapply(function(ellipse, method) {
        if (ellipse$robust) cbind(Method = method, ellipse$data)
        else cbind(Method = method, ellipse$data, Weight = 1)
      }, ellipse = tol_ellipse_list, method = methods,
      SIMPLIFY = FALSE, USE.NAMES = FALSE)
    }
    # combine data sets for scatter plots
    data <- do.call(rbind, data_list)
  } else if (sum_robust == 1L) data <- tol_ellipse_list[[which(robust)]]$data
  else data <- tol_ellipse_list[[1]]$data
  # combine data for ellipses
  ellipse_list <- mapply(function(ellipse, method) {
    cbind(Method = method, ellipse$ellipse)
  }, ellipse = tol_ellipse_list, method = methods,
  SIMPLIFY = FALSE, USE.NAMES = FALSE)
  ellipses <- do.call(rbind, ellipse_list)
  # combine data for lines representing (partial) effects
  line_list <- mapply(function(ellipse, method) {
    line <- ellipse$line
    if (!is.null(line)) cbind(Method = method, line)
  }, ellipse = tol_ellipse_list, method = methods,
  SIMPLIFY = FALSE, USE.NAMES = FALSE)
  lines <- do.call(rbind, line_list)
  # return results
  out <- list(data = data, ellipse = ellipses, line = lines,
              horizontal = tol_ellipse_list[[1]]$horizontal,
              vertical = tol_ellipse_list[[1]]$vertical,
              partial = partial, robust = robust, methods = methods)
  class(out) <- "tol_ellipse"
  out
}


## workhorse function to compute ellipse based on center and covariance matrix
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
