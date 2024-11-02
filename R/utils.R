# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------


## internal functions required for fast-and-robust bootstrap

## get control arguments for psi function as used in a given model fit
get_psi_control <- function(object) object$control[c("tuning.psi", "psi")]

## compute matrix for linear correction
# (see Salibian-Barrera & Van Aelst, 2008)
# The definition of the weigths in Salibian-Barrera & Van Aelst (2008) does not
# include the residual scale, whereas the robustness weights in lmrob() do.
# Hence the residual scale shows up in Equation (16) of Salibian-Barrera & Van
# Aelst (2008), but here the residual scale is already included in the weights.
correction_matrix <- function(X, weights, residuals, scale, control) {
  tmp <- Mpsi(residuals/scale, cc=control$tuning.psi, psi=control$psi, deriv=1)
  solve(crossprod(X, tmp * X)) %*% crossprod(weights * X)
}


## internal functions to compute contrasts

# compute contrasts
get_contrasts <- function(x, combinations = NULL, type = "estimates") {
  # compute combinations if not supplied
  if (is.null(combinations)) {
    indices <- get_contrast_indices(x)
    combinations <- combn(indices, 2, simplify = FALSE)
  }
  # define function to compute contrasts
  if (type == "estimates") fun <- get_original_contrast
  else if (type == "absolute") fun <- get_absolute_contrast
  else stop("type of contrasts not implemented")
  # compute contrasts
  contrasts <- sapply(combinations, fun, x = x)
  if (is.vector(contrasts)) {
    # When computing contrasts of the indirect effect estimates on the original
    # data, add names to the contrasts.  This is not necessary when computing
    # contrasts of bootstrap replicates of the indirect effects, since the
    # "boot" object doesn't use names, and the names of the bootstrap estimates
    # are copied from the estimates on the original data.
    n_contrasts <- length(contrasts)
    names(contrasts) <- get_contrast_names(n_contrasts)
  }
  contrasts
}

# obtain indices to be used for computing contrasts
get_contrast_indices <- function(x) {
  if (is.matrix(x) || is.data.frame(x)) seq_len(ncol(x))
  else seq_along(x)
}

# obtain names for contrasts
get_contrast_names <- function(n) {
  if (n > 1) paste0("Contrast", seq_len(n))
  else "Contrast"
}

# compute contrasts as differences of values
get_original_contrast <- function(j, x) {
  if (is.matrix(x) || is.data.frame(x)) x[, j[1]] - x[, j[2]]
  else x[j[1]] - x[j[2]]
}

# compute contrasts as differences of absolute values
get_absolute_contrast <- function(j, x) {
  if (is.matrix(x) || is.data.frame(x)) abs(x[, j[1]]) - abs(x[, j[2]])
  else abs(x[j[1]]) - abs(x[j[2]])
}

# obtain information on how contrasts are computed
get_contrast_info <- function(names, type = "estimates") {
  # compute combinations of names
  combinations <- combn(names, 2, simplify = FALSE)
  n_contrasts <- length(combinations)
  # obtain labels for contrasts
  labels <- get_contrast_names(n_contrasts)
  # obtain information on contrasts
  if (type == "estimates") {
    fun <- function(names) paste(names, collapse = " - ")
  } else if (type == "absolute") {
    fun <- function(names) {
      paste(paste0("|", names, "|"), collapse = " - ")
    }
  } else stop("type of contrasts not implemented")
  # return information on contrasts
  data.frame(Label = labels, Definition = sapply(combinations, fun))
}


## utility function to get names of coefficients
get_effect_names <- function(..., effects = list(...), sep = "_") {
  # loop over effects and get corresponding names
  effect_names <- mapply(function(effect, name) {
    if (is.null(effect)) effect_name <- NULL
    else {
      # define label
      if (name == "total") label <- "Total"
      else if (name == "direct") label <- "Direct"
      else if (name == "indirect") label <- "Indirect"
      else label <- name
      # set label or add label as prefix
      if (length(effect) == 1L) effect_name <- label
      else effect_name <- paste(label, names(effect), sep = sep)
    }
    # return effect name
    effect_name
  }, effect = effects, name = names(effects),
  SIMPLIFY = FALSE, USE.NAMES = FALSE)
  # return effect names
  unlist(effect_names, use.names = FALSE)
}


## The function for bootstrap replicates is required to return a vector.  This
## means that the columns of the bootstrap replicates contain coefficient
## estimates from different models.  This utility function returns the indices
## that correspond to the respective models, which makes it easier to extract
## the desired coefficients.
get_index_list <- function(p_x, p_m, p_covariates, model = "parallel",
                           fit_yx = TRUE) {
  # numbers of coefficients in different models
  if (p_m > 1L && model == "serial") {
    p_mx <- 1L + seq.int(0L, p_m-1L) + p_x + p_covariates
  } else p_mx <- rep.int(1L + p_x + p_covariates, p_m)
  p_ymx <- 1L + p_m + p_x + p_covariates
  p_yx <- if (fit_yx) 1L + p_x + p_covariates else 0L
  p_all <- sum(p_mx, p_ymx, p_yx)
  # indices of vector for each bootstrap replicate
  indices <- seq_len(p_all)
  # the first columns correspond to models m ~ x + covariates
  first <- 1L
  if (p_m == 1L) {
    indices_mx <- seq.int(from = first, length.out = p_mx)
  } else {
    first_mx <- first + c(0L, cumsum(p_mx[-p_m]))
    indices_mx <- mapply(seq.int, from = first_mx, length.out = p_mx,
                         SIMPLIFY = FALSE, USE.NAMES = FALSE)
  }
  # the next columns correspond to model y ~ m + x + covariates
  first <- first + sum(p_mx)
  indices_ymx <- seq.int(from = first, length.out = p_ymx)
  # if applicable, the last column corresponds to model y ~ x + covariates
  if (fit_yx) {
    first <- first + p_ymx
    indices_yx <- seq.int(from = first, length.out = p_yx)
  } else indices_yx <- integer()
  # return list of indices
  list(fit_mx = indices_mx, fit_ymx = indices_ymx, fit_yx = indices_yx)
}


## internal functions related to indirect effects

# obtain names for indirect effects
get_indirect_names <- function(x, m, model = "parallel", sep = "->") {
  # initializations
  p_x <- length(x)
  p_m <- length(m)
  # construct names for indirect effects
  if (p_m == 1L) {
    # single mediator: only use names in case of multiple independent variables
    names <- if (p_x > 1L) x
  } else {
    # prepare names involving mediators
    if (model == "serial") {
      # currently only implemented for two or three hypothesized mediators
      if (p_m == 2L) {
        # two serial mediators
        names <- c(m, paste(m, collapse = sep))
      } else {
        # three serial mediators
        names <- c(m, paste(m[1L], m[2L], sep = sep),
                   paste(m[1L], m[3L], sep = sep),
                   paste(m[2L], m[3L], sep = sep),
                   paste(m, collapse = sep))
      }
    } else names <- m  # parallel mediators
  }
  # return names
  names
}

# obtain labels for indirect effects in printed output of bootstrap tests
get_indirect_labels <- function(p_m, model = "parallel") {
  # determine number of indirect paths (for a given independent variable)
  if (model == "serial") {
    # currently only implemented for two or three hypothesized mediators
    nr_indirect <- if (p_m == 2L) 3L else 7L
  } else nr_indirect <- p_m
  # return labels
  paste0("Indirect", seq_len(nr_indirect))
}


## check whether an object corresponds to a robust model fit

#' @noRd
is_robust <- function(object) UseMethod("is_robust")

#' @noRd
is_robust.reg_fit_mediation <- function(object) {
  robust <- object$robust
  if (is.character(robust)) TRUE else robust
}

#' @noRd
is_robust.cov_fit_mediation <- function(object) object$robust
