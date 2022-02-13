# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------


## internal functions to extract effects from regression fits

## extract effects of the mediation model from regressions
# TODO: Add argument that allows to select which effects to return.  This can
#       speed things up slightly for retest(), where only the indirect effects
#       need to be extracted in the case of updated contrasts.
extract_effects <- function(x, m, family = "gaussian", model = "simple",
                            contrast = FALSE, fit_mx, fit_ymx, fit_yx = NULL) {
  # initializations
  p_x <- length(x)
  p_m <- length(m)
  have_contrast <- is.character(contrast)
  # extract effects for b path
  b <- extract_b(m, fit_ymx)
  # in case of serial mediators, extract effects for d path
  d <- if (model == "serial") extract_d(m, fit_mx)
  # extract effects for a path and compute indirect effects
  if (p_x == 1L) {
    # extract effects for a path
    a <- extract_a(x, fit_mx)
    # compute indirect effects
    indirect <- compute_indirect(a, b, d)
    if (p_m > 1L) {
      # compute total indirect effect
      indirect_total <- sum(indirect)
    }
  } else {
    # extract effects for a path
    a <- lapply(x, extract_a, fit_mx)
    # compute indirect effects
    indirect <- lapply(a, compute_indirect, b, d)
    if (p_m == 1L) {
      # for models with multiple independent variables but a single mediator,
      # ensure we have a vector of indirect effects for computing contrasts
      indirect <- unlist(indirect, use.names = FALSE)
      names(indirect) <- x
    } else {
      # compute total indirect effect
      indirect_total <- sapply(indirect, sum)
      names(indirect_total) <- x
    }
  }
  # extract direct effect(s)
  direct <- extract_direct(x, fit_ymx)
  # compute or extract total effect(s)
  if (family == "gaussian" && is.null(fit_yx)) {
    # compute the total effect based on the indirect and direct effects
    total <- if (p_m == 1L) indirect + direct else indirect_total + direct
  } else {
    # extract total effect from model fit (NAs if model fit is not available)
    total <- extract_total(x, fit_yx)
  }
  # if requested, compute contrasts
  if (have_contrast) {
    # fit_mediation() ensures that this is only true if there are multiple
    # independent variables or multiple mediators
    if (p_x == 1L || p_m == 1L) {
      # multiple independent variables or multiple mediators, but not both:
      # compute contrasts between all individual indirect effects
      contrasts <- get_contrasts(indirect, type = contrast)
    } else {
      # multiple independent variables and multiple mediators: compute
      # contrasts separately for each independent variable such that the
      # number of comparisons doesn't explode
      contrasts <- lapply(indirect, get_contrasts, type = contrast)
    }
  } else contrasts <- NULL
  # make sure we have vector for a path
  if (p_x > 1L) {
    a <- unlist(a, use.names = FALSE)
    if (p_m == 1L) names(a) <- x
    else names(a) <- sapply(x, paste, m, sep = "->", USE.NAMES = FALSE)
  }
  # finalize vector of indirect effects
  if (p_m == 1L) {
    # For multiple independent variables and a single mediator, the list of
    # indirect effects has already been converted to a vector before computing
    # contrasts.  There is also no total indirect effect to add.
    if (p_x > 1L) indirect <- c(indirect, contrasts)
  } else {
    # multiple mediators: combine total indirect effect, individual indirect
    # effects, and (if requested) contrasts
    if (p_x == 1L) {
      # for a single independent variable, combine vectors
      indirect <- c(Total = indirect_total, indirect, contrasts)
    } else {
      ## for multiple independent variables, lists need to be combined first
      # add name of independent variable to names of total indirect effects
      names(indirect_total) <- paste(names(indirect_total), "Total", sep = "_")
      # add name of independent variable to names of individual indirect effects
      indirect <- mapply(function(current_x, current_indirect) {
        names(current_indirect) <- paste(current_x, names(current_indirect),
                                         sep = "->")
        current_indirect
      }, x, indirect, SIMPLIFY = FALSE, USE.NAMES = FALSE)
      # add name of independent variable to names of contrasts
      if (!is.null(contrasts)) {
        contrasts <- mapply(function(current_x, current_contrasts) {
          names(current_contrasts) <- paste(current_x, names(current_contrasts),
                                            sep = "_")
          current_contrasts
        }, x, contrasts, SIMPLIFY = FALSE, USE.NAMES = FALSE)
      }
      # combine lists
      indirect <- lapply(seq_len(p_x), function(j) {
        c(indirect_total[j], indirect[[j]], contrasts[[j]])
      })
      ## combine everything into one vector
      indirect <- unlist(indirect)
    }
  }
  # return list of effects
  tmp <- list(a = a, b = b)
  if (model == "serial") tmp$d <- d
  c(tmp, list(total = total, direct = direct, indirect = indirect))
}


# extract effect(s) for a path for one independent variable
extract_a <- function(x, fit) {
  if (inherits(fit, "list")) {
    # multiple mediators
    sapply(fit, function(current_fit) unname(coef(current_fit)[x]))
  } else {
    # only one mediator
    unname(coef(fit)[x])
  }
}

# extract effect(s) for b path
extract_b <- function(m, fit) {
  # initializations
  p_m <- length(m)
  # extract effect(s)
  if (p_m == 1L) unname(coef(fit)[m])
  else coef(fit)[m]
}

# extract effect(s) for d path
extract_d <- function(m, fit_list) {
  # initializations
  p_m <- length(m)
  # currently only implemented for two or three hypothesized mediators
  if (p_m == 2L) {
    # two serial mediators
    d <- unname(coef(fit_list[[m[2L]]])[m[1L]])
  } else if (p_m == 3L) {
    # three serial mediators
    d <- c(coef(fit_list[[m[2L]]])[m[1L]], coef(fit_list[[m[3L]]])[m[1L:2L]])
    names(d) <- paste(names(d), m[c(2L, 3L, 3L)], sep = "->")
  } else d <- NULL
  # return effect(s)
  d
}

# compute indirect effect(s) for one independent variable
compute_indirect <- function(a, b, d = NULL) {
  # first-level indirect effects: the computation is the same for
  # simple mediation as well as for parallel or serial mediators
  indirect <- a * b
  # compute higher-level indirect effects for serial mediators
  if (!is.null(d)) {
    # extract names of mediators
    m <- names(indirect)
    p_m <- length(m)
    # currently only implemented for two or three hypothesized mediators
    if (p_m == 2L) {
      # two serial mediators
      tmp <- a[1] * d * b[2]
      names(tmp) <- paste(m, collapse = "->")
    } else {
      # three serial mediators
      tmp <- c(a[1] * d[1] * b[2], a[1] * d[2] * b[3],
               a[2] * d[3] * b[3], a[1] * d[1] * d[3] * b[3])
      names(tmp) <- c(names(d), paste(m, collapse = "->"))
    }
    # add to existing vector for indirect effects
    indirect <- c(indirect, tmp)
  }
  # return indirect effect(s)
  indirect
}

# extract direct effect(s)
extract_direct <- function(x, fit) {
  # initializations
  p_x <- length(x)
  # extract effect(s)
  if (p_x == 1L) unname(coef(fit)[x])
  else coef(fit)[x]
}

# extract total effect(s)
extract_total <- function(x, fit = NULL) {
  # initializations
  p_x <- length(x)
  # if available, extract effect(s)
  if (is.null(fit)) {
    # effect(s) not available
    if (p_x == 1L) total <- NA_real_
    else {
      total <- rep.int(NA_real_, p_x)
      names(total) <- x
    }
  } else {
    # extract effect(s)
    if (p_x == 1L) total <- unname(coef(fit)[x])
    else total <- coef(fit)[x]
  }
  # return effect(s)
  total
}


## internal functions to extract bootstrap replicates of effects

## extract bootstrap replicates for the effects in the mediation model
# fit .......... an object of class "reg_fit_mediation" containing a mediation
#                model fit via regressions
# boot ......... an object of bootstrap replicates for regression fits as
#                returned by function boot()
# index_list ... list of index vectors that indicate which columns of the
#                bootstrap replicates correspond to the respective models
# TODO: Add argument that allows to select for which effects to return the
#       bootstrap replicates.  This can speed things up, for instance for
#       summary(), where only the bootstrap replicates of the total effect
#       need to be extracted (and we don't need to compute, e.g., contrasts).
extract_boot <- function(fit, boot, index_list = NULL) {
  # initializations
  p_x <- length(fit$x)
  p_m <- length(fit$m)
  p_covariates <- length(fit$covariates)
  family <- fit$family
  model <- fit$model
  contrast <- fit$contrast                 # FALSE or type of contrast
  have_contrast <- is.character(contrast)  # indicates whether we have contrasts
  have_yx <- !is.null(fit$fit_yx)
  # if not supplied, get list of index vectors that indicate which columns of
  # the bootstrap replicates correspond to the respective models
  if (is.null(index_list)) {
    index_list <- get_index_list(p_x, p_m, p_covariates, model = model,
                                 fit_yx = family != "gaussian" && have_yx)
  }
  # define useful sequences
  seq_x <- seq_len(p_x)
  seq_m <- seq_len(p_m)
  # extract effects for b path
  boot_b <- extract_boot_b(seq_m, indices = index_list$fit_ymx, boot = boot)
  # in case of serial mediators, extract effects for d path
  if (model == "serial") {
    boot_d <- extract_boot_d(index_list$fit_mx, boot = boot)
  } else boot_d <- NULL
  # extract and compute bootstrap replicates for a path and indirect effects
  if (p_x == 1L) {
    # extract bootstrap replicates of effects for a path
    boot_a <- extract_boot_a(1L, indices = index_list$fit_mx,
                             boot = boot, model = model)
    # compute bootstrap replicates of indirect effects
    boot_indirect <- compute_boot_indirect(boot_a, boot_b, boot_d)
    if (p_m > 1L) {
      # compute bootstrap replicates of total indirect effect
      boot_indirect_total <- rowSums(boot_indirect)
    }
  } else {
    # extract bootstrap replicates of effects for a path
    boot_a <- lapply(seq_x, extract_boot_a, indices = index_list$fit_mx,
                     boot = boot, model = model)
    # compute bootstrap replicates of indirect effects
    boot_indirect <- lapply(boot_a, compute_boot_indirect, boot_b, boot_d)
    if (p_m == 1L) {
      # for models with multiple independent variables but a single mediator,
      # ensure we have a matrix of indirect effects for computing contrasts
      boot_indirect <- do.call(cbind, boot_indirect)
    } else {
      # compute total indirect effect
      boot_indirect_total <- sapply(boot_indirect, rowSums)
    }
  }
  # extract direct effect(s)
  boot_direct <- extract_boot_direct(seq_x, p_m = p_m,
                                     indices = index_list$fit_ymx,
                                     boot = boot)
  # compute or extract total effect(s)
  if (family == "gaussian") {
    # compute the total effect based on the indirect and direct effects
    if (p_m == 1L) boot_total <- boot_indirect + boot_direct
    else boot_total <- boot_indirect_total + boot_direct
  } else if (have_yx) {
    # extract total effect(s)
    boot_total <- extract_boot_total(seq_x, indices = index_list$fit_yx,
                                     boot = boot)
  } else {
    # total effect(s) not available
    boot_total <- matrix(NA_real_, nrow = nrow(boot$t), ncol = p_x)
  }
  # if requested, compute bootstrap estimates of contrasts
  if (have_contrast) {
    # fit_mediation() ensures that this is only true if there are multiple
    # independent variables or multiple mediators
    if (p_x == 1L || p_m == 1L) {
      # multiple independent variables or multiple mediators, but not both:
      # compute contrasts between all individual indirect effects
      boot_contrasts <- get_contrasts(boot_indirect, type = contrast)
    } else {
      # multiple independent variables and multiple mediators: compute
      # contrasts separately for each independent variable such that the
      # number of comparisons doesn't explode
      boot_contrasts <- lapply(boot_indirect, get_contrasts, type = contrast)
    }
  } else boot_contrasts <- NULL
  # make sure we have matrix of bootstrap replicates for a path
  if (p_x > 1L) boot_a <- do.call(cbind, boot_a)
  # compute bootstrap estimates for a path
  a <- colMeans(boot_a, na.rm = TRUE)
  # finalize matrix of bootstrap replicates of indirect effects
  if (p_m == 1L) {
    # For multiple independent variables and a single mediator, the list of
    # bootstrap replicates of the indirect effects has already been converted
    # to a matrix before computing contrasts.  There are also no bootstrap
    # replicates of the total indirect effect to add.
    if (p_x > 1L) boot_indirect <- cbind(boot_indirect, boot_contrasts)
  } else {
    # multiple mediators: combine bootstrap replicates of total indirect
    # effect, individual indirect effects, and (if requested) contrasts
    if (p_x == 1L) {
      # for a single independent variable, combine vectors
      boot_indirect <- cbind(boot_indirect_total, boot_indirect,
                             boot_contrasts, deparse.level = 0)
    } else {
      # combine lists
      boot_indirect <- lapply(seq_x, function(j) {
        cbind(boot_indirect_total[, j, drop = FALSE], boot_indirect[[j]],
              boot_contrasts[[j]])
      })
      ## combine everything into one matrix
      boot_indirect <- do.call(cbind, boot_indirect)
    }
  }
  # copy column names for bootstrap replicates from estimates on original data
  colnames(boot_a) <- names(fit$a)
  colnames(boot_b) <- names(fit$b)
  colnames(boot_total) <- names(fit$total)
  colnames(boot_direct) <- names(fit$direct)
  colnames(boot_indirect) <- names(fit$indirect)
  # return list of bootstrap replicates for the different effects
  tmp <- list(a = boot_a, b = boot_b)
  if (model == "serial") {
    colnames(boot_d) <- names(fit[["d"]])
    tmp$d <- boot_d
  }
  c(tmp, list(total = boot_total, direct = boot_direct,
              indirect = boot_indirect))
}

## extract  bootstrap replicates of the effect(s) for the a path for one
## independent variable
# j_x ....... index of independent variable (relative to number of
#             independent variables)
# indices ... vector of indices of where to find regression coefficients of
#             m ~ x + covariates in bootstrap replicates, or a list of such
#             indices in case of multiple mediators
# boot ...... an object of bootstrap replicates for regression fits as returned
#             by function boot()
# model ..... type of mediation model
extract_boot_a <- function(j_x, indices, boot, model) {
  # obtain indices corresponding to effects for a path
  if (inherits(indices, "list")) {
    # multiple mediators
    if (model == "serial") {
      # serial mediators
      p_m <- length(indices)
      indices_a <- mapply("[", indices, seq_len(p_m) + j_x, USE.NAMES = FALSE)
    } else {
      # parallel mediators
      indices_a <- sapply(indices, "[", 1L + j_x)
    }
  } else {
    # only one mediator
    indices_a <- indices[1L + j_x]
  }
  # return bootstrap replicates
  boot$t[, indices_a, drop = FALSE]
}

## extract bootstrap replicates of the effect(s) for the b path
# j_m ....... indices of mediators for which to extract the effect (relative
#             to number of mediators)
# indices ... vector of indices of where to find regression coefficients of
#             y ~ m + x + covariates in bootstrap replicates
# boot ...... an object of bootstrap replicates for regression fits as returned
#             by function boot()
extract_boot_b <- function(j_m, indices, boot) {
  # obtain indices corresponding to effects for b path
  indices_b <- indices[1L + j_m]
  # return bootstrap replicates
  boot$t[, indices_b, drop = FALSE]
}

## extract bootstrap replicates of the effect(s) for the d path
# index_list ... list of index vectors of where to find regression coefficients
#                of m ~ x + covariates in bootstrap replicates
# boot ......... an object of bootstrap replicates for regression fits as
#                returned by function boot()
extract_boot_d <- function(index_list, boot) {
  # initializations
  p_m <- length(index_list)
  # obtain indices corresponding to effects for d path
  j_list <- lapply(seq_len(p_m-1L), function(j) 1L + seq_len(j))
  indices_d <- unlist(mapply("[", index_list[-1L], j_list, USE.NAMES = FALSE),
                      use.names = FALSE)
  # return bootstrap replicates
  boot$t[, indices_d, drop = FALSE]
}

## compute bootstrap replicates of indirect effect(s) for one independent
## variable
compute_boot_indirect <- function(boot_a, boot_b, boot_d = NULL) {
  # compute bootstrap replicates of first-level indirect effects: the
  # computation is the same for simple mediation as well as for parallel
  # or serial mediators
  boot_indirect <- boot_a * boot_b
  # compute bootstrap replicates of higher-level indirect effects for
  # serial mediators
  if (!is.null(boot_d)) {
    # extract number of mediators
    p_m <- ncol(boot_indirect)
    # currently only implemented for two or three hypothesized mediators
    if (p_m == 2L) {
      # two serial mediators
      tmp <- boot_a[, 1] * boot_d * boot_b[, 2]
    } else {
      # three serial mediators
      tmp <- cbind(boot_a[, 1] * boot_d[, 1] * boot_b[, 2],
                   boot_a[, 1] * boot_d[, 2] * boot_b[, 3],
                   boot_a[, 2] * boot_d[, 3] * boot_b[, 3],
                   boot_a[, 1] * boot_d[, 1] * boot_d[, 3] * boot_b[, 3],
                   deparse.level = 0)
    }
    # add to existing matrix for bootstrap replicates of indirect effects
    boot_indirect <- cbind(boot_indirect, tmp)
  }
  # return bootstrap replicates of indirect effect(s)
  boot_indirect
}

## extract bootstrap replicates of the direct effect(s)
# j_x ....... indices of independent variables for which to extract the
#             effect (relative to number of independent variables)
# p_m ....... number of mediators
# indices ... vector of indices of where to find regression coefficients of
#             y ~ m + x + covariates in bootstrap replicates
# boot ...... an object of bootstrap replicates for regression fits as returned
#             by function boot()
extract_boot_direct <- function(j_x, p_m, indices, boot) {
  # obtain indices corresponding to direct effects
  indices_direct <- indices[1L + p_m + j_x]
  # return bootstrap replicates
  boot$t[, indices_direct, drop = FALSE]
}

## extract bootstrap replicates of the direct effect(s)
# j_x ....... indices of independent variables for which to extract the
#             effect (relative to number of independent variables)
# indices ... vector of indices of where to find regression coefficients of
#             y ~ x + covariates in bootstrap replicates
# boot ...... an object of bootstrap replicates for regression fits as returned
#             by function boot()
extract_boot_total <- function(j_x, indices, boot) {
  # obtain indices corresponding to direct effects
  indices_total <- indices[1L + j_x]
  # return bootstrap replicates
  boot$t[, indices_total, drop = FALSE]
}
