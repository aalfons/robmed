# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------


## internal functions to extract effects from regression fits

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

## extract  bootstrap replicates of the effect(s) for the a path for one
## independent variable
# j_x ......... index of independent variable (relative to number of
#               independent variables)
# indices ..... vector of indices of where to find regression coefficients of
#               m ~ x + covariates in bootstrap replicates, or a list of such
#               indices in case of multiple mediators
# bootstrap ... an object of bootstrap replicates for regression fits as
#               returned by "boot"
# model ....... type of mediation model
extract_boot_a <- function(j_x, indices, bootstrap, model) {
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
  bootstrap$t[, indices_a, drop = FALSE]
}

## extract bootstrap replicates of the effect(s) for the b path
# j_m ......... indices of mediators for which to extract the effect (relative
#               to number of mediators)
# indices ..... vector of indices of where to find regression coefficients of
#               y ~ m + x + covariates in bootstrap replicates
# bootstrap ... an object of bootstrap replicates for regression fits as
#               returned by "boot"
extract_boot_b <- function(j_m, indices, bootstrap) {
  # obtain indices corresponding to effects for b path
  indices_b <- indices[1L + j_m]
  # return bootstrap replicates
  bootstrap$t[, indices_b, drop = FALSE]
}

## extract bootstrap replicates of the effect(s) for the d path
# index_list ... list of index vectors of where to find regression coefficients
#                of m ~ x + covariates in bootstrap replicates
# bootstrap .... an object of bootstrap replicates for regression fits as
#                returned by "boot"
extract_boot_d <- function(index_list, bootstrap) {
  # initializations
  p_m <- length(index_list)
  # obtain indices corresponding to effects for d path
  j_list <- lapply(seq_len(p_m-1L), function(j) 1L + seq_len(j))
  indices_d <- unlist(mapply("[", index_list[-1L], j_list, USE.NAMES = FALSE),
                      use.names = FALSE)
  # return bootstrap replicates
  bootstrap$t[, indices_d, drop = FALSE]
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
# j_x ......... indices of independent variables for which to extract the
#               effect (relative to number of independent variables)
# p_m ......... number of mediators
# indices ..... vector of indices of where to find regression coefficients of
#               y ~ m + x + covariates in bootstrap replicates
# bootstrap ... an object of bootstrap replicates for regression fits as
#               returned by "boot"
extract_boot_direct <- function(j_x, p_m, indices, bootstrap) {
  # obtain indices corresponding to direct effects
  indices_direct <- indices[1L + p_m + j_x]
  # return bootstrap replicates
  bootstrap$t[, indices_direct, drop = FALSE]
}

## extract bootstrap replicates of the direct effect(s)
# j_x ......... indices of independent variables for which to extract the
#               effect (relative to number of independent variables)
# indices ..... vector of indices of where to find regression coefficients of
#               y ~ x + covariates in bootstrap replicates
# bootstrap ... an object of bootstrap replicates for regression fits as
#               returned by "boot"
extract_boot_total <- function(j_x, indices, bootstrap) {
  # obtain indices corresponding to direct effects
  indices_total <- indices[1L + j_x]
  # return bootstrap replicates
  bootstrap$t[, indices_total, drop = FALSE]
}
