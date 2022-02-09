# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

# internal function to extract effect for a path for one independent variable
extract_a <- function(x, fit) {
  if (inherits(fit, "list")) {
    # multiple mediators
    sapply(fit, function(current_fit) unname(coef(current_fit)[x]))
  } else {
    # only one mediator
    unname(coef(fit)[x])
  }
}

# internal function to extract effect(s) for b path
extract_b <- function(m, fit) {
  # initializations
  p_m <- length(m)
  # extract effect(s)
  if (p_m == 1L) unname(coef(fit)[m])
  else coef(fit)[m]
}

# internal function to extract effect(s) for d path
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

# internal function to compute indirect effect(s) for one independent variable
compute_indirect <- function(a, b, d = NULL) {
  # first-level indirect effects: computation is the same for
  # simple mediation and parallel or serial mediators
  indirect <- a * b
  # compute higher-level indirect effects for serial mediators
  if (!is.null(d)) {
    # extract names of mediators
    m <- names(indirect)
    # currently only implemented for two or three hypothesized mediators
    if (length(d) == 1L) {
      # two serial mediators
      indirect_serial <- a[1] * d * b[2]
      names(indirect_serial) <- paste(m, collapse = "->")
    } else {
      # three serial mediators
      indirect_serial <- c(a[1] * d[1] * b[2], a[1] * d[2] * b[3],
                           a[2] * d[3] * b[3], a[1] * d[1] * d[3] * b[3])
      names(indirect_serial) <- c(names(d), paste(m, collapse = "->"))
    }
    # add to existing vector for indirect effects
    indirect <- c(indirect, indirect_serial)
  }
  # return indirect effect(s)
  indirect
}

# internal function to extract direct effect(s)
extract_direct <- function(x, fit) {
  # initializations
  p_x <- length(x)
  # extract effect(s)
  if (p_x == 1L) unname(coef(fit)[x])
  else coef(fit)[x]
}

# internal function to extract total effect(s)
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
