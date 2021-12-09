# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------


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
  else stop(sprintf("%s contrasts not implemented", type))
  # compute contrasts
  sapply(combinations, fun, x = x)
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
get_contrast_info <- function(names, type = "estimates", prefix = FALSE) {
  # compute combinations of names
  combinations <- combn(names, 2, simplify = FALSE)
  n_contrasts <- length(combinations)
  # obtain labels for contrasts
  labels <- get_contrast_names(n_contrasts)
  if (prefix) labels <- paste("Indirect", labels, sep = "_")
  # obtain information on contrasts
  if (type == "estimates") {
    fun <- function(names) {
      paste(paste("Indirect", names, sep = "_"), collapse = " - ")
    }
  } else if (type == "absolute") {
    fun <- function(names) {
      paste(paste0("|Indirect_", names, "|"), collapse = " - ")
    }
  } else stop(sprintf("%s contrasts not implemented", type))
  # return information on contrasts
  data.frame(Label = labels, Definition = sapply(combinations, fun))
}
# get_contrast_info <- function(x, m, type = "estimates", prefix = FALSE) {
#   # initializations
#   p_x <- length(x)
#   p_m <- length(m)
#   # names used for indirect effects
#   if (p_x > 1 && p_m > 1) names <- sapply(m, paste, x, sep = ".")
#   else if (p_x > 1) names <- x
#   else if (p_m > 1) names <- m
#   else {
#     # should not happen
#     stop("contrasts are only applicable in case of multiple indirect effects")
#   }
#   # compute combinations of names
#   combinations <- combn(names, 2, simplify = FALSE)
#   n_contrasts <- length(combinations)
#   # obtain labels for contrasts
#   labels <- get_contrast_names(n_contrasts)
#   if (prefix) labels <- paste("ab", labels, sep = "_")
#   # obtain information on contrasts
#   if (type == "estimates") {
#     fun <- function(names) paste(paste("ab", names, sep = "_"), collapse = " - ")
#   } else if (type == "absolute") {
#     fun <- function(names) paste(paste0("|ab_", names, "|"), collapse = " - ")
#   } else stop(sprintf("%s contrasts not implemented", type))
#   # return information on contrasts
#   data.frame(Label = labels, Definition = sapply(combinations, fun))
# }


## internal functions related to indirect effects

# obtain number of indirect effects
get_nr_indirect <- function(x, m, model = "parallel") {
  p_m <- length(m)
  if (model == "serial") {
    # currently only implemented for a single independent variable and
    # two or three hypothesized mediators
    nr_indirect <- if (p_m == 2L) 3L else 7L
  } else {
    # single mediator or parallel multiple mediators
    p_x <- length(x)
    nr_indirect <- p_x * p_m
  }
}

# obtain names for indirect effects in serial multiple mediator model
get_indirect_names <- function(x, m, model = "parallel") {
  p_m <- length(m)
  if (model == "serial") {
    # currently only implemented for a single independent variable and
    # two or three hypothesized mediators
    if (p_m == 2L) {
      # two serial mediators
      c(m, paste(m, collapse = "."))
    } else {
      # three serial mediators
      c(m, paste(m[1L], m[2L], sep = "."), paste(m[1L], m[3L], sep = "."),
        paste(m[2L], m[3L], sep = "."), paste(m, collapse = "."))
    }
  } else {
    # single mediator or parallel multiple mediators
    p_x <- length(x)
    if (p_m == 1L) names <- if (p_x > 1L) x
    else {
      if (p_x == 1L) names <- m
      else names <- unlist(lapply(m, paste, x, sep = "."), use.names = FALSE)
    }
    # return names
    names
  }
}


## check whether an object corresponds to a robust model fit

is_robust <- function(object) UseMethod("is_robust")

is_robust.reg_fit_mediation <- function(object) {
  robust <- object$robust
  if (is.character(robust)) TRUE else robust
}

is_robust.cov_fit_mediation <- function(object) object$robust
