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
  if (prefix) names <- paste("Indirect", names, sep = "_")
  combinations <- combn(names, 2, simplify = FALSE)
  n_contrasts <- length(combinations)
  # obtain labels for contrasts
  labels <- get_contrast_names(n_contrasts)
  if (prefix) labels <- paste("Indirect", labels, sep = "_")
  # obtain information on contrasts
  # FIXME: add prefix 'Indirect' only if 'prefix = TRUE'
  if (type == "estimates") {
    fun <- function(names) paste(names, collapse = " - ")
  } else if (type == "absolute") {
    fun <- function(names) {
      paste(paste0("|", names, "|"), collapse = " - ")
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


## utility function to get names of coefficients
get_effect_names <- function(a = NULL, b = NULL, d = NULL, total = NULL,
                             direct = NULL, indirect = NULL, sep = "_") {
  # names for effect(s) a
  if (is.null(a)) a_names <- NULL
  else a_names <- if (length(a) == 1L) "a" else paste("a", names(a), sep = sep)
  # names for effect(s) b
  if (is.null(b)) b_names <- NULL
  else b_names <- if (length(b) == 1L) "b" else paste("b", names(b), sep = sep)
  # names for effect(s) d
  if (is.null(d)) d_names <- NULL
  else d_names <- if (length(d) == 1L) "d" else paste("d", names(d), sep = sep)
  # names for total effect(s)
  if (is.null(total)) total_names <- NULL
  else if (length(total) == 1L) total_names <- "Total"
  else total_names <- paste("Total", names(total), sep = sep)
  # names for direct effect(s)
  if (is.null(direct)) direct_names <- NULL
  else if (length(direct) == 1L) direct_names <- "Direct"
  else direct_names <- paste("Direct", names(direct), sep = sep)
  # names for indirect effect(s)
  if (is.null(indirect)) indirect_names <- NULL
  else if (length(indirect) == 1L) indirect_names <- "Indirect"
  else indirect_names <- paste("Indirect", names(indirect), sep = sep)
  # return requested names
  c(a_names, b_names, d_names, total_names, direct_names, indirect_names)
}
# get_effect_names <- function(x, m, sep = "_") {
#   # initializations
#   p_x <- length(x)
#   p_m <- length(m)
#   # construct names
#   if (p_x == 1L) {
#     if (p_m == 1L) c("a", "b", "Direct", "Total")
#     else {
#       c(paste("a", m, sep = sep), paste("b", m, sep = sep), "Direct", "Total")
#     }
#   } else {
#     if (p_m == 1L) {
#       c(paste("a", x, sep = sep), "b", paste("Direct", x, sep = sep),
#         paste("Total", x, sep = sep))
#     } else {
#       c(paste("a", sapply(m, paste, x, sep = "."), sep = sep),
#         paste("b", m, sep = sep),
#         paste("Direct", x, sep = sep),
#         paste("Total", x, sep = sep))
#     }
#   }
# }


## The function for bootstrap replicates is required to return a vector.  This
## means that the columns of the bootstrap replicates contain coefficient
## estimates from different models.  This utility function returns the indices
## that correspond to the respective models, which makes it easier to extract
## the desired coefficients.
get_index_list <- function(p_x, p_m, p_covariates, indirect = TRUE,
                           model = "parallel") {
  # initializations
  include_indirect <- indirect
  nr_indirect <- get_nr_indirect(p_x, p_m, model = model)
  # numbers of coefficients in different models
  if (include_indirect) {
    p_indirect <- if (nr_indirect == 1L) 1L else 1L + nr_indirect
  } else p_indirect <- 0L
  if (!is.null(model) && model == "serial") {
    p_mx <- 1L + seq.int(0L, p_m-1L) + p_x + p_covariates
  } else p_mx <- rep.int(1L + p_x + p_covariates, p_m)
  p_ymx <- 1L + p_m + p_x + p_covariates
  p_total <- p_x
  p_all <- sum(p_indirect, p_mx, p_ymx, p_total)
  # indices of vector for each bootstrap replicate
  indices <- seq_len(p_all)
  # the first columns correspond to indirect effect(s) of x on y
  first <- 1L
  if (include_indirect) {
    indices_indirect <- seq.int(first, length.out = p_indirect)
  } else indices_indirect <- integer()
  # the next columns correspond to models m ~ x + covariates
  first <- first + p_indirect
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
  # the last column corresponds to the total effect of x on y
  first <- first + p_ymx
  index_total <- seq.int(from = first, length.out = p_total)
  # return list of indices
  list(indirect = indices_indirect, fit_mx = indices_mx,
       fit_ymx = indices_ymx, total = index_total)
}


## internal functions related to indirect effects

# obtain names for d path in serial multiple mediator models
get_d_names <- function(m, sep = "->") {
  # initializations
  p_m <- length(m)
  # currently only implemented for a single independent variable and
  # two or three hypothesized mediators
  if (p_m == 2L) {
    # two serial mediators
    names <- NULL
  } else {
    # three serial mediators
    responses <- m[-1L]
    predictors <- list(m[1L], m[-3L])
    names_list <- mapply(function(current_response, current_predictors) {
      paste(current_predictors, current_response, sep = sep)
    }, responses, predictors, SIMPLIFY = FALSE, USE.NAMES = FALSE)
    names <- unlist(names_list, use.names = FALSE)
  }
  # return names
  names
}

# obtain number of indirect effects
get_nr_indirect <- function(p_x, p_m, model = "parallel") {
  if (!is.null(model) && model == "serial") {
    # currently only implemented for a single independent variable and
    # two or three hypothesized mediators
    nr_indirect <- if (p_m == 2L) 3L else 7L
  } else {
    # single mediator or parallel multiple mediators
    nr_indirect <- p_x * p_m
  }
  # return number of indirect effects
  nr_indirect
}

# obtain names for indirect effects
get_indirect_names <- function(x, m, model = "parallel", sep = "->") {
  p_m <- length(m)
  if (!is.null(model) && model == "serial") {
    # currently only implemented for a single independent variable and
    # two or three hypothesized mediators
    if (p_m == 2L) {
      # two serial mediators
      c(m, paste(m, collapse = sep))
    } else {
      # three serial mediators
      c(m, paste(m[1L], m[2L], sep = sep), paste(m[1L], m[3L], sep = sep),
        paste(m[2L], m[3L], sep = sep), paste(m, collapse = sep))
    }
  } else {
    # single mediator or parallel multiple mediators
    p_x <- length(x)
    if (p_m == 1L) {
      names <- if (p_x > 1L) x
    } else {
      if (p_x == 1L) names <- m
      else {
        names <- unlist(lapply(m, function(current_m) {
          paste(x, current_m, sep = sep)
        }), use.names = FALSE)
      }
    }
    # return names
    names
  }
}

# obtain labels for indirect effects in printed output of bootstrap tests: only
# relevant in case of multiple independent variables and multiple mediators, or
# in case of a serial multiple mediator model
get_indirect_labels <- function(n) paste0("Indirect", seq_len(n))

# ## obtain labels for indirect effects in printed output of bootstrap tests
#
# get_indirect_labels <- function(object, ...) UseMethod("get_indirect_labels")
#
# get_indirect_labels.reg_fit_mediation <- function(object, ...) {
#   # initializations
#   x <- object$x
#   p_x <- length(x)
#   m <- object$m
#   p_m <- length(m)
#   model <- object$model
#   nr_indirect <- get_nr_indirect(p_x, p_m, model = model)
#   # obtain labels
#   if (model == "serial") labels <- paste0("Indirect", seq_len(nr_indirect))
#   else if (nr_indirect == 1L) labels <- m
#   else labels <- names(object$indirect)
#   # return labels
#   labels
# }
#
# get_indirect_labels.cov_fit_mediation <- function(object, ...) object$m


## check whether an object corresponds to a robust model fit

is_robust <- function(object) UseMethod("is_robust")

is_robust.reg_fit_mediation <- function(object) {
  robust <- object$robust
  if (is.character(robust)) TRUE else robust
}

is_robust.cov_fit_mediation <- function(object) object$robust
