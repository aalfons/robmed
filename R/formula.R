# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

#' Create an object of hypothesized mediators or control variables
#'
#' \code{m()} creates an object of hypothesized mediators, while
#' \code{covariates()} creates an object of control variables.  Usually,
#' these are used in a formula specifying a mediation model.
#'
#' These are essentially wrappers for \code{\link[base]{cbind}()} with a
#' specific class prepended to the class(es) of the resulting object.
#'
#' @param \dots  variables are supplied as arguments, as usual separated by a
#' comma.
#'
#' @return \code{m()} returns an object of class \code{"parallel_mediators"}
#' and \code{covariates()} returns an object of class \code{"covariates"}.
#' Typically, these inherit from class \code{"matrix"}.
#'
#' @note The formula interface is still experimental and may change in future
#' versions.
#'
#' @author Andreas Alfons
#'
#' @seealso \code{\link{fit_mediation}()}, \code{\link{test_mediation}()}
#'
#' @examples
#' data("BSG2014")
#'
#' # inside formula
#' fit_mediation(TeamCommitment ~ m(TaskConflict) + ValueDiversity,
#'               data = BSG2014)
#'
#' # outside formula
#' mediator <- with(BSG2014, m(TaskConflict))
#' fit_mediation(TeamCommitment ~ mediator + ValueDiversity,
#'               data = BSG2014)
#'
#' @keywords utilities
#' @export

m <- function(...) {
  # This is a bit of a hack: returning a data frame would cause an error
  # with model.frame() in fit_mediation().  Instead, use cbind() to create
  # a matrix and store information on variable types as arguments.
  out <- cbind(...)
  # retrieve information on variable types
  df <- data.frame(..., check.names = FALSE, stringsAsFactors = FALSE)
  classes <- lapply(df, class)
  levels <- lapply(df, function(x) {
    if (is.factor(x)) levels(x)
  })
  # add class
  class(out) <- c("parallel_mediators", class(out))
  # add information on variable types as attributes
  attr(out, "column_types") <- classes
  attr(out, "factor_levels") <- levels
  # return object
  out
}


#' @rdname m
#' @export

covariates <- function(...) {
  # This is a bit of a hack: returning a data frame would cause an error
  # with model.frame() in fit_mediation().  Instead, use cbind() to create
  # a matrix and store information on variable types as arguments.
  out <- cbind(...)
  # retrieve information on variable types
  df <- data.frame(..., check.names = FALSE, stringsAsFactors = FALSE)
  classes <- lapply(df, class)
  levels <- lapply(df, function(x) {
    if (is.factor(x)) levels(x)
  })
  # add class
  class(out) <- c("covariates", class(out))
  # add information on variable types as attributes
  attr(out, "column_types") <- classes
  attr(out, "factor_levels") <- levels
  # return object
  out
}


# workhorse function to convert these objects to a data frame
convert_to_df <- function(x, row.names = NULL, optional = FALSE, ...) {
  # initializations
  p <- ncol(x)
  classes <- attr(x, "column_types")
  levels <- attr(x, "factor_levels")
  # convert variables
  variable_list <- mapply(function(j, class, levels) {
    # extract variable
    xj <- x[, j, drop = TRUE]
    # convert to factors if necessary
    if ("factor" %in% class) {
      if (is.character(xj)) xj <- factor(xj, levels = levels)
      else xj <- factor(xj, labels = levels)
    }
    # make sure to assign the exact class(es): this will also make sure that
    # ordered factors are treated correctly
    class(xj) <- class
    # return variable
    xj
  }, j = seq_len(p), class = classes, levels = levels,
  SIMPLIFY = FALSE, USE.NAMES = FALSE)
  # make sure to keep the variable names
  names(variable_list) <- colnames(x)
  # add additional arguments
  variable_list$row.names <- row.names
  variable_list$check.names <- !optional
  variable_list$stringsAsFactors <- FALSE
  # construct data frame
  do.call(data.frame, variable_list)
}

# convert "parallel_mediators" object to data frame
as.data.frame.parallel_mediators <- convert_to_df

# convert "covariates" object to data frame
as.data.frame.covariates <- convert_to_df
