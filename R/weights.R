# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

#' Robustness weights of Huber M-estimation of location and scatter
#'
#' Extract (relative) robustness weights of a Huber M-estimate of location and
#' scatter.
#'
#' @method weights cov_Huber
#'
#' @param object  an object inheriting from class \code{"\link{cov_Huber}"}
#' containing Huber M-estimates of location and scatter.
#' @param type  a character string specifying the type of robustness weights to
#' be extracted.  Possible values are \code{"consistent"} and
#' \code{"relative"}.  The former can be used for a robust transformation of
#' the data such that the covariance matrix of the transformed data is Fisher
#' consistent.  Observations that are not downweighted in general receive a
#' weight larger than 1.  The latter are useful for interpretation, as
#' observations that are not downweighted receive a relative weight of 1.
#' @param \dots  additional arguments are currently ignored.
#'
#' @return A numeric vetor containing the requested robustness weights.
#'
#' @author Andreas Alfons
#'
#' @references
#' Zu, J. and Yuan, K.-H. (2010) Local Influence and Robust Procedures for
#' Mediation Analysis.  \emph{Multivariate Behavioral Research}, \bold{45}(1),
#' 1--44.  doi:10.1080/00273170903504695.
#'
#' @seealso \code{\link{cov_Huber}()}
#'
#' @examples
#' data("BSG2014")
#'
#' # define variables
#' x <- "ValueDiversity"
#' y <- "TeamCommitment"
#' m <- "TaskConflict"
#'
#' # compute Huber M-estimator
#' S <- cov_Huber(BSG2014[, c(x, y, m)])
#' weights(S, type = "relative")
#'
#' @keywords utilities
#'
#' @export

weights.cov_Huber <- function(object, type = c("consistent", "relative"), ...) {
  # initializations
  type <- match.arg(type)
  # extract weights
  weights <- object$weights
  if(type == "consistent") weights <- weights / sqrt(object$tau)
  weights
}
