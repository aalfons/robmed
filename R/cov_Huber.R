# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

#' Huber M-estimator of location and scatter
#'
#' Compute a Huber M-estimator of location and scatter, which is reasonably
#' robust for a small number of variables.
#'
#' An iterative reweighting algorithm is used to compute the Huber
#' M-estimator.  The Huber weight function thereby corresponds to a
#' convex optimization problem, resulting in a unique solution.
#'
#' @aliases print.cov_Huber
#'
#' @param x  a numeric matrix or data frame.
#' @param control  a list of tuning parameters as generated by
#' \code{\link{cov_control}()}.
#' @param \dots  additional arguments can be used to specify tuning parameters
#' directly instead of via \code{control}.
#'
#' @return An object of class \code{"cov_Huber"} with the following components:
#' \item{center}{a numeric vector containing the location vector estimate.}
#' \item{cov}{a numeric matrix containing the scatter matrix estimate.}
#' \item{prob}{numeric; probability for the quantile of the
#' \eqn{\chi^{2}}{chi-squared} distribution used as cutoff point in the Huber
#' weight function.}
#' \item{weights}{a numeric vector containing the relative robustness
#' weights for the observations.}
#' \item{tau}{numeric; correction for Fisher consistency under
#' multivariate normal distributions.}
#' \item{converged}{a logical indicating whether the iterative
#' reweighting algorithm converged.}
#' \item{iterations}{an integer giving the number of iterations required
#' to obtain the solution.}
#'
#' @author Andreas Alfons
#'
#' @references
#' Huber, P.J. (1981) \emph{Robust Statistics}. John Wiley & Sons.
#'
#' Zu, J. and Yuan, K.-H. (2010) Local Influence and Robust Procedures for
#' Mediation Analysis.  \emph{Multivariate Behavioral Research}, \bold{45}(1),
#' 1--44.  doi:10.1080/00273170903504695.
#'
#' @seealso \code{\link{cov_control}()}, \code{\link{test_mediation}()},
#' \code{\link{fit_mediation}()}
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
#' cov_Huber(BSG2014[, c(x, y, m)])
#'
#' @keywords multivariate
#'
#' @export

cov_Huber <- function(x, control = cov_control(...), ...) {
  ## initializations
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  ## check control arguments
  prob <- control$prob
  max_iterations <- control$max_iterations
  tol <- control$tol
  ## compute starting values for location vector and scatter matrix
  initial <- cov_ML(x)
  mu <- initial$center
  Sigma <- initial$cov
  ## perform iterative reweighting
  i <- 0
  if(prob < 1) {
    # define tuning parameters
    # tau is chosen such that E[chi_p^2 u^2(chi_p^2)] / tau = p
    # (i.e., such that Sigma is Fisher consistent)
    kappa <- 1 - prob
    r <- qchisq(prob, df=p)  # squared cutoff point
    tau <- (-2*dchisq(r, df=p) + kappa) * r / p + prob
    # perform iterative reweighting
    continue <- TRUE
    while(continue && i < max_iterations) {
      old_mu <- mu
      old_Sigma <- Sigma
      # compute weights based on mahalanobis distances
      d <- mahalanobis(x, center=mu, cov=Sigma)  # squared mahalanobis distances
      u <- sqrt(ifelse(d > r, r/d, 1))
      # update location vector and scatter matrix
      mu <- apply(x, 2, weighted.mean, w=u)
      Sigma <- crossprod(u * sweep(x, 2, mu, check.margin=FALSE)) / (n*tau)
      # check convergence
      i <- i + 1
      continue <- max(abs(mu-old_mu)) > tol || max(abs(Sigma-old_Sigma)) > tol
    }
    # check if algorithm converged
    if(i == max_iterations && continue) {
      warning(sprintf("no convergence in %d iterations", max_iterations))
    }
    converged <- !continue
  } else converged <- TRUE
  ## return results
  if(i == 0) {
    tau <- 1
    u <- rep.int(1, n)
  }
  result <- list(center=mu, cov=Sigma, prob=prob, weights=u, tau=tau,
                 converged=converged, iterations=i)
  class(result) <- "cov_Huber"
  result
}
