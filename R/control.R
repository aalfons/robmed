# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

#' Tuning parameters for MM-regression
#'
#' Obtain a list with tuning paramters for \code{\link[robustbase]{lmrob}()}.
#'
#' @param efficiency  a numeric value giving the desired efficiency (defaults
#' to 0.85 for 85\% efficiency).
#' @param max_iterations  an integer giving the maximum number of iterations in
#' various parts of the algorithm.
#' @param tol  a small positive numeric value to be used to determine
#' convergence in various parts of the algorithm.
#' @param seed  optional initial seed for the random number generator (see
#' \code{\link{.Random.seed}}).
#'
#' @return A list of tuning parameters as returned by
#' \code{\link[robustbase]{lmrob.control}()}.
#'
#' @note  This is a simplified wrapper function for
#' \code{\link[robustbase]{lmrob.control}()}, as the latter requires detailed
#' knowledge of the MM-type regression algorithm.  Currently only 95\%, 90\%,
#' 85\% (the default) and 80\% efficiency are supported.  For other values,
#' please specify the corresponding tuning parameters in
#' \code{\link[robustbase]{lmrob.control}()} directly.
#'
#' @author Andreas Alfons
#'
#' @references
#' Salibian-Barrera, M. and Yohai, V.J. (1987) A fast algorithm for
#' S-regression estimates. \emph{Journal of Computational and Graphical
#' Statistics}, \bold{15}(2), 414--427.
#'
#' Yohai, V.J. (1987) High breakdown-point and high efficiency estimates for
#' regression. \emph{The Annals of Statistics}, \bold{15}(20), 642--656.
#'
#' @seealso \code{\link[robustbase]{lmrob}()},
#' \code{\link[robustbase]{lmrob.control}()}
#'
#' @examples
#' data("BSG2014")
#'
#' # run fast-and-robust bootstrap test
#' ctrl <- reg_control(efficiency = 0.95)
#' boot <- test_mediation(BSG2014,
#'                        x = "ValueDiversity",
#'                        y = "TeamCommitment",
#'                        m = "TaskConflict",
#'                        control = ctrl)
#' summary(boot)
#'
#' @keywords regression
#'
#' @export

reg_control <- function(efficiency = 0.85, max_iterations = 200,
                        tol = 1e-07, seed = NULL) {
  # check supplied values
  max_iterations <- rep(as.integer(max_iterations), length.out=1)
  if(!is.finite(max_iterations)) max_iterations <- formals()$max_iterations
  else if(max_iterations < 0) max_iterations <- 0
  tol <- rep(as.numeric(tol), length.out=1)
  if(!is.finite(tol)) tol <- formals()$tol
  else if(tol < 0) tol <- 0
  # efficiency allows only 4 reasonable values, otherwise the robust F-test
  # becomes too much of a mess (to find constant to be used for weights)
  efficiency <- rep(as.numeric(efficiency), length.out=1)
  if(!is.finite(efficiency)) efficiency <- formals()$efficiency
  else {
    choices <- c(0.95, 0.9, 0.85, 0.8)
    differences <- abs(efficiency - choices)
    whichMin <- which.min(differences)
    efficiency <- choices[whichMin]
    if(differences[whichMin] > tol) {
      warning("using ", format(100 * efficiency), "% efficiency")
    }
  }
  # assign tuning parameter for given efficiency
  if(efficiency == 0.8)       tuning <- 3.136909
  else if(efficiency == 0.85) tuning <- 3.443689  # default
  else if(efficiency == 0.9)  tuning <- 3.882646
  else if(efficiency == 0.95) tuning <- 4.685061
  # call lmrob.control() to return list of control parameters
  tmp <- list(efficiency = efficiency)
  out <- lmrob.control(seed=seed, tuning.psi=tuning, max.it=max_iterations,
                       k.max=max_iterations, maxit.scale=max_iterations,
                       refine.tol=tol, rel.tol=tol, solve.tol=tol)
  c(tmp, out)
}


#' Tuning parameters for Huber M-estimation of location and scatter
#'
#' Obtain a list with tuning paramters for \code{\link{cov_Huber}()}.
#'
#' @param prob  numeric; probability for the quantile of the
#' \eqn{\chi^{2}}{chi-squared} distribution to be used as cutoff point in the
#' Huber weight function (defaults to 0.95).
#' @param max_iterations  an integer giving the maximum number of iterations in
#' the iteratively reweighted algorithm.
#' @param tol  a small positive numeric value to be used to determine
#' convergence of the iteratively reweighted algorithm.
#'
#' @return A list with components corresponding to the arguments.
#'
#' @author Andreas Alfons
#'
#' @references
#' Huber, P.J. (1981) \emph{Robust statistics}. John Wiley & Sons.
#'
#' @seealso \code{\link{cov_Huber}()}
#'
#' @examples
#' data("BSG2014")
#'
#' # run bootstrap test after winsorization
#' ctrl <- cov_control(prob = 0.95)
#' boot <- test_mediation(BSG2014,
#'                        x = "ValueDiversity",
#'                        y = "TeamCommitment",
#'                        m = "TaskConflict",
#'                        method = "covariance",
#'                        control = ctrl)
#' summary(boot)
#'
#' @keywords multivariate
#'
#' @export

cov_control <- function(prob = 0.95, max_iterations = 200, tol = 1e-07) {
  # check supplied values
  prob <- rep(as.numeric(prob), length.out=1)
  if(!is.finite(prob)) prob <- formals$prob
  else if(prob < 0) prob <- 0
  else if(prob > 1) prob <- 1
  max_iterations <- rep(as.integer(max_iterations), length.out=1)
  if(!is.finite(max_iterations)) max_iterations <- formals()$max_iterations
  else if(max_iterations < 0) max_iterations <- 0
  tol <- rep(as.numeric(tol), length.out=1)
  if(!is.finite(tol)) tol <- formals()$tol
  else if(tol < 0) tol <- 0
  # return list of control parameters
  list(prob=prob, max_iterations=max_iterations, tol=tol)
}
