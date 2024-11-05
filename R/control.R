# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

#' Tuning parameters for robust regression
#'
#' Obtain a list with tuning paramters for the robust MM-estimator of
#' regression from \code{\link[robustbase]{lmrob}()} or median regression
#' from \code{\link[quantreg]{rq}()}.
#'
#' Prior to version 1.2.0, the MM-estimator of regression was the only type
#' of robust regression in \pkg{robmed} that supported control parameters.
#' Starting with version 1.2.0, control parameters can also be passed to median
#' regression, specifically the type of algorithm to be used. Function
#' \code{reg_control()} is an alias for \code{MM_reg_control()} for backwards
#' compatibility, but it is now recommended to use \code{MM_reg_control()} when
#' performing MM-regression and \code{median_reg_control()} when performing
#' median regression.
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
#' @return
#' For \code{MM_reg_control()} and \code{reg_control()}, a list of control
#' parameters for the MM-estimator of regression as returned by
#' \code{\link[robustbase]{lmrob.control}()}.
#'
#' For \code{median_reg_control()}, a list of control parameters for median
#' regression.
#'
#' @note \code{MM_reg_control()} and its alias \code{reg_control()} is a
#' simplified wrapper function for \code{\link[robustbase]{lmrob.control}()},
#' as the latter requires detailed knowledge of the algorithm for the
#' MM-estimator of regression.  Currently only 95\%, 90\%, 85\% (the default)
#' and 80\% efficiency are supported.  For other values, please specify the
#' corresponding tuning parameters in \code{\link[robustbase]{lmrob.control}()}
#' directly.
#'
#' @author Andreas Alfons
#'
#' @references
#' Salibian-Barrera, M. and Yohai, V.J. (2006) A Fast Algorithm for
#' S-regression Estimates. \emph{Journal of Computational and Graphical
#' Statistics}, \bold{15}(2), 414--427.  doi:10.1198/106186006x113629.
#'
#' Yohai, V.J. (1987) High Breakdown-Point and High Efficiency Estimates for
#' Regression. \emph{The Annals of Statistics}, \bold{15}(20), 642--656.
#' doi:10.1214/aos/1176350366.
#'
#' Koenker, R.W. (2005) \emph{Quantile Regression}. Camebridge University
#' Press.
#'
#' @seealso \code{\link[robustbase]{lmrob}()},
#' \code{\link[robustbase]{lmrob.control}()}
#'
#' @examples
#' data("BSG2014")
#'
#' # run fast-and-robust bootstrap test
#' ctrl <- MM_reg_control(efficiency = 0.95)
#' boot <- test_mediation(BSG2014,
#'                        x = "ValueDiversity",
#'                        y = "TeamCommitment",
#'                        m = "TaskConflict",
#'                        level = 0.9,
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


#' @rdname reg_control
#' @export

MM_reg_control <- reg_control



#' @rdname reg_control
#'
#' @param algorithm  a character string specifying the algorithm for computing
#' the median regression fit.  See argument \code{method} in
#' \code{\link[quantreg]{rq}()} for options. If you experience an infinite
#' loop with the default algorithm (\code{"br"}), you may want to try
#' \code{"fn"}.
#'
#' @export

median_reg_control <- function(algorithm = "br") {
  # note: additional arguments for model fit functions such as rq.fit.br()
  #       may be passed down via the ... argument, but this is currently not
  #       implemented as it's generally not recommended to change those
  #       arguments
  # # remove arguments that are not actual control arguments of the algorithm
  # dots <- list(...)
  # remove <- c("tau", "data", "subset", "weights",
  #             "na.action", "model", "contrasts")
  # dots <- dots[setdiff(names(dots), remove)]
  # # return list of control parameters
  # c(list(method = method), dots)
  list(method = algorithm)
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
#' Huber, P.J. (1981) \emph{Robust Statistics}. John Wiley & Sons.
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
