# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

#' @export

p_value <- function(object, ...) UseMethod("p_value")


#' @rdname p_value
#' @method p_value boot_test_mediation
#' @export

p_value.boot_test_mediation <- function(object, digits = 4L, ...) {
  # number of hypothesized mediators
  p_m <- length(object$fit$m)
  # compute p-value
  if(p_m == 1L) {
    # only one mediator
    # set lower bound of significance level to 0
    lower <- 0
    # loop over the number of digits and determine the corresponding digit after
    # the comma of the p-value
    for (digit in seq_len(digits)) {
      # set step size
      step <- 1 / 10^digit
      # reset the significance level to the lower bound as we continue from there
      # with a smaller stepsize
      alpha <- lower
      # there is no rejection at the lower bound, so increase significance level
      # until there is rejection
      reject <- FALSE
      while(!reject) {
        # update lower bound and significance level
        lower <- alpha
        alpha <- alpha + step
        # retest at current significance level and extract confidence interval
        ci <- confint(object$reps, parm = 1L, level = 1 - alpha,
                      alternative = object$alternative, type = object$type)
        # reject if 0 is not in the confidence interval
        reject <- prod(ci) > 0
      }
    }
  } else {
    # multiple mediators
    rn <- rownames(object$ci)
    alpha <- sapply(seq_along(rn), function(j) {
      # set lower bound of significance level to 0
      lower <- 0
      # loop over the number of digits and determine the corresponding digit after
      # the comma of the p-value
      for (digit in seq_len(digits)) {
        # set step size
        step <- 1 / 10^digit
        # reset the significance level to the lower bound as we continue from there
        # with a smaller stepsize
        alpha_j <- lower
        # there is no rejection at the lower bound, so increase significance level
        # until there is rejection
        reject <- FALSE
        while(!reject) {
          # update lower bound and significance level
          lower <- alpha_j
          alpha_j <- alpha_j + step
          # retest at current significance level and extract confidence interval
          # ci <- retest(object, level = 1 - alpha)$ci[m, ]
          ci <- confint(object$reps, parm = j, level = 1 - alpha_j,
                        alternative = object$alternative, type = object$type)
          # reject if 0 is not in the confidence interval
          reject <- prod(ci) > 0
        }
      }
      # return p-value for current indirect effect
      alpha_j
    })
    names(alpha) <- rn
  }
  # return smallest significance level where 0 is not in the confidence interval
  alpha
}
