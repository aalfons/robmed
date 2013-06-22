# ------------------------------------
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ------------------------------------

#' @S3method summary bootMA
summary.bootMA <- function(object, ...) {
  result <- list(object=object, 
                 summaryMX=summary(object$fitMX), 
                 summaryYX=summary(object$fitYX), 
                 summaryYMX=summary(object$fitYMX))
  class(result) <- "summaryMA"
  result
}

#' @S3method summary sobelMA
summary.sobelMA <- summary.bootMA
