# ------------------------------------
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ------------------------------------

## internal function
summaryMA <- function(object, ...) {
  result <- list(object=object, 
                 summaryMX=summary(object$fitMX), 
                 summaryYX=summary(object$fitYX), 
                 summaryYMX=summary(object$fitYMX))
  class(result) <- "summaryMA"
  result
}


#' @S3method summary bootMA
summary.bootMA <- summaryMA


#' @S3method summary sobelMA
summary.sobelMA <- summaryMA
