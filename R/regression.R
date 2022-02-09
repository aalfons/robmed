# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------


## model fitting functions that make summary() work
# This will allow fit_mediation() to fit the regression models without a
# formula, which is necessary to make the formula interface for mediation
# models work.  Otherwise formulas that contain log(.) or similar terms would
# cause errors when fitting the models in fit_mediation().

# OLS regression
lm_fit <- function(x, y, intercept = TRUE) {
  # if requested, add constant for intercept
  if (intercept) {
    n <- nrow(x)
    x <- cbind("(Intercept)" = rep.int(1, n), x)
  }
  # fit the linear model
  fit <- lm.fit(x, y)
  # Add a dummy formula as a 'terms' component to make summary() method work.
  # This 'terms' component needs to have an attribute that specifies whether
  # the model has an intercept, that's all.  It's not used in any other way.
  f <- as.formula(NULL)
  attr(f, "intercept") <- as.integer(intercept)
  fit$terms <- f
  # add the class and return the model fit
  class(fit) <- "lm"
  fit
}

# MM regression
lmrob_fit <- function(x, y, intercept = TRUE, control = reg_control()) {
  # if requested, add constant for intercept
  if (intercept) {
    n <- nrow(x)
    x <- cbind("(Intercept)" = rep.int(1, n), x)
  }
  # fit the linear model
  fit <- lmrob.fit(x, y, control = control)
  # Add a dummy formula as a 'terms' component to make summary() method work.
  # This 'terms' component needs to have an attribute that specifies whether
  # the model has an intercept, that's all.  It's not used in any other way.
  f <- as.formula(NULL)
  attr(f, "intercept") <- as.integer(intercept)
  fit$terms <- f
  # class is already added in lmrob.fit(), so just return the model fit
  fit
}

# quantile (median) regression
rq_fit <- function(x, y, intercept = TRUE, tau = 0.5) {
  # summary() method requires data frame and formula to extract response and
  # predictor matrix
  data <- data.frame(y, x, check.names = FALSE)
  formula <- if (intercept) y ~ . else y ~ 0 + .
  # if requested, add constant for intercept
  if (intercept) {
    n <- nrow(x)
    x <- cbind("(Intercept)" = rep.int(1, n), x)
  }
  # fit the linear model
  fit <- rq.fit(x, y, tau = tau, method = "br")
  fit$method <- "br"
  # construct terms object from formula that specifies whether there is a
  # response and an intercept
  terms <- formula
  attr(terms, "response") <- 1L
  attr(terms, "intercept") <- as.integer(intercept)
  # summary() method requires that the terms object is also an attribute of the
  # data frame
  attr(data, "terms") <- terms
  # add information to model fit
  fit$formula <- formula
  fit$terms <- terms
  fit$model <- data
  fit$tau <- tau
  # add the class and return the model fit
  class(fit) <- "rq"
  fit
}
