# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

# load package and data
library("MASS")
library("robmed")
data("BSG2014")


# empirical case 2
x <- "ValueDiversity"
y <- "TeamCommitment"
m <- "TaskConflict"

# plot data
labs <- c("Value diversity", "Team commitment", "Task conflict")
names(labs) <- c(x, y, m)
plot(BSG2014[, c(x, y, m)], labels = gsub(" ", "\n", labs),
     pch = 21, bg = "black")

# seed of random number generator
RNGversion("3.5.3")
seed <- 20150601

# perform standard method and proposed robust method
set.seed(seed)
standard_boot <- test_mediation(BSG2014, x = x, y = y, m = m, robust = FALSE)
set.seed(seed)
robust_boot <- test_mediation(BSG2014, x = x, y = y, m = m, robust = TRUE)

# The standard method of Preacher & Hayes (2004, 2008) uses normal theory
# t-tests for the significance of effects other than the indirect effect.
# Our method uses normal approximation bootstrap z-tests instead.
summary(standard_boot, other = "theory")
summary(robust_boot)

# determine the smallest significance level alpha for which the
# (1 - alpha) * 100% confidence interval contains does not contain 0
p_value(standard_boot)
p_value(robust_boot)


# let's have a closer look at the effect of x on m

# function for weighted covariance matrix
weighted_cov <- function(x, w, ...) {
  # sweep out columnwise weighted means
  center <- apply(x, 2, weighted.mean, w = w)
  x <- sweep(x, 2, center, check.margin = FALSE)
  # compute weighted crossproduct
  cn <- colnames(x)
  sapply(cn, function(j) sapply(cn, function(i) {
    sum(x[, i] * x[, j] * w)
  })) / (sum(w) - 1)
}

# function to compute an ellipse based on center and covariance matrix
ellipse <- function (center, cov, level = 0.975, n = 100) {
  # extract scales and correlation
  scale <- sqrt(diag(cov))
  r <- cov[1, 2] / prod(scale)
  # compute ellipse
  d <- acos(r)
  a <- seq(0, 2 * pi, length.out = n)
  q <- sqrt(qchisq(level, df = 2))
  x <- q * scale[1] * cos(a + d/2) + center[1]
  y <- q * scale[2] * cos(a - d/2) + center[2]
  xy <- cbind(x, y)
  # add names of variables and return ellipse
  colnames(xy) <- names(center)
  xy
}

# outliers from robust regression
w <- weights(robust_boot$fit$fit_mx, type = "robustness")
# means and covariance matrices
standard_center <- colMeans(BSG2014[, c(x, m)])
standard_cov <- cov(BSG2014[, c(x, m)])
robust_center <- sapply(BSG2014[, c(x, m)], weighted.mean, w = w)
robust_cov <- weighted_cov(BSG2014[, c(x, m)], w = w)
# compute tolerance ellipses
standard_ellipse <- ellipse(standard_center, standard_cov, level = 0.975)
robust_ellipse <- ellipse(robust_center, robust_cov, level = 0.975)
ellipses <- rbind(data.frame(Method = "Standard", standard_ellipse),
                     data.frame(Method = "Robust", robust_ellipse))
# extract coefficients of regression lines
standard_coef <- coef(standard_boot$fit$fit_mx)
robust_coef <- coef(robust_boot$fit$fit_mx)
# create plot
colors <- c("#F8766D", "#00BFC4")
BSG2014$Weight <- w
ggplot() +
  geom_path(aes_string(x = x, y = m, color = "Method"), data = ellipses) +
  scale_color_manual(values = colors) +
  geom_point(aes_string(x = x, y = m, fill = "Weight"), data = BSG2014,
             shape = 21, size = 3) +
  scale_fill_gradient(low = "white", high = "black") +
  geom_abline(intercept = standard_coef[1], slope = standard_coef[2],
              color = colors[1]) +
  geom_abline(intercept = robust_coef[1], slope = robust_coef[2],
              color = colors[2]) +
  labs(x = labs[x], y = labs[m]) + theme_bw()


# transform the variables and re-apply standard method

# function for Box-Cox transformation
box_cox <- function(data, variable) {
  f <- as.formula(paste(variable, "1", sep = "~"))
  bc <- boxcox(f, data = data, plotit = FALSE)
  power <- bc$x[which.max(bc$y)]
  (data[, variable]^power - 1) / power
}

# Box-Cox transformation of independent variable
bc_x <- paste("bc", x, sep = "_")
BSG2014[, bc_x] <- box_cox(BSG2014, x)
# Box-Cox transformation of dependent variable
bc_y <- paste("bc", y, sep = "_")
BSG2014[, bc_y] <- box_cox(BSG2014, y)
# Box-Cox transformation of hypothesized mediator
bc_m <- paste("bc", m, sep = "_")
BSG2014[, bc_m] <- box_cox(BSG2014, m)

# scatterplot matrix of transformed variables
plot(BSG2014[, c(bc_x, bc_y, bc_m)], labels = gsub(" ", "\n", labs),
     pch = 21, bg = "black")

# perform standard method with transformed variables
set.seed(seed)
bc_boot <- test_mediation(BSG2014, x = bc_x, y = bc_y, m = bc_m, robust = FALSE)
summary(bc_boot, other = "theory")
p_value(bc_boot)
