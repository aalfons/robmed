# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

# load packages and data
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
summary(standard_boot, type = "data")
summary(robust_boot)

# determine the smallest significance level alpha for which the
# (1 - alpha) * 100% confidence interval contains does not contain 0
p_value(standard_boot)
p_value(robust_boot)


# let's have a closer look at the effect of x on m

# extract information for diagnostic plot with tolerance ellipse
boot_list <- list(Standard = standard_boot, Robust = robust_boot)
ellipses <- setup_ellipse_plot(boot_list, horizontal = x, vertical = m)

# create plot
colors <- c("#F8766D", "#00BFC4")
ggplot() +
  geom_path(aes(x = x, y = y, color = Method), data = ellipses$ellipse) +
  geom_point(aes(x = x, y = y, fill = Weight), data = ellipses$data,
             shape = 21, size = 3) +
  geom_abline(aes(intercept = intercept, slope = slope, color = Method),
              data = ellipses$line, show.legend = FALSE) +
  scale_color_manual(values = colors) +
  scale_fill_gradient(limits = 0:1, low = "white", high = "black") +
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
summary(bc_boot, type = "data")
p_value(bc_boot)
