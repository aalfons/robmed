# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

# load package and data
library("robmed")
data("BSG2014")


# empirical case 1
x <- "SharedExperience"
y <- "TeamPerformance"
m <- "TMS"

# plot data
labs <- c("Shared experience", "Team performance", "Transactive memory systems")
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
