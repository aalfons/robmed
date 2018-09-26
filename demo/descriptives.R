# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

# load packages and data
library("dplyr")
library("tidyr")
data("BSG2014")

# compute descriptive statistics
descriptives <- BSG2014 %>%
  gather(key = "Variable", value = "Value") %>%
  group_by(Variable) %>%
  summarize(Mean = mean(Value),
            SD = sd(Value),
            Median = median(Value),
            MAD = mad(Value),
            Min = min(Value),
            Max = max(Value))

# round to 3 digits and print descriptive statistics
as.data.frame(descriptives %>% mutate_if(is.numeric, round, digits = 3))

# Spearman's rank correlation is more robust than the Pearson correlation,
# but it needs to be transformed to be consistent for the Pearson correlation
# (see Croux & Dehon 2010)
corSpearman <- function(...) 2 * sin(pi/6 * cor(..., method = "spearman"))

# compute correlation table and round to 3 digits
R <- corSpearman(BSG2014)
round(R, digits = 3)
