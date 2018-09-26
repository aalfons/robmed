# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

## load package and data
library("robmed", quietly = TRUE)
data("BSG2014")


## define variables
x <- "ValueDiversity"
y <- "TeamCommitment"
m <- c("TaskConflict", "ProcessConflict")

## seed of random number generator
seed <- 20150601


## only one mediator

# robust bootstrap test with regression
set.seed(seed)
robust_boot_reg_simple <- test_mediation(BSG2014, x = x, y = y, m = m[1],
                                         test = "boot", robust = TRUE,
                                         method = "regression")

# standard bootstrap test with regression
set.seed(seed)
standard_boot_reg_simple <- test_mediation(BSG2014, x = x, y = y, m = m[1],
                                           test = "boot", robust = FALSE,
                                           method = "regression")

# robust bootstrap test with covariance matrix
set.seed(seed)
robust_boot_cov_simple <- test_mediation(BSG2014, x = x, y = y, m = m[1],
                                         test = "boot", robust = TRUE,
                                         method = "covariance")

# standard bootstrap test with covariance matrix
set.seed(seed)
standard_boot_cov_simple <- test_mediation(BSG2014, x = x, y = y, m = m[1],
                                           test = "boot", robust = FALSE,
                                           method = "covariance")
# robust Sobel test with regression
set.seed(seed)
robust_sobel_reg_simple <- test_mediation(BSG2014, x = x, y = y, m = m[1],
                                          test = "sobel", robust = TRUE,
                                          method = "regression")

# standard Sobel test with regression
set.seed(seed)
standard_sobel_reg_simple <- test_mediation(BSG2014, x = x, y = y, m = m[1],
                                           test = "sobel", robust = FALSE,
                                           method = "regression")

# robust Sobel test with covariance matrix
set.seed(seed)
robust_sobel_cov_simple <- test_mediation(BSG2014, x = x, y = y, m = m[1],
                                         test = "sobel", robust = TRUE,
                                         method = "covariance")

# standard Sobel test with covariance matrix
set.seed(seed)
standard_sobel_cov_simple <- test_mediation(BSG2014, x = x, y = y, m = m[1],
                                           test = "sobel", robust = FALSE,
                                           method = "covariance")


## multiple mediators

# robust bootstrap test with regression
set.seed(seed)
robust_boot_reg_multi <- test_mediation(BSG2014, x = x, y = y, m = m,
                                        test = "boot", robust = TRUE,
                                        method = "regression")

# standard bootstrap test with regression
set.seed(seed)
standard_boot_reg_multi <- test_mediation(BSG2014, x = x, y = y, m = m,
                                          test = "boot", robust = FALSE,
                                          method = "regression")

## save results to package
save(seed, robust_boot_reg_simple, standard_boot_reg_simple,
     robust_boot_cov_simple, standard_boot_cov_simple,
     robust_sobel_reg_simple, standard_sobel_reg_simple,
     robust_sobel_cov_simple, standard_sobel_cov_simple,
     robust_boot_reg_multi, standard_boot_reg_multi,
     file = "tests/verified_results.RData")
