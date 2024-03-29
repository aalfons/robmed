% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/boot.R
\name{boot_samples}
\alias{boot_samples}
\title{Draw bootstrap samples}
\usage{
boot_samples(n, R)
}
\arguments{
\item{n}{an integer giving the number of observations in the original data
set.}

\item{R}{an integer giving the number of bootstrap samples to be generated.}
}
\value{
An object of class \code{"boot_samples"} with the following
components:
\item{indices}{an integer matrix in which each column contains the indices
of the corresponding bootstrap sample.}
\item{seed}{the state of the random number generator before the bootstrap
samples were drawn}
}
\description{
Draw bootstrap samples to be used for (fast-and-robust) bootstrap tests
for mediation analysis.  Note that this function is intended for use in
simulation studies by experienced users.
}
\examples{
\donttest{
# control parameters
n <- 100
a <- b <- c <- 0.4

# generate data
set.seed(20200309)
x <- rnorm(n)
m <- a * x + rnorm(n)
y <- b * m + c * x + rnorm(n)
simulated_data <- data.frame(x, y, m)

# perform boostrap tests
indices <- boot_samples(n, R = 5000)
robust_boot <- test_mediation(simulated_data,
                              x = "x", y = "y", m = "m",
                              robust = TRUE,
                              indices = indices)
summary(robust_boot)
ols_boot <- test_mediation(simulated_data,
                           x = "x", y = "y", m = "m",
                           robust = FALSE,
                           indices = indices)
summary(ols_boot)
}

}
\seealso{
\code{\link{test_mediation}()}
}
\author{
Andreas Alfons
}
\keyword{utilities}
