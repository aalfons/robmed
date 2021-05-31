# robmed: (Robust) Mediation Analysis

Perform mediation analysis via the robust bootstrap test ROBMED ([Alfons, Ates & Groenen, 2021](https://doi.org/10.1177/1094428121999096)).  In addition to ROBMED, several other (bootstrap) tests for mediation analysis are implemented.


## About ROBMED

The robust bootstrap test ROBMED for mediation analysis is less sensitive to deviations from model assumptions (such as outliers or heavily tailed distributions) than the standard bootstrap test of Preacher & Hayes ([2004](http://dx.doi.org/10.3758/BF03206553), [2008](http://dx.doi.org/10.3758/BRM.40.3.879)).  ROBMED utilizes the robust MM-regression estimator ([Yohai, 1987](https://projecteuclid.org/euclid.aos/1176350366)) instead of the OLS estimator for regression, and runs bootstrap tests with the fast and robust bootstrap methodology ([Salibián-Barrera & Zamar, 2002](https://projecteuclid.org/euclid.aos/1021379865); [Salibián-Barrera & Van Aelst, 2008](https://doi.org/10.1016/j.csda.2008.05.007)).

More information can be found in our article:

[Alfons, A., Ates, N.Y., & Groenen, P.J.F. (2021). A Robust Bootstrap Test for
Mediation Analysis. *Organizational Research Methods*. DOI 10.1177/1094428121999096.](https://doi.org/10.1177/1094428121999096)


## Installation

The `R` package `robmed` is on CRAN (The Comprehensive R Archive Network), hence it can be easily installed from the `R` command line via

```
install.packages("robmed")
```


## Bug reports and feature requests

If you experience any problems with `robmed` or if you have any suggestions for additional features, please submit an issue via the *Issues* tab of this repository.  Please have a look at existing issues first to see if your problem for feature request has already been discussed.
