# robmed: (Robust) Mediation Analysis

[![CRAN](https://www.R-pkg.org/badges/version/robmed)](https://CRAN.R-project.org/package=robmed) 


Perform mediation analysis via the robust bootstrap test ROBMED ([Alfons, Ates & Groenen, 2021](https://doi.org/10.1177/1094428121999096)).  In addition to ROBMED, several other (bootstrap) tests for mediation analysis are implemented.


## About ROBMED

The robust bootstrap test ROBMED for mediation analysis is less sensitive to deviations from model assumptions (such as outliers or heavily tailed distributions) than the standard bootstrap test of Preacher & Hayes ([2004](http://doi.org/10.3758/BF03206553), [2008](http://doi.org/10.3758/BRM.40.3.879)).  ROBMED utilizes the robust MM-regression estimator ([Yohai, 1987](https://doi.org/10.1214/aos/1176350366)) instead of the OLS estimator for regression, and runs bootstrap tests with the fast and robust bootstrap methodology ([Salibián-Barrera & Zamar, 2002](https://doi.org/10.1214/aos/1021379865); [Salibián-Barrera & Van Aelst, 2008](https://doi.org/10.1016/j.csda.2008.05.007)).

More information can be found in our article:

Alfons, A., Ates, N.Y., & Groenen, P.J.F. (2021). A Robust Bootstrap Test for Mediation Analysis. *Organizational Research Methods*. DOI [10.1177/1094428121999096](https://doi.org/10.1177/1094428121999096).


## Installation

Package `robmed` is on CRAN (The Comprehensive R Archive Network), hence the latest release can be easily installed from the `R` command line via

```
install.packages("robmed")
```


## Report issues and request features

If you experience any bugs or issues or if you have any suggestions for additional features, please submit an issue via the [*Issues*](https://github.com/aalfons/robmed/issues) tab of this repository.  Please have a look at existing issues first to see if your problem or feature request has already been discussed.


## Contribute to the package

If you want to contribute to the package, you can fork this repository and create a pull request after implementing the desired functionality.


## Ask for help

If you need help using the package, or if you are interested in collaborations related to this project, please get in touch with the [package maintainer](https://personal.eur.nl/alfons/).
