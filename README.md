# robmed: (Robust) Mediation Analysis

[![CRAN](https://www.R-pkg.org/badges/version/robmed)](https://CRAN.R-project.org/package=robmed) 


Perform mediation analysis via the robust bootstrap test ROBMED ([Alfons, Ates & Groenen, 2022a](https://doi.org/10.1177/1094428121999096)).  In addition to ROBMED, various other methods for mediation analysis are implemented.

To cite the robust bootstrap test ROBMED, please use:

Alfons, A., Ates, N.Y., & Groenen, P.J.F. (2022a). A Robust Bootstrap Test for Mediation Analysis. *Organizational Research Methods*, **25**(3), 591--617. doi: [10.1177/1094428121999096](https://doi.org/10.1177/1094428121999096).

To cite package `robmed` in publications, please use:

Alfons, A., Ates, N.Y., & Groenen, P.J.F. (2022b). Robust Mediation Analysis: The `R` Package `robmed`. *Journal of Statistical Software*, **103**(13), 1--45. doi: [10.18637/jss.v103.i13](https://doi.org/10.18637/jss.v103.i13).

Further discussion on robust mediation analysis can be found in:

A. Alfons and D.R. Schley (2025). *Robust Mediation Analysis: What We Talk About When We Talk About Robustness*. PsyArXiv. doi: [10.31234/osf.io/2hqdy](https://doi.org/10.31234/osf.io/2hqdy).


## About ROBMED

The robust bootstrap test ROBMED for mediation analysis is less sensitive to deviations from model assumptions (such as outliers or heavily tailed distributions) than the standard bootstrap test of Preacher & Hayes ([2004](https://doi.org/10.3758/BF03206553), [2008](https://doi.org/10.3758/BRM.40.3.879)).  ROBMED utilizes the robust MM-regression estimator ([Yohai, 1987](https://doi.org/10.1214/aos/1176350366)) instead of the OLS estimator for regression, and runs bootstrap tests with the fast and robust bootstrap methodology ([Salibián-Barrera & Zamar, 2002](https://doi.org/10.1214/aos/1021379865); [Salibián-Barrera & Van Aelst, 2008](https://doi.org/10.1016/j.csda.2008.05.007)).

More information can be found in our article:

Alfons, A., Ates, N.Y., & Groenen, P.J.F. (2022a). A Robust Bootstrap Test for Mediation Analysis. *Organizational Research Methods*, **25**(3), 591--617. doi: [10.1177/1094428121999096](https://doi.org/10.1177/1094428121999096).

Note that the defaults in `robmed` have since changed to use simple percentile bootstrap confidence intervals for the indirect effects (starting with version 1.1.0).


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
