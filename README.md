# __pumpingtest__: R Package for the analysis and evaluation of pumping tests
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/rves)](https://cran.r-project.org/package=pumpingtest)

The pumpingtest package provides functions to analyze and evaluate aquifer test including pumping tests with constant and variable rate, recovery tests, slug tests, and constant head tests (variable discharge).


## Installation 

The pumpingtest package is not available on CRAN and therefore it must be installed from github using:

```r
devtools::install_github("khaors/pumpingtest")
```

## Usage

The pumpingtest package does not require compilations and once installed it can be directly used by loading it: 

```r
library(pumpingtest)
```

## Package Documentation

The package documentation can be accesed [here](https://khaors.github.io/pumpingtest/) 

## Shiny App

This package includes a shiny app called _pumpingtest\_gui_ designed to help in the interpretation of the different types of aquifer tests. This app can be called using:

```r
pumpingtest_gui()
```

## Contact
 [![contributions welcome](https://img.shields.io/badge/contributions-welcome-brightgreen.svg?style=flat)](https://github.com/khaors/pumpingtest/issues)

Want to help or have questions? Contact me directly, use an issue, fork me or submit a pull request.

## License

GPL
