emclustr
========

[![Build Status](https://travis-ci.org/alexWhitworth/emclustr.svg?branch=master)](https://travis-ci.org/alexWhitworth/emclustr.svg?branch=master)
[![Coverage Status](https://coveralls.io/repos/github/alexWhitworth/emclustr/badge.svg?branch=master)](https://coveralls.io/github/alexWhitworth/emclustr?branch=master)

#### R Package for clustering via EM algorithm
**emclustr** is a package for finite mixture modeling via EM algorithm. Its main extensions to the **mclust** package are (1) mixture modeling non-Gaussian distributed data, and (2) mixture modeling for Gaussian distributed data with missing values

#### Repo Notes
- Oct 2015 - v 0.4 -- wrote and passed unit tests for all functions
- Sept 2015 - v 0.3 -- updated documentation, examples, and added test script to GitHub.
- Sept 2015 - v 0.2 -- re-written to use `library(mvtnorm)`, add log-likelihood calculations, and calculate BIC for all models.
- June 2014 - v 0.1 (alpha version)
