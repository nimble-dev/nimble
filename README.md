# NIMBLE

[![Build Status](https://travis-ci.org/nimble-dev/nimble.svg?branch=devel)](https://travis-ci.org/nimble-dev/nimble)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/nimble-dev/nimble?branch=devel&svg=true)](https://ci.appveyor.com/project/nimble-dev/nimble)
[![CRAN](http://www.r-pkg.org/badges/version/nimble)](https://cran.r-project.org/web/packages/nimble)
[![DOI](https://zenodo.org/badge/20771527.svg)](https://zenodo.org/badge/latestdoi/20771527)

[Website](https://r-nimble.org/) |
[Documentation](https://r-nimble.org/manuals/NimbleUserManual.pdf) |
[Examples](https://r-nimble.org/examples) |
[Developing](https://nimble-dev.github.io/nimble-docs)

NIMBLE is an R package for programming with BUGS models and compiling parts of R.
[BUGS](https://www.mrc-bsu.cam.ac.uk/software/bugs) is a probabilistic programming language that makes it easy to rapidly develop hierarchical Bayesian models for statistical analysis.
NIMBLE provides an implementation of the BUGS language in R together with a DSL for writing custom inference algorithms against BUGS models.
NIMBLE's programming paradigm treats probabilistic graphical models as a basic programming construct.

## Installation

### Install prerequisites

NIMBLE needs a C++ compiler and the GNU `make` utility.
Additionally, Mac users need to have Xcode installed and Windows users need Rtools installed.
See the [User Manual](https://r-nimble.org/manuals/NimbleUserManual.pdf#page=26) for more details.

### Install NIMBLE

The easiest way to install NIMBLE is via CRAN:
```r
install.packages("nimble")
```

To install with extra functionality in `compareMCMCs` and `MCMCsuite`, install through the NIBMLE website:
```r
install.packages("nimble", type = "source", repos = "http://r-nimble.org")
```

## Citation

NIMBLE Development Team. 2017.
NIMBLE: An R Package for Programming with BUGS models, Version 0.6-6.
https://r-nimble.org,
https://zenodo.org/record/820704.

## Licenses

Nimble is released under a mixture of licenses,
and depends on additional third-party libraries with compatible licenses.

- Nimble's non-C++ code (R, bash, Make, etc.) is released under
  [Revised BSD](LICENSE).
- Nimble's C++ code is released under
  [GPL 2](https://www.gnu.org/licenses/gpl-2.0.html).
- Nimble's [User Manual](UserManual) is released under the
  [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) license.
- The [Eigen C++ library](http://eigen.tuxfamily.org) included with Nimble is
  licensed under [MPL 2](https://www.mozilla.org/en-US/MPL/2.0/).

## Acknowledgements

The development of NIMBLE has been funded by:

* an NSF Advances in Biological Informatics grant (DBI-1147230) to P. de Valpine, C. Paciorek, and D. Temple Lang;
* an NSF SI2-SSI grant  (ACI-1550488) to P. de Valpine, C. Paciorek, and D. Temple Lang; and
* an NSF Collaborative Research grant (DMS-1622444) to P. de Valpine, A. Rodriguez, and C. Paciorek.

with additional support provided by postdoctoral funding for D. Turek from the Berkeley Institute for Data Science and Google Summer of Code fellowships for N. Michaud (2015) and C. Lewis-Beck (2017).


