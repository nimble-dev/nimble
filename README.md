# NIMBLE

[![Build Status](https://travis-ci.org/nimble-dev/nimble.svg?branch=devel)](https://travis-ci.org/nimble-dev/nimble)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/nimble-dev/nimble?branch=devel&svg=true)](https://ci.appveyor.com/project/nimble-dev/nimble)
[![CRAN](http://www.r-pkg.org/badges/version/nimble)](https://cran.r-project.org/web/packages/nimble)
[![DOI](https://zenodo.org/badge/20771527.svg)](https://zenodo.org/badge/latestdoi/20771527)

[Website](http://r-nimble.org/) |
[Documentation](http://r-nimble.org/manuals/NimbleUserManual.pdf) |
[Examples](https://r-nimble.org/examples) |
[Developing](https://nimble-dev.github.io/nimble-docs)

NIMBLE is an R package for programming with BUGS models and compiling parts of R.
[BUGS](https://www.mrc-bsu.cam.ac.uk/software/bugs) is a probabilistic programming language that makes it easy to rapidly develop hierarchical Bayesian models for statistical analysis.
NIMBLE provides an implementation of the BUGS language in R together with a DSL for writing custom inference algorithms against BUGS models.
NIMBLE's programming paradigm treats probabilistic graphical models as a basic programming construct.

## Installation

NIMBLE can be installed in the usual fashion from CRAN, but Mac users need to have Xcode installed and Windows users need Rtools installed. Please see the User Manual for more details.

One can also install from the package version on our website; in general there is no need to do this, but some additional functionality in `compareMCMCs` and `MCMCsuite` is only available through the version on our website. To install from the website, invoke:
```
install.packages("nimble", type = "source", repos = "http://r-nimble.org")
```

Or use `R CMD INSTALL` after downloading from our [website](http://r-nimble.org/download-nimble).

## Citation

NIMBLE Development Team. 2017.
NIMBLE: An R Package for Programming with BUGS models, Version 0.6-5.
http://r-nimble.org.

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


