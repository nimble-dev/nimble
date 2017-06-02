NIMBLE
======

[![Build Status](https://travis-ci.org/nimble-dev/nimble.svg?branch=devel)](https://travis-ci.org/nimble-dev/nimble)
[![Coverage](https://codecov.io/github/nimble-dev/nimble/branch/devel/graphs/badge.svg)](https://codecov.io/github/nimble-dev/nimble) 
[![CRAN](http://www.r-pkg.org/badges/version/nimble)](https://cran.r-project.org/web/packages/nimble)

NIMBLE is an R package for programming with BUGS models and compiling parts of R.

*   [Website](http://r-nimble.org/)
*   [Examples](https://r-nimble.org/examples)
*   [User manual](http://r-nimble.org/manuals/NimbleUserManual.pdf)
*   [Developer documentation](https://nimble-dev.github.io/nimble-docs)

## Installation

NIMBLE can be installed in the usual fashion from CRAN, but Mac users need to have Xcode installed and Windows users need Rtools installed. Please see the User Manual for more details.

One can also install from the package version on our website; in general there is no need to do this, but some additional functionality in `compareMCMCs` and `MCMCsuite` is only available through the version on our website. To install from the website, invoke:
```
install.packages("nimble", type = "source", repos = "http://r-nimble.org")
```

Or use `R CMD INSTALL` after downloading from our [website](http://r-nimble.org/download-nimble).

## Citation

NIMBLE Development Team. 2017. NIMBLE: An R Package for Programming with BUGS models, Version 0.6-5.   http://r-nimble.org.

## Acknowledgements

The development of NIMBLE has been funded by:

* an NSF Advances in Biological Informatics grant (DBI-1147230) to P. de Valpine, C. Paciorek, and D. Temple Lang;
* an NSF SI2-SSI grant  (ACI-1550488) to P. de Valpine, C. Paciorek, and D. Temple Lang; and
* an NSF Collaborative Research grant (DMS-1622444) to P. de Valpine, A. Rodriguez, and C. Paciorek.

with additional support provided by postdoctoral funding for D. Turek from the Berkeley Institute for Data Science and Google Summer of Code fellowships for N. Michaud (2015) and C. Lewis-Beck (2017).


