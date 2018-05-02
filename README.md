# NIMBLE

[![Build Status](https://travis-ci.org/nimble-dev/nimble.svg?branch=devel)](https://travis-ci.org/nimble-dev/nimble)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/nimble-dev/nimble?branch=devel&svg=true)](https://ci.appveyor.com/project/nimble-dev/nimble)
[![CRAN](http://www.r-pkg.org/badges/version/nimble)](https://cran.r-project.org/web/packages/nimble)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1211191.svg)](https://zenodo.org/record/1211191)
[![Google Group](https://img.shields.io/badge/google-group-blue.svg)](https://groups.google.com/forum/#!forum/nimble-users)

[Website](https://r-nimble.org/) |
[Documentation](https://r-nimble.org/manuals/NimbleUserManual.pdf) |
[Examples](https://r-nimble.org/examples) |
[Developing](https://nimble-dev.github.io/nimble-docs)

NIMBLE is an R package for hierarchical statistical modeling (aka
graphical modeling).  It enables writing general models along with
methods such as Markov chain Monte Carlo (MCMC), particle filtering
(aka sequential Monte Carlo), and other general methods.

For writing statistical models, NIMBLE adopts and extends the BUGS
language, making it largely compatible with
[BUGS](https://www.mrc-bsu.cam.ac.uk/software/bugs) and
[JAGS](http://mcmc-jags.sourceforge.net/).  NIMBLE makes BUGS
extensible, allowing users to add new functions and new distributions.

For writing algorithms (aka analysis methods), NIMBLE provides a
model-generic programming system embedded within R.  This provides
control over models as generic objects and mathematical manipulation
of model variables. In this way, NIMBLE's programming paradigm treats
probabilistic graphical models as a basic programming construct.

Both models and algorithms are compiled via generating customized C++
and providing seamless interfaces to compiled C++ from R.

NIMBLE's most developed methods are for MCMC.  Users can easily
customize sampler configurations from R and write new samplers in
NIMBLE's algorithm programming system.

Developers of new computational statistical methods can build them in
NIMBLE to gain the benefits of its graphical modeling language,
compilation, and distribution via [CRAN](https://cran.r-project.org/).

## Installation

### Install prerequisites

NIMBLE needs a C++ compiler and the GNU `make` utility.
Typically, Mac users can obtain these by installing Xcode, including
command line utilities, while Windows users can obtain them by
installing [Rtools](https://cran.r-project.org/bin/windows/Rtools/).
See the [User Manual](https://r-nimble.org/manuals/NimbleUserManual.pdf#page=26) for more details.

### Install NIMBLE

The easiest way to install NIMBLE is via CRAN:
```r
install.packages("nimble")
```

To install with extra functionality in `compareMCMCs` and `MCMCsuite`, install through the NIMBLE website:
```r
install.packages("nimble", type = "source", repos = "http://r-nimble.org")
```

## Citation

In published work that uses or mentions NIMBLE, please cite:

de Valpine, P., D. Turek, C.J. Paciorek, C. Anderson-Bergman,
D. Temple Lang, and R. Bodik. 2017. Programming with models: writing
statistical algorithms for general model structures with
NIMBLE. Journal of Computational and Graphical Statistics 26:403-413. [https://doi.org/10.1080/10618600.2016.1172487.](https://doi.org/10.1080/10618600.2016.1172487)

To cite a version of the package, please cite:

NIMBLE Development Team. 2018.
NIMBLE: An R Package for Programming with BUGS models, Version 0.6-10.
https://r-nimble.org,
https://zenodo.org/record/1174525.

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


