#!/usr/bin/env Rscript

requirements <- c(
    'igraph',
    'coda',
    'testthat',
    'R6',
    'mvtnorm',  ## needed for test-distributions.R
    'abind',    ## needed for test-compareMCMCs.R
    'ggplot2',  ## needed for test-compareMCMCs.R
    'covr',     ## needed for code coverage reports
    'pracma',   ## for AD
    'numDeriv',  ## for AD
    'mcmcse'    ## for MCEM
    # 'lme4'     ## for test-ADlaplace.R
    )     

for (package in requirements) {
    if (!suppressPackageStartupMessages(require(package,
                                                character.only = TRUE))) {
        install.packages(package, repos = 'https://cran.r-project.org')
    }
}

## Apparently a bug in Matrix (as of early 2024) is causing an issue (https://bioconductor.org/packages/devel/bioc/vignettes/dreamlet/inst/doc/errors.html) that is causing test-ADlaplace.R failures when fitting a model with lmer.

install.packages("lme4", type = "source")
