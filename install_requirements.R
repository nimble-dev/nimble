#!/usr/bin/env Rscript

requirements <- c(
    'igraph',
    'coda',
    'testthat',
    'R6',
    'mvtnorm',  ## needed for test-distributions.R
    'abind',    ## needed for test-compareMCMCs.R
    'ggplot2',  ## needed for test-compareMCMCs.R
    'covr'      ## needed for code coverage reports
    )     

for (package in requirements) {
    if (!suppressPackageStartupMessages(require(package,
                                                character.only = TRUE))) {
        install.packages(package, repos = 'http://cran.us.r-project.org')
    }
}

