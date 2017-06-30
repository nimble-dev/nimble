#!/usr/bin/env Rscript

requirements <- c(
    'igraph',
    'coda',
    'testthat',
    'mvtnorm',  ## needed for test-distributions.R
    'abind',    ## needed for test-compareMCMCs.R
    'covr')     ## needed for code coverage reports

for (package in requirements) {
    if (!suppressPackageStartupMessages(require(package,
                                                character.only = TRUE))) {
        install.packages(package, repos = 'http://cran.us.r-project.org')
    }
}

## Tensorflow requires custom installation.
if(!require(tensorflow)) {
    install.packages('devtools')
    library(devtools)
    devtools::install_github('nimble-dev/r-tensorflow@upgrade-pip')
}
library(tensorflow)
tryCatch({
    ## Calling tf$ triggers loading of the python tensorflow library. 
    cat('Found Tensorflow version', tf$`__version__`, '\n')
}, error = function(e) {
    install_tensorflow()
})
