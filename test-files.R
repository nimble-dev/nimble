#!/usr/bin/env Rscript
# This runs a file through the testthat package, and provides a little more
# detail than running a test script directly.

library(methods)
library(testthat)
suppressPackageStartupMessages(library(nimble))

nimbleOptions(verboseErrors = TRUE)
nimbleOptions(showCompilerOutput = TRUE)

filenames <- commandArgs(trailingOnly = TRUE)
if(length(filenames) == 0) {
    cat('Usage: test-files.R FILENAMES\n')
    quit('no', 1)
}

for(filename in filenames) {
    test_file(filename, reporter = c('tap', 'summary'))
}
