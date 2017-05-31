#!/usr/bin/env Rscript
# This script runs most tests in nimble/inst/tests/ prioritized by duration.

# Avoid running these blacklisted tests, since they take too long.
blacklist <- c('test-Math2.R', 'test-Mcmc2.R', 'test-Mcmc3.R', 'test-Filtering2.R')

allTests <- list.files('packages/nimble/inst/tests')
allTests <- allTests[grepl('test-.*\\.R', allTests)]

# Sort tests by duration, running the shortest tests first.
testTimes <- read.csv('test_times.csv', sep = '\t', header = TRUE, row.names = 'filename')
for (test in allTests) {
    if (!(test %in% row.names(testTimes)) && !(test %in% blacklist)) {
        # Tests without timing data are probably new, so run them first.
        testTimes[test, 'time'] <- 0
    }
}
testTimes <- testTimes[order(testTimes),, drop = FALSE]
allTests <- row.names(testTimes)

# Run each test in a separate process to avoid dll garbage overload.
for (test in allTests) {
    runViaTestthat <- TRUE  # TODO Fix test-size.R and set this to TRUE.
    if (runViaTestthat) {
        name <- gsub('test-(.*)\\.R', '\\1', test)
        script = paste0('library(methods); library(testthat); library(nimble); test_package("nimble", "', name, '")')
        command <- c('Rscript', '-e', shQuote(script))
    } else {
        command <- c('Rscript', file.path('packages', 'nimble', 'inst', 'tests', test))
    }
    if(system2('/usr/bin/time', c('-v', command))) {
        stop(paste('Failed', test))
    }
}
