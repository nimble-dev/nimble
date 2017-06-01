#!/usr/bin/env Rscript
# This script runs most tests in nimble/inst/tests/ prioritized by duration.

# Avoid running these blacklisted tests, since they take too long.
blacklist <- c('test-Math2.R', 'test-Mcmc2.R', 'test-Mcmc3.R', 'test-Filtering2.R')
if(Sys.getenv('SKIP_EXPENSIVE_TESTS') == 'TRUE') {
    blacklist <- c(blacklist, 'test-mcmc.R', 'test-numericTypes.R')
}
cat('SKIPPING', blacklist, sep = '\n  ')

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

# Run under /usr/bin/time -v if possible, to gather timing information.
if (system2('/usr/bin/time', c('-v', 'echo', 'Running under /usr/bin/time -v'))) {
    runner <- 'Rscript'
} else {
    runner <- c('/usr/bin/time', '-v', 'Rscript')
}

# Run each test in a separate process to avoid dll garbage overload.
for (test in allTests) {
    cat('--------------------------------------------------------------------------------\n')
    cat('TESING', test, '\n')
    runViaTestthat <- FALSE  # TODO Fix test-size.R and set this to TRUE.
    if (runViaTestthat) {
        name <- gsub('test-(.*)\\.R', '\\1', test)
        script = paste0('library(methods); library(testthat); library(nimble); test_package("nimble", "', name, '")')
        command <- c(runner, '-e', shQuote(script))
    } else {
        command <- c(runner, file.path('packages', 'nimble', 'inst', 'tests', test))
    }
    if(system2(command[1], tail(command, -1))) {
        stop(paste('FAILED', test))
    }
    cat('PASSED', test, '\n')
}
