#!/usr/bin/env Rscript
# This script runs most tests in nimble/inst/tests/ prioritized by duration.
# To run only some of the tests, define the environment variable NIMBLE_TEST_BATCH=N
# where N is a number in 1,...,5.
# To see which tests will be run, pass the argument --dry-run, for example
#   ./run_tests.R --dry-run                        # Run all tests.                        
#   $ NIMBLE_TEST_BATCH=1 ./run_tests.R --dry-run  # Run one batch of tests.

# Avoid running these blacklisted tests, since they take too long.
blacklist <- c('test-Math2.R', 'test-Mcmc2.R', 'test-Mcmc3.R', 'test-Filtering2.R')
# Avoid running these tests since they test experimental features.
blacklist <- c(blacklist, 'test-ADfunctions.R', 'test-ADmodels.R', 'test-externalCalls.R')
cat('SKIPPING', blacklist, sep = '\n  ')

allTests <- list.files('packages/nimble/inst/tests')
allTests <- allTests[grepl('test-.*\\.R', allTests)]
allTests <- setdiff(allTests, blacklist)

# Sort tests by duration, running the shortest tests first.
testTimes <- read.csv('test_times.csv', sep = '\t', header = TRUE, row.names = 'filename')
for (test in allTests) {
    if (!(test %in% row.names(testTimes))) {
        # Tests without timing data are probably new, so run them first.
        testTimes[test, 'time'] <- 0.1  # Bogus very short duration.
    }
}
testTimes <- testTimes[order(testTimes),, drop = FALSE]
allTests <- intersect(row.names(testTimes), allTests)

# Parallelize tests by splitting them up into batches.
# We use the Best Fit Decreasing heuristic to approximatly solve this 1-D bin packing problem.
totalBatches = 5
testBatch = as.integer(Sys.getenv('NIMBLE_TEST_BATCH'))
if (!is.na(testBatch)) {
    if (testBatch < 1 || testBatch > totalBatches) {
        stop(paste0('Invalid NIMBLE_TEST_BATCH=', Sys.getenv('NIMBLE_TEST_BATCH')))
    }
    decreasingTests <- rev(allTests)
    binTimes <- -1e-6 * seq(totalBatches)  # Bias big tests towards later batches.
    allTests <- character()
    for (test in decreasingTests) {
        bin <- which.min(binTimes)
        binTimes[bin] <- binTimes[bin] + testTimes[test, 'time']
        if (bin == testBatch) {
            allTests = c(test, allTests)  # Test shorter tests first.
        }
    }
}

cat('PLANNING TO TEST', allTests, sep = '\n  ')
cat('PREDICTED DURATION =', sum(testTimes[allTests, 'time']), 'sec\n')
if ('--dry-run' %in% commandArgs(trailingOnly = TRUE)) quit()

# Run under /usr/bin/time -v if possible, to gather timing information.
if (system2('/usr/bin/time', c('-v', 'echo', 'Running tests under /usr/bin/time -v'))) {
    cat('Warning: Unable to run tests under /usr/bin/time -v\n')
    runner <- 'Rscript'
} else {
    runner <- c('/usr/bin/time', '-v', 'Rscript')
}

# Run each test in a separate process to avoid dll garbage overload.
for (test in allTests) {
    cat('--------------------------------------------------------------------------------\n')
    cat('TESTING', test, '\n')
    runViaTestthat <- TRUE
    if (runViaTestthat) {
        name <- gsub('test-(.*)\\.R', '\\1', test)
        script = paste0('library(methods); library(testthat); library(nimble); test_package("nimble", "^', name, '$")')
        command <- c(runner, '-e', shQuote(script))
    } else {
        command <- c(runner, file.path('packages', 'nimble', 'inst', 'tests', test))
    }
    if(system2(command[1], tail(command, -1))) {
        stop(paste('FAILED', test))
    }
    cat('PASSED', test, '\n')
}
