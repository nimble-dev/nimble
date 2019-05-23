#!/usr/bin/env Rscript

help_message <-
"Run tests in nimble/inst/tests/ prioritized by duration.
Usage:
  ./run_tests.R       [OPTIONS]   # Run the default set of tests.
  ./run_tests.R NAMES [OPTIONS]   # Run a custom set of tests, e.g. 'math'

Options:
  --parallel  Runs test in parallel.
  --stop      Stop after the first failure.
  --dry-run   Prints which tests would be run, then exits. This is useful to
              see what tests are available.

Environment Variables:
  NIMBLE_TEST_BATCH=N
    This allows distributed testing across machines. To run only one batch, set
    the NIMBLE_TEST_BATCH=N where N is a number in 1,...,5. For example

      NIMBLE_TEST_BATCH=1 ./run_tests.R --dry-run   # Run one batch of tests.
"

# Parse command line options.
argv <- commandArgs(trailingOnly = TRUE)
optionDryRun <- ('--dry-run' %in% argv)
optionParallel <- ('--parallel' %in% argv)
reporter <- if ('--stop' %in% argv) 'c("stop", "summary")' else '"summary"'
if ('-h' %in% argv || '--help' %in% argv) {
    cat(help_message)
    quit()
}

# Determine which tests to run.
if (length(grep('^-', argv, invert = TRUE))) {
    # Run only tests specified on commmand line.
    allTests <- paste0('test-', argv[!grepl('^-', argv)], '.R')
} else {
    # Run a default set of tests.
    allTests <- list.files('packages/nimble/inst/tests')
    allTests <- allTests[grepl('test-.*\\.R', allTests)]

    # Avoid running these blacklisted tests, since they take too long.
    blacklist <- c(
        'test-Math2.R',
        'test-Mcmc2.R',
        'test-Mcmc3.R',
        'test-Filtering2.R',
        'test-benchmark-building-steps.R')
    # Avoid running these tests since they test experimental features.
    blacklist <- c(
        blacklist,
        'test-ADfunctions.R',
        'test-ADmodels.R')
        ## 'test-benchmarks.R')  # some issue with version conflicts causing tensorflow to fail on Travis with errors such as 'nimble-tensorflow_11_20_18_17_45.so: undefined symbol: TF_DeleteImportGraphDefOptions'
    cat('SKIPPING', blacklist, sep = '\n  ')
    allTests <- setdiff(allTests, blacklist)
}

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
if (optionDryRun) quit()

# Run under /usr/bin/time -v if possible, to gather timing information.
runner <- 'Rscript'
if (optionParallel || system2('/usr/bin/time', c('-v', 'echo'), stderr = FALSE)) {
    cat('Not running tests under /usr/bin/time -v\n')
} else {
    cat('Running tests under /usr/bin/time -v\n')
    runner <- c('/usr/bin/time', '-v', 'Rscript')
}

# Run under exec_wait if sys package is installed, to support CTRL+C interrupts.
if (require(sys)) {
    custom_system2 <- sys::exec_wait
    custom_shQuote <- function(x) x  # exec_wait doesn't like shQuote.
} else {
    cat('Missing suggested package sys, falling back to system2\n')
    custom_system2 <- system2
    custom_shQuote <- shQuote
}

# Run each test in a separate process to avoid dll garbage overload.
runTest <- function(test, logToFile = FALSE, runViaTestthat = TRUE) {
    if (!logToFile) cat('--------------------------------------------------------------------------------\n')
    cat('TESTING', test, '\n')
    if (runViaTestthat) {
        name <- gsub('test-(.*)\\.R', '\\1', test)
        script <- paste0('library(methods);',
                         'library(testthat);',
                         'library(nimble);',
                         'tryCatch(test_package("nimble", "^', name, '$",',
                         '                      reporter = ', reporter, '),',
                         '  error = function(e) quit(status = 1))')
        command <- c(runner, '-e', custom_shQuote(script))
    } else {
        command <- c(runner, file.path('packages', 'nimble', 'inst', 'tests', test))
    }
    Sys.setenv(MAKEFLAGS = '-j1')  # Work around broken job pipe when GNU make is run under mclapply.
    if (logToFile) {
        logDir <- '/tmp/log/nimble'
        dir.create(logDir, recursive = TRUE, showWarnings = FALSE)
        stderr.log <- file.path(logDir, paste0('test-', name, '.stderr'))
        stdout.log <- file.path(logDir, paste0('test-', name, '.stdout'))
        if (custom_system2(command[1], tail(command, -1), stderr.log, stdout.log)) {
            cat('\x1b[31mFAILED\x1b[0m', test, 'See', stderr.log, stdout.log, '\n')
            return(TRUE)
        }
    } else {
        if (custom_system2(command[1], tail(command, -1))) {
            stop(paste('\x1b[31mFAILED\x1b[0m', test))
        }
    }
    cat('\x1b[32mPASSED\x1b[0m', test, '\n')
    return(FALSE)
}

if (optionParallel) {
    if (!require(parallel)) stop('Missing parallel package, required for --parallel')
    cores <- detectCores()
    cat('PARALLELIZING OVER', cores, 'CORES\n')
    failed <- mclapply(allTests, runTest, logToFile = TRUE,
                       mc.cores = cores, mc.preschedule = FALSE, mc.cleanup = TRUE)
    numFailed <- sum(unlist(failed))
    if (numFailed == 0) {
        cat('PASSED all tests\n')
    } else {
        stop(paste('FAILED', numFailed, 'tests'))
    }
} else {
    for (test in allTests) {
        runTest(test)
    }
}
