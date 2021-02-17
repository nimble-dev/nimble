#!/usr/bin/env Rscript

help_message <-
"Run tests in nimbleSMC/tests/testthat
Usage:
  ./run_tests.R       [OPTIONS]   # Run the default set of tests.
  ./run_tests.R NAMES [OPTIONS]   # Run a custom set of tests, e.g. 'filtering'
Options:
  --parallel  Runs test in parallel.
  --stop      Stop after the first failure.
  --dry-run   Prints which tests would be run, then exits. This is useful to
              see what tests are available.
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

# Parallelize tests by splitting them up into batches.
# We use the Best Fit Decreasing heuristic to approximatly solve this 1-D bin packing problem.

cat('PLANNING TO TEST ALL SMC TESTS', sep = '\n  ')
if (optionDryRun) quit()

# Run under /usr/bin/time -v if possible, to gather timing information.
runner <- 'Rscript'
if (system2('/usr/bin/time', c('-v', 'echo'), stderr = FALSE)) {
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
runTests <- function(logToFile = FALSE, runViaTestthat = TRUE) {
    if (!logToFile) cat('--------------------------------------------------------------------------------\n')
    cat('TESTING ALL TESTS', '\n')
    if (runViaTestthat) {
        script <- paste0('library(methods);',
                         'library(testthat);',
                         'library(nimble);',
                         'library(nimbleSMC);',
                         'tryCatch(test_package("nimbleSMC", ',
                         '                      reporter = ', reporter, '),',
                         '  error = function(e) quit(status = 1))')
        command <- c(runner, '-e', custom_shQuote(script))
    } else {
        stop("Not set up to run outside of testthat")
    }
    Sys.setenv(MAKEFLAGS = '-j1')  # Work around broken job pipe when GNU make is run under mclapply.
    if (logToFile) {
        logDir <- 'log'
        dir.create(logDir, recursive = TRUE, showWarnings = FALSE)
        stderr.log <- file.path(logDir, paste0('test-', name, '.stderr'))
        stdout.log <- file.path(logDir, paste0('test-', name, '.stdout'))
        if (custom_system2(command[1], tail(command, -1), stderr.log, stdout.log)) {
            cat('\x1b[31mFAILED\x1b[0m', ' See', stderr.log, stdout.log, '\n')
            return(TRUE)
        }
    } else {
        if (custom_system2(command[1], tail(command, -1))) {
            stop(paste('\x1b[31mFAILED\x1b[0m'))
        }
    }
    cat('\x1b[32mPASSED\x1b[0m', '\n')
    return(FALSE)
}

runTests()

print("here")
testthat::test_package('nimbleSMC', reporter='summary')
dir=system.file(file.path('tests', 'testthat'), package = 'nimbleSMC')
print(dir)
print(list.files(dir))
print('done')
