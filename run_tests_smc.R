
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

# Determine which tests to run.
if (length(grep('^-', argv, invert = TRUE))) {
    # Run only tests specified on commmand line.
    allTests <- paste0('test-', argv[!grepl('^-', argv)], '.R')
} else {
    # Run a default set of tests.
    allTests <- list.files('packages/nimbleSMC/tests/testthat')
    allTests <- allTests[grepl('test-.*\\.R', allTests)]
}

# Parallelize tests by splitting them up into batches.
# We use the Best Fit Decreasing heuristic to approximatly solve this 1-D bin packing problem.

cat('PLANNING TO TEST', allTests, sep = '\n  ')
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
runTest <- function(test, logToFile = FALSE, runViaTestthat = TRUE) {
    if (!logToFile) cat('--------------------------------------------------------------------------------\n')
    cat('TESTING', test, '\n')
    if (runViaTestthat) {
        name <- gsub('test-(.*)\\.R', '\\1', test)
        script <- paste0('library(methods);',
                         'library(testthat);',
                         'library(nimble);',
                         'library(nimbleSMC);',
                         'tryCatch(test_dir("packages/nimbleSMC/tests/testthat", "^', name, '$",',
                         '                      reporter = ', reporter, '),',
                         '  error = function(e) quit(status = 1))')
        command <- c(runner, '-e', custom_shQuote(script))
    } else {
        command <- c(runner, file.path('packages', 'nimbleSMC', 'testthat', 'tests', test))
    }
    Sys.setenv(MAKEFLAGS = '-j1')  # Work around broken job pipe when GNU make is run under mclapply.
    if (logToFile) {
        logDir <- 'log'
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


for (test in allTests) {
  runTest(test)
}
