#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = TRUE)
if(length(args) != 1 || '-h' %in% args || '--help' %in% args) {
    cat(
        'Usage: Rscript sort-tests.R TEST_LOG_FILE',
        'Sort the test-devel.sh script by running the cheapest tests first.',
        'The TEST_LOG_FILE should be either a file that you manually create',
        'using someting like',
        '  test-devel.sh 2>&1 > test_log.txt',
        'or a travis log file that you download from',
        '  https://travis-ci.org/nimble-dev/nimble/builds',
        sep = '\n'
    )
    quit('no', 1)
}

# Find list of test files sorted by name.
lines = readLines(args[1])
col_usr <- lines[grepl('\tUser time .seconds.:', lines)]
col_sys <- lines[grepl('\tSystem time .seconds.:', lines)]
col_time <- as.numeric(gsub('.*: ', '', col_usr)) + as.numeric(gsub('.*: ', '', col_sys))
col_name <- lines[grepl('\tCommand being timed:', lines)]
col_name <- gsub('.* ', '', gsub('"', '', col_name))
df <- data.frame(name = col_name, time = col_time)
df <- df[order(df$time),]
new_tests <- df$name

# Create new script file that is sorted by name.
lines <- readLines('test-devel.sh')
header_lines <- lines[!grepl('\\$RUNNER', lines)]
old_runner_lines <- lines[grepl('\\$RUNNER', lines)]
old_tests <- gsub('\\s+#.*', '', gsub('^.RUNNER ', '', old_runner_lines))
sorted_new_tests <- new_tests[order(new_tests)]
sorted_old_tests <- old_tests[order(old_tests)]
if (!all(sorted_new_tests == sorted_old_tests)) {
    stop(paste('Set of tests changed since profiling. Please rerun profile and try again.',
               'OLD TESTS:', paste0(sorted_old_tests, collapse = '\n'),
               'NEW TESTS:', paste0(sorted_new_tests, collapse = '\n'),
               sep = '\n'))
}
new_runner_lines <- paste0('$RUNNER ', df$name, '  # takes ', df$time, ' seconds')
lines <- c(header_lines, new_runner_lines)
writeLines(lines, con = file('test-devel.sh', 'w'))
