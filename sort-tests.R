#!/usr/bin/env Rscript
# This script parses a travis test log file to update test_times.csv.

args = commandArgs(trailingOnly = TRUE)
if(length(args) != 1 || '-h' %in% args || '--help' %in% args) {
    cat('Usage: Rscript sort-tests.R TEST_LOG_FILE',
        'Sort the test-devel.sh script by running the cheapest tests first.',
        'The TEST_LOG_FILE should be either a file that you manually create',
        'using someting like',
        '  run_tests.sh 2>&1 > test_log.txt',
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
col_name <- gsub('.*(\\<test-.*\\.R).*', '\\1', col_name)
df <- data.frame(time = col_time, filename = col_name)
df <- df[order(df$time),]
write.table(df, 'test_times.csv', sep = '\t', row.names = FALSE, quote = FALSE)
