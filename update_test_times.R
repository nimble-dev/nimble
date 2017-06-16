#!/usr/bin/env Rscript
# This script parses a travis test log file to update test_times.csv.

args = commandArgs(trailingOnly = TRUE)
if (length(args) != 1 || '-h' %in% args || '--help' %in% args) {
    cat('Usage: Rscript update_test_times.R TEST_LOG_FILE',
        'Sort the test-devel.sh script by running the cheapest tests first.',
        'The TEST_LOG_FILE should be either a file that you manually create',
        'using someting like',
        '  run_tests.sh 2>&1 > test_log.txt',
        'or a travis log file that you download from',
        '  https://travis-ci.org/nimble-dev/nimble/builds',
        'It is fine if the test log has run only a subset of tests;',
        'this script will only update the times of the tests that were run.',
        sep = '\n'
    )
    quit('no', 1)
}

# Read old times in case some tests were missing from the log.
df <- read.table('test_times.csv', sep = '\t', header = TRUE)

#  Search log file for tests.
for (logFileName in args) {
    lines = readLines(logFileName)
    col_name <- lines[grepl('\tCommand being timed:.*nimble', lines)]
    col_name <- gsub('.*\\^(.*)\\$.*', 'test-\\1.R', col_name)
    col_usr <- lines[grepl('\tUser time .seconds.:', lines)]
    col_sys <- lines[grepl('\tSystem time .seconds.:', lines)]
    col_time <- as.numeric(gsub('.*: ', '', col_usr)) + as.numeric(gsub('.*: ', '', col_sys))
    col_time <- col_time[(1 + length(col_time) - length(col_name)):length(col_time)]
    df <- rbind(df, data.frame(time = col_time, filename = col_name))
    df <- df[!duplicated(df$filename, fromLast = TRUE),]
}

# Write table sorted by time.
df <- df[order(df$time),]
write.table(df, 'test_times.csv', sep = '\t', row.names = FALSE, quote = FALSE)
