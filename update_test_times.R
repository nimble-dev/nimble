#!/usr/bin/env Rscript
# This script parses a travis test log file to update test_times.csv.

args = commandArgs(trailingOnly = TRUE)
if (length(args) != 1 || '-h' %in% args || '--help' %in% args) {
    cat('Usage: Rscript update_test_times.R [TEST_LOG_FILES | BUILD_NUMBER]',
        'Sort the test-devel.sh script by running the cheapest tests first.',
        'The TEST_LOG_FILES should be either file that you manually create',
        'using someting like',
        '  run_tests.sh 2>&1 > test_log.txt',
        'or travis logs that you download from',
        '  https://travis-ci.org/nimble-dev/nimble/builds/<BUILD_NUMBER>',
        'If a build number is specified, this script will downlod the logs.',
        'It is fine if the test logs have run only a subset of tests;',
        'this script will only update the times of the tests that were run.',
        sep = '\n'
    )
    quit('no', 1)
}

# Read old times in case some tests were missing from the log.
df <- read.table('test_times.csv', sep = '\t', header = TRUE)

# Parse a log file.
addLines <- function(df, lines) {
    col_name <- lines[grepl('\tCommand being timed:.*nimble', lines)]
    col_name <- gsub('.*\\^(.*)\\$.*', 'test-\\1.R', col_name)
    col_usr <- lines[grepl('\tUser time .seconds.:', lines)]
    col_sys <- lines[grepl('\tSystem time .seconds.:', lines)]
    col_time <- as.numeric(gsub('.*: ', '', col_usr)) + as.numeric(gsub('.*: ', '', col_sys))
    col_time <- col_time[(1 + length(col_time) - length(col_name)):length(col_time)]
    cat('Found', length(col_name), 'test results\n')
    df <- rbind(df, data.frame(time = col_time, filename = col_name))
    df <- df[!duplicated(df$filename, fromLast = TRUE),]
    return(df)
}

# Find log files locally or from travis.
for (arg in args) {
    if (is.na(suppressWarnings(as.integer(arg)))) {
        # Assume arg is a filename.
        cat('Reading local travis log', arg, '\n')
        lines <- readLines(arg)
        df <- addLines(df, lines)
    } else {
        # Assume arg is a travis build id.
        if (!require('httr')) stop('Missing required package httr')
        if (!require('jsonlite')) stop('Missing required package jsonlite')
        url <- paste0('https://api.travis-ci.org/builds/', arg)
        response <- GET(url)
        if (response$status_code != 200) stop(paste('Failed to GET', url))
        build <- fromJSON(content(response, type = 'text', encoding = 'UTF-8'))
        jobs <- build$matrix$id[build$matrix$config$os == 'linux']
        for (job in jobs) {
            url <- paste0('https://api.travis-ci.org/jobs/', job, '/log.txt?deansi=true')
            cat('Downloading travis log from', url, '\n')
            response <- GET(url)
            if (response$status_code != 200) stop(paste('Failed to GET', url))
            lines <- unlist(strsplit(content(response, type = 'text', encoding = 'UTF-8'), '\r'))
            df <- addLines(df, lines)
        }
    }
}

# Write table sorted by time.
df <- df[order(df$time),]
write.table(df, 'test_times.csv', sep = '\t', row.names = FALSE, quote = FALSE)
