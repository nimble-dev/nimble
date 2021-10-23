#!/usr/bin/env Rscript

revDeps <- c(
    'nimbleSMC'
    )     

for (package in revDeps) {
    if (!suppressPackageStartupMessages(require(package,
                                                character.only = TRUE))) {
        install.packages(package, repos = 'http://cran.us.r-project.org')
    }
}

