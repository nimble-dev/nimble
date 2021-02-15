#!/usr/bin/env Rscript

requirements <- c(
    'nimbleSMC'
)     

for (package in requirements) {
    if (!suppressPackageStartupMessages(require(package,
                                                character.only = TRUE))) {
        install.packages(package, repos = 'http://cran.us.r-project.org')
    }
}
