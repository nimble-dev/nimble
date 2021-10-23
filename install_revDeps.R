#!/usr/bin/env Rscript

revDeps <- c(
    'nimbleSMC'
    )     


## Note that if don't specify repos, the RStudio repo seems to install from
## Linux binary (that is what the .tar.gz is?) and doesn't install tests.
for (package in revDeps) {
    if (!suppressPackageStartupMessages(require(package,
                                                character.only = TRUE))) {
        install.packages(package, type = 'source', INSTALL_opts = '--install-tests', repos = 'http://cran.us.r-project.org')
    }
}

