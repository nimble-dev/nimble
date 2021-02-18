#!/usr/bin/env Rscript

requirements <- c(
    'nimbleSMC'
)     

for (package in requirements) {
    if (suppressPackageStartupMessages(require(package,
                                               character.only = TRUE))) {
        remove.packages(package)  ## need to make sure install tests so this will force re-installation
    }
    install.packages(package, repos = 'http://cran.us.r-project.org', INSTALL_opts = '--install-tests')
}
