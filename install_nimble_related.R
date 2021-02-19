#!/usr/bin/env Rscript

requirements <- c(
    'nimbleSMC'
)     

for (package in requirements) {
    if (suppressPackageStartupMessages(require(package,
                                               character.only = TRUE))) {
        ## We need to make sure to install tests.
        ## This will force re-installation, dealing with fact that Travis seems to cache installed packages.
        remove.packages(package)  
    }
    install.packages(package, repos = 'http://cran.us.r-project.org', INSTALL_opts = '--install-tests')
}
