#!/usr/bin/env Rscript

revDeps <- c(
    'nimbleSMC'
    )     


## Note that if don't specify repos, the RStudio repo seems to install from
## Linux binary (that is what the .tar.gz is?) and doesn't install tests.
for (package in revDeps) {
        install.packages(package, type = 'source', INSTALL_opts = '--install-tests')
    }
    print(list.files(system.file('tests','testthat', package='nimbleSMC')))
    print(list.files(system.file('', package='nimbleSMC')))
}

