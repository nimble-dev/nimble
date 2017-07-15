try({
    packageMakeconfFile <- file.path(R.home('etc'), .Platform$r_arch, 'Makeconf')
    if(file.exists(packageMakeconfFile)) {
        file.copy(packageMakeconfFile, 'Makeconf')
        lines <- readLines('Makeconf')
        lines <- gsub('-g', '', lines)
        writeLines(lines, con = 'Makeconf')
    }
    })
