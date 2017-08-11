#!/usr/bin/env Rscript

header <- readLines('GPL_header.txt')
penultimate_line <- ' * https://www.R-project.org/Licenses/'

cppFiles = c(Sys.glob('nimble/inst/include/nimble/*.h'),
             Sys.glob('nimble/inst/CppCode/*.cpp'))
for (filename in cppFiles) {
    lines <- readLines(filename)
    old <- grep(penultimate_line, lines)
    if (length(old)) {
        lines <- lines[(old[1] + 3):length(lines)]
    }
    lines <- c(header, '', lines)
    writeLines(lines, con = filename)
}
