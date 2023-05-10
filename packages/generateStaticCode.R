#!/usr/bin/env Rscript

# This script is used to generate the C++ predefinedNimbleList.* files in nimble/inst.
# You should run this script if you make changes to the internal representation of nimbleLists
# or if you add a new predefined nimbleList type (i.e. with predefined = TRUE).
 
directions <- 'Directions:
1. Install clang-format.
   Ubuntu: sudo apt-get install clang-format
   OS X: brew install clang-format
   Windows: Download from http://llvm.org/builds
2. Temporarily set GENERATE_STATIC_CODE <- TRUE in the nimbleList function
   in nimbleList_core.R, and reinstall nimble.
3. Run this script.
4. Review changes this script has made to the predefinedNimbleLists.* files.
5. Revert the temporary change to GENERATE_STATIC_CODE from step 2,
   and reinstall nimble.
6. If desired (default = yes), manually update
   inst/include/nimble/dynamicRegistrations.h to have entries for the three
   new SEXP-returning functions for each predefined nimbleList.
   See examples in that file.
'

how_to_add_a_new_predefined_nimbleList <- 'How to:
1. Enter the nimbleList call in nimbleList_core, with predefined = TRUE.  See examples there.
2. If the new new predefined nimbleList should be returned from any DSL functions (not nimbleFunctions), add an entry to nimbleListReturningFunctionList.  Study examples there, looking C++ code as well.
3. Go through the above Directions.
4. The result will be that the C++ for the new nimbleList has been added to the
 package source code and will be used from there.
'

# NB: This script did not work when adding AGHQuad_params and AGHQuad_summary,
#  so those were done by hand-coding.

library(methods)
suppressPackageStartupMessages(library(nimble))

# Warn user if script is being used incorrectly.
(function(){
    con <- file(file.path('nimble', 'R', 'nimbleList_core.R'))
    lines <- readLines(con)
    close(con)
    if (sum(grep('GENERATE_STATIC_CODE <- TRUE', lines)) == 0) {
        stop(directions, call. = FALSE)
    }
})()

# This finds the latest generated .h and .cpp files.
findGeneratedSources <- function() {
    root <- file.path(tempdir(), 'nimble_generatedCode')
    files <- file.info(list.files(root, full.names = TRUE))
    files <- files[order(files$mtime, decreasing = TRUE),]
    paths <- rownames(files)
    paths_h <- paths[grep('\\bP_\\w+\\.h$', paths)]
    paths_cpp <- paths[grep('\\bP_\\w+\\.cpp$', paths)]
    if(length(paths_h) == 0 || length(paths_cpp) == 0) stop('Failed to generate sources')
    ret <- system2('clang-format', list('-i', paths_h[1], paths_cpp[1]))
    if(ret != 0) stop('Please install clang-format and try again')
    return(list(.h = paths_h[1], .cpp = paths_cpp[1]))
}

readRelevantLines <- function(filename) {
    con <- file(filename, open = 'r')
    lines = readLines(con)
    close(con)

    # Parse file using a little state machine.
    step <- 'HEADER'
    for(i in seq_along(lines)) {
        if (step == 'HEADER') {
            if (grepl('P_1_REMOVE_THIS_CODE', lines[i])
                & grepl('include', lines[i])) {
                beginBody <- i+1
                step <- 'BODY'
            }
        } else if (step == 'BODY') {
            if (grepl('REMOVE_THIS_CODE', lines[i])) {
                endBody <- i - 1
                break
            }
        }
    }
    return(lines[beginBody:endBody])
}

writeFile <- function(filename, lines) {
    cat('Generating', filename, '\n')
    con <- file(filename, open = 'w')
    writeLines(lines, con)
    close(con)
    ret <- system2('clang-format', c('--sort-includes=0', c('-i', "-style=file", filename)))
    if(ret != 0) stop('Please install clang-format and try again')
}

main <- function() {
    # Generate code with some unwanted garbage.
    nimFun <- nimbleFunction(
        name = 'REMOVE_THIS_CODE',
        setup = TRUE,
        run = function() {
            return(0)
            returnType(double(0))
        },
        methods = list(
            eigenStub = function() {
                return(eigenNimbleList$new())
                returnType(eigenNimbleList())
            },
            svdStub = function() {
                return(svdNimbleList$new())
                returnType(svdNimbleList())
            },
            optimResultStub = function() {
                return(optimResultNimbleList$new())
                returnType(optimResultNimbleList())
            },
            optimControlStub = function() {
                return(optimControlNimbleList$new())
                returnType(optimControlNimbleList())
            },
            ADStub = function() {
              return(ADNimbleList$new())
              returnType(ADNimbleList())
            },
            WAICStub = function() {
              return(waicList$new())
              returnType(waicList())
            },
            WAICdetailsStub = function() {
              return(waicDetailsList$new())
              returnType(waicDetailsList())
            }
        )
    )()

    # Compilation should with a SHLIBCreationError due to duplicate symbols,
    # but any other error is unexpected.    
    tryCatch(compileNimble(nimFun), error = function(e){
        if(!inherits(e, 'SHLIBCreationError')) stop(e)
    })

    # Read the relevant parts of both files.
    files <- findGeneratedSources()
    lines_h <- readRelevantLines(files$.h)
    lines_cpp <- readRelevantLines(files$.cpp)

    # Patch header to avoid #including nimOptim.h.
    badIncludeLine <- grep('#include <nimble/nimOptim.h>', lines_h)[1]
    lines_h <- lines_h[-badIncludeLine]

    lastIncludeLine <- max(which(grepl('#include', lines_h)))
    lines_h <- append(lines_h, '#undef eval', lastIncludeLine)
    
    lastIncludeLine <- max(which(grepl('#include', lines_cpp)))
    lines_cpp <- append(lines_cpp, '#undef eval', lastIncludeLine)

  # Write the header.
    provenance <- c(
        '// DO NOT EDIT BY HAND.',
        '// This file was automatically generated by nimble/packages/generateStaticCode.R',
        ''
    )
    writeFile(
        file.path('nimble', 'inst', 'include', 'nimble', 'predefinedNimbleLists.h'),
        c(
            provenance,
            '#ifndef __NIMBLE_PREDEFINEDNIMBLELISTS_H',
            '#define __NIMBLE_PREDEFINEDNIMBLELISTS_H',
            '',
            '#include <nimble/smartPtrs.h>',
            lines_h,
            '#endif  // __NIMBLE_PREDEFINEDNIMBLELISTS_H'
        )
    )
    writeFile(
        file.path('nimble', 'inst', 'CppCode', 'predefinedNimbleLists.cpp'),
        c(
            provenance,
            '#include <nimble/predefinedNimbleLists.h>',
            lines_cpp
        )
    )
}

main()
