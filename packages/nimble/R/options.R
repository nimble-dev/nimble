# environment that holds user-provided information
nimbleUserNamespace <- as.environment(list()) 
# new.env() here fails with: Error in as.environment(pos) : using 'as.environment(NULL)' is defunct when testing package loading during INSTALL

#' options used for NIMBLE package
#'
#' @details These options are for development use at this point.
.nimbleOptions <- as.environment(
    list(
        convertSingleVectorsToScalarsInSetupArgs = TRUE,
        messagesWhenBuildingOrFinalizingCppObjects = FALSE,
        indexDrop = TRUE,
        debugRCfunProcessing = FALSE,
        debugCppLineByLine = FALSE,
        debugNFProcessing = FALSE,
        compileAltParamFunctions = TRUE,
        verifyConjugatePosteriors = FALSE,        ## verifies the correct posterior is created for any conjugate samplers, at run-time
                                                 ## The sampler functions are created at loading time, so this only works at installation
        includeCPPdists = TRUE,    ## includes dists.cpp and nimDists.cpp in the compilation.  Momentarily we have a problem on Windows.
        processBackwardsModelIndexRanges = FALSE,    ## if FALSE (default), for(i in 9:7) in model code becomes for(i in 7).  if TRUE, becomes for(i in c(9, 8, 7))
        prioritizeColonLikeBUGS = TRUE, ## if FALSE, 1:2 + 1 evaluates to 2:3, consistent with R.  If TRUE, it evalutes to 1:3, consistent with BUGS 
        useNewNimCopy = TRUE, ## for development purposes.  FALSE will give 0.3-1 behavior

        ## default settings for MCMC samplers
        ## (formerly controlDefaultList appearing in MCMCspec.R)
        MCMCcontrolDefault_adaptive = TRUE,
        MCMCcontrolDefault_adaptScaleOnly = FALSE,
        MCMCcontrolDefault_adaptInterval = 200,
        MCMCcontrolDefault_scale = 1,
        MCMCcontrolDefault_propCov = 'identity',
        MCMCcontrolDefault_sliceWidth = 1,
        MCMCcontrolDefault_sliceMaxSteps = 100
    )
    )

# sets a single option
setNimbleOption <- function(name, value) {
    assign(name, value, envir = .nimbleOptions)
    invisible(value)
}

#' Get NIMBLE Option
#'
#' Allow the user to get the value of a global _option_
#' that affects the way in which NIMBLE operates
#' 
#' @param x a character string holding an option name 
#' @author Christopher Paciorek
#' @export
#' @return The value of the option.
#' @examples
#' getNimbleOption('verifyConjugatePosteriors')
getNimbleOption <- function(x) {
    get(x, envir = .nimbleOptions)
}

#' NIMBLE Options Settings
#'
#' Allow the user to set and examine a variety of global _options_
#' that affect the way in which NIMBLE operates
#' 
#' @param ... any options to be defined as one or more 'name = value' pairs.
#' Options can also be passed by giving a single unnamed argument that is a named list.
#' @author Christopher Paciorek
#' @export
#' @details \code{nimbleOptions} mimics \code{options}. Invoking
#' \code{nimbleOptions()} with no arguments returns a list with the
#'   current values of the options.  To access the value of a single option,
#'    one should use \code{getNimbleOption()}.
#' @return When invoked with no arguments, a list with the current values of all options. 
#' @examples
#' nimbleOptions(verifyConjugatePosteriors = FALSE)
nimbleOptions <- function(...) {
    if(!is.list(args)) args <- list(...)
    if(length(args)) {
        for(i in seq_along(args))
            setNimbleOption(names(args)[[i]], args[[i]])
        invisible(args)
    } else return(as.list(.nimbleOptions))
}

