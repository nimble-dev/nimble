# environment that holds user-provided information

nimbleUserNamespace <- as.environment(list(sessionSpecificDll = NULL)) 
# new.env() here fails with: Error in as.environment(pos) : using 'as.environment(NULL)' is defunct when testing package loading during INSTALL

# options used for NIMBLE package
# These options are for development use at this point.
.nimbleOptions <- as.environment(
    list(
        useCppADoptimize = TRUE,
        useADcholAtomic = TRUE, # If TRUE, use nimble's CppAD atomic for cholesky decomposition
        useADsolveAtomic = TRUE, # If TRUE, use nimble's CppAD atomic for matrix inverse
        useADmatMultAtomic = TRUE, # If TRUE, use nimble's CppAD atomic class for %*%
        useADmatInverseAtomic = TRUE, # If TRUE, use nimble's CppAD atomic class for inverse
      #  groupDetermWithGivenInCondIndSets = TRUE, # used in getConditionallyIndependentSets
        use_C_getParents = FALSE,
        useSafeDeparse = TRUE,
        useNewConfigureMCMC = TRUE,
        oldConjugacyChecking = FALSE,
        disallow_multivariate_argument_expressions = TRUE,
        stop_after_processing_model_code = FALSE,
        enableModelMacros = FALSE,
        enableMacroComments = FALSE,
        codeInMacroComments = FALSE,
        allowDynamicIndexing = TRUE,
        nimbleProjectForTesting = NULL,  ## only used by withTempProject and compileNimble in testing code.
        stopCompilationBeforeLinking = NULL,
        experimentalNewSizeProcessing = FALSE,
        experimentalSelfLiftStage = FALSE,
        enableSpecialHandling = FALSE,
        pauseAfterWritingFiles = FALSE,
        CppAD_directory = NA,
        enableDerivs = TRUE,
        buildModelDerivs = FALSE,
        doADerrorTraps = TRUE,
        useADreconfigure = TRUE,
        determinePredictiveNodesInModel = TRUE,
        getDependenciesIncludesPredictiveNodes = TRUE,    ## may be toggled off during MCMC sampler building
        convertSingleVectorsToScalarsInSetupArgs = TRUE,
        messagesWhenBuildingOrFinalizingCppObjects = FALSE,
        indexDrop = TRUE,
        debugRCfunProcessing = FALSE,
        debugCppLineByLine = FALSE,
        debugNFProcessing = FALSE,
        debugSizeProcessing = FALSE,
        compileAltParamFunctions = TRUE,
        includeCPPdists = TRUE,    ## includes dists.cpp and nimDists.cpp in the compilation.  Momentarily we have a problem on Windows.
        processBackwardsModelIndexRanges = FALSE,    ## if FALSE (default), for(i in 9:7) in model code becomes for(i in 7).  if TRUE, becomes for(i in c(9, 8, 7))
        prioritizeColonLikeBUGS = TRUE, ## if FALSE, 1:2 + 1 evaluates to 2:3, consistent with R.  If TRUE, it evalutes to 1:3, consistent with BUGS 
        useNewNimCopy = TRUE, ## for development purposes.  FALSE will give 0.3-1 behavior
        compileOnly = FALSE,
        buildInterfacesForCompiledNestedNimbleFunctions = FALSE,   ## provides interfaces, i.e. named access in R, to all variables in nested compiled nimbleFunctions
        clearNimbleFunctionsAfterCompiling = FALSE,
        checkModel = FALSE,
        checkNimbleFunction = TRUE,
        checkDuplicateNodeDefinitions = TRUE,
        precleanCompilation = TRUE,
        verbose = TRUE,
        verboseErrors = FALSE,

        ## verifies the correct posterior is created for any conjugate samplers, at run-time.
        ## if this option is changed, then congugate sampler functions can be rebuilt using:
        ## buildConjugateSamplerFunctions()
        verifyConjugatePosteriors = FALSE,

        showCompilerOutput = FALSE,

        MCMCprogressBar = TRUE,
        MCMCsaveHistory = FALSE,
        MCMCmultivariateNodesAsScalars = FALSE,
        MCMCmonitorAllSampledNodes = FALSE,
        MCMCuseConjugacy = TRUE,
        MCMCorderPriorSamplesSamplersFirst = TRUE,
        MCMCorderPosteriorPredictiveSamplersLast = TRUE,
        MCMCusePredictiveDependenciesInCalculations = FALSE,
        MCMCusePosteriorPredictiveSampler = TRUE,
        MCMCwarnUnsampledStochasticNodes = TRUE,
        MCMCRJcheckHyperparam = TRUE,
        MCMCenableWAIC = FALSE,
        useClearCompiledInADTesting = TRUE,
        unsupportedDerivativeHandling = 'error', # default is error, other options are 'warn' and 'ignore'. Handled in updateADproxyModelMethods in cppDefs_nimbleFunction.R
        errorIfMissingNFVariable = TRUE,
        stopOnSizeErrors = TRUE,
        useOldcWiseRule = FALSE, # This is a safety toggle for one change in sizeBinaryCwise, 1/24/23. After a while we can remove this.
        stripUnusedTypeDefs = TRUE
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
    option <- try(get(x, envir = .nimbleOptions), silent = TRUE)
    if(inherits(option, 'try-error'))
        return(NULL)
    return(option)
}

#' NIMBLE Options Settings
#'
#' Allow the user to set and examine a variety of global _options_
#' that affect the way in which NIMBLE operates. Call \code{nimbleOptions()}
#' with no arguments to see a list of available opions.
#' 
#' @param ... any options to be defined as one or more \code{name = value} pairs
#' or as a single \code{list} of \code{name=value} pairs.
#' @author Christopher Paciorek
#' @export
#'
#' @details \code{nimbleOptions} mimics \code{options}. Invoking
#' \code{nimbleOptions()} with no arguments returns a list with the
#'   current values of the options.  To access the value of a single option,
#'    one should use \code{getNimbleOption()}.
#'
#' @return
#' When invoked with no arguments, returns a list with the current values of all options.
#' When invoked with one or more arguments, returns a list of the the updated options with their updated values.
#'
#' @examples
#' # Set one option:
#' nimbleOptions(verifyConjugatePosteriors = FALSE)
#'
#' # Compactly print all options:
#' str(nimbleOptions(), max.level = 1)
#'
#' # Save-and-restore options:
#' old <- nimbleOptions()                    # Saves old options.
#' nimbleOptions(showCompilerOutput = TRUE,
#'               verboseErrors = TRUE)       # Sets temporary options.
#' # ...do stuff...
#' nimbleOptions(old)                        # Restores old options.
nimbleOptions <- function(...) {
    invisibleReturn <- FALSE
    args <- list(...)
    if (!length(args)) {
        # Get all nimble options.
        return(as.list(.nimbleOptions))
    }
    if (length(args) == 1 && is.null(names(args)) && is.list(args[[1]])) {
        # Unpack a single list of many args.
        args <- args[[1]]
    }
    if (is.null(names(args))) {
        # Get some nimble options.
        args <- unlist(args)
    } else {
        # Set some nimble options.
        for(i in seq_along(args)) {
            setNimbleOption(names(args)[[i]], args[[i]])
        }
        args <- names(args)
        invisibleReturn <- TRUE
    }
    out <- as.list(.nimbleOptions)[args]
    if(length(out) == 1) out <- out[[1]]
    if(invisibleReturn) return(invisible(out)) else return(out)
}

#' Temporarily set some NIMBLE options.
#'
#' @param options a list of options suitable for \code{nimbleOptions}.
#' @param expr an expression or statement to evaluate.
#' @return expr as evaluated with given options.
#' @export
#'
#' @examples
#' \dontrun{
#' if (!(getNimbleOption('showCompilerOutput') == FALSE)) stop()
#' nf <- nimbleFunction(run = function(){ return(0); returnType(double()) })
#' cnf <- withNimbleOptions(list(showCompilerOutput = TRUE), {
#'     if (!(getNimbleOption('showCompilerOutput') == TRUE)) stop()
#'     compileNimble(nf)
#' })
#' if (!(getNimbleOption('showCompilerOutput') == FALSE)) stop()
#' }
withNimbleOptions <- function(options, expr) {
    old <- nimbleOptions()
    cleanup <- substitute(do.call(nimbleOptions, old))
    do.call(on.exit, list(cleanup, add = TRUE))
    nimbleOptions(options)
    return(expr)
}

