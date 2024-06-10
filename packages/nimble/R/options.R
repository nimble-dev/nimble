# environment that holds user-provided information

nimbleUserNamespace <- as.environment(list(sessionSpecificDll = NULL)) 
# new.env() here fails with: Error in as.environment(pos) : using 'as.environment(NULL)' is defunct when testing package loading during INSTALL

nimbleUserNamespace$.optimizers <- as.environment(list())

#' Set or get an optimization function to be used by \code{nimOptim}
#'
#' Add, check, or remove an R optimization function to/from NIMBLE's set of
#' registered optimization functions that can be called from \code{nimOptim}.
#'
#' @param name character string, giving the name of optimization method that
#'   can be referred to by the \code{method} argument of \code{nimOptim}
#'   (aka \code{optim} in a nimbleFunction).
#'
#' @param value An optimization function with specifications described below and
#'   in \code{\link{nimOptim}}. If \code{value} is \code{NULL}, then \code{name}
#'   will be found in NIMBLE's set of registered optimizer names. If
#'   \code{value} is missing, the registered optimizer for \code{name} will be
#'   returned.
#'
#' @details
#'
#' When programming in nimbleFunctions, \code{optim}, which is converted
#' automatically to \code{\link{nimOptim}}, provides a generalization of R's
#' \code{optim} methods for optimization. If one of the supported original
#' \code{optim} methods is not chosen with the \code{method} argument to
#' \code{nimOptim}, an arbitrary method name can be given. If that name has been
#' registered as a \code{name} by a call to \code{nimOptimMethod}, then the
#' corresponding function (\code{value}) will be called for optimization.
#'
#' The function \code{value} must perform minimization. If the call to
#' \code{nimOptim} includes a control list with either \code{fnscale} (which, if
#' negative, turns the problem into a maximization) or \code{parscale}, these
#' will be managed by \code{nimOptim} outside of the optimizer such that the
#' optimization should be minimization.
#'
#' The function \code{value} must take named arguments \code{par} (initial
#' parameter vector), \code{fn} (objective function), \code{gr} (optional
#' gradient function), \code{he} (optional Hessian function), \code{lower}
#' (vector of lower bounds), \code{upper} (vector of upper bounds),
#' \code{control} (arbitrary control list), and \code{hessian} (logical
#' indicating whether a Hessian at the optimum is requested). It must return a
#' list with elements elements \code{par} (parameter values of the optimium,
#' i.e., "arg min"), \code{value} (function value at the minimum),
#' \code{convergence} (should be 0 if convergence occurred), \code{counts}
#' (optional vector of counts of calls to \code{fn}, \code{gr}, and \code{he}),
#' \code{evaluations} (optional total function evaluations), \code{message}
#' (optional character message), and \code{hessian} (optional Hessian matrix,
#' which may be NULL).
#'
#' If the call to \code{nimOptim} has \code{hessian=TRUE}, that will be passed
#' as \code{hessian=TRUE} to the optimizer. However, if the optimizer returns a
#' \code{NULL} in the \code{hessian} element of the return list, then
#' \code{nimOptim} will calculate the Hessian by finite element differences.
#' Hence, an optimizer need not provide a Hessian capability.
#'
#' The \code{control} list passed from \code{nimOptim} to the optimizer will
#' have only a limited set of the \code{optim} control list options. These will
#' include \code{abstol}, \code{reltol}, \code{maxit}, and \code{trace}. The
#' optimizer may use these as it wishes. Other control options for a particular
#' optimizer must be managed in some other way.
#'
#' Note that it is possible to use \code{browser()} inside of \code{value}, or
#' to set \code{debug(value)}, to enter a browser when the optimizer
#' (\code{value}) is called and then inspect its arguments to make sense of the
#' situation.
#'
#' This whole feature is particularly helpful when the nimbleFunction using
#' \code{nimOptim} has been compiled by \code{compileNimble}. Many optimizers
#' are available through R, so \code{nimOptim} arranges to call a named
#' (registered) optimizer in R, while providing \code{fn} and optionally
#' \code{gr} or \code{he} as functions that will call the compiled (by nimble)
#' versions of the corresponding functions provided in the call to
#' \code{nimOptim}.
#'
#' R's optimizer \code{nlminb} is automatically registered under the name
#' \code{"nlminb"}.
#'
#' @export
nimOptimMethod <- function(name, value) {
  if(missing(value))
    nimbleUserNamespace$.optimizers[[name]]
  else
    nimbleUserNamespace$.optimizers[[name]] <- value
}

nimOptimMethod("nlminb",
               function(par, fn, gr, he, lower, upper, control, hessian) {
                 control_nlminb <- list(
                   abs.tol = control$abstol,
                   rel.tol = control$reltol,
                   iter.max = control$maxit,
                   trace = control$trace
                 )
                 invalid <- function(x) is.null(x) || is.na(x) || is.infinite(x)
                 if(invalid(control_nlminb$abs.tol)) control_nlminb$abs.tol <- 0
                 if(invalid(control_nlminb$rel.tol)) control_nlminb$rel.tol <- 1e-10
                 if(invalid(control_nlminb$iter.max)) control_nlminb$iter.max <- 150
                 if(invalid(control_nlminb$trace)) control_nlminb$trace <- 0
                 # NB: control$parscale and control$fnscale are applied internally
                 result <- nlminb(par, objective=fn, gradient=gr, hessian=he,
                                  lower = lower, upper = upper, control=control_nlminb)
                 result$value <- result$objective
                 result$objective <- NULL
                 result$counts <- result$evaluations
                 result$evaluations <- NULL
                 # We could do hessian here like this, but we will do it in C++
                 # if not returned from here, so we ignore it here.
                 ## if(isTRUE(hessian)) {
                 ##   # do we need to worry if control has a parscale element?
                 ##   hessian_result <- optimHess(result$par, fn=fn, gr=gr, control=control)
                 ## }
                 result
               }
               )

# options used for NIMBLE package
# These options are for development use at this point.
.nimbleOptions <- as.environment(
    list(
        newAssignSizeChecks = TRUE, # If TRUE, do additional checking when [matrix|vector] <- [matrix|vector]. Eventually if successful this option can be removed. It is a failsafe to disable new checks if they give spurious traps.
        allowNFobjInModel = TRUE, # If TRUE, allow use of nimbleFunctions with setup code as model dist or fxn.
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
        stripUnusedTypeDefs = TRUE,
        digits = NULL,
        enableVirtualNodeFunctionDefs = FALSE
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

