## creates unique labels ('nfRefClass1') for the reference class names for nimbleFunctions
nf_refClassLabelMaker <- labelFunctionCreator('nfRefClass')

#' create a virtual nimbleFunction, a base class for other nimbleFunctions
#'
#' define argument types and returnType for the \code{run} function and any \code{methods}, to be used in the \code{contains} argument of \code{nimbleFunction}
#'
#' @param contains Not yet functional
#' @param run      A NIMBLE function that will only be used to inspect its argument types and returnType.
#' @param methods  An optional named list of NIMBLE functions that will also only be used for inspecting argument types and returnTypes.
#' @param name     An optional name used internally by the NIMBLE compiler.  This is usually omitted and NIMBLE provides one.
#' @param methodControl An optional list that allows specification of methods with defaults. 
#'
#' @author NIMBLE development team
#'
#' @export
#'
#' @details See the NIMBLE \href{https://r-nimble.org/html_manual/cha-welcome-nimble.html}{User Manual} section on nimbleFunctionLists for explanation of how to use a virtual nimbleFunction.
#'
#' @seealso \code{\link{nimbleFunction}}
#' 
#' @return An object that can be passed as the \code{contains} argument to \code{nimbleFunction} or as the argument to \code{nimbleFunctionList}
nimbleFunctionVirtual <- function(contains = NULL,
                                  run = function() { },
                                  methods     = list(),
                                  name        = NA,
                                  methodControl = list()) {
    virtual <- TRUE
    if(is.na(name)) name <- nf_refClassLabelMaker()
    className <- name
    ## We make this look like a nimbleFunction in relevant ways for compilation
    methodList <- c(list(run = run), methods)   # create a list of the run function, and all other methods
    methodList <- lapply(methodList, nfMethodRC, check = FALSE)
    generatorFunction <- function() {}
    force(contains)
    nfRefClassDef <- nfRefClass <- NULL ## Existence of these makes this treated like a nfGenerator
    environment(generatorFunction) <- GFenv <- new.env()
    parent.env(GFenv) <- parent.frame()
    for(var in c('generatorFunction','nfRefClassDef','nfRefClass','run','methods','methodList','name', 'className', 'contains', 'virtual', 'methodControl')) {
        GFenv[[var]] <- get(var)
    }
    return(generatorFunction)
}

#' create a nimbleFunction
#'
#' create a nimbleFunction from a setup function, run function, possibly other methods, and possibly inheritance via \code{contains}
#'
#' @param setup An optional R function definition for setup processing.
#' @param run An optional NIMBLE function definition that executes the primary job of the nimbleFunction
#' @param methods An optional named list of NIMBLE function definitions for other class methods.
#' @param globalSetup For internal use only
#' @param contains An optional object returned from \code{\link{nimbleFunctionVirtual}} that defines arguments and returnTypes for \code{run} and/or methods, to which the current nimbleFunction must conform
#' @param buildDerivs A list of names of function methods for which to build derivatives capabilities.
#' @param name An optional name used internally, for example in generated C++ code.  Usually this is left blank and NIMBLE provides a name.
#' @param check Boolean indicating whether to check the run code for function calls that NIMBLE cannot compile. Checking can be turned off for all calls to \code{nimbleFunction} using \code{nimbleOptions(checkNimbleFunction = FALSE)}.
#' @param where An optional \code{where} argument passed to \code{setRefClass} for where the reference class definition generated for this nimbleFunction will be stored.  This is needed due to R package namespace issues but should never need to be provided by a user.
#'
#' @author NIMBLE development team
#'
#' @export
#'
#' @details
#' This is the main function for defining nimbleFunctions.  A lot of information is provided in the NIMBLE \href{https://r-nimble.org/html_manual/cha-welcome-nimble.html}{User Manual}, so only a brief summary will be given here.
#'
#' If a \code{setup} function is provided, then \code{nimbleFunction} returns a generator: a function that when called with arguments for the setup function will execute that function and return a specialized nimbleFunction.   The \code{run} and other methods can be called using \code{$} like in other R classes, e.g. \code{nf$run()}. The methods can use objects that were created in or passed to the \code{setup} function.
#'
#' If no \code{setup} function is provided, then \code{nimbleFunction} returns a function that executes the \code{run} function.  It is not a generator in this case, and no other \code{methods} can be provided.
#'
#' If one wants a generator but does not need any setup arguments or code, \code{setup = TRUE} can be used.
#'
#' See the NIMBLE \href{https://r-nimble.org/html_manual/cha-welcome-nimble.html}{User Manual} for examples.
#'
#' For more information about the \code{contains} argument, see the section on nimbleFunctionLists.
nimbleFunction <- function(setup         = NULL,
                           run           = function() { },
                           methods       = list(),
                           globalSetup   = NULL,
                           contains      = NULL,
                           buildDerivs  = list(),
                           name          = NA,
                           check         = getNimbleOption('checkNimbleFunction'),
                           where         = getNimbleFunctionEnvironment()
                           ) {
    force(where) # so that we can get to namespace where a nf is defined by using topenv(parent.frame(2)) in getNimbleFunctionEnvironment()
    if(is.logical(setup)) if(setup) setup <- function() {} else setup <- NULL

    ## Ceck for correct entries in `buildDerivs` separately from `nfMethodRC$new()` because
    ## that only has access to `thisBuildDerivs`, and we need to check if `buildDerivs` is set
    ## for the method on which `nimDerivs` is called.
    tmp <- sapply(c(list(run = run), methods), nf_checkDSLcode_buildDerivs, buildDerivs)
    
    if(is.null(setup)) {
        if(length(methods) > 0) stop('Cannot provide multiple methods if there is no setup function.  Use "setup = function(){}" or "setup = TRUE" if you need a setup function that does not do anything', call. = FALSE)
        if(!is.null(contains)) stop('Cannot provide a contains argument if there is no setup function.  Use "setup = function(){}" or "setup = TRUE" if you need a setup function that does not do anything', call. = FALSE)
        thisBuildDerivs <- FALSE
        if(isTRUE(getNimbleOption("enableDerivs"))) {
            if(isTRUE(buildDerivs)) buildDerivs <- list(run = list()) ## empty list means TRUE with no configuration information
            if(isFALSE(buildDerivs)) buildDerivs <- list()
            if(identical(buildDerivs, 'run')) buildDerivs <- list(run = list())
            thisBuildDerivs <- buildDerivs[['run']]
            if(is.null(thisBuildDerivs)) thisBuildDerivs <- FALSE
        }
        return(RCfunction(run, name = name, check = check, buildDerivs = thisBuildDerivs, where = where))
    }

    if(isTRUE(getNimbleOption("enableDerivs")) && isTRUE(buildDerivs))
        stop("'buildDerivs' cannot be 'TRUE' when a setup function is provided. Please specify the specific method(s) for which 'buildDerivs' should be set.")

    virtual <- FALSE
    # we now include the namespace in the name of the RefClass to avoid two nfs having RefClass of same name but existing in different namespaces
    if(is.na(name)) name <- nf_refClassLabelMaker(envName = environmentName(where))
    className <- name
    methodList <- c(list(run = run), methods)   # create a list of the run function, and all other methods
    # simply pass in names of vars in setup code so that those can be used in nf_checkDSLcode; to be more sophisticated we would only pass vars that are the result of nimbleListDefs or nimbleFunctions
    if(isTRUE(getNimbleOption('enableDerivs'))
       && length(buildDerivs)>0) {
        ## convert buildDerivs to a format of name = list(controls...)
        if(is.character(buildDerivs)) {
            buildDerivs <- structure(
                lapply(buildDerivs, function(x) list()),
                names = buildDerivs)
        }
    } else if(!isTRUE(getNimbleOption('enableDerivs'))
              && length(buildDerivs) > 0)
        stop('To build nimbleFunction derivatives, you must first set "nimbleOptions(enableDerivs = TRUE)".')
    origMethodList <- methodList
    methodList <- list()
    setupVarNames = c(all.vars(body(setup)), names(formals(setup)))

    for(iM in seq_along(origMethodList)) {
        thisBuildDerivs <- FALSE
        if(getNimbleOption('enableDerivs')
           && length(buildDerivs)>0) {
            thisBuildDerivs <-  !is.null(buildDerivs[[ names(origMethodList)[iM] ]])
        }
        methodList[[iM]] <- nfMethodRC(origMethodList[[iM]],
                                       check = check,
                                       methodNames = names(origMethodList),
                                       setupVarNames = setupVarNames,
                                       buildDerivs = thisBuildDerivs,
                                       where = where)
    }
    names(methodList) <- names(origMethodList)
    
    ## methodList <- lapply(methodList,
    ##                      nfMethodRC,
    ##                      check = check,
    ##                      methodNames = names(methodList),
    ##                      setupVarNames = c(all.vars(body(setup)), names(formals(setup))),
    ##                      where = where)
    ## if(nimbleOptions('enableDerivs')
    ##    && length(buildDerivs)>0) {
    ##     for(iM in seq_along(methodList)) {
    ##         thisBuildDerivs <-  buildDerivs[[ names(methodList)[iM] ]]
    ##         if(!is.null(thisBuildDerivs))
    ##             methodList[[iM]]$buildDerivs <- thisBuildDerivs
    ##     }
    ## }
    
    ## record any setupOutputs declared by setupOutput()
    setupOutputsDeclaration <- nf_processSetupFunctionBody(setup, returnSetupOutputDeclaration = TRUE)
    declaredSetupOutputNames <- nf_getNamesFromSetupOutputDeclaration(setupOutputsDeclaration)
    rm(setupOutputsDeclaration)
    ## create the reference class definition
    
    nfRefClassDef <- nf_createRefClassDef(setup, methodList, className, globalSetup, declaredSetupOutputNames, contains)
    nfRefClass    <- eval(nfRefClassDef)
    .namesToCopy <- nf_namesNotHidden(names(nfRefClass$fields()))
    .namesToCopyFromGlobalSetup <- intersect(.namesToCopy,
                                             if(!is.null(globalSetup)) nf_assignmentLHSvars(body(globalSetup)) else character(0))
    .namesToCopyFromSetup <- setdiff(.namesToCopy, .namesToCopyFromGlobalSetup)
    ## create a list to hold all specializations (instances) of this nimble function.  The following objects are accessed in environment(generatorFunction) in the future
    ## create the generator function, which is returned from nimbleFunction()
    generatorFunction <- eval(nf_createGeneratorFunctionDef(setup))
    force(contains) ## eval the contains so it is in this environment
    formals(generatorFunction) <- nf_createGeneratorFunctionArgs(setup, parent.frame())
    environment(generatorFunction) <- GFenv <- new.env()
    parent.env(GFenv) <- parent.frame()

        .globalSetupEnv <- new.env()
    if(!is.null(globalSetup)) {
        if(!is.function(globalSetup)) stop('If globalSetup is not NULL, it must be a function', call. = FALSE)
        if(!length(formals(globalSetup))==0) stop('globalSetup cannot take input arguments', call. = FALSE)
        eval(body(globalSetup), envir = .globalSetupEnv)
    }

    for(var in c('generatorFunction','nfRefClassDef','nfRefClass','setup','run','methods','methodList','name', 'className', 'contains', 'buildDerivs', 'virtual', '.globalSetupEnv', '.namesToCopy', '.namesToCopyFromGlobalSetup', '.namesToCopyFromSetup','declaredSetupOutputNames','.globalSetupEnv')) {
        GFenv[[var]] <- get(var)
    }
    return(generatorFunction)
}

## See https://github.com/nimble-dev/nimble/wiki/Developer-backdoor-to-manually-replace-generated-C
`specialHandling<-` <- function(x, value) {
    if(is.rcf(x)) {
        environment(x)[['.specialHandling']] <- value
        return(x)
    }
    if(inherits(x, 'modelBaseClass')) {
        assign('.specialHandling', value, x)
        return(x)
    }
    stop('special handling only works for RCfunctions and models')
}

specialHandling <- function(nf) {
    if(is.rcf(nf)) {
        return(environment(nf)[['.specialHandling']])
    }
    if(inherits(nf, 'modelBaseClass')) {
        return(nf$.specialHandling)
    }
    stop('special handling only works for RCfunctions and models')
}

filenameFromSpecialHandling <- function(x, jnk) {
    if(is.rcf(x)) {
        SH <- specialHandling(x)
        if(is.null(SH)) return(NULL)
        return(SH$filename)
    }
    if(inherits(x, 'modelBaseClass')) {
        if(exists('.specialHandling', envir = x))
            return(x$.specialHandling)
        else
            return(NULL)
    }
    NULL
}

# for now export this as R<3.1.2 give warnings if don't

#' Class \code{nimbleFunctionBase}
#' @aliases nimbleFunctionBase
#' @export
#' @description
#' Classes used internally in NIMBLE and not expected to be called directly by users.
nimbleFunctionBase <- setRefClass(Class = 'nimbleFunctionBase', 
                                  fields = list(
                                      .generatorFunction = 'ANY',
                                      .CobjectInterface = 'ANY', 
                                      .newSetupLinesProcessed = 'ANY'
                                  ),
                                  methods = list(
                                      initialize = function(...)
                                          callSuper(...),
                                      getDefinition = function()
                                          nimble:::getDefinition(.self)
                                  ))	#	$runRelated


## template for the reference class internal to all nimble functions
nf_createRefClassDef <- function(setup, methodList, className = nf_refClassLabelMaker(), globalSetup, declaredSetupOutputNames,
                                 contains = NULL) {
    finalMethodList <- lapply(methodList, function(nfMethodRCobject) nfMethodRCobject$generateFunctionObject())
    finalMethodList[['show']] <- eval(substitute(function() writeLines(paste0('reference class object for nimble function class ', className)),
                                                 list(className = className)))
    if(!is.null(contains)) {
      finalMethodList <- c(finalMethodList, nf_getBaseClassMethods(methodList, contains))
    }
    substitute(
        setRefClass(Class   = NFREFCLASS_CLASSNAME,
                    fields  = NFREFCLASS_FIELDS,
                    methods = NFREFCLASS_METHODS, 
                    contains = 'nimbleFunctionBase',	#	$runRelated
                    where   = where),
        list(NFREFCLASS_CLASSNAME = className,
             NFREFCLASS_FIELDS    = nf_createRefClassDef_fields(setup, methodList, globalSetup, declaredSetupOutputNames),
             NFREFCLASS_METHODS   = finalMethodList
            )
        )
}

nf_getBaseClassMethods <- function(methodList, contains) {
  contains_env <- environment(contains)
  baseClassMethodNames <- names(contains_env$methods)
  ## Including run seems to make this complete, although in fact missing run args to
  ## nimbleFunction will result in an empty function, so it won't be missing.
  if(is.function(contains_env$run)) baseClassMethodNames <- c("run", baseClassMethodNames)
  for(mn in baseClassMethodNames) {
    reqd <- TRUE # reqd FALSE means the method might be taken from contains
    if(is.logical(contains_env$methodControl[[mn]][['required']]))
      reqd <- contains_env$methodControl[[mn]][['required']][1]
    provided <- mn %in% names(methodList)
    if(!reqd) {
      if(!provided)
        methodList[[mn]] <- contains_env$methods[[mn]]
    }
    if(reqd) {
      if(!provided)
        messageIfVerbose("  [Warning] method ", mn, " is required from the contains (base) class, but was not provided.")
    }
  }
  methodList
}

## creates a list of the fields (setupOutputs) for a nimble function reference class
nf_createRefClassDef_fields <- function(setup, methodList, globalSetup, declaredSetupOutputNames) {
    setupOutputNames <- nf_createSetupOutputNames(setup, methodList, declaredSetupOutputNames, globalSetup)
    if(FALSE) print(setupOutputNames)
    fields <- as.list(rep('ANY', length(setupOutputNames)))
    names(fields) <- setupOutputNames
    return(fields)
}

nf_createSetupOutputNames <- function(setup, methodList, declaredSetupOutputNames, globalSetup) {
    setupOutputNames <- character(0)
    setupOutputNames <- c(setupOutputNames, names(formals(setup)))   # add all setupArgs to potential setupOutputs
    setupOutputNames <- c(setupOutputNames, nf_assignmentLHSvars(body(setup)), if(!is.null(globalSetup)) nf_assignmentLHSvars(body(globalSetup)) else character())  # add all variables on LHS of <- in setup to potential setupOutputs
    setupOutputNames <- intersect(setupOutputNames, nf_createAllNamesFromMethodList(methodList))
    setupOutputNames <- c(setupOutputNames, declaredSetupOutputNames)
    setupOutputNames <- unique(setupOutputNames)
    return(setupOutputNames)
}

nf_assignmentLHSvars <- function(code) {
  if(!is.call(code))     return(character(0))
  isAssign <- code[[1]] == '<-' | code[[1]] == '='
  if(!isAssign)  return(unique(unlist(lapply(as.list(code), nf_assignmentLHSvars))))
  if(isAssign)  return(c(nf_getVarFromAssignmentLHScode(code[[2]]), nf_assignmentLHSvars(code[[3]])))
}

## determines the name of the target variable, from the LHS code of an `<-` assignment statement
nf_getVarFromAssignmentLHScode <- function(code) {
    if(is.name(code)) return(deparse(code))
    return(nf_getVarFromAssignmentLHScode(code[[2]]))
}

## creates a list of all the names of all variables and functions in the code of methodList functions
nf_createAllNamesFromMethodList <- function(methodList, onlyArgsAndReturn = F) {
    methodBodyListCode <- list()
    if(!onlyArgsAndReturn)
      methodBodyListCode <- lapply(methodList, function(f) f$code)
    methodReturnListCode <- lapply(methodList, function(f) f$returnType) ##might need changing
    methodArgListCode <- lapply(methodList, function(f) f$argInfo$argList[[1]])
    methodListCode <- c(methodBodyListCode, methodArgListCode, methodReturnListCode)
    if(length(methodListCode) > 0)
    return(unique(unlist(lapply(methodListCode, function(code) all.names(code)))))
}

nf_getNamesFromSetupOutputDeclaration <- function(setupOutputsDeclaration) {
    if(setupOutputsDeclaration[[1]] != 'setupOutputs') stop('something went wrong')
    return(unlist(lapply(setupOutputsDeclaration[-1], function(so) { if(is.call(so)) stop('cannot have a call inside setupOutputs() declaration') else deparse(so) } )))
}

## processing of all objects to become NF member data
## needs to be exported as otherwise use of nimble::: in `nf_createGeneratorFunctionDef()` gives R CMD check NOTE
#' @export
nf_preProcessMemberDataObject <- function(obj) {
    if(is(obj, 'CmodelBaseClass')) {
        warning('This nimbleFunction was passed a *compiled* model object.\nInstead, the corresponding *uncompiled* model object was used.', call. = FALSE)
        return(obj$Rmodel)
    }
    return(obj)
}

## definition for the nimble function generator (specializer)
nf_createGeneratorFunctionDef <- function(setup) {
  generatorFunctionDef <- substitute(
    function() {
      SETUPCODE                    # execute setupCode
      nfRefClassObject <- nfRefClass()   # create an object of the reference class
      nfRefClassObject$.generatorFunction <- generatorFunction   # link upwards to get the generating function of this nf
      ## assign setupOutputs into reference class object
      if(!getNimbleOption('compileOnly'))
        for(.var_unique_name_1415927 in .namesToCopyFromGlobalSetup)    { 
          nfRefClassObject[[.var_unique_name_1415927]] <- nf_preProcessMemberDataObject(get(.var_unique_name_1415927, envir = .globalSetupEnv)) 
        }
      for(.var_unique_name_1415927 in .namesToCopyFromSetup)    {
        nfRefClassObject[[.var_unique_name_1415927]] <- nf_preProcessMemberDataObject(get(.var_unique_name_1415927)) 
      }
      return(nfRefClassObject)
    },
    list(SETUPCODE = nf_processSetupFunctionBody(setup, returnCode = TRUE)))
  generatorFunctionDef[[4]] <- NULL
  return(generatorFunctionDef)
}

nf_processSetupFunctionBody <- function(setup, returnCode = FALSE, returnSetupOutputDeclaration = FALSE) {
    code <- body(setup)
    returnLineNum <- 0
    for(i in seq_along(code)) {
        if(is.call(code[[i]])) {
            if(is.name(code[[i]][[1]])) {
                if(code[[i]][[1]] == 'setupOutputs') {
                    returnLineNum <- i
                    break;
                }
            }
        }
    }
    if(sum(all.names(code) == 'setupOutputs') > 1) stop('multiple setupOutputs() declarations in nimbleFunction setup argument; only one allowed')
    if(returnLineNum == 0) {
        ## no setupOutputs() declaration found; default behavior
        setupOutputDeclaration <- quote(setupOutputs())
    } else {
        ## setupOutputs() declaration was found
        setupOutputDeclaration <- code[[returnLineNum]]
        code[returnLineNum] <- NULL
    }
    if('list' %in% all.names(setupOutputDeclaration))  stop('setupOutputs(...) declaration should not include \'list()\'')
    if(returnCode)                   return(code)
    if(returnSetupOutputDeclaration) return(setupOutputDeclaration)
    stop('must specify either returnCode=TRUE or returnSetupOutputDeclaration=TRUE')
}

## generates the argument list for the generator function
nf_createGeneratorFunctionArgs <- function(setup, pf) {
    generatorFunctionArgs <- lapply(formals(setup), function(arg) { if(is.blank(arg)) arg else eval(arg, pf) })
    return(generatorFunctionArgs)
}

## helper function for creating argument lists
nf_createAList <- function(argNames) {
    aListText <- if(length(argNames) > 0) {
        argNames <- paste0('`', argNames, '`')
        paste0('alist(', paste(argNames,'=',collapse=','), ')')
    } else { 'alist()' }
    eval(parse(text = aListText, keep.source = FALSE))
}

## returns names which don't begin with '.'
nf_namesNotHidden <- function(names) {
    names[!grepl('^\\.', names)]
}
