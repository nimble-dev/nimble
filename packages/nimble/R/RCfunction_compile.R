## Small class for information on compilation of each nf method

RCfunctionCompileClass <- setRefClass(
    'RCfunctionCompileClass',
    fields = list(
        origRcode = 'ANY',
        origLocalSymTab = 'ANY',
        nimExpr = 'ANY',
        newLocalSymTab = 'ANY',
        returnSymbol = 'ANY',
        newRcode = 'ANY',
        typeEnv = 'ANY'	#environment
    ),
    methods = list(
        initialize = function(...){
            typeEnv <<- new.env()
            callSuper(...)
        }
    ))

RCvirtualFunProcessing <- setRefClass(
    'RCvirtualFunProcessing',
    fields = list(
        name = 'ANY',		##character
        RCfun = 'ANY',          ##nfMethodRC
        nameSubList = 'ANY',
        compileInfo = 'ANY',    ## RCfunctionCompileClass``
        const = 'ANY',
        neededRCfuns = 'ANY',
        initialTypeInferenceDone = 'ANY'
    ),
    methods = list(
        initialize = function(f = NULL, funName, const = FALSE) {
            const <<- const
            if(!is.null(f)) {
                if(missing(funName)) {
                    sf <- substitute(f)
                    name <<- Rname2CppName(deparse(sf))
                } else {
                    name <<- funName
                }
                if(is.function(f)) {
                    RCfun <<- nfMethodRC$new(f)
                } else {
                    if(!inherits(f, 'nfMethodRC')) {
                        stop('Error: f must be a function or an RCfunClass')
                    }
                    RCfun <<- f
                }
                compileInfo <<- RCfunctionCompileClass$new(
                    origRcode = RCfun$code,
                    newRcode = RCfun$code)
            }
            neededRCfuns <<- list()
            initialTypeInferenceDone <<- FALSE
        },
        showCpp = function() {
            writeCode(
                nimGenerateCpp(compileInfo$nimExpr, compileInfo$newLocalSymTab)
            )
        },
        setupSymbolTables = function(parentST = NULL,
                                     neededTypes = NULL,
                                     nimbleProject = NULL) {
            argInfoWithMangledNames <- RCfun$argInfo
            numArgs <- length(argInfoWithMangledNames)
            if(numArgs > 0) {
                argIsBlank <- unlist(
                    lapply(RCfun$argInfo,
                           identical,
                           formals(function(a) {})[[1]])
                )
                ## N.B: It seems to be impossible to store the value of a blank argument, formals(function(a) {})[[1]], in a variable
                if(any(argIsBlank)) {
                    stop(paste0("Type declaration missing for argument(s) ",
                                paste(names(RCfun$argInfo)[argIsBlank],
                                      collapse = ", ")),
                         call. = FALSE)
                }
            }
            argInfoWithMangledNames <- RCfun$argInfo
            numArgs <- length(argInfoWithMangledNames)
            if(numArgs>0)
                names(argInfoWithMangledNames) <-
                    paste0("ARG",
                           1:numArgs,
                           "_",
                           Rname2CppName(names(argInfoWithMangledNames)),"_")
            nameSubList <<- lapply(names(argInfoWithMangledNames),
                                   as.name)
            names(nameSubList) <<- names(RCfun$argInfo)
            
            ## This will only handle basic types.
            ## nimbleLists will be added below
            compileInfo$origLocalSymTab <<-
                argTypeList2symbolTable(argInfoWithMangledNames,
                                        neededTypes,
                                        names(RCfun$argInfo)) 
            compileInfo$newLocalSymTab <<-
                argTypeList2symbolTable(argInfoWithMangledNames,
                                        neededTypes,
                                        names(RCfun$argInfo))

            ## This modifies the symTab and returns types needed for
            ## input or return.

            ## Currently the only non-basic type resolved in this step
            ## would be nimbleLists.
            if(is.null(nimbleProject))
                nimbleProject <- get('nimbleProject', envir = RCfun)
            neededRCfuns <<-
                resolveUnknownTypes(compileInfo$origLocalSymTab,
                                    neededTypes,
                                    nimbleProject)
            resolveUnknownTypes(compileInfo$newLocalSymTab,
                                c(neededRCfuns, neededTypes),
                                nimbleProject)
            
            if(!is.null(parentST)) {
                compileInfo$origLocalSymTab$setParentST(parentST)
                compileInfo$newLocalSymTab$setParentST(parentST)
            }
            compileInfo$returnSymbol <<-
                argType2symbol(RCfun$returnType,
                               neededTypes,
                               "return",
                               "returnType")
            updatedReturn <-
                resolveOneUnknownType(compileInfo$returnSymbol,
                                      c(neededRCfuns, neededTypes),
                                      nimbleProject)
            compileInfo$returnSymbol <<- updatedReturn[[1]]
            neededRCfuns <<- c(neededRCfuns, updatedReturn[[2]])
        },
        process = function(...) {
            if(inherits(compileInfo$origLocalSymTab, 'uninitializedField')) {
                setupSymbolTables()
            }
        },
        printCode = function() {
            writeCode(nimDeparse(compileInfo$nimExpr))
        }
    )
)


RCfunction <- function(f, name = NA, returnCallable = TRUE, check, enableDerivs = FALSE, where = NULL) {
    if(is.na(name))
        name <- rcFunLabelMaker(envName = environmentName(where))
    nfm <- nfMethodRC$new(f, name, check = check, enableDerivs = enableDerivs)
    if(returnCallable)
        nfm$generateFunctionObject(keep.nfMethodRC = TRUE, where = where)
    else
        nfm
}

is.rcf <- function(x, inputIsName = FALSE) {
    if(inputIsName)
        x <- get(x)
    if(inherits(x, 'nfMethodRC'))
        return(TRUE)
    if(is.function(x)) {
        if(is.null(environment(x)))
            return(FALSE)
        if(exists('nfMethodRCobject', envir = environment(x), inherits = FALSE))
            return(TRUE)
    }
    FALSE
}

rcFunLabelMaker <- labelFunctionCreator('rcFun')

nf_substituteExceptFunctionsAndDollarSigns <- function(code, subList) {
    ## this is almost like doing substitute with code and subList, but it doesn't traverse the RHS of a '$' operator
    ## and it doesn't replace and function names
    if(is.character(code)) return(code)
    if(is.numeric(code)) return(code)
    if(is.logical(code)) return(code)
    if(is.null(code)) return(code)
    if(is.name(code)) {
        maybeAns <- subList[[as.character(code)]]
        return(
            if(is.null(maybeAns))
                code
            else
                maybeAns
        )
    }
    if(is.call(code)) {
        if(length(code) == 1)
            return(code)
        else {
            if(is.call(code[[1]]))
                indexRange <- 1:length(code)
            else {
                if(as.character(code[[1]])=='$')
                    indexRange <- 2
                else 
                    indexRange <- 2:length(code)
            }
        }
        for(i in indexRange)
            code[[i]] <-
                nf_substituteExceptFunctionsAndDollarSigns(code[[i]], subList)
        return(code)
    }
    if(is.list(code)) {
        ## Keyword processing stage of compilation may have stuck
        ## lists into the argument list of a call (for maps)
        code <- lapply(code,
                       nf_substituteExceptFunctionsAndDollarSigns,
                       subList)
        return(code)
    }
    stop(paste("Error doing replacement for code ", deparse(code)))
}

RCfunProcessing <- setRefClass(
    'RCfunProcessing',
    contains = 'RCvirtualFunProcessing',
    fields = list(),
    methods = list(
        process = function(debug = FALSE,
                           debugCpp = FALSE,
                           debugCppLabel = character(),
                           doKeywords = TRUE,
                           nimbleProject = NULL,
                           initialTypeInferenceOnly = FALSE) {
            if(!is.null(nimbleOptions()$debugRCfunProcessing)) {
                if(nimbleOptions()$debugRCfunProcessing) {
                    debug <- TRUE
                    writeLines('Debugging RCfunProcessing (nimbleOptions()$debugRCfunProcessing is set to TRUE)') 
                }
            }
          
            if(debug) {
                writeLines('**** READY to start RCfunProcessing::process *****')
                browser()
            }
            
            if(is.null(nimbleProject))
                nimbleProject <- get('nimbleProject', envir = RCfun)

            if(!initialTypeInferenceDone) {
                
                if(!is.null(nimbleOptions()$debugCppLineByLine)) {
                    if(nimbleOptions()$debugCppLineByLine) {
                        debugCpp <- TRUE
                        if(length(debugCppLabel) == 0) debugCppLabel <- name
                    }
                }
                
                if(doKeywords) {
                    matchKeywords()
                    processKeywords()
                }
                
                if(inherits(compileInfo$origLocalSymTab,
                            'uninitializedField')) {

                    setupSymbolTables()
                }
                initialTypeInferenceDone <<- TRUE
            }

            if(initialTypeInferenceOnly)
                return(NULL);
            
            if(debug) {
                writeLines('**** READY FOR makeExprClassObjects *****')
                browser()
            }
            if(length(nameSubList) > 0)
                compileInfo$newRcode <<-
                    nf_substituteExceptFunctionsAndDollarSigns(compileInfo$newRcode,
                                                               nameSubList)
            ## set up exprClass object
            compileInfo$nimExpr <<- RparseTree2ExprClasses(compileInfo$newRcode)
            
            if(debug) {
                print('nimDeparse(compileInfo$nimExpr)')
                writeCode(nimDeparse(compileInfo$nimExpr))
                writeLines('***** READY FOR processSpecificCalls *****')
                browser()
            }
            
            exprClasses_processSpecificCalls(compileInfo$nimExpr,
                                             compileInfo$newLocalSymTab)

            if(debug) {
                print('nimDeparse(compileInfo$nimExpr)')
                writeCode(nimDeparse(compileInfo$nimExpr))
                writeLines('***** READY FOR buildInterms *****')
                browser()
            }
            
            ## build intermediate variables
            exprClasses_buildInterms(compileInfo$nimExpr)
            
            if(debug) {
                print('nimDeparse(compileInfo$nimExpr)')
                writeCode(nimDeparse(compileInfo$nimExpr))
                print('compileInfo$newLocalSymTab')
                print(compileInfo$newLocalSymTab)
                writeLines('***** READY FOR initSizes *****')
                browser()
            }
            compileInfo$typeEnv <<-
                exprClasses_initSizes(compileInfo$nimExpr,
                                      compileInfo$newLocalSymTab,
                                      returnSymbol = compileInfo$returnSymbol)
            if(debug) {
                print('ls(compileInfo$typeEnv)')
                print(ls(compileInfo$typeEnv))
                print('lapply(compileInfo$typeEnv, function(x) x$show())')
                lapply(compileInfo$typeEnv, function(x) x$show())
                writeLines('***** READY FOR setSizes *****')
                browser()
            }
            
            compileInfo$typeEnv[['neededRCfuns']] <<- list()
            compileInfo$typeEnv[['.AllowUnknowns']] <<- TRUE ## will be FALSE for RHS recursion in setSizes
            compileInfo$typeEnv[['.ensureNimbleBlocks']] <<- FALSE ## will be TRUE for LHS recursion after RHS sees rmnorm and other vector dist "r" calls.
            compileInfo$typeEnv[['.allowFunctionAsArgument']] <<- FALSE ## will be TRUE when recursing on optim. See sizeOptim.
            compileInfo$typeEnv[['.nimbleProject']] <<- nimbleProject

            passedArgNames <-
                as.list(compileInfo$origLocalSymTab$getSymbolNames()) 
            names(passedArgNames) <-
                compileInfo$origLocalSymTab$getSymbolNames() 
            compileInfo$typeEnv[['passedArgumentNames']] <<- passedArgNames ## only the names are used.
            compileInfo$typeEnv[['nameSubList']] <<- nameSubList

            ## exprClasses_setSizes contains lots of reticulated error trapping.
            ## If options('error') is not NULL (typically options(error = recover)), let that take precedent for error trapping of exprClasses_setSizes.
            ## Otherwise, wrap the setSizes call in our own error-trapping that outputs: any local message from a stop() inside a size handler,
            ##    a message from here showing the larger block of code in which the error occurred, and, if nimbleOptions('verboseErrors') is TRUE, the call stack. 
            if(!is.null(options('error')))
                exprClasses_setSizes(compileInfo$nimExpr,
                                     compileInfo$newLocalSymTab,
                                     compileInfo$typeEnv)
            else tryCatch(
                     withCallingHandlers(
                         exprClasses_setSizes(compileInfo$nimExpr,
                                              compileInfo$newLocalSymTab,
                                              compileInfo$typeEnv),
                         error = function(e) {
                             stack <- sapply(sys.calls(), deparse)
                             .GlobalEnv$.nimble.traceback <-
                                 capture.output(traceback(stack))
                         }),
                     error = function(e) {
                         eMessage <- if(!is.null(e$message))
                                         paste0(as.character(e$message),"\n")
                                     else ""
                         message <-
                             paste(eMessage,
                                   'There was some problem in the the setSizes processing step for this code:',
                                   paste(deparse(compileInfo$origRcode), collapse = '\n'))
                         if(getNimbleOption('verboseErrors')) {
                             message <- paste(message,
                                              'Internal Error:',
                                              e,
                                              'Traceback:',
                                              paste0(.GlobalEnv$.nimble.traceback, collapse = '\n'),
                                              sep = '\n')
                         }
                         stop(message, call. = FALSE)
                     })
            neededRCfuns <<- c(neededRCfuns,
                               compileInfo$typeEnv[['neededRCfuns']])
            if(debug) {
                print('compileInfo$nimExpr$print(showType = TRUE) -- broken')
                print('compileInfo$nimExpr$print(showAssertions = TRUE) -- possible broken')
                writeLines('***** READY FOR insertAssertions *****')
                browser()
            }
            
            if(nimbleOptions('experimentalNewSizeProcessing')) {
                exprClasses_setToEigenize(compileInfo$nimExpr,
                                          compileInfo$newLocalSymTab,
                                          compileInfo$typeEnv)
            }
            
            tryResult <- try(exprClasses_insertAssertions(compileInfo$nimExpr))
            if(inherits(tryResult, 'try-error')) {
                stop(
                    paste('There is some problem at the insertAdditions processing step for this code:\n',
                          paste(deparse(compileInfo$origRcode),
                                collapse = '\n'),
                          collapse = '\n'),
                    call. = FALSE)
            }
            
            if(debug) {
                print('compileInfo$nimExpr$print(showAssertions = TRUE)')
                compileInfo$nimExpr$print(showAssertions = TRUE)
                print('compileInfo$nimExpr$print(showToEigenize = TRUE)')
                compileInfo$nimExpr$print(showToEigenize = TRUE)
                print('nimDeparse(compileInfo$nimExpr)')
                writeCode(nimDeparse(compileInfo$nimExpr))
                writeLines('***** READY FOR labelForEigenization *****')
                browser()
            }
            
            tryResult <- try(
                exprClasses_labelForEigenization(compileInfo$nimExpr)
            )
            if(inherits(tryResult, 'try-error')) {
                stop(
                    paste('There is some problem at the Eigen labeling processing step for this code:\n',
                          paste(deparse(compileInfo$origRcode),
                                collapse = '\n'),
                          collapse = '\n'),
                    call. = FALSE)
            }
            if(debug) {
                print('nimDeparse(compileInfo$nimExpr)')
                writeCode(nimDeparse(compileInfo$nimExpr))
                writeLines('***** READY FOR eigenize *****')
                browser()
            }
            tryResult <- try(
                exprClasses_eigenize(compileInfo$nimExpr,
                                     compileInfo$newLocalSymTab,
                                     compileInfo$typeEnv))
            if(inherits(tryResult, 'try-error')) {
                stop(
                    paste('There is some problem at the Eigen processing step for this code:\n',
                          paste(deparse(compileInfo$origRcode),
                                collapse = '\n'),
                          collapse = '\n'),
                    call. = FALSE)
            }
            if(debug) {
                print('nimDeparse(compileInfo$nimExpr)')
                writeCode(nimDeparse(compileInfo$nimExpr))
                print('compileInfo$newLocalSymTab')
                print(compileInfo$newLocalSymTab)
                print('ls(compileInfo$typeEnv)')
                print(ls(compileInfo$typeEnv))
                
                writeLines('***** READY FOR liftMaps*****')
                browser()
            }
            
            exprClasses_liftMaps(compileInfo$nimExpr,
                                 compileInfo$newLocalSymTab,
                                 compileInfo$typeEnv)
            if(debug) {
                print('nimDeparse(compileInfo$nimExpr)')
                writeCode(nimDeparse(compileInfo$nimExpr))
                writeLines('***** READY FOR cppOutput*****')
                browser()
            }
            
            if(debugCpp) {
                if(debug) writeLines('*** Inserting debugging')
                exprClasses_addDebugMarks(compileInfo$nimExpr,
                                          paste(debugCppLabel, name))
                if(debug) {
                    print('nimDeparse(compileInfo$nimExpr)')
                    writeCode(nimDeparse(compileInfo$nimExpr))
                    writeLines('***** READY FOR cppOutput*****')
                    browser()
                }
            }
            if(debug & debugCpp) {
                print('writeCode(nimGenerateCpp(compileInfo$nimExpr, newMethods$run$newLocalSymTab))')
                writeCode(
                    nimGenerateCpp(compileInfo$nimExpr,
                                   compileInfo$newLocalSymTab)
                )
            }
        },
        processKeywords = function(nfProc = NULL) {
            compileInfo$newRcode <<-
                processKeywords_recurse(compileInfo$origRcode, nfProc)
        },
        matchKeywords = function(nfProc = NULL) {
            compileInfo$origRcode <<-
                matchKeywords_recurse(compileInfo$origRcode, nfProc) ## nfProc needed for member functions of nf objects
        }
    )
)
