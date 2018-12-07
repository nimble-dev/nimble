## A small class for information deduced about a variable in a BUGS model
varInfoClass <- setRefClass('varInfoClass',
                            fields = list(
                                varName = 'ANY',
                                mins = 'ANY',
                                maxs = 'ANY',
                                nDim = 'ANY',
                                anyStoch = 'ANY',
                                anyDynamicallyIndexed = 'ANY'))

## A small class for information about a node in the igraph
graphNode <- setRefClass(
    Class = 'graphNode',
    fields = list(
        nodeName         = 'ANY',
        graphID          = 'ANY',
        type             = 'ANY',
        originNodeName   = 'ANY',
        nodeFunctionName = 'ANY'
    )
)

#' Class for NIMBLE model definition
#'
#' Class for NIMBLE model definition that is not usually needed directly by a user.
#'
#' @details See \code{?modelBaseClass} for information about creating NIMBLE BUGS models.
modelDefClass <- setRefClass('modelDefClass',
                             fields = list(
                                 ## set in the call modelDefClass$new(name)
                                 name = 'ANY',
                                 
                                 ## the following are all set in setupModel()
                                 BUGScode = 'ANY',  ## original BUGS code, set in assignBUGScode()
                                 constantsEnv = 'ANY', ## environment with constants, set in assignConstants()
                                 constantsList = 'ANY',  ## named list with constants, set in assignConstants()
                                 constantsNamesList = 'ANY', ## list of constants name objects, set in assignConstants()
                                 constantsScalarNamesList = 'ANY', ## could eventually replace constantsNamesList. added for newNodeFxns
                                 dimensionsList = 'ANY',		#list		   ## list of provided dimension information, set in assignDimensions()
                                 contexts = 'ANY',				#list 			 ## list of BUGScontextClass objects
                                 declInfo = 'ANY',				#list				 ## list of BUGSdeclInfo objects
                                 varInfo = 'ANY',			#list	  ## list of varInfoClass objects, set in genVarInfo()
                                 logProbVarInfo = 'ANY',	#list	  ## list of varInfoClass objects, set in genLogProbVarInfo()
                                 isDataVarInfo = 'ANY', 	#list		## list of varInfoClass objects, set in genIsDataVarInfo()
                                 varNames = 'ANY',  ## vector of all model variable names, set in genVarNames()
                                 unknownIndexNames = 'ANY', ## vector of unknown index variable names
                                 symTab = 'ANY',  ## symbolTable object, set in buildSymbolTable()
                                 graph = 'ANY',     ## igraph object, set in buildIgraph()
                                 graphNodesList = 'ANY',   ## list of graphNode objects, set in genGraphNodesList()
                                 maps = 'ANY',   ## object of mapsClass, set in buildMaps()
                                 numNodeFunctions = 'ANY',  ## FIXME: obsolete as only used in buildNodeFunctions_old()
                                 
                                 modelClass = 'ANY',   ## custom model class
                                 modelValuesClassName = 'ANY',    ## set in setModelValuesClassName()
                                 modelValuesClass = 'ANY', ## custom model values class
                                 classEnvironment = 'ANY'	#environment		 # environment in which the reference classes will be defined.
                             ),
                             
                             methods = list(
                                 initialize = function(...){
                                 	name <<- character()
                                 	dimensionsList <<- list() 	
                                 	contexts <<- list()
                                 	declInfo <<- list()
                                 	varInfo <<- list()
                                 	logProbVarInfo <<- list()
                                 	isDataVarInfo <<- list()
                                 	classEnvironment <<- new.env()
                                 	callSuper(...)
                                 },
                                 setupModel = function(code, constants, dimensions, inits, data, debug) {},
                                 
                                 ## the following are all run, in this order, by setupModel():
                                 setModelValuesClassName        = function() {},
                                 assignBUGScode                 = function() {},
                                 assignConstants                = function() {},
                                 assignDimensions               = function() {},
                                 initializeContexts             = function() {},
                                 processBUGScode                = function() {},
                                 splitConstantsAndData          = function() {},
                                 addMissingIndexing             = function() {},
                                 processBoundsAndTruncation     = function() {},
                                 expandDistributions            = function() {},
                                 checkMultivarExpr              = function() {},
                                 processLinks                   = function() {},
                                 reparameterizeDists            = function() {},
                                 replaceAllConstants            = function() {},
                                 liftExpressionArgs             = function() {},
                                 addRemainingDotParams          = function() {},
                                 addIndexVarsToDeclInfo         = function() {},
                                 genSymbolicParentNodes         = function() {},
                                 genUnknownIndexDeclarations    = function() {},
                                 genReplacementsAndCodeReplaced = function() {},
                                 genAltParamsModifyCodeReplaced = function() {},
                                 genBounds                      = function() {},
                                 genReplacedTargetValueAndParentInfo = function() {},
                                 removeEmptyBUGSdeclarations    = function() {},
                                 genIsDataVarInfo               = function() {},
                                 genVarNames                    = function() {},
                                 buildSymbolTable               = function() {},
                                 genGraphNodesList              = function() {},
                                                                  
                                 newModel                       = function() {},
                                 fixRStudioHanging              = function() {},
                                 printDI                        = function() {},
                                 
                                 genNodeInfo3                   = function() {},
                                 genVarInfo3                    = function() {},
                                 # put check of var dims here?
                                 addUnknownIndexVars            = function() {},
                                 findDynamicIndexParticipants   = function() {},
                                 addFullDimExtentToUnknownIndexDeclarations = function() {},
                                 genExpandedNodeAndParentNames3 = function() {},
                                 stripUnknownIndexInfo          = function() {},
                                 #These functions are NOT run inside of setupModel
                                 nodeName2GraphIDs = function(){},
                                 graphIDs2indexedNodeInfo = function(){},
                                 nodeName2LogProbName = function(){}
                             ))



## This is the master entry function
##
## NOTES:
## (1) for the moment, we bail on 'lifting' any expression with vectorized indexing, e.g. x[1:10], UNLESS it's a call to chol(x[..., ...]).
##     this occurs in isExprLiftable().
## (2) in generating replacement code, we do *not* replace ':' expressions.
##     this takes place due to a single line, near the end of genReplacementsAndCodeRecurse() in nimbleBUGS_class_BUGSdeclClass.R
##     further, nameMashupFromExpr(expr) in nimbleBUGS_utils.R throws an error if expr contains a ':'
##
modelDefClass$methods(setupModel = function(code, constants, dimensions, inits, data, userEnv, debug = FALSE) {
    if(debug) browser()
    code <- codeProcessIfThenElse(code, constants, userEnv) ## evaluate definition-time if-then-else
    if(nimbleOptions("enableModelMacros")) code <- codeProcessModelMacros(code)
    setModelValuesClassName()         ## uses 'name' field to set field: modelValuesClassName
    assignBUGScode(code)              ## uses 'code' argument, assigns field: BUGScode.  puts codes through nf_changeNimKeywords
    assignConstants(constants)        ## uses 'constants' argument, sets fields: constantsEnv, constantsList, constantsNamesList
    assignDimensions(dimensions, inits, data)      ## uses 'dimensions' argument, sets field: dimensionList
    initializeContexts()              ## initializes the field: contexts
    processBUGScode(userEnv = userEnv)                 ## uses BUGScode, sets fields: contexts, declInfo$code, declInfo$contextID
    if(nimbleOptions("stop_after_processing_model_code")) {
        print(code)
        stop(paste0('Stopped after processing model code because\n',
                    'nimbleOptions("stop_after_processing_model_code") is TRUE\n'),
             call.=FALSE)
    }

    ## We will try to infer sizes later
    ##addMissingIndexing()              ## overwrites declInfo, using dimensionsList, fills in any missing indexing
    splitConstantsAndData()           ## deals with case when data is passed in as constants
    addMissingIndexing()              ## overwrites declInfo, using dimensionsList, fills in any missing indexing
    processBoundsAndTruncation()      ## puts bound expressions into declInfo, including transforming T(ddist(),lower,upper); need to do this before expandDistributions(), which is not set up to handle T() wrapping; need to save bound info for later use in reparameterizeDists() -- hence temporarily stored in boundExprs (can't put in code because it would be stripped out in expandDistributions, though alternative is to modify expandDistributions to add lower,upper back into code)
    expandDistributions()             ## overwrites declInfo for stochastic nodes: calls match.call() on RHS      (uses distributions$matchCallEnv)
    if(getNimbleOption('disallow_multivariate_argument_expressions'))
        checkMultivarExpr()               ## checks that multivariate params are not expressions
    processLinks()                    ## overwrites declInfo (*and adds*) for nodes with link functions           (uses linkInverses)
    reparameterizeDists()             ## overwrites declInfo when distribution reparameterization is needed       (uses distributions), keeps track of orig parameter in .paramName; also processes bound info to evaluate in context of model
    replaceAllConstants()
    liftExpressionArgs()              ## overwrites declInfo (*and adds*), lifts expressions in distribution arguments to new nodes.  does NOT lift '.param' names or 'lower' or 'upper'
    addRemainingDotParams()           ## overwrites declInfo, adds any additional .paramNames which aren't there  (uses distributions)
    replaceAllConstants()             ## overwrites declInfo with constants replaced; only replaces scalar constants
    addIndexVarsToDeclInfo()          ## sets field declInfo[[i]]$indexVariableExprs from contexts.  must be after overwrites of declInfo
    genSymbolicParentNodes()          ## sets field declInfo[[i]]$symbolicParentNodes. must be after overwrites of declInfo
    genUnknownIndexDeclarations()     ## add 'lifted' declarations for unknownIndex entities; needs symbolicParentNodes to exist in order to know what declarations have unknown indices
    genReplacementsAndCodeReplaced()  ## sets fields: declInfo[[i]]$replacements, $codeReplaced, $replacementNameExprs, $logProbNodeExpr
    genAltParamsModifyCodeReplaced()  ## sets field declInfo[[i]]$altParams, and modifies $codeReplaced to not include .param arguments (if stochastic)
    genBounds()                       ## sets field declInfo[[i]]$boundExprs, and (if not truncated) modifies $codeReplaced to omit lower and upper arguments (if stochastic)
    ## From here down is the new "version 3" processing
    genReplacedTargetValueAndParentInfo() ## In each declInfo[[i]], symbolicParentNodesReplaced, rhsVars, targetIndexNamePieces, and parentIndexNamePieces set
    genNodeInfo3(debug = debug)           ## In each contexts[[i]], replacementsEnv set. In each declInfo[[i]], replacementsEnv, unrolledIndicesMatrix, outputSize, and numUnrolledNodes set
    genVarInfo3()                         ## Sets varInfo[[nodeNames]] and logProbVarInfo[[nodeNames]] with varInfoClass objects (varName mins, maxs, nDim, anyStoch)
    addUnknownIndexVars(debug = debug)    ## adds elements to varInfo for the unknownIndex entities and creates unknownIndexNames to store names of these lifted variables; needs to occur after genVarInfo3 as it uses the varInfo of relevant parent variables
    findDynamicIndexParticipants() ## strip out USED_IN_INDEX() wrapping; if we move to dynamically updating the graph, this will be augmented to find variable elements involved in dynamic indexing (this may need to be split in two pieces as we may need this info at the vertex level, in which case some processing is needed after genExpandedNodeAndParentNames3)
    addFullDimExtentToUnknownIndexDeclarations() ## update unknownIndex declarations with full extent of relevant parent variable; this splits parent variable and does edge determination (from parent variable to unknownIndex variable) but without splitting based on unknown indices
    genExpandedNodeAndParentNames3(debug = debug) ## heavy processing: all graphIDs, maps, graph, nodeNames etc. built here
    stripUnknownIndexInfo()               ## removes unknownIndex declarations and vars
    checkForSelfParents()                 ## Checks to see if any nodes are their own parents, and errors out if so
    maps$setPositions3()                  ## Determine top, latent and end nodes
    buildSymbolTable()                    ## 
    genIsDataVarInfo()                    ## only the maxs is ever used, in newModel
    genVarNames()                         ## sets varNames <<- c(names(varInfo), names(logProbVarInfo))
    return(NULL)        
})

codeProcessIfThenElse <- function(code, constants, envir = parent.frame()) {
    codeLength <- length(code)
    if(code[[1]] == '{') {
        if(codeLength > 1) for(i in 2:codeLength) code[[i]] <- codeProcessIfThenElse(code[[i]], constants, envir)
        return(code)
    } 
    if(code[[1]] == 'for') {
        code[[4]] <- codeProcessIfThenElse(code[[4]], constants, envir)
        return(code)
    }
    if(codeLength > 1)
        if(code[[1]] == 'if') {
            constantsEnv <- as.environment(constants)
            parent.env(constantsEnv) <- envir
            evaluatedCondition <- try(eval(code[[2]], constantsEnv))
            if(inherits(evaluatedCondition, "try-error")) 
                stop("Cannot evaluate condition of 'if' statement: ", deparse(code[[2]]), ".\nCondition must be able to be evaluated based on values in 'constants'.")
            if(evaluatedCondition) return(codeProcessIfThenElse(code[[3]], constants, envir))
            else {
                if(length(code) == 4) return(codeProcessIfThenElse(code[[4]], constants, envir))
                else return(NULL)
            }
        } else
            return(code)
    else
        return(code)
}

## This function recurses through a block of code and expands any submodels
codeProcessModelMacros <- function(code,
                                   recursionLabels = character()) {
    expandRecursionLabels <- function(possibleMacroName,
                                      labels = character()) {
        paste0(possibleMacroName,
               if(length(labels) > 0)
                   paste0('(expanded from ',
                          paste(labels, collapse = '-->'),
                          ')')                          
               else
                   character()
               )
    }
    codeLength <- length(code)
    ## First check if this is the start of a curly-bracketed block
    if(code[[1]] == '{') {
        if(codeLength > 1)
            ## Recurse on each line
            for(i in 2:codeLength)
                code[[i]] <- codeProcessModelMacros(code[[i]])
        return(code)
    }
    ## If this is a for loop, recurse on the body of the loop
    if(code[[1]] == 'for') {
        code[[4]] <- codeProcessModelMacros(code[[4]])
        return(code)
    }
    ## Check if this line invokes a submodel.
    ## This can be done in two ways:
    ## (i) node1 [<- | ~] <macro name>(...)
    ## or
    ## (ii) <macro name>(args)
    ##
    ## The first version is more BUGS-like.
    ## The second version allows more full control.

    ## Initialize possibleMacroName assuming version (ii):
    possibleMacroName <- deparse(code[[1]])
    ## If it is really version (i), possibleMacroName will be
    ## ~ or <- and should be updated to the call on the right-hand side:
    if(possibleMacroName %in% c('<-', '~')) {
        possibleMacroName <- deparse(code[[3]][[1]])
    }
    if(exists(possibleMacroName)) { ## may need to provide an envir argument
        possibleMacro <- get(possibleMacroName) ## ditto
        if(inherits(possibleMacro, "model_macro")) {
            expandedInfo <- try(possibleMacro$process(code))
            if(inherits(expandedInfo, 'try-error'))
                stop(paste0("Model macro ",
                            expandRecursionLabels(
                                possibleMacroName,
                                recursionLabels
                            ),
                            " failed."),
                     call. = FALSE)
            if(!is.list(expandedInfo))
                stop(paste0("Model macro ",
                            expandRecursionLabels(
                                possibleMacroName,
                                recursionLabels
                            ),
                            " should return a list with an element named ",
                            "'code'.  It did not return a list."),
                     call. = FALSE)
            if(!is.call(expandedInfo[['code']]))
                stop(paste0("Model macro ",
                            expandRecursionLabels(
                                possibleMacroName,
                                recursionLabels
                            ),
                            " should return a list with an element named ",
                            "'code' that is a call."),
                     call. = FALSE)
            ## Return object is a list so we can ossibly extract other
            ## content in the future.  We recurse on the returned code
            ## to expand macros that it might contain.
            code <- codeProcessModelMacros(expandedInfo$code,
                                           c(recursionLabels, possibleMacroName)
                                           )
        }
    }
    code
}

modelDefClass$methods(setModelValuesClassName = function() {
    modelValuesClassName <<- paste0(Rname2CppName(name), '_MV_', nimbleUniqueID())
})
modelDefClass$methods(assignBUGScode = function(code) {
    ## uses 'code' argument, assigns field: BUGScode
    BUGScode <<- nf_changeNimKeywords(code)
})
modelDefClass$methods(assignConstants = function(constants) {
    ## uses 'constants' argument, sets fields: constantsEnv, constantsList, constantsNamesList
    constantsEnv <<- new.env()
    if(length(constants) > 0) {
        if(!is.list(constants) || is.null(names(constants)))   stop('constants argument must be a named list')
        list2env(constants, constantsEnv)
        constantsList <<- constants
        constantsNamesList <<- lapply(ls(constants), as.name)
        constantLengths <- unlist(lapply(constants, length))
        if(any(constantLengths > 1)) {
            iLong <- which(constantLengths > 1)
            ## message(paste0('Constant(s) ', paste0(names(constants)[iLong], sep=" ", collapse = " "), ' are non-scalar and may be handled as data if necessary.'))
            ## note some of the processing behind this message occurs in BUGSmodel between making the model def and the model
            constantsScalarNamesList <<- constantsNamesList[-iLong]
        } else
            constantsScalarNamesList <<- constantsNamesList 
    } else {
        constantsList <<- list()
        names(constantsList) <<- character(0)
        constantsNamesList <<- list()
        constantsScalarNamesList <<- list()
    }
})
modelDefClass$methods(assignDimensions = function(dimensions, initsList, dataList) {
    ## uses 'dimensions' argument, sets field: dimensionList
    
    # first, add the provided dimensions
    dL <- dimensions
    if(is.null(dL))
        dL <- list()
    
    # add dimensions of any *non-scalar* constants to dimensionsList
    # we'll try to be smart about this: check for duplicate names in constants and dimensions, and make sure they agree
    for(i in seq_along(constantsList)) {
        constName <- names(constantsList)[i]
        ## constDim <- if(is.null(dim(constantsList[[i]]))) length(constantsList[[i]]) else dim(constantsList[[i]])
#       constDim <- nimbleInternalFunctions$dimOrLength(constantsList[[i]], scalarize = TRUE)
        constDim <- nimbleInternalFunctions$dimOrLength(constantsList[[i]], scalarize = FALSE)  # don't scalarize as want to preserve dims as provided by user, e.g. for 1x1 matrices
        if(length(constDim) == 1 && constDim == 1)
            constDim <- numeric(0)  # but for 1-length vectors treat as scalars as that is how handled in system
        if(constName %in% names(dL)) {
            if(!identical(as.numeric(dL[[constName]]), as.numeric(constDim))) {
                stop('inconsistent dimensions between constants and dimensions arguments: ', constName)
            }
        } else {
            dL[[constName]] <- constDim
        }
    }

    # add dimensions of any *non-scalar* inits to dimensionsList
    # we'll try to be smart about this: check for duplicate names in inits and dimensions, and make sure they agree
    for(i in seq_along(initsList)) {
        initName <- names(initsList)[i]
        initDim <- nimbleInternalFunctions$dimOrLength(initsList[[i]], scalarize = FALSE)  # don't scalarize as want to preserve dims as provided by user, e.g. for 1x1 matrices
        if(!(length(initDim) == 1 && initDim == 1)) {  # i.e., non-scalar inits; 1-length vectors treated as scalars and not passed along as dimension info to avoid conflicts between scalars and one-length vectors/matrices/arrays in various places
            if(initName %in% names(dL)) {
                if(!identical(as.numeric(dL[[initName]]), as.numeric(initDim))) {
                    warning('inconsistent dimensions between inits and dimensions arguments: ', initName, '; ignoring dimensions in inits.')
                }
            } else {
                dL[[initName]] <- initDim
            }
        }
    }

    # add dimensions of any *non-scalar* data to dimensionsList
    # we'll try to be smart about this: check for duplicate names in data and dimensions, and make sure they agree
    # main use case here is when user provides RHS only variable as data
    for(i in seq_along(dataList)) {
        dataName <- names(dataList)[i]
        dataDim <- nimbleInternalFunctions$dimOrLength(dataList[[i]], scalarize = FALSE)  # don't scalarize as want to preserve dims as provided by user, e.g. for 1x1 matrices
        if(!(length(dataDim) == 1 && dataDim == 1)) {  # i.e., non-scalar data; 1-length vectors treated as scalars and not passed along as dimension info to avoid conflicts between scalars and one-length vectors/matrices/arrays in various places
            if(dataName %in% names(dL)) {
                if(!identical(as.numeric(dL[[dataName]]), as.numeric(dataDim))) {
                    warning('inconsistent dimensions between data and dimensions arguments: ', dataName, '; ignoring dimensions in data.')
                }
            } else {
                dL[[dataName]] <- dataDim
            }
        }
    }
    
    dimensionsList <<- dL
})

modelDefClass$methods(initializeContexts = function() {
    ## initializes the field: contexts
    ## there is always a context #1 that is the empty context. this sets it up.
    BUGScontextClassObject <- BUGScontextClass$new()
    BUGScontextClassObject$setup(singleContexts = list())
    contexts[[1]] <<- BUGScontextClassObject
})

reprioritizeColonOperator <- function(code) {
    split.code <- strsplit(deparse(code), ":")
    if(length(split.code[[1]]) == 2) return(parse(text = paste0("(", split.code[[1]][1], "):(", split.code[[1]][2], ")"), keep.source = FALSE)[[1]])
    if(length(split.code[[1]]) > 2) stop(paste0('Error with this code: ', deparse(code)))
    return(code)
}

modelDefClass$methods(processBUGScode = function(code = NULL, contextID = 1, lineNumber = 0, userEnv) {
    ## uses BUGScode, sets fields: contexts, declInfo$code, declInfo$contextID.
    ## all processing of code is done by BUGSdeclClass$setup(code, contextID).
    ## all processing of contexts is done by BUGScontextClass$setup()
    if(is.null(code)) {
        code <- BUGScode 
        declInfo <<- list()
    }
    for(i in 1:length(code)) {
        if(code[[i]] == '{') if(length(code[[i]])==1) next  ## skip { lines
        lineNumber <- lineNumber + 1
        if(code[[i]][[1]] == '~' || code[[i]][[1]] == '<-') {  ## a BUGS declaration
            iAns <- length(declInfo) + 1
            BUGSdeclClassObject <- BUGSdeclClass$new() ## record the line number (a running count of non-`{` lines) for use in naming nodeFunctions later
            if(code[[i]][[1]] == '~') {
                code[[i]] <- replaceDistributionAliases(code[[i]])
                checkUserDefinedDistribution(code[[i]], userEnv)
            }
            if(code[[i]][[1]] == '<-')
                checkForDeterministicDorR(code[[i]])

            BUGSdeclClassObject$setup(code[[i]], contextID, lineNumber)
            declInfo[[iAns]] <<- BUGSdeclClassObject
        }
        if(code[[i]][[1]] == 'for') {        ## e.g. (for i in 1:N).  New context (for-loop info) needed
            indexVarExpr <- code[[i]][[2]]   ## This is the `i`
            if(length(contexts) > 0) {
                if(as.character(indexVarExpr) %in%
                   contexts[[contextID]]$indexVarNames)
                    stop(paste0(
                        "Variable ",
                        as.character(indexVarExpr),
                        " used multiple times as for loop index in nested\n",
                        "loops.\n",
                        "If your model has macros or if-then-else blocks\n",
                        "you can inspect the processed model code by doing\n",
                        "nimbleOptions(stop_after_processing_model_code = TRUE)\n",
                        "before calling nimbleModel.\n"
                    ),
                    call. = FALSE)
            }
            indexRangeExpr <- code[[i]][[3]] ## This is the `1:N`
            if(nimbleOptions()$prioritizeColonLikeBUGS)
                indexRangeExpr <- reprioritizeColonOperator(indexRangeExpr)
            nextContextID <- length(contexts) + 1
            forCode <- code[[i]][1:3]        ## This is the (for i in 1:N) without the code block
            forCode[[3]] <- indexRangeExpr
            singleContexts <- c(
                if(contextID == 1) NULL
                else contexts[[contextID]]$singleContexts, ## concatenate any current contexts
                list(BUGSsingleContextClass$new(indexVarExpr = indexVarExpr,       ## Add the new context
                                                indexRangeExpr = indexRangeExpr,
                                                forCode = forCode))
            )
            BUGScontextClassObject <- BUGScontextClass$new()
            BUGScontextClassObject$setup(singleContexts = singleContexts)
            contexts[[nextContextID]] <<- BUGScontextClassObject
            if(length(code[[i]][[4]])==1) {
                stop(paste0('Error, not sure what to do with ', deparse(code[[i]])))
            }
            recurseCode <- if(code[[i]][[4]][[1]] == '{') {
                code[[i]][[4]]
            } else {
                substitute( {ONELINE}, list(ONELINE = code[[i]][[4]]))
            }
            lineNumber <- processBUGScode(recurseCode, nextContextID, lineNumber = lineNumber, userEnv = userEnv)  ## Recursive call to process the contents of the for loop
        }
        if(code[[i]][[1]] == '{') {  ## recursive call to a block contained in a {}, perhaps as a result of processCodeIfThenElse
            lineNumber <- processBUGScode(code[[i]], contextID, lineNumber = lineNumber, userEnv = userEnv)
        }
        if(!deparse(code[[i]][[1]]) %in% c('~', '<-', 'for', '{')) 
            stop("Error: ", deparse(code[[i]][[1]]), " not allowed in BUGS code in ", deparse(code[[i]]))
    }
    lineNumber
})

# check if distribution is defined and if not, attempt to register it
checkUserDefinedDistribution <- function(code, userEnv) {
    dist <- as.character(code[[3]][[1]])
    if(dist %in% c("T", "I")) 
        dist <- as.character(code[[3]][[2]][[1]])
    if(!dist %in% distributions$namesVector)
        if(!exists('distributions', nimbleUserNamespace, inherits = FALSE) || !dist %in% nimbleUserNamespace$distributions$namesVector) {
            registerDistributions(dist, userEnv)
            cat("NIMBLE has registered ", dist, " as a distribution based on its use in BUGS code. Note that if you make changes to the nimbleFunctions for the distribution, you must call 'deregisterDistributions' before using the distribution in BUGS code for those changes to take effect.\n", sep = "") 
        }
}
        

replaceDistributionAliases <- function(code) {
    dist <- as.character(code[[3]][[1]])
    trunc <- FALSE
    if(dist %in% c("T", "I")) {
        dist <- as.character(code[[3]][[2]][[1]])
        trunc <- TRUE
    }
    if(dist %in% names(distributionAliases)) {
        dist <- as.name(distributionAliases[dist])
        if(trunc) code[[3]][[2]][[1]] <- dist else code[[3]][[1]] <- dist
    }
    return(code)
}

checkForDeterministicDorR <- function(code) {
    if(is.call(code[[3]])) {
        drFuns <- c(distribution_dFuns, distribution_rFuns)
        if(exists("distributions", nimbleUserNamespace, inherits = FALSE)) {
            dFunsUser <- get('namesVector', nimbleUserNamespace$distributions)
            drFuns <- c(drFuns, dFunsUser, paste0("r", stripPrefix(dFunsUser)))
        }
        if(as.character(code[[3]][[1]]) %in% drFuns)
            warning("Model includes deterministic assignment using '<-' of the result of a density ('d') or simulation ('r') calculation. This is likely not what you intended in: ", deparse(code), ".")
    }
    return(NULL)
}

modelDefClass$methods(splitConstantsAndData = function() {
    # removes items from constantsNamesList that appear as variables in declInfo
    # also, move detected data to 'data'
    # this deals with case when 'data' are passed in as 'constants'
    if(length(constantsNamesList)) {
        vars <- sapply(declInfo, function(x) x$targetVarName)
        constantsNames <- as.character(constantsNamesList)
        newDataVars <- constantsNames[constantsNames %in% vars]
        if(length(newDataVars)) {
            if(nimbleOptions('verbose')) cat("Detected", paste(newDataVars, collapse = ','), "as data within 'constants'.\n")
            constantsNamesList <<- constantsNamesList[!constantsNames %in% vars]
            constantsScalarNamesList <<- constantsScalarNamesList[ !(as.character(constantsScalarNamesList) %in% newDataVars) ]
            constantsList[newDataVars] <<- NULL
            for(varName in newDataVars) eval(substitute(rm(varName, envir = constantsEnv), list(varName = varName)))
        }
    }
})
 


modelDefClass$methods(addMissingIndexing = function() {
    ## overwrites declInfo, using dimensionsList, fills in any missing indexing
    for(i in seq_along(declInfo)) {
        BUGSdecl <- declInfo[[i]]
        newCode <- addMissingIndexingRecurse(BUGSdecl$code, dimensionsList)
        
        BUGSdeclClassObject <- BUGSdeclClass$new()
        BUGSdeclClassObject$setup(newCode, BUGSdecl$contextID, BUGSdecl$sourceLineNumber)
        declInfo[[i]] <<- BUGSdeclClassObject
    }
})

addMissingIndexingRecurse <- function(code, dimensionsList) {
    if(!is.call(code)) return(code)   # returns names, numbers
    if(code[[1]] != '[') {
        for(i in seq_along(code))     code[[i]] <- addMissingIndexingRecurse(code[[i]], dimensionsList)
        return(code)
    }
    if(code[[1]] != '[')   stop('something went wrong: expecting a [')
    ## code must be an indexing call, e.g. x[.....]
    if(length(code[[2]]) > 1 && code[[2]][[1]] == '$'){
      code[[2]][[2]] <- addMissingIndexingRecurse(code[[2]][[2]], dimensionsList)
      return(code)
    }
    if(!any(code[[2]] == names(dimensionsList))) {
      ## dimension information was NOT provided for this variable
      ## let's check to make sure all indexes are present
      if(any(unlist(lapply(as.list(code), is.blank)))) {
        stop(paste0('Error: This part of NIMBLE is still under development.', '\n',
                    'The model definition included the expression \'', deparse(code), '\', which contains missing indices.', '\n',
                    'There are two options to resolve this:', '\n',
                    '(1) Explicitly provide the missing indices in the model definition (e.g., \'', deparse(example_fillInMissingIndices(code)), '\'), or', '\n',
                    '(2) Provide the dimensions of variable \'', code[[2]], '\' via the \'dimensions\' argument to nimbleModel(), e.g.,', '\n',
                    '    nimbleModel(code, dimensions = list(', code[[2]], ' = ', deparse(example_getMissingDimensions(code)), '))', '\n',
                    'Thanks for bearing with us.'), call. = FALSE)
      }
      ## and to recurse on all elements
      for(i in seq_along(code))     code[[i]] <- addMissingIndexingRecurse(code[[i]], dimensionsList)
      return(code)
    }
    if(any(code[[2]] == names(dimensionsList))) {
      dimensions <- dimensionsList[[as.character(code[[2]])]]
      ## dimension information WAS provided for this variable
      ## first, just check that the dimensionality of the node is consistent
      if(length(code) != length(dimensions)+2)   stop(paste0('inconsistent dimensionality provided for node \'', code[[2]], '\''))
      ## then, fill in any missing indicies, and recurse on all other elements
      for(i in seq_along(code)) {
        if(is.blank(code[[i]])) {
          code[[i]] <- substitute(1:TOP, list(TOP = as.numeric(dimensions[i-2])))
        } else {
          code[[i]] <- addMissingIndexingRecurse(code[[i]], dimensionsList)
        }
      }
      return(code)
    }
    stop('something went wrong')
}

example_fillInMissingIndices <- function(code) {
    as.call(lapply(as.list(code), function(el) if(is.blank(el)) quote(1:10) else el))
}
example_getMissingDimensions <- function(code) {
    cCall <- quote(c())
    for(i in seq_along(code)[-c(1,2)]) {
        cCall[[i-1]] <- parse(text = paste0('dim', i-2, '_max'))[[1]]
    }
    return(cCall)
}

modelDefClass$methods(processBoundsAndTruncation = function() {
    ## for non-truncated declarations, extracts range info from distribution; for truncated declarations, pulls bounds out of T() syntax
    for(i in seq_along(declInfo)) {
        
        BUGSdecl <- declInfo[[i]]

        if(BUGSdecl$type != 'stoch') next
        callName <- BUGSdecl$distributionName ## replaces deparse(BUGSdecl$valueExpr[[1]])
        if(!(callName %in% c("T", "I"))) {
            truncated <- FALSE
            boundExprs <- getDistributionInfo(callName)$range
        } else {
            truncated <- TRUE
            if(callName == "I")
                warning(paste0("Interpreting I(,) as truncation (equivalent to T(,)) in ", deparse(BUGSdecl$code), "; this is only valid when ", deparse(BUGSdecl$targetExpr), " has no unobserved (stochastic) parents."))
                
            newCode <- BUGSdecl$code
            newCode[[3]] <- BUGSdecl$valueExpr[[2]]  # insert the core density function call

            distName <- as.character(newCode[[3]][[1]])
            if(!getAllDistributionsInfo('pqAvail')[distName]) 
                stop("Cannot implement truncation for ", distName, "; 'p' and 'q' functions not available.")

            distRange <- getDistributionInfo(distName)$range
            boundExprs <- distRange

        
            if(length(BUGSdecl$valueExpr) >= 3 && BUGSdecl$valueExpr[[3]] != "") 
                boundExprs$lower <- BUGSdecl$valueExpr[[3]]
            if(length(BUGSdecl$valueExpr) >= 4 && BUGSdecl$valueExpr[[4]] != "") 
                boundExprs$upper <- BUGSdecl$valueExpr[[4]]
            if(length(BUGSdecl$valueExpr) != 4)
                warning(paste0("Lower and upper bounds not supplied for T(); proceeding with bounds: (",
                               paste(boundExprs, collapse = ','), ")."))
        
            BUGSdecl$code <- newCode
        }
        BUGSdeclClassObject <- BUGSdeclClass$new()
        BUGSdeclClassObject$setup(BUGSdecl$code, BUGSdecl$contextID, BUGSdecl$sourceLineNumber, truncated, boundExprs)
        declInfo[[i]] <<- BUGSdeclClassObject
    }
})


modelDefClass$methods(expandDistributions = function() {
    ## overwrites declInfo for stochastic nodes: calls match.call() on RHS (uses distributions$matchCallEnv)
    for(i in seq_along(declInfo)) {
        
        BUGSdecl <- declInfo[[i]]
        if(BUGSdecl$type != 'stoch') next
        
        newCode <- BUGSdecl$code
        newCode[[3]] <- evalInDistsMatchCallEnv(BUGSdecl$valueExpr)
        
        BUGSdeclClassObject <- BUGSdeclClass$new()
        BUGSdeclClassObject$setup(newCode, BUGSdecl$contextID, BUGSdecl$sourceLineNumber, BUGSdecl$truncated, BUGSdecl$boundExprs)
        declInfo[[i]] <<- BUGSdeclClassObject
    }
})

modelDefClass$methods(checkMultivarExpr = function() {
    checkForExpr <- function(expr) {
        ##output <- FALSE
        if(length(expr) == 1 && class(expr) %in% c("name", "numeric")) return(FALSE)
        if(!deparse(expr[[1]]) == '[') return(TRUE)
        ## recurse only on the first argument of the `[`
        return(checkForExpr(expr[[2]]))
        ## Previously we recursed more completely.  Now we stop because expressions
        ## inside `[` are allowed.
        ## if(!deparse(expr[[1]]) %in% c('[', ':')) return(TRUE)
        ## for(i in 2:length(expr)) 
        ##     if(checkForExpr(expr[[i]])) output <- TRUE
        ## return(output)
    }

    for(i in seq_along(declInfo)) {
        BUGSdecl <- declInfo[[i]]
        if(BUGSdecl$type != 'stoch') next
        ## dist <- deparse(BUGSdecl$valueExpr[[1]])
        dist <- BUGSdecl$distributionName
        ## The following line is a one-time insertion to break testing in any case
        ## where the condition fails.
        ## If the condition is always met, we can use BUGSdecl$distributionName in place of deparse(BUGSdecl$valueExpr[[1]])
        ## if(dist != BUGSdecl$distributionName)
        ##     stop(paste0("dist (", dist,") != BUGSdecl$distributionName (",BUGSdecl$distributionName,")"))
        types <- nimble:::distributions[[dist]]$types
        if(is.null(types)) next
        ## tmp <- strsplit(types, " = ")
        ## nms <- sapply(tmp, `[[`, 1)
        # originally was only checking for expr in multivar dist:
        ## if('value' %in% nms) next
        ## distDim <- parse(text = tmp[[which(nms == 'value')]])[[2]][[2]]
        ## if(distDim < 1) next
        if(length(BUGSdecl$valueExpr) > 1) {
            for(k in 2:length(BUGSdecl$valueExpr)) {
                paramName <- names(BUGSdecl$valueExpr)[k]
                nDim <- types[[paramName]][['nDim']]
                if(is.numeric(nDim))
                    if(nDim == 0) next
                if(checkForExpr(BUGSdecl$valueExpr[[k]])) {
                    ## Draft gentler warning for possible future adoption: message("Warning about parameter '", names(BUGSdecl$valueExpr)[k], "' of distribution '", dist, "': This multivariate parameter is provided as an expression.  If this is a costly calculation, try making it a separate model declaration for it to improve efficiency.")
                    stop("Error with parameter '", names(BUGSdecl$valueExpr)[k], "' of distribution '", dist, "': multivariate parameters cannot be expressions; please define the expression as a separate deterministic variable and use that variable as the parameter.")  
                }
            }
        }
    }
})

modelDefClass$methods(processLinks = function() {
    ## overwrites declInfo (*and adds*) for nodes with link functions (uses linkInverses)
    newDeclInfo <- list()
    for(i in seq_along(declInfo)) {
        BUGSdecl <- declInfo[[i]]
        nextNewDeclInfoIndex <- length(newDeclInfo) + 1
        if(is.null(BUGSdecl$transExpr))     { newDeclInfo[[nextNewDeclInfoIndex]] <- BUGSdecl; next }
        linkText <- deparse(BUGSdecl$transExpr)
        if(!(linkText %in% names(linkInverses)))    stop(paste('Error, unknown link function:',linkText))
        
        if(BUGSdecl$type == 'stoch') {   # stochastic node
            code <- BUGSdecl$code
            code[[2]] <- parse(text = paste0(linkText, '_', BUGSdecl$targetNodeName), keep.source = FALSE)[[1]]  
            
            newRHS <- linkInverses[[linkText]]
            newRHS[[2]] <- code[[2]]
            newCode <- substitute(A <- B, list(A = BUGSdecl$targetNodeExpr, B = newRHS))
            
            BUGSdeclClassObject <- BUGSdeclClass$new()
            BUGSdeclClassObject$setup(code, BUGSdecl$contextID, BUGSdecl$sourceLineNumber, BUGSdecl$truncated, BUGSdecl$boundExprs)
            newDeclInfo[[nextNewDeclInfoIndex]]     <- BUGSdeclClassObject
            
            BUGSdeclClassObject <- BUGSdeclClass$new()
            BUGSdeclClassObject$setup(newCode, BUGSdecl$contextID, BUGSdecl$sourceLineNumber, BUGSdecl$truncated, BUGSdecl$boundExprs)
            newDeclInfo[[nextNewDeclInfoIndex + 1]] <- BUGSdeclClassObject
            
        } else {    # deterministic node
            newRHS <- linkInverses[[linkText]]
            newRHS[[2]] <- BUGSdecl$code[[3]]
            newLHS <- BUGSdecl$targetNodeExpr
            newCode <- substitute(A <- B, list(A = newLHS, B = newRHS))
            
            BUGSdeclClassObject <- BUGSdeclClass$new()
            BUGSdeclClassObject$setup(newCode, BUGSdecl$contextID, BUGSdecl$sourceLineNumber, BUGSdecl$truncated, BUGSdecl$boundExprs)
            newDeclInfo[[nextNewDeclInfoIndex]] <- BUGSdeclClassObject
        }
    }  # close loop over declInfo
    declInfo <<- newDeclInfo
})

modelDefClass$methods(reparameterizeDists = function() {
    for(i in seq_along(declInfo)) {
        BUGSdecl <- declInfo[[i]]     ## grab this current BUGS declation info object
        if(BUGSdecl$type == 'determ')  next  ## skip deterministic nodes
        code <- BUGSdecl$code   ## grab the original code
        valueExpr <- BUGSdecl$valueExpr   ## grab the RHS (distribution)
        distName <- as.character(valueExpr[[1]])
        if(!(distName %in% getAllDistributionsInfo('namesVector')))    stop('unknown distribution name: ', distName)      ## error if the distribution isn't something we recognize
        distRule <- getDistributionInfo(distName)
        numArgs <- length(distRule$reqdArgs)
        newValueExpr <- quote(dist())       ## set up a parse tree for the new value expression
        newValueExpr[[1]] <- as.name(distName)     ## add in the distribution name
        if(numArgs==0) { ## for dflat, or a user-defined distribution might have 0 arguments
          nonReqdArgExprs <- NULL
          boundExprs <- BUGSdecl$boundExprs
        } else {   
          newValueExpr[1 + (1:numArgs)] <- rep(NA, numArgs)      ## fill in the new parse tree with required arguments
          names(newValueExpr)[1 + (1:numArgs)] <- distRule$reqdArgs    ## add names for the arguments
          
          params <- if(length(valueExpr) > 1) as.list(valueExpr[-1]) else structure(list(), names = character()) ## extract the original distribution parameters
          
          if(identical(sort(names(params)), sort(distRule$reqdArgs))) {
            matchedAlt <- 0
          } else {
            matchedAlt <- NULL; count <- 0
            while(is.null(matchedAlt) && count < distRule$numAlts) {
              count <- count + 1
              if(identical(sort(unique(distRule$alts[[count]])), sort(unique(names(params)))))
                matchedAlt <- count
            }
            if(is.null(matchedAlt)) stop(paste0('bad parameters for distribution ', deparse(valueExpr), '. (No available re-parameterization found.)'), call. = FALSE)
          }
          nonReqdArgs <- names(params)[!(names(params) %in% distRule$reqdArgs)]
          for(iArg in 1:numArgs) {   ## loop over the required arguments
            reqdArgName <- distRule$reqdArgs[iArg]
            ## if it was supplied, copy the supplied expression "as is"
            if(reqdArgName %in% names(params)) {
              newValueExpr[[iArg + 1]] <- params[[reqdArgName]];
              next
            }
            if(!matchedAlt) error("Something wrong - looking for alternative parameterization but supplied args are same as required args: ", deparse(valueExpr))
            if(!reqdArgName %in% names(distRule$exprs[[matchedAlt]]))
              stop('Error: could not find ', reqdArgName, ' in alternative parameterization number ', matchedAlt, ' for: ', deparse(valueExpr), '.')
            transformedParameterPT <- distRule$exprs[[matchedAlt]][[reqdArgName]]
            ## fixing issue of pathological-case model variable names, e.g.,
            ## y ~ dnorm(0, tau = sd)
            ## DT, Feb 2016.
            ##for(nm in c(nonReqdArgs, distRule$reqdArgs))
            namesToSubstitute <- intersect(c(nonReqdArgs, distRule$reqdArgs), all.vars(transformedParameterPT))
            for(nm in namesToSubstitute) {
              ## loop thru possible non-canonical parameters in the expression for the canonical parameter
              if(is.null(params[[nm]])) stop('this shouldn\'t happen -- something wrong with my understanding of parameter transformations')
              transformedParameterPT <- parseTreeSubstitute(pt = transformedParameterPT, pattern = as.name(nm), replacement = params[[nm]])
            }
            newValueExpr[[iArg + 1]] <- transformedParameterPT
          }
          
                                        # evaluate boundExprs in context of model
          boundExprs <- BUGSdecl$boundExprs
          reqdParams <- as.list(newValueExpr[-1])
          for(iBound in 1:2) {
            if(!is.numeric(boundExprs[[iBound]])) {
                                        # only expecting boundExprs to be functions of reqdArgs
              if(length(intersect(nonReqdArgs, all.vars(boundExprs[[iBound]]))))
                stop("Expecting expressions for distribution range for ", distName, " to be functions only of required arguments, namely the parameters used in the 'Rdist' element.")
              namesToSubstitute <- intersect(c(distRule$reqdArgs), all.vars(boundExprs[[iBound]]))
              for(nm in namesToSubstitute) {
                if(is.null(params[[nm]])) stop('this shouldn\'t happen -- something wrong with my understanding of parameter transformations')
                boundExprs[[iBound]] <- parseTreeSubstitute(pt = boundExprs[[iBound]], pattern = as.name(nm), replacement = params[[nm]])
              }
            }
          }
          
          ## hold onto the expressions for non-required args
          nonReqdArgExprs <- params[nonReqdArgs]    ## grab the non-required args from the original params list
          names(nonReqdArgExprs) <- if(length(nonReqdArgExprs) > 0) paste0('.', names(nonReqdArgExprs)) else character(0)  ## append '.' to the front of all the old (reparameterized away) param names
          
                                        # insert altParams and bounds into code
        }
        names(boundExprs)[names(boundExprs) %in% c('lower', 'upper')] <- paste0(names(boundExprs)[names(boundExprs) %in% c('lower', 'upper')], '_')
        newValueExpr <- as.call(c(as.list(newValueExpr), nonReqdArgExprs, boundExprs))
        newCode <- BUGSdecl$code
        newCode[[3]] <- newValueExpr
        
        BUGSdeclClassObject <- BUGSdeclClass$new()
                                        # note at this point boundExprs set back to NULL as all info in lower,upper in valueExpr
        BUGSdeclClassObject$setup(newCode, BUGSdecl$contextID, BUGSdecl$sourceLineNumber, BUGSdecl$truncated, NULL)
        declInfo[[i]] <<- BUGSdeclClassObject
      }  # close loop over declInfo
  })
    
modelDefClass$methods(addRemainingDotParams = function() {
    for(iDecl in seq_along(declInfo)) {
        BUGSdecl <- declInfo[[iDecl]]     ## grab this current BUGS declation info object
        if(BUGSdecl$type == 'determ')  next  ## skip deterministic nodes
        valueExpr <- BUGSdecl$valueExpr   ## grab the RHS (distribution)
        newValueExpr <- valueExpr
        defaultParamExprs <- getDistributionInfo(as.character(newValueExpr[[1]]))$altParams
        if(length(defaultParamExprs) == 0)   next   ## skip if there are no altParams defined in distributions
        
        defaultParamNames <- names(defaultParamExprs)
        defaultDotParamNames <- paste0('.', defaultParamNames)
        for(iParam in seq_along(defaultDotParamNames)) {
            dotParamName <- defaultDotParamNames[iParam]
            if(!(dotParamName %in% names(newValueExpr))) {
                defaultParamExpr <- defaultParamExprs[[iParam]]
                subParamExpr <- eval(substitute(substitute(EXPR, as.list(valueExpr)[-1]), list(EXPR=defaultParamExpr)))
                newValueExpr[[dotParamName]] <- subParamExpr
            }
        }
        newCode <- BUGSdecl$code
        newCode[[3]] <- newValueExpr
        BUGSdeclClassObject <- BUGSdeclClass$new()
        BUGSdeclClassObject$setup(newCode, BUGSdecl$contextID, BUGSdecl$sourceLineNumber, BUGSdecl$truncated, BUGSdecl$boundExprs)
        declInfo[[iDecl]] <<- BUGSdeclClassObject
    }
})

modelDefClass$methods(replaceAllConstants = function() {
    ## overwrites declInfo with constants replaced; only replaces scalar constants
    ## does both LHS and RHS of each BUGSdecl code
    for(i in seq_along(declInfo)) {
        newCode <- replaceConstantsRecurse(declInfo[[i]]$code, constantsEnv, constantsNamesList)$code
        
        BUGSdeclClassObject <- BUGSdeclClass$new()
        BUGSdeclClassObject$setup(newCode, declInfo[[i]]$contextID, declInfo[[i]]$sourceLineNumber, declInfo[[i]]$truncated, declInfo[[i]]$boundExprs)
        declInfo[[i]] <<- BUGSdeclClassObject
    }
})

neverReplaceable <- list(
    ## only the names matter, any non-null value will do
    chol = TRUE,
    inverse = TRUE,
    CAR_calcNumIslands = TRUE,
    CAR_calcC = TRUE,
    CAR_calcM = TRUE,
    CAR_calcEVs2 = TRUE,
    CAR_calcEVs3 = TRUE
)

replaceConstantsRecurse <- function(code, constEnv, constNames, do.eval = TRUE) {
    ## This takes as input a call and an environment or list of constants (only names matter)
    ## It replaces constants that involve no indexing
    ## E.g. dnorm(x[N], sd) , where N is a constant, gets N replaced
    ## but dnorm(x[blockID[i]], sd), where i is a for-loop index, does not get replaced at this step
    cLength <- length(code)
    if(cLength == 1) {
        if(is.name(code)) {
            if( any(code == constNames)) {                
                if(do.eval) {
                    origCode <- code
                    code <- as.numeric(eval(code, constEnv))
                    if(length(code) != 1) warning(paste('Code', deparse(origCode),'was given as known but evaluates to a non-scalar.  This is probably not what you want.'))
                }
                return(list(code = code,
                            replaceable = TRUE))
            }
            return(list(code = code,
                        replaceable = FALSE))
        }
        if( is.numeric(code)) {
            return(list(code = code,
                        replaceable = TRUE))
        }
    }
    if(is.call(code)) {
        if(code[[1]] == '[') {
            replacements <- lapply(code[-c(1,2)],
                                   function(x) replaceConstantsRecurse(x, constEnv, constNames))
            for(i in 1:length(replacements)) {
                code[[i+2]] <- replacements[[i]]$code
            }
            replaceables <- unlist(lapply(replacements, function(x) x$replaceable))
            allReplaceable <- all(replaceables) & do.eval
            repVar <- replaceConstantsRecurse(code[[2]], constEnv, constNames, FALSE)
            code[[2]] <- repVar$code
            if(allReplaceable & repVar$replaceable) {
                testcode <- as.numeric(eval(code, constEnv))
                if(length(testcode) == 1) code <- testcode
            }
            return(list(code = code,
                        replaceable = allReplaceable & repVar$replaceable))
        }
        ## call that is not '['
        if(cLength > 1) {
            if(as.character(code[[1]]) %in% c('<-', '~')) {
                replacements <- c(list(replaceConstantsRecurse(code[[2]], constEnv, constNames, FALSE)),
                                  lapply(code[-c(1,2)], function(x) replaceConstantsRecurse(x, constEnv, constNames) ) )
                replacements[[1]]$replaceable <- FALSE
            } else {
                replacements <- lapply(code[-1], function(x) replaceConstantsRecurse(x, constEnv, constNames))
            }
            for(i in 1:length(replacements)) {
                code[[i+1]] <- replacements[[i]]$code
            }
            replaceables <- unlist(lapply(replacements, function(x) x$replaceable))
            allReplaceable <- all(replaceables)
        } else {
            allReplaceable <- TRUE
        }
        if(allReplaceable) {
            if(!any(code[[1]] == getAllDistributionsInfo('namesVector'))) {
                callChar <- as.character(code[[1]])
                if(exists(callChar, constEnv)) {
                    # if(callChar != ':') {
                    if(!is.vectorized(code)) {
                        if(is.null(neverReplaceable[[callChar]])) {
                            if(class(get(callChar, constEnv)) == 'function') {
                                testcode <- as.numeric(eval(code, constEnv))
                                if(length(testcode) == 1) code <- testcode
                            }
                        }
                    }
                }
            }
        }
        return(list(code = code, replaceable = allReplaceable))
    }
    stop('Error, hit end')
}

liftedCallsDoNotAddIndexing <- c(
    'CAR_calcNumIslands'
)

liftedCallsGetIndexingFromArgumentNumbers <- list(
    CAR_calcC = c(1),
    CAR_calcM = c(1),
    CAR_calcEVs2 = c(2),
    CAR_calcEVs3 = c(3)
)

modelDefClass$methods(liftExpressionArgs = function() {
    ## overwrites declInfo (*and adds*), lifts any expressions in distribution arguments to new nodes
    newDeclInfo <- list()
    for(i in seq_along(declInfo)) {
        BUGSdecl <- declInfo[[i]]        ## grab this BUGS declaration info
        valueExpr <- BUGSdecl$valueExpr  ## extract original valueExpr
        newValueExpr <- valueExpr        ## newValueExpr is initially a copy of the old one
        
        nextNewDeclInfoIndex <- length(newDeclInfo) + 1
        
        if(BUGSdecl$type == 'stoch') {
            params <- as.list(valueExpr[-1])   ## extract the original distribution parameters
            paramNames <- names(valueExpr)[-1]
            types <- nimble:::distributions[[BUGSdecl$distributionName]]$types
            ## types may be NULL if all are scalar
            
            for(iParam in seq_along(params)) {
                if(grepl('^\\.', names(params)[iParam]) || names(params)[iParam] %in% c('lower_', 'upper_'))   next        ## skips '.param' names, 'lower', and 'upper'; we do NOT lift these
                paramExpr <- params[[iParam]]
                paramName <- paramNames[iParam]
                if(!isExprLiftable(paramExpr, types[[paramName]]))    next     ## if this param isn't an expression, go ahead to next parameter
                requireNewAndUniqueDecl <- any(contexts[[BUGSdecl$contextID]]$indexVarNames %in% all.vars(paramExpr))
                uniquePiece <- if(requireNewAndUniqueDecl) paste0("_L", BUGSdecl$sourceLineNumber) else ""
                newNodeNameExpr <- as.name(paste0('lifted_', Rname2CppName(paramExpr, colonsOK = TRUE), uniquePiece))   ## create the name of the new node ##nameMashup
                if(deparse(paramExpr[[1]]) %in% liftedCallsDoNotAddIndexing) {   ## skip adding indexing to mixed-size calls
                    newNodeNameExprIndexed <- newNodeNameExpr
                } else {
                    newNodeNameExprIndexed <- addNecessaryIndexingToNewNode(newNodeNameExpr, paramExpr, contexts[[BUGSdecl$contextID]]$indexVarExprs)  ## add indexing if necessary
                }
                
                newValueExpr[[iParam + 1]] <- newNodeNameExprIndexed  ## update the newValueExpr
                
                newNodeCode <- substitute(LHS <- RHS, list(LHS = newNodeNameExprIndexed, RHS = paramExpr))     ## create code line for declaration of new node
                ## if requireNewAndUniqueDecl is TRUE, the _L# is appended to the newNodeNameExpr and it should be impossible for this to be TRUE:
                identicalNewDecl <- checkForDuplicateNodeDeclaration(newNodeCode, newNodeNameExprIndexed, newDeclInfo)
                
                if(!identicalNewDecl) {
                    BUGSdeclClassObject <- BUGSdeclClass$new()
                    BUGSdeclClassObject$setup(newNodeCode, BUGSdecl$contextID, BUGSdecl$sourceLineNumber, FALSE, NULL)   ## keep new declaration in the same context, regardless of presence/absence of indexing
                    newDeclInfo[[nextNewDeclInfoIndex]] <- BUGSdeclClassObject
                    
                    nextNewDeclInfoIndex <- nextNewDeclInfoIndex + 1     ## update for lifting other nodes, and re-adding BUGSdecl at the end
                }
            }    # closes loop over params
        }        
        newCode <- BUGSdecl$code
        newCode[[3]] <- newValueExpr
        
        BUGSdeclClassObject <- BUGSdeclClass$new()
        BUGSdeclClassObject$setup(newCode, BUGSdecl$contextID, BUGSdecl$sourceLineNumber, BUGSdecl$truncated, BUGSdecl$boundExprs)
        newDeclInfo[[nextNewDeclInfoIndex]] <- BUGSdeclClassObject    ## regardless of anything, add BUGSdecl itself in
    }    # closes loop over declInfo
    declInfo <<- newDeclInfo
})
isExprLiftable <- function(paramExpr, type = NULL) {
    ## determines whether a parameter expression is worthy of lifiting up to a new node
    if(is.name(paramExpr))       return(FALSE)
    if(is.numeric(paramExpr))    return(FALSE)
    if(is.call(paramExpr)) {
        callText <- getCallText(paramExpr)
        if(callText == 'chol')         return(TRUE)    ## do lift calls to chol(...)
        if(callText == 'inverse')      return(TRUE)    ## do lift calls to inverse(...)
        if(callText == 'CAR_calcNumIslands') return(TRUE)    ## do lift calls to CAR_calcNumIslands(...)
        if(callText == 'CAR_calcC')    return(TRUE)    ## do lift calls to CAR_calcC(...)
        if(callText == 'CAR_calcM'  )  return(TRUE)    ## do lift calls to CAR_calcM(...)
        if(callText == 'CAR_calcEVs2') return(TRUE)    ## do lift calls to CAR_calcEVs2(...)
        if(callText == 'CAR_calcEVs3') return(TRUE)    ## do lift calls to CAR_calcEVs3(...)
        if(length(paramExpr) == 1)     return(FALSE)   ## don't lift function calls with no arguments
        if(callText == '[')            return(FALSE)   ## don't lift simply indexed expressions:  x[...]
        nDim <- type[['nDim']]
        if(is.numeric(nDim))
            if(nDim > 0)               return(FALSE)  ## beyond above cases, don't lift non-scalar arguments
                                                       ## This case comes after '[' to avoid using '[' as regexp in grepl
        ## if(getCallText(paramExpr) == '[') { ## these lines are for future handling of foo()[]
        ##     if(is.name(paramExpr))          return(FALSE)   ## don't lift simply indexed expressions:  x[...]
        ##                                     return(TRUE)    ## do lift foo(x)[...]
        ## }
        if(is.vectorized(paramExpr))        return(FALSE)   ## don't lift any expression with vectorized indexing,  funName(x[1:10])
        return(TRUE)
    }
    stop(paste0('error, I need to figure out how to process this parameter expression: ', deparse(paramExpr)))
}
addNecessaryIndexingToNewNode <- function(newNodeNameExpr, paramExpr, indexVarExprs) {
    if(is.call(paramExpr) && deparse(paramExpr[[1]]) %in% names(liftedCallsGetIndexingFromArgumentNumbers))
        return(addNecessaryIndexingFromArgumentNumbers(newNodeNameExpr, paramExpr, indexVarExprs))
    usedIndexVarsList <- indexVarExprs[indexVarExprs %in% all.vars(paramExpr)]    # this extracts any index variables which appear in 'paramExpr'
    vectorizedIndexExprsList <- extractAnyVectorizedIndexExprs(paramExpr)    # creates a list of any vectorized (:) indexing expressions appearing in 'paramExpr'
    neededIndexExprsList <- c(usedIndexVarsList, vectorizedIndexExprsList)
    if(length(neededIndexExprsList) == 0)  return(newNodeNameExpr)  # no index variables, or vectorized indexing, return the (un-indexed) name expression
    newNodeNameExprIndexed <- substitute(NAME[], list(NAME = newNodeNameExpr))
    newNodeNameExprIndexed[3:(2+length(neededIndexExprsList))] <- neededIndexExprsList
    return(newNodeNameExprIndexed)
}
addNecessaryIndexingFromArgumentNumbers <- function(newNodeNameExpr, paramExpr, indexVarExprs) {
    paramExprCallName <-  as.character(paramExpr[[1]])
    argNumbers <- liftedCallsGetIndexingFromArgumentNumbers[[paramExprCallName]]
    argList <- as.list(paramExpr[argNumbers + 1])    ## +1 to skip past the function name (first element)
    neededIndexExprsList <-  lapply(argList, function(x) x[[3]])
    newNodeNameExprIndexed <- substitute(NAME[], list(NAME = newNodeNameExpr))
    newNodeNameExprIndexed[3:(2+length(neededIndexExprsList))] <- neededIndexExprsList
    return(newNodeNameExprIndexed)
}
extractAnyVectorizedIndexExprs <- function(expr) {
    if(!(':' %in% all.names(expr)))    return(list())
    if(!is.call(expr))     return(list())
    if(expr[[1]] == ':')     return(expr)
  ##  if(expr[[1]] == '[')    return(as.list(expr[-c(1,2)])) ## 
    ret <- unlist(lapply(expr[-1], function(i) extractAnyVectorizedIndexExprs(i)))
    if(is.null(ret)) return(list()) else return(ret)
}
checkForDuplicateNodeDeclaration <- function(newNodeCode, newNodeNameExprIndexed, newDeclInfo) {
    for(i in seq_along(newDeclInfo)) {
        if(identical(newNodeNameExprIndexed, newDeclInfo[[i]]$targetExpr)) {
            ## we've found a node declaration with exactly the same LHS, which is a mangling of the RHS during lifting
            if(!identical(newNodeCode, newDeclInfo[[i]]$code))   { stop('something fishy going on with our new node declarations.....') }
            return(TRUE)   ## indicate that we found a matching node declaration
        }
    }
    return(FALSE)     ## a duplicate node entry was *not* found
}
modelDefClass$methods(addIndexVarsToDeclInfo = function() {
    ## sets field declInfo[[i]]$indexVariableExprs from contexts.  must be after overwrites of declInfo
    for(i in seq_along(declInfo)) {
        declInfo[[i]]$setIndexVariableExprs(contexts[[declInfo[[i]]$contextID]]$indexVarExprs)
    }
})
modelDefClass$methods(genSymbolicParentNodes = function() {
    ## sets field declInfo[[i]]$symbolicParentNodes. must be after overwrites of declInfo
    
    nimFunNames <- getAllDistributionsInfo('namesExprList')

    for(i in seq_along(declInfo)){
        declInfo[[i]]$genSymbolicParentNodes(constantsNamesList, contexts[[declInfo[[i]]$contextID]], nimFunNames, contextID = declInfo[[i]]$contextID)
    }
})


modelDefClass$methods(genReplacementsAndCodeReplaced = function() {
    ## sets fields declInfo[[i]]$replacements, $codeReplaced, and $replacementNameExprs
    
    nimFunNames <- getAllDistributionsInfo('namesExprList')
    
    for(i in seq_along(declInfo)) {
        declInfo[[i]]$genReplacementsAndCodeReplaced(constantsNamesList, contexts[[declInfo[[i]]$contextID]], nimFunNames)
    }
})

modelDefClass$methods(genReplacedTargetValueAndParentInfo = function() {
    
    nimFunNames <- distributions$namesExprList
    
    for(i in seq_along(declInfo)) {
        declInfo[[i]]$genReplacedTargetValueAndParentInfo(constantsNamesList, contexts[[declInfo[[i]]$contextID]],
                                                          nimFunNames, contextID = declInfo[[i]]$contextID)
    }
    NULL
})

modelDefClass$methods(genAltParamsModifyCodeReplaced = function() {
    ## sets field declInfo[[i]]$altParams, and modifies $codeReplaced to not include .param arguments (if stochastic)
    for(i in seq_along(declInfo)) {
        declInfo[[i]]$genAltParamsModifyCodeReplaced()
    }
})

modelDefClass$methods(genBounds = function() {
    ## sets field declInfo[[i]]$boundExprs and if not truncated modifies $codeReplaced to remove lower,upper
    for(i in seq_along(declInfo)) {
        declInfo[[i]]$genBounds()
    }
})

## genNodeInfo3
modelDefClass$methods(genNodeInfo3 = function(debug = FALSE) {
    ## This uses the contexts (for loops) to create an environment called replacementsEnv that has unrolled indices and replacements
    ## First it iterates through contexts, working with all lines of BUGS code in the same context.
    ## Then it iterates through each line of BUGS code, refining the results needed by it.
    if(debug) browser()

    ## 1. Iterate over context (for loops), where 1st context will always be no-for-loop
    for(i in seq_along(contexts)) {
        boolContext <- unlist(lapply(declInfo, function(x) x$contextID == i))                    ## TRUE for BUGS lines (declInfo elements) that use this context
        allReplacements <- do.call('c', lapply(declInfo[boolContext], `[[`, 'replacements'))     ## Collect replacement expressions from all lines
        if(length(allReplacements) > 0) {
            allReplacementNameExprs <- do.call('c', lapply(declInfo[boolContext], `[[`, 'replacementNameExprs')) ## names of allReplacements as expressions
            boolNotDup <- !duplicated(allReplacements) ## remove duplicates, e.g. if i+1 appears in two lines in the same expression, we only it once as a replacement

            ## This makes an environment with a vector of each replacement and for-loop index from executing the for-loops
            unrolledContextAndReplacementsEnv <- expandContextAndReplacements(allReplacements[boolNotDup], allReplacementNameExprs[boolNotDup], contexts[[i]], constantsEnv)
            ## record the environment in the context
            contexts[[i]]$replacementsEnv <<- unrolledContextAndReplacementsEnv
        } else {
            ## if there were no replacements:
            contexts[[i]]$replacementsEnv <<- NULL
        }
    }

    ## There is a tricky disinction of cases that come out of previous step, from expandContextAndReplacements
    ## If there is a context with NO replacements (must mean that all non-for-loop lines have no replacements)
    ##      then expandContextAndReplacements returns NULL
    ## If there is a context with NO indices but >0 replacements, then a valid result comes back from expandContextAndReplacements
    ##      with outputSize = 1
    ## If there is a context with >0 indices but none of them yield nodes (e.g. they descend, for(i in 2:1), which by default option results in no iteration)
    ##      then there WILL be an environment from expandContextAndReplacements but it will have outputSize == 0
    ## In later processing we need to know when such declarations (they are in BUGS code, but due to indices they do nothing) occur
    ##      so we record numUnrolledNodes.
    
    ## 2. iterate over declInfo (one entry for each BUGS declaration)
    for(i in seq_along(declInfo)) {
        ## extract information needed for each declInfo and each parentExpr
        ## note that sometimes lifted nodes don't need some or any of the indices from the context
        ## e.g. there might be a lifted node involving no indices, but it is still embedded in the same context
        BUGSdecl <- declInfo[[i]]
        context <- contexts[[BUGSdecl$contextID]]

        ## set up the BUGSdecl$replacementsEnv and related bits.
        ## These may get changed again below
        if(is.null(context$replacementsEnv)) { ## This would occur if there were no for loops and no replacements
            BUGSdecl$replacementsEnv <- NULL
            BUGSdecl$unrolledIndicesMatrix <- matrix(nrow = 0, ncol = 0)
            BUGSdecl$outputSize <- 0
            BUGSdecl$numUnrolledNodes <- 1
            next
        } else { ## there was at least for loop and/or at least one replacement
            ## copy (by reference) the replacementsEnv to this BUGSdecl
            BUGSdecl$replacementsEnv <- context$replacementsEnv
            BUGSdecl$outputSize <- BUGSdecl$replacementsEnv$outputSize
            BUGSdecl$numUnrolledNodes <- BUGSdecl$outputSize ## will only be 0 if the for loops all had numeric(0) index ranges, like for(i in 2:1)
        }

        ## Pick out which parts of the context (which for loop indices) are used in this BUGSdecl
        if(nimbleOptions()$allowDynamicIndexing) {
            ## We need NA as indexExpr for useContext to correctly determine not to use context for dynamic indexes.
            indexExprWithNA <- lapply(BUGSdecl$indexExpr, function(x) if(isDynamicIndex(x)) as.numeric(NA) else x)
            useContext <- unlist(lapply(context$singleContexts, function(x) isNameInExprList(x$indexVarExpr, indexExprWithNA)))
        } else useContext <- unlist(lapply(context$singleContexts, function(x) isNameInExprList(x$indexVarExpr, BUGSdecl$indexExpr)))
        ## We want to do something like cbind(i, i_plus_1, j) to make a matrix of unrolled indices
        ## To do that we need to construct the cbind expression with the name expressions needed and then eval it in replacementsEnv

        ## We will include anything that is not a list
        rNEtoInclude <- unlist(lapply(names(BUGSdecl$replacementNameExprs), function(x) !is.list( BUGSdecl$replacementsEnv[[x]])))
        
        BUGSdecl$replacementsEnv$cbindArgExprs <- unique(c(context$indexVarExprs[useContext], BUGSdecl$replacementNameExprs[rNEtoInclude]))
        BUGSdecl$unrolledIndicesMatrix <- with(BUGSdecl$replacementsEnv, do.call('cbind', cbindArgExprs))
        rm(list = 'cbindArgExprs', envir = BUGSdecl$replacementsEnv)

        ## 
        if(!all(useContext)) {
            if(!any(useContext)) {
                ## A line that is in contexts but doesn't use any of them can arise from lifting.  In such a case, keep only the first row.
                boolUse <- c(TRUE, rep(FALSE, BUGSdecl$outputSize-1))
            } else {
                ## if not every context was used, some cleanup is needed: the unrolledIndicesMatrix may have duplicates to remove
                boolUse <- !duplicated(BUGSdecl$unrolledIndicesMatrix[, context$indexVarNames[useContext] ] )
            }
            if(!is.null(BUGSdecl$unrolledIndicesMatrix)) BUGSdecl$unrolledIndicesMatrix <- BUGSdecl$unrolledIndicesMatrix[boolUse, , drop = FALSE]

            ## And then give a new replacementsEnv to the BUGSdecl from the no-duplicates matrix
            ## BUGSdecl$replacementsEnv <- list2env(as.data.frame(BUGSdecl$unrolledIndicesMatrix)) ## used to work until there were lists involved
            BUGSdecl$replacementsEnv <- new.env()
            for(repName in names(BUGSdecl$replacementNameExprs)) BUGSdecl$replacementsEnv[[repName]] <- context$replacementsEnv[[repName]][boolUse]
            
            BUGSdecl$replacementsEnv$outputSize <- sum(boolUse)
            BUGSdecl$outputSize <- BUGSdecl$replacementsEnv$outputSize
            BUGSdecl$numUnrolledNodes <- BUGSdecl$outputSize
        }
        if(is.null(BUGSdecl$unrolledIndicesMatrix)) BUGSdecl$unrolledIndicesMatrix <- matrix(nrow = 0, ncol = 0)
    }
})

isNameInExprList <- function(target, codeList) {
    for(i in seq_along(codeList)) {
        if(isNameInExpr(target, codeList[[i]])) return(TRUE)
    }
    FALSE
}

isNameInExpr <- function(target, code) {
    if(length(code) == 1) return(identical(target, code))
    for(i in seq_along(code)) {
        if(isNameInExpr(target, code[[i]])) return(TRUE)
    }
    FALSE
}

nm_seq_noDecrease <- function(a, b) {
    if(a > b) {
        numeric(0)
    } else {
        a:b
    }
}

expandContextAndReplacements <- function(allReplacements, allReplacementNameExprs, context, constantsEnv) {
##    browser()
    ## allReplacements is a list like
    ## list(i = i, i_plus_1 = i+1, mean_x_1to5 = mean(x[1:5]))
    ## context is a BUGScontextClass object
    ## constantsEnv is an environment with constants that can be used to permanently replace values in the allReplacements code

    numContexts <- length(context$singleContexts)
    if(numContexts == 0) { ## it has no indices or known indices
        if(length(allReplacements)==0) {
            context$replacementsEnv <<- NULL
            return(NULL)
        }
    }
    
    ## when done, we will have created a new environment and want to remove the constants from it
    namesToRemoveAtEnd <- ls(constantsEnv)
    constantsEnvCopy <- list2env(as.list(constantsEnv))
    ## some replacements like min(j:100) should no longer be needed but are still there

    ## If this all works, useContext can be removed
    useContext <- rep(TRUE, numContexts)
    
    valueVarNames <- if(numContexts > 0) paste0("INDEXVALUE_", 1:numContexts, "_") else character(0)
    ## indexRecordingCode gives lines of code like "INDEXVALUE_1_[iAns] <- i". This will later have its name changed to "i"
    indexRecordingCode <- vector('list', length = numContexts)
    for(i in seq_along(context$singleContexts)) {
        if(useContext[i])
            indexRecordingCode[[i]] <- substitute(V[iAns] <- index, list(V = as.name(valueVarNames[i]), index = context$singleContexts[[i]]$indexVarExpr))
    }

    numReplacements <- length(allReplacements)
    useReplacement <- unlist(lapply(allReplacementNameExprs, function(x) { ## do not use replacements that are identical to indexVars
        for(i in seq_along(context$singleContexts)) {
            if( identical(context$singleContexts[[i]]$indexVarExpr, x) ) return(FALSE)
        }
        return(TRUE)
    }))
    ## replacementRecordingCode gives lines of code like "i_plus_1[iAns] <- i+1"
    replacementRecordingCode <- vector('list', length = numReplacements)
    for(i in seq_along(replacementRecordingCode)) {
        if(useReplacement[i])
            replacementRecordingCode[[i]] <- substitute(A[[iAns]] <- B, list(A = allReplacementNameExprs[[i]], B = allReplacements[[i]])) 
    }

    ## From here through the while loop combines the for loops from the contexts, with the replacementRecordingCode and indexRecordingCode in the innermost
    innerLoopCode <- as.call(c(list(quote(`{`)), replacementRecordingCode, indexRecordingCode, quote(iAns <- iAns + 1)))

    innerLoopCode <- context$embedCodeInForLoop(innerLoopCode, useContext)
    ## at this point "innerLoopCode" has the full loop  ## determineContextSize does something similar -- creates and executes nested for loops -- only for the purpose of counting how big the result will be
    outputSize <- determineContextSize(context, useContext, constantsEnvCopy)
    for(i in seq_along(context$singleContexts)) {
        if(useContext[i])
            assign(valueVarNames[i], numeric(outputSize), constantsEnvCopy)
    }
    for(i in seq_along(replacementRecordingCode)) {
        if(useReplacement[i])
            assign(names(allReplacements)[i], vector('list', length = outputSize), constantsEnvCopy)
    }
    assign("iAns", 1, constantsEnvCopy)
    eval(innerLoopCode, constantsEnvCopy)
    for(i in seq_along(context$singleContexts)){
        if(useContext[i]) {
            constantsEnvCopy[[ as.character(context$singleContexts[[i]]$indexVarExpr) ]] <- constantsEnvCopy[[ valueVarNames[i] ]]
            rm(list = valueVarNames[i], envir = constantsEnvCopy)
        }
    }
    ## Turn lists into vectors when all elements are scalars.  When not, ensure all list elements are numeric, not integer, to avoid compiler mix-ups.
    for(i in seq_along(allReplacementNameExprs)) {
        if(useReplacement[i]) {
            unlistScalarCode <- substitute( {
                FOO_allScalar <- all(unlist(lapply(VARNAME, function(x) length(x) == 1)))
                if(FOO_allScalar) VARNAME <- unlist(VARNAME) ## Ok to have integers here
                else {
                    for(FOO_i in seq_along(VARNAME)) storage.mode(VARNAME[[FOO_i]]) <- 'double' ## but not here
                    rm(FOO_i)
                }
                rm(FOO_allScalar)
            }, list(VARNAME = allReplacementNameExprs[[i]]) )
            eval(unlistScalarCode, envir = constantsEnvCopy)
            ##rm(list = 'FOO_allScalar', envir = constantsEnvCopy)
        }
    }

    rm(list = c(namesToRemoveAtEnd, 'iAns'), envir = constantsEnvCopy)
    assign("outputSize", outputSize, constantsEnvCopy)
    return(constantsEnvCopy) ## becomes replacementsEnv
}

determineContextSize <- function(context, useContext = rep(TRUE, length(context$singleContexts)), evalEnv = new.env()) {
    ## could improve this by checking for nested loops that don't use indices from outer loops
    innerLoopCode <- quote(iAns <- iAns + 1)
    innerLoopCode <- context$embedCodeInForLoop(innerLoopCode, useContext)

    assign("iAns", 0L, evalEnv)
    test <- try(eval(innerLoopCode, evalEnv))
    if(is(test, 'try-error'))
        stop("Could not evaluate loop syntax: is indexing information provided via 'constants'?")
    ans <- evalEnv$iAns
    rm(list = c('iAns', context$indexVarNames[useContext]), envir = evalEnv)
    return(ans)
}

#### wrapAsNumeric <- function(code) substitute(as.numeric(X), list(X = code))     ## as.numeric() flattens everything to a vector
wrapAsNumeric <- function(code) substitute({ value <- CODE;   storage.mode(value) <- 'numeric';   value }, list(CODE = code))

removeNonScalarElementsFromList <- function(lst) {
    scalarFlags <- unlist(lapply(lst, function(el) { if(!is.null(dim(el))) return(FALSE);   if(length(el)>1) return(FALSE);   return(TRUE) }))
    return(lst[scalarFlags])
}
modelDefClass$methods(removeEmptyBUGSdeclarations = function() {
    numberIndexedNodes <- unlist(lapply(declInfo, function(di) length(di$indexedNodeInfo)))
    exptyDeclInfoIndexes <- (numberIndexedNodes == 0)
    declInfo[exptyDeclInfoIndexes] <<- NULL
})


## turn x[ a, b:c ] into x[a[iAns], b[iAns], c[iAns] ]
insertSubIndexExpr <- function(code, indexCode) {
    if(length(code) == 1) {
        if(is.name(code)) 
            return(substitute(a[i], list(a = code, i = indexCode)))
        return(code)
    }
    if(code[[1]] == '[') {
        if(length(code)==2) return(code) ## not sure if this can occur
        for(i in 3:length(code)) code[[i]] <- insertSubIndexExpr(code[[i]], indexCode)
        return(code)
    }
    if(code[[1]] == ':') {
        for(i in 2:length(code)) code[[i]] <- insertSubIndexExpr(code[[i]], indexCode)
        return(code)
    }
}

removeColonOperator <- function(code) {
     if(length(code) == 1) {
        return(code)
    }
    if(code[[1]] == '[') {
        if(length(code)==2) return(code) ## not sure if this can occur
        for(i in 3:length(code)) code[[i]] <- removeColonOperator(code[[i]])
        return(code)
    }
    if(code[[1]] == ':') {
        return(code[[2]])
    }
}

## utility function used by makeVertexNamesFromIndexArray2
makeSplitInfo <- function(splitIndices) {
    r <- range(splitIndices)     ## [min, max]
    vec <- r[1]!=r[2]            ## logical: is it a vector
    contig <- if(vec)            ## logical: is it contiguous
                  ## first condition checks for gaps like [1 3 4], in which case (4-1+1) != length(c(1,3,4))
                  ## but this won't catch cases like [1 3 4; 1 2 3 4]
                  diff(r)+1 == length(unique(splitIndices)) &
                      ## second condition checks for unequal counts of some indices: for [1 3 4; 1 2 3 4], there will be unequal counts
                      length(unique(table(splitIndices))) == 1
              else
                  TRUE   ## doesn't really matter below if we define a scalar as contiguous or not contiguous
    c(r, as.numeric(vec), as.numeric(contig)) ## make all numeric for future rbind into matrix
}


## This function takes a array of vertex indices and returns a set of node names based on shared indices
## e.g. c(1, 1, 2, 2, 2) for varName 'x' would return 'x[1:2]' and 'x[3:5]'
## when there are non-contiguous indices, `%.s%` is used instead of `:`.  `%.s%`  can be parsed but is not intended to be evaluated
## e.g. c(1, 2, 2, 1, 1) would return x[1 %.s% 5] and x[2:3]
## Interwoven non-contiguous indices such as c(1, 2, 1, 2, 2) should never occur based on how splitVertices works
makeVertexNamesFromIndexArray2 <- function(indArr, minInd = 1, varName) {
    dims <- dim(indArr)
    nDim <- length(dims)
    indArr[indArr < minInd] <- NA
    if(all(is.na(indArr))) return(list(indices = integer(), names = character()))
    ## arrayWithIndices is a list of arrays of the same size as x
    ## In the first, elements give the row index, e.g. [1 1 1; 2 2 2; 3 3 3]
    ## In the second, elements give the column index, e.g. [1 2 3; 1 2 3; 1 2 3]
    ## etc.
    arrayWithIndices <- vector('list', length = nDim)
    arrayWithIndices[[1]] <- array( rep(1:dims[1], prod(dims[-1])), dims)
    ## We could reduce memory footprint by doing all steps on each dimension before building arrayWithIndices for new dimension
    if(nDim > 1) {
        permutation <- 1:nDim
        for(iD in 2:nDim) {
            usePerm <- permutation
            usePerm[1] <- iD
            usePerm[iD] <- 1
            ##            arrayWithIndices[[iD]] <- aperm(arrayWithIndices[[1]], usePerm)
            arrayWithIndices[[iD]] <- aperm(array( rep(1:dims[iD], prod(dims[-1])), dims[usePerm]), usePerm)
        }
    }
    ## splits is a list of vectors of the indices that share the same indArr value
    ## e.g. if indArr is [1 1 1; 2 2 2; 2 2 2]
    ## Then the first element of splits will be `1`: [1 1 1] and the second will be `2`: [2 3 2 3 2 3].  These are vectors of row numbers
    ## And the second element of splits will be `1`: [1 2 3] and the second will be `2`: [1 1 2 2 3 3].  These are vectors of column numbers
    ## etc.
    splits <- lapply(arrayWithIndices, split, indArr)

    ## info is a list of summaries from the splits
    ## each entry is (min, max, 0/1 for vector, 0/1 for contiguous)
    info <- lapply(splits, lapply, makeSplitInfo)

    ## From here on is the construction of string labels from the info
    dimStrings <- lapply(info, function(x) {
        all <- do.call('rbind', x)   ## This makes a table with a row for each unique element of indArr
        ## and columns following the order of info
        seps <- rep(':', nrow(all))  ## initialize seps and modify later if needed 
        scal <- all[,3]==0           ## which rows are for scalar elements
        seps[scal] <- ''             ## set the sep for scalars to ''
        seps[all[,4]==0] <- '%.s%'    ## for rows that are not contiguous, use i %.s% j. Any actual call to %.s% results in an error.
        maxStrs <- as.character(all[,2]) ## maximums
        maxStrs[scal] <- ''              ## clear maximums for scalars 
        paste0(all[,1], seps, maxStrs)   ## paste minimum-separator-maximum
    })
    dimStrings[['sep']] <- ', '
    newNames <- paste0(varName, '[',  do.call('paste', dimStrings), ']') ## paste together pieces from different dimensions
    list(indices = as.integer(names(splits[[1]])), names = newNames)
}

## This converts i%.s%j (indicating a split vertex) to i:j
## This will only be called for RHSonly nodes.
## These are correct, even after splitting, to use in eval(parse(...)) in the "vars2..." environments
## But split vertices that are LHSinferred keep their %.s% and should never end up evaluated in a "vars2..." environment
convertSplitVertexNamesToEvalSafeNames <- function(origNames, ok = TRUE) {
    boolConvert <- grepl("%.s%", origNames)
    convertOneName <- function(oneName) {
        parsedName <- parse(text = oneName, keep.source = FALSE)[[1]]
        boolSplitByIndex <- rep(FALSE, length(parsedName)-2)
        for(i in 3:length(parsedName)) {
            if(!is.name(parsedName[[i]])) {
                if(deparse(parsedName[[i]][[1]]) == '%.s%') {
                    parsedName[[i]][[1]] <- quote(`:`)
                    boolSplitByIndex[i-2] <- TRUE
                }
            }
        }
        deparse(parsedName)
    }
    origNames[boolConvert] <- unlist(lapply(origNames[boolConvert], convertOneName))
    origNames
}

splitVertexIDsToElementIDs <- function(var2vertexID, nextVertexID) {
    vertexIDtable <- tabulate(var2vertexID) ## this fills in 0s for all elements
    boolMultiple <- vertexIDtable > 1
    newElementID_2_vertexID <- numeric()
    if(any(boolMultiple)) {
        IDsToFix <- which(boolMultiple)
        for(ID in IDsToFix) {
            newElementIDs <- nextVertexID-1 + 1:vertexIDtable[ID]
            var2vertexID[which(var2vertexID == ID)] <- newElementIDs
            newElementID_2_vertexID <- c(newElementID_2_vertexID, rep(ID, vertexIDtable[ID]))
            nextVertexID <- nextVertexID + vertexIDtable[ID]
        }
    }
    list(var2vertexID = var2vertexID, nextVertexID = nextVertexID, newElementID_2_vertexID = newElementID_2_vertexID)
}

splitCompletionForOrigNodes <- function(var2nodeOrigID, var2vertexID, maxOrigNodeID, nextVertexID) {
    ## could potentially be combined with collecting edges from inferred vertices, but let's wait
    origIDtable <- tabulate(var2nodeOrigID, maxOrigNodeID)
    vertexIDtable <- tabulate(var2vertexID, maxOrigNodeID)
    boolInOrig <- origIDtable > 0
    ok <- vertexIDtable[boolInOrig] == origIDtable[boolInOrig] | vertexIDtable[boolInOrig] == 0
    if(any(!ok)) {
        IDsToFix <- which(boolInOrig)[!ok]
        for(ID in IDsToFix) {
            var2vertexID[var2vertexID == ID] <- nextVertexID
            nextVertexID <- nextVertexID + 1
        }
    }
    list(var2vertexID = var2vertexID, nextVertexID = nextVertexID) 
}

splitVertices <- function(var2vertexID, unrolledBUGSindices, indexExprs = NULL, indexNames = NULL, parentExpr, parentExprReplaced = NULL, parentIndexNamePieces, replacementNameExprs, nextVertexID, maxVertexID, rhsVarInfo = NULL, debug = FALSE) {
    if(debug) browser()
    ## parentIndexNamePieces: when there is an NA for a dynamic index, we should assume a split on the full range of values.
    ## This is the same case as anyContext = TRUE and all(useContent) = TRUE, i.e. boolUseUnrolledRow <- rep(TRUE, nrow(unrolledBUGSindices))
    ## The "context" label may sometimes apply to x[dynamicI] even with no context.
    ## actually (6/12/17) note that no splitting based on dynamic indices will occur here because we don't enter splitVertices with any NA in parentExpr indices

    ## 1. Determine which indexExprs are in parentExpr
    useContext <- unlist(lapply(indexExprs, isNameInExprList, parentExpr))
    if(nimbleOptions()$allowDynamicIndexing) {
        dynamicIndices <- detectDynamicIndexes(parentExpr)
        numDynamicIndices <- sum(dynamicIndices)
        if(numDynamicIndices) {
            stop("splitVertices: no unknown (NA) indices should be detected here; please contact the NIMBLE development team.")
         ##   if(length(parentExpr) > 3 || length(indexExprs) > 1) stop("splitVertices: dynamic indexing with multiple indices not yet allowed in: ", deparse(parentExpr))  ## allowing this to pass for now to handle mu[3,NA,i]
            useDynamicIndices <- TRUE
            ranges <- data.frame(rbind(rhsVarInfo$mins[dynamicIndices], rhsVarInfo$maxs[dynamicIndices]))
            unrolledVarIndices <- as.matrix(do.call(expand.grid, lapply(ranges, function(x) x[1]:x[2])))
            dimnames(unrolledVarIndices)[[2]] <- NULL
        } else useDynamicIndices <- FALSE
    }
    
    anyContext <- any(useContext)
    ## 2. Use unique or duplicated on unrolledBUGSindices to get a needed set
    if(anyContext) {
        if(!(all(useContext)))
            boolUseUnrolledRow <- !duplicated(unrolledBUGSindices[, indexNames[useContext] ]) ## relies on column ordering
        else
            boolUseUnrolledRow <- rep(TRUE, nrow(unrolledBUGSindices))
    }
    ## 3. Determine if indices are all scalar
        allScalar <- TRUE
        ## parentIndexNamePieces NA: need to decide if these count for scalar processing. I think they would count as scalar, using the 
    vectorIndices <- lapply(parentIndexNamePieces, function(x) {
        if(is.list(x)) {
            allScalar <<- FALSE;
            return(TRUE)
        }
        FALSE
    })

    ## step 4 evaporated
    if(all(is.na(var2vertexID)))
        currentVertexCounts <- rep(0, maxVertexID)
    else
        currentVertexCounts <- tabulate(var2vertexID, max(max(var2vertexID, na.rm = TRUE), maxVertexID))
    ## 5. Set up initial table of vertexIDcounts

    ## 6. All scalar case: iterate or vectorize via cbind and put new vertexIDs over -1s
    if(allScalar) {
        if(anyContext) { ## parentIndexNamePieces NA: This block is where we want to go.  Goal of next step is to set varIndicesToUse. For an NA, we want the full extent.
            boolIndexNamePiecesExprs <- !unlist(lapply(parentIndexNamePieces, is.numeric)) 
            if(all(boolIndexNamePiecesExprs)) {
                varIndicesToUse <- unrolledBUGSindices[ boolUseUnrolledRow, unlist(parentIndexNamePieces) ]
            } else {
                if(nimbleOptions()$allowDynamicIndexing) {
                    boolIndexNamePiecesExprs <- boolIndexNamePiecesExprs[!dynamicIndices]
                    parentIndexNamePieces <- parentIndexNamePieces[!dynamicIndices]
                }
                varIndicesToUse <- matrix(nrow = sum(boolUseUnrolledRow), ncol = length(parentIndexNamePieces))
                varIndicesToUse[, boolIndexNamePiecesExprs] <- unrolledBUGSindices[ boolUseUnrolledRow, unlist(parentIndexNamePieces)[boolIndexNamePiecesExprs]]
                indexPieceNumericInds <- which(!boolIndexNamePiecesExprs)
                for(iii in seq_along(indexPieceNumericInds))
                    varIndicesToUse[, indexPieceNumericInds[iii] ] <-
                        parentIndexNamePieces[[ indexPieceNumericInds[iii] ]]
            }
        }    
        else { ## is it hard-coded like x ~ dnorm(mu[1,2], 1)
            if(is.null(parentIndexNamePieces))
                stop("Error in splitVertices: you may have omitted indexing for a multivariate variable: ", as.character(parentExprReplaced), ".")
            if(nimbleOptions()$allowDynamicIndexing && numDynamicIndices == length(dynamicIndices)) {
                varIndicesToUse <- NULL  # only dynamic indexing
            } else {
                if(nimbleOptions()$allowDynamicIndexing) {
                    parentIndexNamePieces <- parentIndexNamePieces[!dynamicIndices]
                    parentExprReplaced <- parentExprReplaced[c(1:2,which(!dynamicIndices)+2)]
                }
                if(length(parentIndexNamePieces)==1) varIndicesToUse <- as.numeric(parentExprReplaced[[3]])
                else {
                    varIndicesToUse <- matrix(0, nrow = 1, ncol = length(parentIndexNamePieces))
                    for(iI in 1:ncol(varIndicesToUse)) varIndicesToUse[1, iI] <- as.numeric(parentExprReplaced[[iI+2]])
                }
            }
        }

        ## now expand.grid (actually Cartesian product) with necessary varInfo
        if(nimbleOptions()$allowDynamicIndexing) 
            if(useDynamicIndices) {
                if(is.null(varIndicesToUse)) {
                    varIndicesToUse <- unrolledVarIndices
                } else {
                    tmp <- as.matrix(merge(unrolledVarIndices, varIndicesToUse, by.x = NULL, by.y = NULL)) # incorrect column order
                    varIndicesToUse <- 0; length(varIndicesToUse) <- length(tmp)
                    dim(varIndicesToUse) <- dim(tmp)
                    # put back in correct order 
                    varIndicesToUse[ , dynamicIndices] <- tmp[ , 1:numDynamicIndices]
                    if(numDynamicIndices < length(dynamicIndices))
                        varIndicesToUse[ , !dynamicIndices] <- tmp[ , (numDynamicIndices+1):ncol(varIndicesToUse)]
                }
            }
        
        ## parentIndexNamePieces Should there be a unique in one of the next lines? Or varIndicesToUse <- unique(varIndicesToUse).
        ## OR use a !duplicated construction in boolUseUnrolledRow <- rep(TRUE, nrow(unrolledBUGSindices)) above
        currentVertexIDs <- var2vertexID[varIndicesToUse]
        needsVertexID <- is.na(currentVertexIDs)
        numNewVertexIDs <- sum(needsVertexID)
        if(numNewVertexIDs > 0) {
            var2vertexID[varIndicesToUse][needsVertexID] <- nextVertexID - 1 + 1:numNewVertexIDs
            nextVertexID <- nextVertexID + numNewVertexIDs
        }
        ## Still need to look for splits on other existing vertexIDs
        oldIndices <- !needsVertexID
        existingIndices <- currentVertexIDs[oldIndices]
##        uniqueExistingIndices <- unique(existingIndices)
        currentCountsOldIndices <- currentVertexCounts[ existingIndices ]
        needsSplit <- currentCountsOldIndices > 1

        numNeedSplit <- sum(needsSplit)
        if(numNeedSplit > 0) {
            var2vertexID[varIndicesToUse][ oldIndices][needsSplit] <- nextVertexID - 1 + 1:numNeedSplit
            nextVertexID <- nextVertexID + numNeedSplit
        }
        ## This can result in skips, e.g. if a previous BUGSdecl labeled a vector with a new vectorID, and that now gets split in scalar vectorID labels
        ## Then the earlier vectorID will be gone forever
        ## later we will clean up skips and splits of original nodes
    } else {

        if(anyContext) {
            colNums <- 1:ncol(unrolledBUGSindices)
            names(colNums) <- dimnames(unrolledBUGSindices)[[2]]
            newIndexExprs <- lapply(replacementNameExprs, function(x) substitute(unrolledBUGSindices[iRow, iCol], list(X = x, iCol = colNums[as.character(x)])))
            accessExpr <- eval( substitute( substitute(AE, newIndexExprs), list(AE = parentExprReplaced) ) )
            iRowRange <- (1:nrow(unrolledBUGSindices))[boolUseUnrolledRow]
        } else {
            accessExpr <- parentExprReplaced
            iRowRange <- 1
        }
        accessExpr[[2]] <- quote(var2vertexID)
        assignExprNextVID <- substitute( A <- nextVertexID, list(A = accessExpr))
        assignExprCVIDB <- assignExprNextVID
        assignExprCVIDB[[3]] <- quote(currentVertexIDblock)
        ## construct argList
        for(iRow in iRowRange) {
            currentVertexIDblock <- eval(accessExpr)
            uniqueCurrentVertexIDs <- unique(c(currentVertexIDblock))
            if(length(uniqueCurrentVertexIDs)==1) { ## current block has only 1 ID
                if( is.na(uniqueCurrentVertexIDs[1]) |      ## It's all unassigned OR
                   currentVertexCounts[ uniqueCurrentVertexIDs[1] ] != length(currentVertexIDblock) ) { ## It does not fully cover  existing vertexID
                    if(!is.na(uniqueCurrentVertexIDs[1])) currentVertexCounts[ uniqueCurrentVertexIDs[1] ] <- currentVertexCounts[ uniqueCurrentVertexIDs[1] ] - sum(currentVertexIDblock == uniqueCurrentVertexIDs[1])
                    eval(assignExprNextVID) ## var2vertexID[ all indexing stuff ] <- nextVertexID
                    nextVertexID <- nextVertexID + 1
                }
            } else { ## need to iterate through IDs
                for(VID in uniqueCurrentVertexIDs) {
                    if(is.na(VID)) {
                        boolIsNA <- is.na(currentVertexIDblock)
                        currentVertexIDblock[ boolIsNA ] <- nextVertexID
                        currentVertexCounts[nextVertexID] <- currentVertexCounts[nextVertexID] + sum(boolIsNA)
                        nextVertexID <- nextVertexID + 1
                    } else {
                        boolWithinBlock <- currentVertexIDblock == VID
                        numWithinBlock <- sum(boolWithinBlock, na.rm = TRUE)
                        if(currentVertexCounts[ VID ] != numWithinBlock) {
                            currentVertexIDblock[boolWithinBlock] <- nextVertexID
                            currentVertexCounts[VID] <- currentVertexCounts[VID] - numWithinBlock
                            currentVertexCounts[nextVertexID] <- currentVertexCounts[nextVertexID] + numWithinBlock
                            nextVertexID <- nextVertexID + 1
                        }
                    }
                }
                eval(assignExprCVIDB) ## var2vertexIDs[ all indexing stuff] <- currentVertexIDblock
            }
        }
    }
    list(var2vertexID = var2vertexID, nextVertexID = nextVertexID)
    
}

collectInferredVertexEdges <- function(var2nodeID, var2vertexID) {
    splitVerts <- lapply(split(var2vertexID, var2nodeID), unique)
    nodeIDs <- as.numeric(names(splitVerts))
    edges <- mapply(function(from , to ) {
        if(length(to) == 1)
            if(from == to)
                return(NULL)
        matrix(c( rep(from, length(to)), to), nrow = length(to) )
    }, nodeIDs, splitVerts, SIMPLIFY = FALSE)
    edges <- do.call('rbind', edges)
    list(edgesFrom = edges[,1], edgesTo = edges[,2])
}

collectEdges <- function(var2vertexID, unrolledBUGSindices, targetIDs, indexExprs = NULL, parentExprReplaced = NULL, parentIndexNamePieces, replacementNameExprs, debug = FALSE) {
    ## This collects edges from extracted information of BUGS declarations
    ## var2vertexID is the shape of a rhsVar (right hand side variable, e.g. "mu") giving vertexIDs
    ## targetIDs are the IDs of all LHS nodes (i.e. from unrolling any for loops)
    ## unrolledBUGSindices is the unrolled index matrix
    ## indexExprs gives the expressions of the for-loop indices, e.g. i and j for a pair of nested loops using those
    ## parentExprReplaced gives the list of parent node expressions after replacements, e.g. "mu[i_plus_1, 2]" if the original code had "mu[i+1,2]" on the RHS
    ## parentIndexNamePieces gives a corresponding list of the index expressions of each parentExprReplaced, e.g. i_plus_i, 2.
    ## replacementNameExpressions has the names of things that are replacements, e.g. i_plus_1
    if(debug) browser()
    anyContext <- ncol(unrolledBUGSindices) > 0 
    if(length(anyContext)==0) stop('collectEdges: problem with anyContext')

    if(nimbleOptions()$allowDynamicIndexing) {
        ## replace NA with 1 to index into first element, since for unknownIndex vars there should be only one vertex
        for(iii in seq_along(parentIndexNamePieces))
            if(length(parentIndexNamePieces[[iii]]) == 1 && isDynamicIndex(parentIndexNamePieces[[iii]]))
                parentIndexNamePieces[[iii]] <- 1
        if(length(parentExprReplaced) >= 3) {
            indexExprs <- parentExprReplaced[3:length(parentExprReplaced)]
            for(iii in seq_along(indexExprs)) 
            if(length(indexExprs[[iii]]) == 1 && is.numeric(indexExprs[[iii]]) && is.na(indexExprs[[iii]]))  ## is.numeric check avoids warning when check k[i] type situation
                parentExprReplaced[[iii+2]] <- 1
        }
    }
    
    allScalar <- TRUE
    vectorIndices <- lapply(parentIndexNamePieces, function(x) {if(is.list(x)) {allScalar <<- FALSE; return(TRUE)}; FALSE})
        
    if(allScalar) {
        if(is.null(parentIndexNamePieces)) varIndicesToUse <- rep(1, length(targetIDs))
        else {
            if(anyContext) {
                boolIndexNamePiecesExprs <- !unlist(lapply(parentIndexNamePieces, is.numeric)) ##!is.numeric(parentIndexNamePieces)
                if(all(boolIndexNamePiecesExprs)) {
                    test <- try(varIndicesToUse <- unrolledBUGSindices[ , unlist(parentIndexNamePieces), drop = FALSE])
                    if(inherits(test, 'try-error')) stop('collectEdges: problem with unrolledBUGSindices')
                } else {
                    varIndicesToUse <- matrix(nrow = nrow(unrolledBUGSindices), ncol = length(parentIndexNamePieces))
                    varIndicesToUse[, boolIndexNamePiecesExprs] <- unrolledBUGSindices[ , unlist(parentIndexNamePieces)[boolIndexNamePiecesExprs], drop = FALSE]
                    indexPieceNumericInds <- which(!boolIndexNamePiecesExprs)
                    for(iii in seq_along(indexPieceNumericInds)) varIndicesToUse[, indexPieceNumericInds[iii] ] <- parentIndexNamePieces[[ indexPieceNumericInds[iii] ]]
                }
            } else {
                if(length(parentIndexNamePieces)==1) varIndicesToUse <- rep(as.numeric(parentExprReplaced[[3]]), length(targetIDs))
                else {
                    varIndicesToUse <- matrix(0, nrow = 1, ncol = length(parentIndexNamePieces))
                    for(iI in 1:ncol(varIndicesToUse)) varIndicesToUse[1, iI] <- as.numeric(parentExprReplaced[[iI+2]])
                }
            }
        }

        edgesFrom <- var2vertexID[varIndicesToUse]
        edgesTo <- targetIDs
        
    } else {
        if(anyContext) {
           if(length(ncol(unrolledBUGSindices))==0) stop('collectEdges: problem with unrolledBUGSindices')
            colNums <- 1:ncol(unrolledBUGSindices)
            names(colNums) <- dimnames(unrolledBUGSindices)[[2]]
            newIndexExprs <- lapply(replacementNameExprs, function(x) substitute(unrolledBUGSindices[iRow, iCol], list(X = x, iCol = colNums[as.character(x)])))
            accessExpr <- eval( substitute( substitute(AE, newIndexExprs), list(AE = parentExprReplaced) ) )
            iRowRange <- (1:nrow(unrolledBUGSindices))
        } else {
            accessExpr <- parentExprReplaced
            iRowRange <- 1
        }
        accessExpr[[2]] <- quote(var2vertexID)

        uniqueCurrentVertexIDsList <- lapply(iRowRange, function(iRow) unique(as.numeric(eval(accessExpr))))
        edgesFrom <- do.call('c', uniqueCurrentVertexIDsList)
        edgesToList <- lapply(iRowRange, function(iRow) rep(targetIDs[iRow], length(unique(as.numeric(eval(accessExpr))))))
        edgesTo <- do.call('c', edgesToList)
    }
    list(edgesFrom = edgesFrom, edgesTo = edgesTo)    
}


## pull out dynamic indexing info for use in constraining range in nodeFunctions and then strip out USED_IN_INDEX tagging and replace .DYN_INDEX tagged indexing code with NA
## original plan was for some code here (if based on variables) or later (if based on vertices) to find the elements used in dynamic indexing for when we planned to dynamically update the graph
modelDefClass$methods(findDynamicIndexParticipants = function() {
    if(nimbleOptions()$allowDynamicIndexing) {
        for(iDI in seq_along(declInfo)) {
            if(declInfo[[iDI]]$type == "unknownIndex") next
            declInfo[[iDI]]$dynamicIndexInfo <<- list()
            for(iSPN in seq_along(declInfo[[iDI]]$symbolicParentNodesReplaced)) {
                symbolicParent <- declInfo[[iDI]]$symbolicParentNodesReplaced[[iSPN]]
                dynamicIndexes <- detectDynamicIndexes(symbolicParent)
                ## We do not yet check bounds of inner indexes in nested indexing. To do so we need to
                ## find dynamic indexing within a USED_IN_INDEX() and add to dynamicIndexInfo;
                ## then in nodeFunctions we need nested if statements so inner index is checked first.
                ## That being said, compiled execution will error out with appropriate out of bounds error
                ## because C++ will put an out-of-bound value in for 'k' in k[d[0]] or k[d[1342134]].
                if(any(dynamicIndexes)) {
                    indexedVar <- stripUnknownIndexFromVarName(deparse(symbolicParent[[2]]))
                    numSPNR <- length(declInfo[[iDI]]$symbolicParentNodesReplaced)
                    for(iIndex in which(dynamicIndexes)) {
                        declInfo[[iDI]]$dynamicIndexInfo[[length(declInfo[[iDI]]$dynamicIndexInfo) + 1]] <<-
                            list(indexCode = stripDynamicallyIndexedWrapping(symbolicParent[[2+iIndex]]),
                                 lower = varInfo[[indexedVar]]$mins[iIndex],
                                 upper = varInfo[[indexedVar]]$maxs[iIndex])
                        declInfo[[iDI]]$symbolicParentNodes[[iSPN]][[2+iIndex]] <<- as.numeric(NA) ## Indexing code is not needed anymore.
                        if(iSPN <= numSPNR)
                            declInfo[[iDI]]$symbolicParentNodesReplaced[[iSPN]][[2+iIndex]] <<- as.numeric(NA) ## Indexing code is not needed anymore.
                    }
                }
            }
            declInfo[[iDI]]$symbolicParentNodes <<- lapply(declInfo[[iDI]]$symbolicParentNodes, stripIndexWrapping)
            declInfo[[iDI]]$symbolicParentNodesReplaced <<- lapply(declInfo[[iDI]]$symbolicParentNodesReplaced, stripIndexWrapping)
        }
    }
})

modelDefClass$methods(addFullDimExtentToUnknownIndexDeclarations = function() {
    if(nimbleOptions()$allowDynamicIndexing) {
        for(iDI in seq_along(declInfo)) {
            if(declInfo[[iDI]]$type == "unknownIndex") {
                parentExpr <- declInfo[[iDI]]$symbolicParentNodes[[1]]
                parentExprReplaced <- declInfo[[iDI]]$symbolicParentNodesReplaced[[1]]
                targetExpr <- declInfo[[iDI]]$targetExpr
                targetExprReplaced <- declInfo[[iDI]]$targetExprReplaced
                varName <- declInfo[[iDI]]$rhsVars[1] # deparse(parentExpr[[2]])
                dynamicIndices <- detectDynamicIndexes(parentExpr)
                ranges <- data.frame(rbind(varInfo[[varName]]$mins[dynamicIndices], varInfo[[varName]]$maxs[dynamicIndices]))
                if(any(ranges[2, ] == 1))
                    stop("Variable ", varName, " is dynamically-indexed but has at least one dimension of length one, probably because NIMBLE could not automatically determine its dimensionality. Please provide dimensions via the 'dimensions' argument.")
                fullExtent <- lapply(ranges, function(x) 
                    substitute(X:Y, list(X = x[1], Y = x[2])))
                parentExpr[which(dynamicIndices)+2] <- fullExtent
                parentExprReplaced[which(dynamicIndices)+2] <- fullExtent
                targetExpr[which(dynamicIndices)+2] <- fullExtent
                targetExprReplaced[which(dynamicIndices)+2] <- fullExtent
                declInfo[[iDI]]$symbolicParentNodes[[1]] <<- parentExpr
                declInfo[[iDI]]$symbolicParentNodesReplaced[[1]] <<- parentExprReplaced
                declInfo[[iDI]]$targetExpr <<- targetExpr
                declInfo[[iDI]]$targetExprReplaced <<- targetExprReplaced
                count <- 1
                for(p in seq_along(declInfo[[iDI]]$parentIndexNamePieces[[1]])) # only one parent by construction
                    if(dynamicIndices[p]) {
                        declInfo[[iDI]]$parentIndexNamePieces[[1]][[p]] <<- as.list(ranges[[count]])
                        count <- count + 1
                    }
                ## count <- 1
                ## for(p in seq_along(declInfo[[iDI]]$targetIndexNamePieces))
                ##      if(dynamicIndices[p]) {
                ##          declInfo[[iDI]]$targetIndexNamePieces[[p]] <<- as.list(ranges[[count]])
                ##          count <- count + 1
                ##      }
            } 
        }
    }
})

modelDefClass$methods(genExpandedNodeAndParentNames3 = function(debug = FALSE) {
## This is a (ridiculously long) function that works through the BUGS lines (declInfo elements) and varInfo to set up the graph and maps
    if(debug) browser()
    ## 1. initialize origMaps:
    ##     "orig" in the label refers to these being the IDs that will be initially assigned.
    ##            Later the IDs will be changed when the graph is topologically sorted
    vars_2_nodeOrigID <- new.env()     ## IDs for node function labels, e.g. for x[1:4] might be (NA, NA, NA, 2) if x[1:3] is RHS-only and x[4] is LHS
    vars_2_vertexOrigID <- new.env()   ## IDs for node labels, e.g. x[1:4] might be: 1, 1, 1, 2
    vars2LogProbName <- new.env()      ## e.g. "x" might be "logProb_x" if it has any LHS
    ##vars2LogProbID <- new.env()        ## yiels LogProbIDs, which are not sorted in any way and we might move away from.
    ## 1b. set up variables in all the environments
    for(iV in seq_along(varInfo)) {
        varName <- varInfo[[iV]]$varName
        if(varInfo[[iV]]$nDim > 0) {
            vars_2_nodeOrigID[[varName]] <- array(as.numeric(NA), dim = varInfo[[iV]]$maxs)
            if(nimbleOptions()$allowDynamicIndexing)
                if(varName %in% unknownIndexNames) next ## unknownIndex objects have no logProb
            vars2LogProbName[[varName]] <- array(dim = varInfo[[iV]]$maxs)
           ## vars2LogProbID[[varName]] <- array(dim = varInfo[[iV]]$maxs)
            storage.mode(vars2LogProbName[[varName]]) <- 'character'
        } else {
            vars_2_nodeOrigID[[varName]] <- as.numeric(NA)
            if(nimbleOptions()$allowDynamicIndexing)
                if(varName %in% unknownIndexNames) next ## unknownIndex objects have no logProb
            vars2LogProbName[[varName]] <- as.character(NA)
            ##vars2LogProbID[[varName]] <- as.numeric(NA)
        }
    }

    ## 2. collect names and do eval to create variables in vars_2_nodeOrigID
    next_origID <- 1              ## this is a counter for the next origID to be assigned
    allNodeNames <- character()   ## vector of all node names
    types <- character()          ## parallel vector of types such as "stoch", "determ", "RHSonly", and "LHSinferred"
    ## 2b. Use LHS declarations create nodes:
    ## note we need to create nodes for unknownIndex vars because origID is used in creating edges from parent variable to unknownIndex variable
    if(debug) browser()
    for(iDI in seq_along(declInfo)) {
        BUGSdecl <- declInfo[[iDI]]
        if(BUGSdecl$numUnrolledNodes == 0) next
        lhsVar <- BUGSdecl$targetVarName
        nDim <- varInfo[[lhsVar]]$nDim
        if(nDim > 0) {  ## pieces is a list of index text to construct node names, e.g. list("1", c("1:2", "1:3", "1:4"), c("3", "4", "5"))
            pieces <- vector('list', nDim)
            for(i in 1:nDim) {
                indexNamePieces <- BUGSdecl$targetIndexNamePieces[[i]]
                if(is.list(indexNamePieces)) {
                    pieces[[i]] <- paste0( if(is.character(indexNamePieces[[1]]))
                                               BUGSdecl$replacementsEnv[[ indexNamePieces[[1]] ]]
                                           else indexNamePieces[[1]],
                                          ':',
                                          if(is.character(indexNamePieces[[2]]))
                                              BUGSdecl$replacementsEnv[[ indexNamePieces[[2]] ]]
                                          else indexNamePieces[[2]] )
                } else {
                    pieces[[i]] <- if(is.character(indexNamePieces)) BUGSdecl$replacementsEnv[[ indexNamePieces ]] else indexNamePieces
                }
            }
            pieces[['sep']] <- ', '
            BUGSdecl$nodeFunctionNames <- paste0(lhsVar, '[', do.call('paste', pieces), ']')  ## create the names.  These are LHS so they are nodeFunctions
            allNodeNames <- c(allNodeNames, BUGSdecl$nodeFunctionNames)
            types <- c(types, rep(BUGSdecl$type, length(BUGSdecl$nodeFunctionNames) ) )       ## append vector of "stoch" or "determ" to types vector
            BUGSdecl$origIDs <- next_origID -1 + (1:length(BUGSdecl$nodeFunctionNames))       ## record the original IDs used here
            next_origID <- next_origID + length(BUGSdecl$nodeFunctionNames)
            
            ## Fill in the vars_2_nodeOrigID elements
            if(is.environment(BUGSdecl$replacementsEnv)) { ## this means there was some replacement involved in this BUGS line
                BUGSdecl$replacementsEnv[['origIDs']] <- BUGSdecl$origIDs
                IDassignCode <- insertSubIndexExpr(BUGSdecl$targetExprReplaced, quote(iAns))
                ## make a line of code to evaluate in the replacementsEnv
                forCode <- substitute( for(iAns in 1:OUTPUTSIZE) ASSIGNCODE <- origIDs[iAns], list(OUTPUTSIZE = BUGSdecl$outputSize, ASSIGNCODE = IDassignCode) )
                BUGSdecl$replacementsEnv[[lhsVar]] <- vars_2_nodeOrigID[[lhsVar]]
                eval(forCode, envir = BUGSdecl$replacementsEnv)
                vars_2_nodeOrigID[[lhsVar]] <- BUGSdecl$replacementsEnv[[lhsVar]]
                rm(list = lhsVar, envir = BUGSdecl$replacementsEnv)
            } else {
                ## If no replacementsEnv was set up, then there were no index variables (only numerics)
                eval(substitute(A <- B, list(A = BUGSdecl$targetExprReplaced, B = BUGSdecl$origIDs)), envir = vars_2_nodeOrigID)
            }
        } else { ## nDim == 0
            BUGSdecl$nodeFunctionNames <- lhsVar ## name simply for the variable
            allNodeNames <- c(allNodeNames, BUGSdecl$nodeFunctionNames)
            types <- c(types, BUGSdecl$type)
            BUGSdecl$origIDs <- next_origID
            vars_2_nodeOrigID[[lhsVar]] <- next_origID
            next_origID <- next_origID + 1
        }
    }

    maxOrigNodeID <- next_origID - 1
    if(nimbleOptions()$allowDynamicIndexing) {
        nodeNamesLHSall <- allNodeNames[types != "unknownIndex"]
        numNodeFunctions <<- maxOrigNodeID - sum(types == "unknownIndex")
    } else {
        nodeNamesLHSall <- allNodeNames
        numNodeFunctions <<- maxOrigNodeID
    }
    allNewVertexNames <- character(0)
    allNewVertexIDs <- integer(0)

    ## 2b. Collect logProbNames and do eval to create variables in vars2LogProbName and vars2LogProbID
    ## This section is very similar to above, except that index ranges are collapsed to their beginning.
    ## E.g. a LHS "x[1:5]" must be a multivariate node, but it only needs "logProb[1]", not "logProb[1:5]"
    if(debug) browser()
    nextLogProbID <- 1
    logProbNames <- character()
    logProbIDs <- integer()
    for(iDI in seq_along(declInfo)) {
        BUGSdecl <- declInfo[[iDI]]
        if(BUGSdecl$type == 'determ' || BUGSdecl$type == "unknownIndex") next
        if(BUGSdecl$numUnrolledNodes == 0) next
        lhsVar <- BUGSdecl$targetVarName
        logProbVarName <- makeLogProbName(lhsVar)
        nDim <- logProbVarInfo[[logProbVarName]]$nDim
        if(nDim > 0) {
            pieces <- vector('list', nDim)
            for(i in 1:nDim) {
                indexNamePieces <- BUGSdecl$targetIndexNamePieces[[i]]
                if(is.list(indexNamePieces)) { 
                    pieces[[i]] <- if(is.character(indexNamePieces[[1]]))
                                       BUGSdecl$replacementsEnv[[ indexNamePieces[[1]] ]]
                                   else indexNamePieces[[1]] ## for anything with a colon, this takes only the start value (cf step 2. above)
                } else {
                    pieces[[i]] <- if(is.character(indexNamePieces)) BUGSdecl$replacementsEnv[[ indexNamePieces ]] else indexNamePieces
                }
            }
            pieces[['sep']] <- ', '
            newLogProbNames <- paste0(logProbVarName, '[', do.call('paste', pieces), ']') 
            logProbNames <- c(logProbNames, newLogProbNames)
            newLogProbIDs <- nextLogProbID - 1 + 1:length(newLogProbNames)
            logProbIDs <- c(logProbIDs, newLogProbIDs)
            nextLogProbID <- nextLogProbID + length(newLogProbNames)

            targetExprWithMins <- removeColonOperator(BUGSdecl$targetExprReplaced) ## this strips colons, leaving only start value
            if(is.environment(BUGSdecl$replacementsEnv)) {
                IDassignCode <- insertSubIndexExpr(targetExprWithMins, quote(iAns))
                BUGSdecl$replacementsEnv[['logProbIDs']] <- newLogProbIDs
                forCode <- substitute( for(iAns in 1:OUTPUTSIZE) ASSIGNCODE <- logProbIDs[iAns], list(OUTPUTSIZE = BUGSdecl$outputSize, ASSIGNCODE = IDassignCode) )
##                BUGSdecl$replacementsEnv[[lhsVar]] <- vars2LogProbID[[lhsVar]]
##                eval(forCode, envir = BUGSdecl$replacementsEnv)
##                vars2LogProbID[[lhsVar]] <- BUGSdecl$replacementsEnv[[lhsVar]]
##                rm(list = c(lhsVar, 'logProbIDs'), envir = BUGSdecl$replacementsEnv)

                BUGSdecl$replacementsEnv[['logProbIDs']] <- newLogProbNames
                BUGSdecl$replacementsEnv[[lhsVar]] <- vars2LogProbName[[lhsVar]]
                eval(forCode, envir = BUGSdecl$replacementsEnv)
                vars2LogProbName[[lhsVar]] <- BUGSdecl$replacementsEnv[[lhsVar]]
                rm(list = c(lhsVar, 'logProbIDs'), envir = BUGSdecl$replacementsEnv)
            } else {
                ## If no replacementsEnv was set up, then there were no index variables (only numerics)
##                eval(substitute(A <- B, list(A = targetExprWithMins, B = newLogProbIDs)), envir = vars2LogProbID)
                eval(substitute(A <- B, list(A = targetExprWithMins, B = newLogProbNames)), envir = vars2LogProbName)
            }
        } else { ## nDim == 0
            logProbNames <- c(logProbNames, logProbVarName)
            logProbIDs <- c(logProbIDs, nextLogProbID)
            ##vars2LogProbID[[lhsVar]] <- nextLogProbID
            vars2LogProbName[[lhsVar]] <- logProbVarName
            nextLogProbID <- nextLogProbID + 1
        }
    }
    
    ## 3. determine total model size from all objects
    ##    and initialize vars_2_vertexOrigID from vars_2_nodeOrigID
    ## A "vertex" is any vertex that will go in the graph.  This can include nodeFunctions, fractured nodeFunctions (e.g. if x[1:2] is declared but also used as x[1] and x[2], and RHSonly parts
    ## for now treat unknownIndex vars as part of totModelSize (otherwise would need to avoid having elements for unknownIndex vars later)
    totModelSize <- 0
    for(iV in seq_along(varInfo)) {
        totModelSize <- totModelSize + if(varInfo[[iV]]$nDim == 0) 1 else prod(varInfo[[iV]]$maxs)
        varName <- varInfo[[iV]]$varName
        vars_2_vertexOrigID[[ varName ]] <- vars_2_nodeOrigID[[ varName ]]
    }
    

    ## 4. Use RHS pieces to split vertices in vars_2_vertexOrigID
    ##    E.g. say x[1:2] ~ dmnorm(...)
    ##             for(i in 1:2) z[i] <- x[i]
    ##    From above, there is a nodeOrigID for "x[1:2]"
    ##    But now, based on the RHS usage of x[1] and x[2], there needs to be a separate vertexID for each of them
    ##
    if(debug) browser()
    nextVertexID <- maxOrigNodeID+1 ## origVertexIDs start after origNodeIDs.  These will be sorted later.
    for(iDI in seq_along(declInfo)) {  ## Iterate over BUGS declarations (lines)
        BUGSdecl <- declInfo[[iDI]]
        if(BUGSdecl$numUnrolledNodes == 0) next
        rhsVars <- BUGSdecl$rhsVars
        for(iV in seq_along(rhsVars)) {  ## Iterate over the RHS variables in a BUGS line
            rhsVar <- rhsVars[iV]
            if(rhsVar %in% unknownIndexNames) next ## vertex splitting occurs on lifted unknownIndex declaration, not on unknownIndex usage on RHS
            nDim <- varInfo[[rhsVar]]$nDim
            if(nDim != length(BUGSdecl$parentIndexNamePieces[[iV]]))  # this check should be redundant with equivalent check in genVarInfo3
               stop("Dimension of ", rhsVar, " is ", nDim, ", which does not match its usage in '", deparse(BUGSdecl$code), "'.")
            if(nDim > 0) { ## Make the split.  This function is complicated.
                splitAns <- splitVertices(vars_2_vertexOrigID[[rhsVar]], BUGSdecl$unrolledIndicesMatrix,
                                          contexts[[BUGSdecl$contextID]]$indexVarExprs, contexts[[BUGSdecl$contextID]]$indexVarNames,
                                          BUGSdecl$symbolicParentNodes[[iV]], BUGSdecl$symbolicParentNodesReplaced[[iV]],
                                          BUGSdecl$parentIndexNamePieces[[iV]], BUGSdecl$replacementNameExprs, nextVertexID, totModelSize, varInfo[[rhsVar]])
                
                vars_2_vertexOrigID[[rhsVar]] <- splitAns[[1]]
                nextVertexID <- splitAns[[2]]
            } else { ## The RHS var is scalar
                if(is.na(vars_2_vertexOrigID[[rhsVar]])) { ## And it wasn't on the LHS, so it still has NA
                    vars_2_vertexOrigID[[rhsVar]] <- nextVertexID  ## so assign a vertexOrigID
                    nextVertexID <- nextVertexID + 1
                }
            }
        }
    }

    ## 5. In cases where a vertex was split leaving an incomplete part of an original node, complete the split by creating an inferred vertex for the remaining part
    ##    E.g. consider x[1:4] ~ dmnorm(...)
    ##                  for(i in 1:2) z[i] ~ x[i]
    ##    So x[1] and x[2] were split by the above section, but x[3:4] still has the same origID as x[1:4].
    ##    In this step a piece like x[3:4] is identified and gets a new unique vertexID
    if(debug) browser()
    for(iV in seq_along(varInfo)) {
        if(nimbleOptions()$allowDynamicIndexing)
            if(varInfo[[iV]]$varName %in% unknownIndexNames) next 
        if(varInfo[[iV]]$nDim > 0) {
            varName <- varInfo[[iV]]$varName
            if(varInfo[[iV]]$nDim > 0) {
                splitFix <- splitCompletionForOrigNodes(vars_2_nodeOrigID[[varName]], vars_2_vertexOrigID[[varName]], maxOrigNodeID, nextVertexID)
                vars_2_vertexOrigID[[varName]] <- splitFix[[1]]
                nextVertexID <- splitFix[[2]]
            }
        }
    }
    if(debug) browser()
    
    ## 6. Make vertex names
    ##    E.g. from the results of the previous step, we may now need vertex names "x[1]", "x[2]" and "x[3:4]"
    ##    Those are constructed here from the vars_2_vertexOrigID elements
    numVertices <- nextVertexID - 1
    vertexID_2_nodeID <- c(1:maxOrigNodeID, integer(numVertices - maxOrigNodeID) )
    for(iV in seq_along(varInfo)) {
        varName <- varInfo[[iV]]$varName
        
        if(varInfo[[iV]]$nDim > 0) {
            newVertexNames <- makeVertexNamesFromIndexArray2(vars_2_vertexOrigID[[varName]], maxOrigNodeID + 1, varName = varName)
            allNewVertexNames <- c(allNewVertexNames, newVertexNames$names)
            allNewVertexIDs <- c(allNewVertexIDs, newVertexNames$indices)
        } else {
            if(vars_2_vertexOrigID[[varName]] > maxOrigNodeID) {
                allNewVertexNames <- c(allNewVertexNames, varName)
                allNewVertexIDs <- c(allNewVertexIDs, vars_2_vertexOrigID[[varName]])
            }
        }
    }
    allVertexNames <- c(allNodeNames, character(length(allNewVertexNames)))
    allVertexNames[allNewVertexIDs] <- allNewVertexNames ## indexed by vertexID

    ## 7. Re-order vertexIDs to make them contiguous
    ##    The above steps can result in gaps in the ID sequences, so here we re-order to plug those gaps
    ##    Sorting comes later
    ##
    ## Make a vector to map contiguous IDs to the original vertex IDs
    ##   The originalNodeIDs, which are at the beginning of the vertexIDs, are already guaranteed to be continuous
    if(debug) browser()
    contigID_2_origVertexID <- c(1:maxOrigNodeID, sort(allNewVertexIDs))
    numContigVertices <- length(contigID_2_origVertexID)
    ## Now make a vector to map the other way
    origVertexID_2_contigID <- numeric(nextVertexID-1)
    origVertexID_2_contigID[contigID_2_origVertexID] <- 1:length(contigID_2_origVertexID)
    ## Re-order the vertexNames accordingly (gaps would have "" entries, I think)
    allVertexNames <- allVertexNames[contigID_2_origVertexID]

    ## 7b. check for existence of non-orig nodes and add them to types vector
    ## for now these are all labeled as "RHSonly" and later we'll use the vectorID_2_nodeID to relabel some as "LHSinferred"
    if(debug) browser()
    if(length(allVertexNames) > maxOrigNodeID) {
        nodeNamesRHSonly <- allVertexNames[ (maxOrigNodeID + 1) : length(allVertexNames) ]
        types <- c(types, rep('RHSonly', length(allVertexNames) - maxOrigNodeID))
    } else
        nodeNamesRHSonly <- character()
    
    ## 7c. Re-label the IDs in the vars_2_vertexOrigID with the contiguous version
    for(iV in seq_along(varInfo)) {
        temp <- vars_2_vertexOrigID[[varInfo[[iV]]$varName]]
        if(!all(is.na(temp))) vars_2_vertexOrigID[[varInfo[[iV]]$varName]][] <- origVertexID_2_contigID[temp] 
    }

    ## 7d. Treat unknownIndex variables like RHSonly with NAs in vars_2_nodeOrigID
    if(nimbleOptions()$allowDynamicIndexing) {
        for(iV in seq_along(varInfo)) {
            varName <- varInfo[[iV]]$varName
            if(varName %in% unknownIndexNames) 
                vars_2_nodeOrigID[[ varName ]][!is.na(vars_2_nodeOrigID[[ varName ]])] <- NA
        }
    }
        
    ## 8. Collect edges from RHS vars to LHS nodes
    if(debug) browser()
    edgesFrom <- numeric(0)          ## Set up aligned vectors of edgeFrom, edgesTo, and the parentExprID of an edge, i.e. which part of a pulled apart expression is an edge from.  E.g x ~ dnorm(a, b) would make an edgeFrom entry with the ID of a, edgeTo with ID of x, and parentExprID would be 1 since the a would be the first piece of the RHS as it is pulled apart.  The next edge would from from b, to x, with parentExprID of 2.
    edgesTo <- numeric(0)
    edgesParentExprID <- numeric(0)
    for(iDI in seq_along(declInfo)) {    ## Iterate over BUGS lines
        BUGSdecl <- declInfo[[iDI]]
        if(BUGSdecl$numUnrolledNodes == 0) next
        rhsVars <- BUGSdecl$rhsVars
        for(iV in seq_along(rhsVars)) {  ## Iterate over RHS vars
            
            rhsVar <- rhsVars[iV]
            nDim <- varInfo[[rhsVar]]$nDim  ## Collect edges (noting in a case like x ~ dnorm(a, a) there could be 2 edges from a to x)
            ## Note also that all edges from unrolling any for loops are collected at once
            ## e.g. for(i in 1:10) x[i] ~ dnorm(mu, sigma) will collect all 10 edges from mu to x[1]...x[10] at once
            newEdges <- collectEdges(vars_2_vertexOrigID[[rhsVar]], BUGSdecl$unrolledIndicesMatrix, BUGSdecl$origIDs, contexts[[BUGSdecl$contextID]]$indexVarExprs,
                                      BUGSdecl$symbolicParentNodesReplaced[[iV]],
                                      BUGSdecl$parentIndexNamePieces[[iV]], BUGSdecl$replacementNameExprs)
            
            edgesFrom <- c(edgesFrom, newEdges[[1]])
            edgesTo <- c(edgesTo, newEdges[[2]])
            edgesParentExprID <- c(edgesParentExprID, rep(iV, length(newEdges[[1]])))
        }
    }

   
    ## 9. Collect edges from original nodes to inferred vertices
    ##    e.g. When x[1:2] has been fractured into x[1] and x[2], there are edges from x[1:2] to x[1] and x[2]
    ##         Note that in getDependencies("x[1]"), the result will include "x[1:2]" as the nodeFunction of "x[1]"
    if(debug) browser()
    for(iV in seq_along(varInfo)) {
        varName <- varInfo[[iV]]$varName
        if(nimbleOptions()$allowDynamicIndexing)
            if(varName %in% unknownIndexNames)  ## inferred vertices are from parent variable not unknownIndex variable
                next
        newEdges <- collectInferredVertexEdges(vars_2_nodeOrigID[[varName]], vars_2_vertexOrigID[[varName]])
        edgesFrom <- c(edgesFrom, newEdges[[1]])
        edgesTo <- c(edgesTo, newEdges[[2]])
        edgesParentExprID <- c(edgesParentExprID, rep(NA, length(newEdges[[1]])))
        if(nimbleOptions()$allowDynamicIndexing) {
            if(!varName %in% unknownIndexNames)
                vertexID_2_nodeID[newEdges[[2]]] <- newEdges[[1]]
        } else vertexID_2_nodeID[newEdges[[2]]] <- newEdges[[1]]
    }

    ## 9b. truncate vertexID_2_nodeID and relabel some vertices as LHSinferred
   if(debug) browser()
    vertexID_2_nodeID <- vertexID_2_nodeID[1:numContigVertices]
    types[ types == 'RHSonly' & vertexID_2_nodeID != 0] <- 'LHSinferred' ## The types == 'RHSonly' could be obtained more easily, since it will be a single set of FALSES followed by a single set of TRUES, but anyway, this works
    unknownIndexNodes <- types == 'unknownIndex'
    types[ types == 'unknownIndex' ] <- 'LHSinferred' ## treat unknownIndex nodes as LHSinferred so graph dependency calculations simply pass through them

    ## 9c. for RHSonly names that have a splitVertex indication (i %.s% j), convert back to colon (i:j)
    ##     because these can then be correctly used in eval(parse(...)) in one of the vars2... environments
    allVertexNames[types == 'RHSonly'] <- convertSplitVertexNamesToEvalSafeNames(allVertexNames[types == 'RHSonly'])
    
    ## 10. Build the graph for topological sorting
    if(debug) browser()
    graph <<- graph.empty()
    graph <<- add.vertices(graph, length(allVertexNames), name = allVertexNames) ## add all vertices at once
    allEdges <- as.numeric(t(cbind(edgesFrom, edgesTo)))
    graph <<- add.edges(graph, allEdges)                                         ## add all edges at once

    ## 11. Topologically sort and re-index all objects with vertex IDs
    if(debug) browser()
    newGraphID_2_oldGraphID <- as.numeric(topological.sort(graph, mode = 'out'))
    oldGraphID_2_newGraphID <- sort(newGraphID_2_oldGraphID, index = TRUE)$ix
    graph <<- permute.vertices(graph, oldGraphID_2_newGraphID)  # re-label vertices in the graph

    ## 11b. make new maps that use the sorted IDS
    if(debug) browser()
    vars_2_nodeID <- new.env()       ## this will become maps$vars2graphID_functions
    vars_2_vertexID <- new.env()     ## this will become maps$vars2graphID_values
    for(iV in seq_along(varInfo)) {  ## for each variable, populate vars_2_nodeID and vars_2_vettexID
        temp <- vars_2_nodeOrigID[[varInfo[[iV]]$varName]]
        if(!all(is.na(temp))) vars_2_nodeID[[varInfo[[iV]]$varName]] <- oldGraphID_2_newGraphID[temp] else vars_2_nodeID[[varInfo[[iV]]$varName]] <- temp
        temp <- vars_2_vertexOrigID[[varInfo[[iV]]$varName]]
        if(all(is.na(temp))) cat(paste('Something weird: all vertex IDs NA for variable', varInfo[[iV]]$varName))
        temp[temp==-1] <- NA
        vars_2_vertexID[[varInfo[[iV]]$varName]] <- oldGraphID_2_newGraphID[temp]

        if(varInfo[[iV]]$nDim > 0) dim(vars_2_nodeID[[varInfo[[iV]]$varName]]) <- dim(vars_2_vertexID[[varInfo[[iV]]$varName]]) <- dim(vars_2_vertexOrigID[[varInfo[[iV]]$varName]])        
    }

    ## 11c. re-index the graphIDs in the BUGSdecl objects and record graphID_2_declID (mapping from a graphID to its BUGS declaration ID, corresponding to declInfo order)
    graphID_2_declID <- numeric(numVertices)
    for(iDI in seq_along(declInfo)) {
        BUGSdecl <- declInfo[[iDI]]
        if(BUGSdecl$numUnrolledNodes == 0) next
        BUGSdecl$graphIDs <- oldGraphID_2_newGraphID[ BUGSdecl$origIDs ]
        graphID_2_declID[ BUGSdecl$graphIDs ] <- iDI
    }
    if(nimbleOptions()$allowDynamicIndexing) {
        ## inserted declInfo elements for unknownIndex vars will be removed later
        unknownIndexDecl <- which(sapply(declInfo, function(x) x$type == "unknownIndex"))
        graphID_2_declID[graphID_2_declID %in% unknownIndexDecl] <- 0 
    }
    
    ## 11d. set up the elementIDs
    ## These provide an ID for each scalar, but any new IDs created here are not in the graph
    ## for now we have elementIDs for all unknownIndex vars, as it's easier to generate them than not
    ## given that initial element vector is based on vertices and unknownIndex vars need to be in vertices;
    ## they presumably won't cause any issues as they shouldn't ever be used
    if(debug) browser()
    vars_2_elementID <- new.env()
    maxVertexID <- numContigVertices
    nextElementID <- numContigVertices+1
    origElementID_2_vertexID <- 1:maxVertexID
    elementNames <- character(totModelSize)
    checkingTotModelSize <- 0
    for(iV in seq_along(varInfo)) {
        vN <- varInfo[[iV]]$varName
        newSplit <- splitVertexIDsToElementIDs(vars_2_vertexID[[vN]], nextElementID)
        vars_2_elementID[[vN]] <- newSplit[[1]]
        nextElementID <- newSplit[[2]]
        origElementID_2_vertexID <-  c(origElementID_2_vertexID, newSplit[[3]])
        origElementID_2_vertexID[unique(newSplit[[3]])] <- 0 ## remove IDs that were replaced
    }
    ##  make the element IDs contiguous and match the sorting
    boolStillInUse <- origElementID_2_vertexID != 0
    sortOrder <- order(origElementID_2_vertexID[boolStillInUse])
    sortOrder2 <- order(sortOrder)
    origElemID_2_contigSortedElementID <- numeric(length(origElementID_2_vertexID))
    origElemID_2_contigSortedElementID[boolStillInUse] <- sortOrder2

    for(iV in seq_along(varInfo)) {
        vN <- varInfo[[iV]]$varName
        vars_2_elementID[[vN]] <- origElemID_2_contigSortedElementID[ vars_2_elementID[[vN]] ]
        if(varInfo[[iV]]$nDim > 0)  dim(vars_2_elementID[[vN]]) <- dim(vars_2_vertexID[[vN]]) 
    }
    elementID_2_vertexID <- 1:length(sortOrder2)
    elementID_2_vertexID[sortOrder2] <- origElementID_2_vertexID[boolStillInUse]   

    for(iV in seq_along(varInfo)) {
        ## construct element names
        vN <- varInfo[[iV]]$varName
        dims <- dim(vars_2_elementID[[vN]])
        if(is.null(dims)) dims <- length(vars_2_elementID[[vN]])
        if(length(dims)==1 & dims[1] == 1) {
            newElementNames <- vN
            checkingTotModelSize <- checkingTotModelSize + 1
            if(!is.na(vars_2_elementID[[vN]]))
            elementNames[ vars_2_elementID[[vN]] ] <- newElementNames
        } else {
            fullyIndexedText <- paste0(vN, '[', paste0(1, ':', dims, collapse = ',') , ']')
            newElementNames <- nl_expandNodeIndex(fullyIndexedText)
            checkingTotModelSize <- checkingTotModelSize + length(newElementNames)
            boolNotNA <- as.logical(!is.na(vars_2_elementID[[vN]]))
            elementNames[ as.numeric(vars_2_elementID[[vN]])[boolNotNA] ] <- newElementNames[boolNotNA]
        }

        
    }
    if(checkingTotModelSize != totModelSize) {
        stop(paste0("Something is inconsistent in the model.\n",
                    "Check for conflicting declarations.\n",
                    "If your model has macros or if-then-else blocks\n",
                    "you can inspect the processed model code by doing\n",
                    "nimbleOptions(stop_after_processing_model_code = TRUE)\n",
                    "before calling nimbleModel.\n"
                    ),
             call. = FALSE)
    }

    ## set up the vars2graphID_functions_and_RHSonly 
    ## This will be the same as vars_2_nodeID but with NAs filled in from vars_2_vertexID
    if(debug) browser()
    vars_2_nodeID_noNAs <- new.env()
    for(iV in seq_along(varInfo)) {
        vN <- varInfo[[iV]]$varName
        vars_2_nodeID_noNAs[[vN]] <- vars_2_nodeID[[vN]]
        boolNA <- is.na(vars_2_nodeID_noNAs[[vN]])
        if(any(boolNA))
            vars_2_nodeID_noNAs[[vN]][boolNA] <-  vars_2_vertexID[[vN]][boolNA] 
    }

    ## 12. Set up things needed for maps.
    if(debug) browser()
    maps <<- mapsClass$new()
    maps$elementID_2_vertexID <<- elementID_2_vertexID
    maps$vars2ID_elements <<- vars_2_elementID
    maps$elementNames <<- elementNames
    maps$vars2GraphID_functions_and_RHSonly <<- vars_2_nodeID_noNAs
    maps$graphID_2_declID <<- graphID_2_declID
    maps$graphID_2_nodeName <<- allVertexNames[newGraphID_2_oldGraphID]
    maps$types <<- types[newGraphID_2_oldGraphID]
    maps$notStoch <<- maps$types != 'stoch'
    maps$nodeNamesLHSall <<- nodeNamesLHSall
    maps$nodeNamesRHSonly <<- maps$graphID_2_nodeName[maps$types == 'RHSonly'] ##nodeNamesRHSonly
    maps$nodeNames <<- maps$graphID_2_nodeName

    if(any(duplicated(maps$nodeNames[!unknownIndexNodes[newGraphID_2_oldGraphID]]))) {  ## x[k[i],block[i]] can lead to duplicated nodeNames for unknownIndex declarations; this should be ok, though there is inefficiency in having a vertex in the graph for each element of second index instead of collapsing into one vertex per unique value.
        stop(
            paste0("There are multiple definitions for nodes:",
                   paste(maps$nodeNames[duplicated(maps$nodeNames[!unknownIndexNodes[newGraphID_2_oldGraphID]])],
                         collapse = ','), "\n",
                   "If your model has macros or if-then-else blocks\n",
                   "you can inspect the processed model code by doing\n",
                   "nimbleOptions(stop_after_processing_model_code = TRUE)\n",
                   "before calling nimbleModel.\n"
                   ),
            call. = FALSE)
    }
    if(debug) browser()
    
    newVertexID_2_nodeID <- vertexID_2_nodeID [ newGraphID_2_oldGraphID ]
    bool <- newVertexID_2_nodeID != 0
    newVertexID_2_nodeID[bool] <- oldGraphID_2_newGraphID[newVertexID_2_nodeID]
    maps$vertexID_2_nodeID <<- newVertexID_2_nodeID
    
    maps$graphID_2_nodeFunctionName <<- maps$graphID_2_nodeName
    maps$graphID_2_nodeFunctionName[bool] <<- maps$graphID_2_nodeName[ newVertexID_2_nodeID ]

    if(debug) browser()
    maps$vars2GraphID_values <<- vars_2_vertexID
    maps$vars2GraphID_functions <<- vars_2_nodeID

    if(debug) browser()

    maps$vars2LogProbName <<- vars2LogProbName
##    maps$vars2LogProbID <<- vars2LogProbID
##    maps$logProbIDs_2_LogProbName <<- logProbNames

    graphID_2_logProbName <- paste0("logProb_", gsub(":[0123456789]+", "", maps$graphID_2_nodeName))
    graphID_2_logProbName[ maps$types != "stoch"] <- NA
    maps$graphID_2_logProbName <<- graphID_2_logProbName
    
    maps$edgesFrom <<- oldGraphID_2_newGraphID[edgesFrom]
    maps$edgesTo <<- oldGraphID_2_newGraphID[edgesTo]
    maps$edgesParentExprID <<- edgesParentExprID
    edgesLevels <- if(length(maps$vertexID_2_nodeID) > 0) 1:length(maps$vertexID_2_nodeID) else numeric(0)
    ##edgesLevels <- if(length(maps$edgesFrom) > 0) 1:max(max(maps$edgesFrom), max(maps$edgesTo)) else numeric(0)
    fedgesFrom <- factor(maps$edgesFrom, levels = edgesLevels) ## setting levels ensures blanks inserted into the splits correctly
    maps$edgesFrom2To <<- split(maps$edgesTo, fedgesFrom)
    maps$edgesFrom2ParentExprID <<- split(maps$edgesParentExprID, fedgesFrom)
    maps$graphIDs <<- 1:length(maps$graphID_2_nodeName)

    maps$nimbleGraph <<- nimbleGraphClass()
    maps$nimbleGraph$setGraph(maps$edgesFrom, maps$edgesTo, maps$edgesParentExprID, maps$vertexID_2_nodeID, maps$types, maps$graphID_2_nodeName, length(maps$graphID_2_nodeName))
##    maps$nodeName_2_graphID <<- list2env( nodeName2GraphIDs(maps$nodeNames) )
##    maps$nodeName_2_logProbName <<- list2env( nodeName2LogProbName(maps$nodeNames) )

    ## A new need for new node function system:
    graphID_2_unrolledIndicesMatrixRow <- rep(-1L, (length(maps$graphIDs)))
    for(iDI in seq_along(declInfo)) {
        BUGSdecl <- declInfo[[iDI]]
        unrolledRows <- nrow(BUGSdecl$unrolledIndicesMatrix)
        if(unrolledRows == 0) {
            if(BUGSdecl$numUnrolledNodes == 1) ## a singleton declaration
                graphID_2_unrolledIndicesMatrixRow[BUGSdecl$graphIDs] <- 0
            else
                stop(paste('confused assigning unrolledIndicesMatrixRow in case with no unrolledRows by numUnrolledNodes != 1 for code', deparse(BUGSdecl$code)))
        } else {
            theseGraphIDs <- BUGSdecl$graphIDs
            graphID_2_unrolledIndicesMatrixRow[theseGraphIDs] <- 1:length(theseGraphIDs)
        }
    }
    graphID_2_unrolledIndicesMatrixRow[graphID_2_unrolledIndicesMatrixRow==-1] <- NA
    maps$graphID_2_unrolledIndicesMatrixRow <<- graphID_2_unrolledIndicesMatrixRow
    NULL
})


modelDefClass$methods(genVarInfo3 = function() {
    ## First set up varInfo's for all LHS variables and collect anyStoch.
    ## That allows determination of when logProb information needs to be collected
    unknownIndexNames <<- NULL

    for(iDI in seq_along(declInfo)) {
        BUGSdecl <- declInfo[[iDI]]
        if(nimbleOptions()$allowDynamicIndexing)
            if(BUGSdecl$type == "unknownIndex") next  # handled in addUnknownIndexVars to avoid all the processing here
        if(BUGSdecl$numUnrolledNodes == 0) next
        ## LHS:
        lhsVar <- BUGSdecl$targetVarName
        if(!(lhsVar %in% names(varInfo))) {
            nDim <- if(length(BUGSdecl$targetNodeExpr)==1) 0 else length(BUGSdecl$targetNodeExpr)-2
            varInfo[[lhsVar]] <<- varInfoClass$new(varName = lhsVar,
                                                   mins = rep(10000000, nDim),
                                                   maxs = rep(0, nDim),
                                                   nDim = nDim,
                                                   anyStoch = FALSE,
                                                   anyDynamicallyIndexed = FALSE)
        }
        varInfo[[lhsVar]]$anyStoch <<- varInfo[[lhsVar]]$anyStoch | (BUGSdecl$type == 'stoch')
    }

    anyStoch = unlist(lapply(varInfo, `[[`, 'anyStoch'))
    logProbVarInfo <<- lapply(varInfo[anyStoch], function(x)
        varInfoClass$new(varName = makeLogProbName(x$varName),
                         mins = rep(10000000, x$nDim),
                         maxs = rep(0,      x$nDim),
                         nDim = x$nDim,
                         anyStoch = FALSE))
    names(logProbVarInfo) <<- lapply(logProbVarInfo, `[[`, 'varName')
    
    for(iDI in seq_along(declInfo)) {
        BUGSdecl <- declInfo[[iDI]]
        if(nimbleOptions()$allowDynamicIndexing)
            if(BUGSdecl$type == "unknownIndex") next  # handled in addUnknownIndexVars to avoid all the processing here
        if(BUGSdecl$numUnrolledNodes == 0) next
        ## LHS:
        lhsVar <- BUGSdecl$targetVarName
        anyStoch <- varInfo[[lhsVar]]$anyStoch
        if(anyStoch) lhsLogProbVar <- makeLogProbName(lhsVar)
        if(varInfo[[lhsVar]]$nDim > 0) {
            for(iDim in 1:varInfo[[lhsVar]]$nDim) {
                indexNamePieces <- BUGSdecl$targetIndexNamePieces[[iDim]] 
                if(is.list(indexNamePieces)) { ## a list would be made if there is a ':' operator in the index expression
                    indsLow <- if(is.numeric(indexNamePieces[[1]])) indexNamePieces[[1]] else BUGSdecl$replacementsEnv[[ indexNamePieces[[1]] ]]
                    indsHigh <- if(is.numeric(indexNamePieces[[2]])) indexNamePieces[[2]] else BUGSdecl$replacementsEnv[[ indexNamePieces[[2]] ]]
                    rangeIndsLow <- range(indsLow)
                    varInfo[[lhsVar]]$mins[iDim] <<- min(varInfo[[lhsVar]]$mins[iDim], rangeIndsLow[1])
                    varInfo[[lhsVar]]$maxs[iDim] <<- max(varInfo[[lhsVar]]$maxs[iDim], max(indsHigh))

                    if(anyStoch) {
                        logProbVarInfo[[lhsLogProbVar]]$mins[iDim] <<- min(logProbVarInfo[[lhsLogProbVar]]$mins[iDim], rangeIndsLow[1])
                        logProbVarInfo[[lhsLogProbVar]]$maxs[iDim] <<- max(logProbVarInfo[[lhsLogProbVar]]$maxs[iDim], rangeIndsLow[2]) ## This collapses i:j to be i for logProb purposes because a multivariate node needs only a single logProb value
                    }
                } else { 
                    inds <- if(is.numeric(indexNamePieces)) indexNamePieces else BUGSdecl$replacementsEnv[[ indexNamePieces ]]
                    rangeInds <- range(inds)
                    varInfo[[lhsVar]]$mins[iDim] <<- min(varInfo[[lhsVar]]$mins[iDim], rangeInds[1])
                    varInfo[[lhsVar]]$maxs[iDim] <<- max(varInfo[[lhsVar]]$maxs[iDim], rangeInds[2])

                    if(anyStoch) {
                        ## Do it independently because RHS declarations can pick up different min/max info for the vars that is not needed for the logProbVars
                        logProbVarInfo[[lhsLogProbVar]]$mins[iDim] <<- min(logProbVarInfo[[lhsLogProbVar]]$mins[iDim], rangeInds[1])
                        logProbVarInfo[[lhsLogProbVar]]$maxs[iDim] <<- max(logProbVarInfo[[lhsLogProbVar]]$maxs[iDim], rangeInds[2])
                    }
                }
            }
        }

        ## RHS:
        rhsVars <- BUGSdecl$rhsVars

        for(iV in seq_along(rhsVars)) {
            rhsVar <- rhsVars[iV]
            if(nimbleOptions()$allowDynamicIndexing)
                rhsVar <- stripUnknownIndexFromVarName(rhsVar) # symbolicParentNodes already expressed in terms of unknownIndex entity but we want the inferred variable to be the original one
            if(!(rhsVar %in% names(varInfo))) {
                if(!nimbleOptions()$allowDynamicIndexing) {
                    nDim <- if(length(BUGSdecl$symbolicParentNodes[[iV]])==1)
                                0
                            else 
                                length(BUGSdecl$symbolicParentNodes[[iV]])-2
                } else {
                    tmp <- stripIndexWrapping(BUGSdecl$symbolicParentNodes[[iV]])
                    nDim <- if(length(tmp)==1) 0 else length(tmp)-2
                    
                }
                varInfo[[rhsVar]] <<- varInfoClass$new(varName = rhsVar,
                                                       mins = rep(100000, nDim),
                                                       maxs = rep(0, nDim),
                                                       nDim = nDim,
                                                       anyStoch = FALSE,
                                                       anyDynamicallyIndexed = FALSE)
            }
            if(varInfo[[rhsVar]]$nDim != length(BUGSdecl$parentIndexNamePieces[[iV]]))
                stop("Dimension of ", rhsVar, " is ", varInfo[[rhsVar]]$nDim, ", which does not match its usage in '", deparse(BUGSdecl$code), "'.")
            if(varInfo[[rhsVar]]$nDim > 0) {
                for(iDim in 1:varInfo[[rhsVar]]$nDim) {
                    indexNamePieces <- BUGSdecl$parentIndexNamePieces[[iV]][[iDim]]
                    if(is.null(indexNamePieces)) stop(paste0('There is a problem with some indexing in this line: ', deparse(BUGSdecl$codeReplaced), '.\nOne way this can happen is if you wanted to provide a vector of indices but did not include it in constants.'))
                    if(is.list(indexNamePieces)) { ## a list would be made if there is a ':' operator in the index expression
                        indsLow <- if(is.numeric(indexNamePieces[[1]])) indexNamePieces[[1]] else BUGSdecl$replacementsEnv[[ indexNamePieces[[1]] ]]
                        indsHigh <- if(is.numeric(indexNamePieces[[2]])) indexNamePieces[[2]] else BUGSdecl$replacementsEnv[[ indexNamePieces[[2]] ]]
                        varInfo[[rhsVar]]$mins[iDim] <<- min(varInfo[[rhsVar]]$mins[iDim], min(indsLow))
                        varInfo[[rhsVar]]$maxs[iDim] <<- max(varInfo[[rhsVar]]$maxs[iDim], max(indsHigh))
                    } else {
                        ## If the index is dynamic (marked by NA), there is nothing to learn about index range of the variable.
                        if(nimbleOptions()$allowDynamicIndexing)
                            if(isDynamicIndex(indexNamePieces)) {
                                varInfo[[rhsVar]]$mins[iDim] <<- min(varInfo[[rhsVar]]$mins[iDim], 1) ## o.w., never changed from 1e5 if only on RHS and in 'dimensions' input
                                varInfo[[rhsVar]]$maxs[iDim] <<- max(varInfo[[rhsVar]]$maxs[iDim], 1) ## o.w., can end up with (1,0) as (min,max) before 'dimensions' are used
                                next
                            }
                        ## Otherwise extend the range of known mins and maxs based on this expression
                        inds <- if(is.numeric(indexNamePieces)) indexNamePieces else BUGSdecl$replacementsEnv[[ indexNamePieces ]]
                        rangeInds <- range(inds)
                        varInfo[[rhsVar]]$mins[iDim] <<- min(varInfo[[rhsVar]]$mins[iDim], rangeInds[1])
                        varInfo[[rhsVar]]$maxs[iDim] <<- max(varInfo[[rhsVar]]$maxs[iDim], rangeInds[2])
                    }
                }
            }
        }
    }
    
    ## now use dimensionsList, to check / update varInfo
    for(i in seq_along(dimensionsList)) {
        dimVarName <- names(dimensionsList)[i]
        if(!(dimVarName %in% names(varInfo))) next
        if(length(dimensionsList[[dimVarName]]) != varInfo[[dimVarName]]$nDim)   stop('inconsistent dimensions for variable ', dimVarName)
        if(any(dimensionsList[[dimVarName]] < varInfo[[dimVarName]]$maxs))  stop(paste0('dimensions specified are smaller than model specification for variable \'', dimVarName, '\''))
        varInfo[[dimVarName]]$maxs <<- dimensionsList[[dimVarName]]
    }

    ## check for maxs < mins; this would generally be from a BUGS syntax error,
    ## e.g., for(i in 1:4) y[k] ~ dnorm(0,1);
    ## in some cases these would be caught by the check for mins or maxs zero or less
    ## but this error message is more informative
    if(any(sapply(varInfo, function(x) length(x$mins) && length(x$maxs) &&
                                       any(x$mins > x$maxs)))) {
        problemVars <- which(sapply(varInfo, function(x) any(x$mins > x$maxs)))
        stop("genVarInfo3: indexing error found for model variable(s): ",
             paste(names(varInfo)[problemVars], "; please check that variables used for indexing are properly defined in the relevant for loop(s)", collapse = ' '))
    }

    ## check for mins or maxs zero or less (these trigger various errors including R crashes)
    if(any(sapply(varInfo, function(x) length(x$mins) && min(x$mins) < 1)) ||
       any(sapply(varInfo, function(x) length(x$maxs) && min(x$maxs) < 1))) {
           problemVars <- c(which(sapply(varInfo, function(x) min(x$mins)) < 1),
                         which(sapply(varInfo, function(x) min(x$maxs)) < 1))
           stop("genVarInfo3: index value of zero or less found for model variable(s): ",
                paste(names(varInfo)[problemVars], collapse = ' '))
    }
})

modelDefClass$methods(addUnknownIndexVars = function(debug = FALSE) {
    if(debug) browser()
    unknownIndexNames <<- NULL
    if(nimbleOptions()$allowDynamicIndexing) 
        for(iDI in seq_along(declInfo)) {
            BUGSdecl <- declInfo[[iDI]]
            if(BUGSdecl$type != 'unknownIndex') next
            lhsVar <- BUGSdecl$targetVarName
            if(!(lhsVar %in% names(varInfo))) {
                if(length(BUGSdecl$rhsVars) > 1)
                    stop("addUnknownIndexVars: more than one right-hand side variable in unknownIndex declaration: ",
                         deparse(BUGSdecl$code))
                varInfo[[lhsVar]] <<- varInfo[[BUGSdecl$rhsVars]]$copy()
                varInfo[[lhsVar]]$varName <<- lhsVar
                varInfo[[lhsVar]]$anyStoch <<- FALSE
                unknownIndexNames <<- c(unknownIndexNames, lhsVar)
                varInfo[[BUGSdecl$rhsVars]]$anyDynamicallyIndexed <<- TRUE
            } else stop("addUnknownIndexVars: ", lhsVar, " already present in varInfo. This code should not have been triggered.")
        }
})


## This removes temporary declarations and vars created because of dynamic indexing.
modelDefClass$methods(stripUnknownIndexInfo = function() {
    if(nimbleOptions()$allowDynamicIndexing) {
        declInfo[sapply(declInfo, function(x) x$type == 'unknownIndex')] <<- NULL
        sapply(unknownIndexNames, function(x) varInfo[[x]] <<- NULL) # At one point, this was needed for isData stuff since we have unknownIndex vars as part of graph, but doesn't seem to be needed anymore.
    }
})

modelDefClass$methods(genUnknownIndexDeclarations = function() {
    if(nimbleOptions()$allowDynamicIndexing) {
        nimFunNames <- getAllDistributionsInfo('namesExprList')
        for(i in seq_along(declInfo)){
            for(p in seq_along(declInfo[[i]]$symbolicParentNodes)) {
                parentExpr <- stripIndexWrapping(declInfo[[i]]$symbolicParentNodes[[p]])
                dynamicIndices <- detectDynamicIndexes(parentExpr)
                if(sum(dynamicIndices) && !any(sapply(declInfo, function(x) identical(x$targetExpr, parentExpr)))) {
                    ## don't create declaration if there is already one for the exact same unknownIndex target
                    BUGSdeclClassObject <- BUGSdeclClass$new()
                    lhsCode <- parentExpr
                    rhsCode <- lhsCode
                    rhsCode <- stripUnknownIndexFromVarNameInBracketExpr(rhsCode)
                    newCode <- substitute(LHS <- RHS, list(LHS = lhsCode, RHS = rhsCode))
                    BUGSdeclClassObject$setup(newCode, declInfo[[i]]$contextID, declInfo[[i]]$sourceLineNumber)
                    BUGSdeclClassObject$setIndexVariableExprs(contexts[[declInfo[[i]]$contextID]]$indexVarExprs)
                    BUGSdeclClassObject$genSymbolicParentNodes(constantsNamesList, contexts[[declInfo[[i]]$contextID]], nimFunNames,
                                                              contextID = declInfo[[i]]$contextID)
                    ##                    BUGSdeclClassObject$genReplacementsAndCodeReplaced(constantsNamesList, contexts[[declInfo[[i]]$contextID]], nimFunNames)
                    ##                    BUGSdeclClassObject$genReplacedTargetValueAndParentInfo(constantsNamesList, contexts[[declInfo[[i]]$contextID]], nimFunNames)
                    ## insert info that is usually done in genNodeInfo3; passing through genNodeInfo3 would be complicated as that operates by context
                    ##BUGSdeclClassObject$unrolledIndicesMatrix <- declInfo[[i]]$unrolledIndicesMatrix
                    ##BUGSdeclClassObject$replacementsEnv <- declInfo[[i]]$replacementsEnv
                                        #BUGSdeclClassObject$outputSize <- 1 # not sure if needed
                    ##BUGSdeclClassObject$numUnrolledNodes <- 1 # might be more, but I think the only relevant distinction in terms of how this is used is 0 vs. 1
                    BUGSdeclClassObject$type <- "unknownIndex"
                    declInfo[[length(declInfo)+1]] <<- BUGSdeclClassObject
                }
            }
        }
    }
})

modelDefClass$methods(genIsDataVarInfo = function() {
    ## uses varInfo to set field: isDataVarInfo
    isDataVarInfo <<- lapply(varInfo, function(x)
        varInfoClass$new(
            varName = x$varName,
            mins = x$mins,
            maxs = x$maxs,
            nDim = x$nDim,
            anyStoch = FALSE))
    names(isDataVarInfo) <<- lapply(isDataVarInfo, `[[`, 'varName')
})

modelDefClass$methods(genVarNames = function() {
    ## uses varInfo and logProbVarInfo to set field: varNames
    varNames <<- c(names(varInfo), names(logProbVarInfo))
})

modelDefClass$methods(buildSymbolTable = function() {
    ## uses varInfo and logProbVarInfo to set field: symTab
    
    st <- symbolTable()
    
    for(vI in c(varInfo)) {
        st$addSymbol(symbolDouble(name = vI$varName, nDim = vI$nDim, size = vI$maxs))
    }

    for(vI in c(logProbVarInfo)) {
        st$addSymbol(symbolDouble(name = vI$varName, nDim = vI$nDim, size = vI$maxs))
    }
    
    symTab <<- st
})


modelDefClass$methods(newModel = function(data = list(), inits = list(), where = globalenv(), modelName = character(), check = getNimbleOption('checkModel'), calculate = TRUE, debug = FALSE) {
    if(debug) browser()
    if(inherits(modelClass, 'uninitializedField')) {
        vars <- lapply(varInfo, `[[`, 'maxs')
        logProbVars <- lapply(logProbVarInfo, `[[`, 'maxs')
        isDataVars <- lapply(isDataVarInfo, `[[`, 'maxs')
        modelClass <<- RMakeCustomModelClass(vars = c(vars, logProbVars),
                                             className = paste0(name,'_modelClass_', nimbleUniqueID()),
                                             isDataVars = isDataVars,
                                             modelDef = .self,
                                             where = where)
        classEnvironment <<- where
    }
    if(inherits(modelValuesClass, 'uninitializedField')) {
        modelValuesClass <<- makeCustomModelValuesClass(symTab, modelValuesClassName, where = where, addUniqueID = FALSE, modelDef = .self)
    }
    ## would be better to build this into the class
    ## if(missing(modelName)) modelName <- name
    ## if(modelName == character(0)) stop('Error, empty name for a new model', call. = FALSE)
    model <- modelClass(name = modelName)
    model$setGraph(graph)
    model$buildNodeFunctions(where = where, debug = debug)
    model$buildNodesList() ## This step makes RStudio choke, we think from circular reference classes -- fixed, by not displaying Global Environment in RStudio

    ## handling for JAGS style inits (a list of lists)
    ## added Oct 2015, DT
    if(length(inits) > 0 && is.list(inits[[1]])) {
        message('detected JAGS style initial values, provided as a list of lists...  using the first set of initial values')
        inits <- inits[[1]]
    }
    
    if(length(data) + length(inits) > 0)
        if(nimbleOptions('verbose')) message("setting data and initial values...")
    model$setData(data)
    # prevent overwriting of data values by inits
    if(FALSE) {  # should now be handled by checking if setInits tries to overwrite data nodes
        for(varName in intersect(names(inits), model$getVarNames())) {
            dataVars <- model$isData(varName)
            if(sum(dataVars) && !identical(data[[varName]][dataVars],
                                           inits[[varName]][dataVars])) {
                                        # only warn if user passed conflicting actual values
                nonNAinits <- !is.na(inits[[varName]][dataVars])
                if(!identical(data[[varName]][dataVars][nonNAinits],
                              inits[[varName]][dataVars][nonNAinits]))
                    warning("newModel: Conflict between 'data' and 'inits' for ", varName, "; using values from 'data'.\n")
                inits[[varName]][dataVars] <- data[[varName]][dataVars]
            }
        }
    }
    nonVarIndices <- !names(inits) %in% model$getVarNames()
    if(sum(nonVarIndices))
        warning("newModel: ", paste(names(inits)[nonVarIndices], collapse = ','),
                " ", ifelse(sum(nonVarIndices) > 1, "are", "is"), " not ", ifelse(sum(nonVarIndices) > 1, "variables", "a variable"), " in the model; initial value ignored.")
    model$setInits(inits[!nonVarIndices])
    ## basic size/dimension, NA checking
    if(calculate) {
        if(nimbleOptions('verbose')) message("running calculate on model (any error reports that follow may simply reflect missing values in model variables) ... ", appendLF = FALSE)
        result <- try(model$calculate(), silent = TRUE)
        if(nimbleOptions('verbose')) 
            if(is(result, 'try-error')) 
                message(geterrmessage()) else message("")  # this ensures a single newline is included
    }
    if(nimbleOptions('verbose')) message("checking model sizes and dimensions...", appendLF = FALSE)
    model$checkBasics()
    if(nimbleOptions('verbose')) message("")  # appends newline   
    ## extended model checking via calculate; disabled by default as of July 2016
    if(check) {
        if(nimbleOptions('verbose')) message("checking model calculations...")
        model$check()
    }
    fixRStudioHanging(model)
    return(model)
})

modelDefClass$methods(fixRStudioHanging = function(model) {
    ## hopefully the *final* work needed towards fixing the RStudio "hanging" problem. . . .
    ## added by DT, Nov 2015
    nullStrMethod <- function(object, ...) str(NULL)
    classNames <- c(as.character(class(model)),
                    as.character(class(model$defaultModelValues)),
                    unique(sapply(model$nodeFunctions, function(nf) as.character(class(nf)))))
    for(name in c(classNames, "modelDefClass", "igraph")) {
        eval(substitute(NAME <- METHOD,
                        list(NAME   = as.name(paste0('str.', name)),
                             METHOD = nullStrMethod)),
             envir = globalenv())
    }
})

modelDefClass$methods(show = function() {
    writeLines("modelDefClass object with fields and methods:")
    writeLines(ls(.self))
})

modelDefClass$methods(printDI = function() {
    for(i in seq_along(declInfo)) {
        BUGSdecl <- declInfo[[i]]
        cat(paste0('[[', i, ']]  '))
        lapply(contexts[[BUGSdecl$contextID]]$singleContexts, function(x) cat(paste0('for(', x$indexVarExpr, ' in ', deparse(x$indexRangeExpr), ') {{{   ')))
        cat(paste0(deparse(BUGSdecl$code)))
        cat(paste0(rep('   }}}', length(contexts[[BUGSdecl$contextID]]$singleContexts)), collapse=''))
        cat('\n')
    }
})

modelDefClass$methods(graphIDs2indexedNodeInfo = function(graphIDs) {
    declIDs <- maps$graphID_2_declID[graphIDs]
    rowIndices <- maps$graphID_2_unrolledIndicesMatrixRow[graphIDs]
    ## populateNodeFxnVectorNew_copyFromRobject relies on the following order (not names)
    list(declIDs = as.integer(declIDs), unrolledIndicesMatrixRows = as.integer(rowIndices))
})

modelDefClass$methods(nodeName2GraphIDs = function(nodeName, nodeFunctionID = TRUE, unique = TRUE){
    if(length(nodeName) == 0)
        return(NULL)
    ## If unique is FALSE, we still use unique for each element of nodeName
    ## but we allow non-uniqueness across elements in the result
    if(nodeFunctionID) {
        if(unique)
            output2 <- unique(parseEvalNumericMany(nodeName, env = maps$vars2GraphID_functions_and_RHSonly))
        else
            output2 <- unlist(lapply(parseEvalNumericManyList(nodeName, env = maps$vars2GraphID_functions_and_RHSonly), unique))
    } else {
        output2 <- unique(parseEvalNumericMany(nodeName, env = maps$vars2ID_elements))
    }
    output <- output2
    return(output[!is.na(output)])
})

## next two functions work for properly formed nodeNames.
modelDefClass$methods(nodeName2LogProbName = function(nodeName){ ## used in 3 places: MCMC_build, valuesAccessorVector, and cppInterfaces_models
    ## This function needs better processing.
    if(length(nodeName) == 0)
        return(NULL)
    
##     ## 1. so this needs to first get to a nodeFunctionID
##     graphIDs <- unique(unlist(sapply(nodeName, parseEvalNumeric, env = maps$vars2GraphID_functions, USE.NAMES = FALSE)))
##     ## 2. get node function names
##     fullNodeNames <- maps$graphID_2_nodeName[graphIDs]
##     ## 3 get corresponding logProbNames
##     output <- unique(unlist(sapply(fullNodeNames, parseEvalCharacter, env = maps$vars2LogProbName, USE.NAMES = FALSE)))

##     ##graphIDs2 <- unique(parseEvalNumericMany(nodeName, env = maps$vars2GraphID_functions))
##     ##fullNodeNames2 <- maps$graphID_2_nodeName[graphIDs2]
##     ##output2 <- unique(parseEvalCharacterMany(fullNodeNames2, env = maps$vars2LogProbName))
## ##    output2 <- unique(parseEvalCharacterMany(nodeName, env = maps$vars2LogProbName)) ## output2
## ##    if(!identical(output[!is.na(output)], as.character(output2[!is.na(output2)]))) browser()
##     output <- output[!is.na(output)]

    graphIDs2 <- unique(parseEvalNumericMany(nodeName, env =  maps$vars2GraphID_functions))
    output2 <- maps$graphID_2_logProbName[graphIDs2]
    output2 <- output2[!is.na(output2)]

    ## fullNodeNames2 <- maps$graphID_2_nodeName[graphIDs [maps$types[graphIDs] == "stoch"] ]
    ## output2 <- gsub(":[0123456789]+", "", fullNodeNames2 )
    ## output2 <- output2[!is.na(output2)]
    ## output2 <- paste0("logProb_", output2)
##    if(!identical(output, output2)) browser()

##    output <- output2
    return(output2)
##    return(output[!is.na(output)])
})

## modelDefClass$methods(nodeName2LogProbID = function(nodeName){ ## used only in cppInterfaces_models
##     ## I think this will only work if nodeName is already ensured to be a node function name
## 	if(length(nodeName) == 0)
## 		return(NULL)
## 	output <- unique(unlist(sapply(nodeName, parseEvalNumeric, env = maps$vars2LogProbID, USE.NAMES = FALSE) ) ) 
## 	return(output[!is.na(output)])
## })

parseEvalNumeric <- function(x, env){
    ans <- eval(parse(text = x, keep.source = FALSE)[[1]], envir = env)
    as.numeric(ans)
}

parseEvalNumericManyFindErrors <- function(x, env) {
    problems <- list()
    for(thisx in x) {
        oneResult <- try(parseEvalNumeric(thisx, env), silent = TRUE)
        if(inherits(oneResult, 'try-error')) {
            problems[[ length(problems) + 1]] <- oneResult[1]
            if(length(problems) >= 10)
                return(problems)
        }
    }
    return(problems)
}

parseEvalNumericManyHandleError <- function(cond, x, env) {
    problems <- parseEvalNumericManyFindErrors(x, env)
    if(length(problems)==0) message(paste0('There an unknown problem looking for variables ', paste0(x, collapse=','), ' in the model.\n'))
    else {
        message(paste0('One or more errors occurred looking for variables in a model (first 10 shown below).\n',
                       'These messages may be cryptic, but generally the variable or expression somewhere in each message was not valid in a model:\n',
                       paste0(unlist(problems), collapse = ''))) 
    }
    invokeRestart('abort')
}

parseEvalNumericMany <- function(x, env) {
    withCallingHandlers(
        if(length(x) > 1) {
            as.numeric(eval(parse(text = paste0('c(', paste0(x, collapse=','),')'), keep.source = FALSE)[[1]], envir = env))
        } else 
            as.numeric(eval(parse(text = x, keep.source = FALSE)[[1]], envir = env))
       ,
        error = function(cond) {
           parseEvalNumericManyHandleError(cond, x, env)
        }
    )
}


parseEvalNumericManyList <- function(x, env) {
    withCallingHandlers(
        eval(.Call(makeParsedVarList, x), envir = env)
        ## Above line replaces:
        ## if(length(x) > 1) {
        ##     eval(parse(text = paste0('list(', paste0("as.numeric(",x,")", collapse=','),')'), keep.source = FALSE)[[1]], envir = env)
        ## } else 
        ##     eval(parse(text = paste0('list(as.numeric(',x,'))'), keep.source = FALSE)[[1]], envir = env)
       ,
        error = function(cond) {
            parseEvalNumericManyHandleError(cond, x, env)
        }
    )
}

parseEvalCharacter <- function(x, env){
    ans <- eval(parse(text = x, keep.source = FALSE)[[1]], envir = env)
    as.character(ans)
}

parseEvalCharacterMany <- function(x, env){
    if(length(x) > 1) {
        return(as.character(eval(parse(text = paste0('c(', paste0(x, collapse=','),')'), keep.source = FALSE)[[1]], envir = env)))
    } else 
        as.character(eval(parse(text = x, keep.source = FALSE)[[1]], envir = env))
}

getDependencyPaths <- function(nodeID, maps, nodeIDrow = NULL) {
    newNodes <- maps$edgesFrom2To[[nodeID]]
    newPEIDs <- maps$edgesFrom2ParentExprID[[nodeID]]
    if(length(newNodes) > 0) {
        if(is.null(nodeIDrow)) nodeIDrow <- c(nodeID, NA)
        nodeAndPEID_list <- split(cbind(newNodes, newPEIDs, deparse.level = 0), seq_along(newNodes))
        ans <- do.call('c', lapply(nodeAndPEID_list, 
               function(x) {
                   if(maps$notStoch[x[1]]) ## not stochastic so recurse
                       ans2 <- lapply(getDependencyPaths(x[1], maps = maps, nodeIDrow = x),
                                      function(z) rbind(nodeIDrow, z, deparse.level = 0))
                   else  ## It is stochastic so append x and terminate
                       ans2 <- list(rbind(nodeIDrow, x, deparse.level = 0))
                   ans2
               }))
        ans
    } else
        NULL
}

stripUnknownIndexFromVarName <- function(varName) {
    if(length(grep("^\\..+_unknownIndex.*", varName))) {
        tmp <- sub("_unknownIndex.*", "", varName)
        return(sub("^\\.", "", tmp))
    } else return(varName)
}

stripUnknownIndexFromVarNameInBracketExpr <- function(parentExpr) {
    parentExpr[[2]] <- as.name(stripUnknownIndexFromVarName(parentExpr[[2]]))
    return(parentExpr)
}

addUnknownIndexToVarName <- function(varName, extraText, contextID = NULL)
    return(paste0(".", varName, "_unknownIndex_", extraText,
                  ifelse(is.null(contextID), "", paste0("_context", contextID))))

addUnknownIndexToVarNameInBracketExpr <- function(parentExpr, contextID = NULL) {
    parentExpr[[2]] <- as.name(addUnknownIndexToVarName(parentExpr[[2]], Rname2CppName(parentExpr), contextID))
    return(parentExpr)
}

detectDynamicIndexes <- function(expr) {
    if(length(expr) == 1 || expr[[1]] != "[") return(FALSE) # stop("whichDynamicIndices: 'expr' should be a bracket expression")
    return(sapply(expr[3:length(expr)], isDynamicIndex)) 
}

modelDefClass$methods(checkForSelfParents = function(){
    if(any(maps$edgesFrom == maps$edgesTo)) {
        problemNodes <- maps$edgesFrom[maps$edgesFrom == maps$edgesTo]
        stop(paste("In building model, each of the following nodes",
                   "has itself as a parent node:",
                   paste(maps$graphID_2_nodeName[problemNodes], collapse = ", ")),
             call. = FALSE)
    }
})

