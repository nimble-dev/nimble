## A small class for information deduced about a variable in a BUGS model
varInfoClass <- setRefClass('varInfoClass',
                            fields = list(
                                varName = 'ANY',
                                mins = 'ANY',
                                maxs = 'ANY',
                                nDim = 'ANY',
                                anyStoch = 'ANY'))

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
                                 dimensionsList = 'ANY',		#list		   ## list of provided dimension information, set in assignDimensions()
                                 contexts = 'ANY',				#list 			 ## list of BUGScontextClass objects
                                 declInfo = 'ANY',				#list				 ## list of BUGSdeclInfo objects
                                 nodeInfo = 'ANY', 			#list				## list of nodeInfoClass objects, set in genNodeInfo()
                                 varInfo = 'ANY',			#list	  ## list of varInfoClass objects, set in genVarInfo()
                                 logProbVarInfo = 'ANY',	#list	  ## list of varInfoClass objects, set in genLogProbVarInfo()
                                 isDataVarInfo = 'ANY', 	#list		## list of varInfoClass objects, set in genIsDataVarInfo()
                                 varNames = 'ANY',  ## vector of all model variable names, set in genVarNames()
                                 symTab = 'ANY',  ## symbolTable object, set in buildSymbolTable()
                                 graph = 'ANY',     ## igraph object, set in buildIgraph()
                                 graphNodesList = 'ANY',   ## list of graphNode objects, set in genGraphNodesList()
                                 maps = 'ANY',   ## object of mapsClass, set in buildMaps()
                                 numNodeFunctions = 'ANY',
                                 
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
                                 	nodeInfo <<- list()
                                 	varInfo <<- list()
                                 	logProbVarInfo <<- list()
                                 	isDataVarInfo <<- list()
                                 	classEnvironment <<- new.env()
                                 	callSuper(...)
                                 },
                                 setupModel = function(code, constants, dimensions, debug) {},
                                 
                                 ## the following are all run, in this order, by setupModel():
                                 setModelValuesClassName        = function() {},
                                 assignBUGScode                 = function() {},
                                 assignConstants                = function() {},
                                 assignDimensions               = function() {},
                                 initializeContexts             = function() {},
                                 processBUGScode                = function() {},
                                 splitConstantsAndData          = function() {},
                                 addMissingIndexing             = function() {},
                                 removeTruncationWrapping       = function() {},
                                 expandDistributions            = function() {},
                                 processLinks                   = function() {},
                                 reparameterizeDists            = function() {},
                                 addRemainingDotParams          = function() {},
                                 insertDistributionBounds       = function() {},
                                 replaceAllConstants            = function() {},
                                 liftExpressionArgs             = function() {},
                                 addIndexVarsToDeclInfo         = function() {},
                                 genSymbolicParentNodes         = function() {},
                                 genReplacementsAndCodeReplaced = function() {},
                                 genAltParamsModifyCodeReplaced = function() {},
                                 genReplacedTargetValueAndParentInfo = function() {},
                                 genNodeInfo                    = function() {},
                                 removeEmptyBUGSdeclarations    = function() {},
                                 genVarInfo                     = function() {},
                                 genLogProbVarInfo              = function() {},
                                 genIsDataVarInfo               = function() {},
                                 genVarNames                    = function() {},
                                 buildSymbolTable               = function() {},
                                 buildIgraph                    = function() {},
                                 genGraphNodesList              = function() {},
                                 buildMaps2                     = function() {},
                                 
                                 newModel   = function() {},
                                 printDI    = function() {},
                                 
                                 genNodeInfo3                    = function() {},
                                 genVarInfo3                     = function() {},
                                 genExpandedNodeAndParentNames3  = function() {},
                                 
                                 #These functions are NOT run inside of setupModel
                                 nodeName2GraphIDs = function(){},
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
## set v3 = FALSE to use old processing
modelDefClass$methods(setupModel = function(code, constants, dimensions, debug = FALSE, v3 = TRUE, debugV3 = FALSE) {
    if(debug) browser()
    setModelValuesClassName()         ## uses 'name' field to set field: modelValuesClassName
    assignBUGScode(code)              ## uses 'code' argument, assigns field: BUGScode.  puts codes through nf_changeNimKeywords
    assignConstants(constants)        ## uses 'constants' argument, sets fields: constantsEnv, constantsList, constantsNamesList
    assignDimensions(dimensions)      ## uses 'dimensions' argument, sets field: dimensionList
    initializeContexts()              ## initializes the field: contexts
    processBUGScode()                 ## uses BUGScode, sets fields: contexts, declInfo$code, declInfo$contextID
    ## We will try to infer sizes later
    ##addMissingIndexing()              ## overwrites declInfo, using dimensionsList, fills in any missing indexing
    splitConstantsAndData()           ## deals with case when data is passed in as constants
    addMissingIndexing()              ## overwrites declInfo, using dimensionsList, fills in any missing indexing
    removeTruncationWrapping()        ## transforms T(ddist(),lower,upper) to put bounds into declInfo
    expandDistributions()             ## overwrites declInfo for stochastic nodes: calls match.call() on RHS      (uses distributions$matchCallEnv)
    processLinks()                    ## overwrites declInfo (*and adds*) for nodes with link functions           (uses linkInverses)
    reparameterizeDists()             ## overwrites declInfo when distribution reparameterization is needed       (uses distributions), keeps track of orig parameter in .paramName
    addRemainingDotParams()           ## overwrites declInfo, adds any additional .paramNames which aren't there  (uses distributions)
    insertDistributionBounds()        # put lower and upper into code expression
    replaceAllConstants()             ## overwrites declInfo with constants replaced; only replaces scalar constants
    liftExpressionArgs()              ## overwrites declInfo (*and adds*), lifts expressions in distribution arguments to new nodes.  does NOT lift '.param' names
    addIndexVarsToDeclInfo()          ## sets field declInfo[[i]]$indexVariableExprs from contexts.  must be after overwrites of declInfo
    genSymbolicParentNodes()          ## sets field declInfo[[i]]$symbolicParentNodes. must be after overwrites of declInfo
    genReplacementsAndCodeReplaced()  ## sets fields: declInfo[[i]]$replacements, $codeReplaced, $replacementNameExprs, $logProbNodeExpr
    genAltParamsModifyCodeReplaced()  ## sets field declInfo[[i]]$altParams, and modifies $codeReplaced to not include .param arguments (if stochastic)
    if(v3) {
        ## This is the new "version 3" processing
        if(debug) browser()
        genReplacedTargetValueAndParentInfo()
        genNodeInfo3(debug = debugV3)
        genVarInfo3()
        genExpandedNodeAndParentNames3(debug = debugV3)
        maps$setPositions3()
        buildSymbolTable()
        genIsDataVarInfo() ## only the maxs is ever used, in newModel
        genVarNames()
        return(NULL)
    }
    genNodeInfo()    #*               ## sets fields: declInfo[[i]]$indexedNodeInfo, and nodeInfo. these never change.
    removeEmptyBUGSdeclarations()     ## removes any declInfo[[i]] BUGSdecl objects for which length(indexedNodeInfo) == 0
    genVarInfo()                      ## uses nodeInfo and dimensionsList to set field: varInfo
    genLogProbVarInfo()               ## uses varInfo to set field: logProbVarInfo
    genIsDataVarInfo()                ## uses varInfo to set field: isDataVarInfo
    genVarNames()                     ## uses varInfo and logProbVarInfo to set field: varNames
    buildSymbolTable()                ## uses varInfo and logProbVarInfo to set field: symTab
    buildIgraph()      #*             ## uses declInfo and declInfo[[i]]$indexedNodeInfo to set field: graph
     ## old: use next two lines
    ##    genGraphNodesList()  #*           ## uses graphIDs, nodeInfo and varInfo to create graphNodesList, which is accessed by modeDef$graphNodes()
    ##    buildMaps()
    ## uses everything above to set field: maps
    ## new: use next one line
    buildMaps2()
})

modelDefClass$methods(setModelValuesClassName = function() {
    ## use this line to always ensure a unique internal refClass name
    modelValuesClassName <<- paste0(Rname2CppName(name), '_MV_', nimbleUniqueID())
    ## or this line if not.
    ##modelValuesClassName <<- paste0(name, '_MV')
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
    } else {
        constantsList <<- list()
        names(constantsList) <<- character(0)
        constantsNamesList <<- list()
    }
})
modelDefClass$methods(assignDimensions = function(dimensions) {
    ## uses 'dimensions' argument, sets field: dimensionList
    
    # first, add the provided dimensions
    dL <- dimensions
    
    # add dimensions of any *non-scalar* constants to dimensionsList
    # we'll try to be smart about this: check for duplicate names in constants and dimensions, and make sure they agree
    for(i in seq_along(constantsList)) {
        constName <- names(constantsList)[i]
        ## constDim <- if(is.null(dim(constantsList[[i]]))) length(constantsList[[i]]) else dim(constantsList[[i]])
        constDim <- dimOrLength(constantsList[[i]])
        if(constName %in% names(dL)) {
            if(!identical(as.numeric(dL[[constName]]), as.numeric(constDim))) {
                stop('inconsistent dimensions between constants and dimensions arguments')
            }
        } else {
            dL[[constName]] <- constDim
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

modelDefClass$methods(processBUGScode = function(code = NULL, contextID = 1, lineNumber = 0) {
    ## uses BUGScode, sets fields: contexts, declInfo$code, declInfo$contextID.
    ## all processing of code is done by BUGSdeclClass$setup(code, contextID).
    ## all processing of contexts is done by BUGScontextClass$setup()
    if(is.null(code)) {
        code <- BUGScode 
        declInfo <<- list()
    }
    for(i in 1:length(code)) {
        if(code[[i]] == '{') next  ## skip { lines
        lineNumber <- lineNumber + 1
        if(code[[i]][[1]] == '~' || code[[i]][[1]] == '<-') {  ## a BUGS declaration
            iAns <- length(declInfo) + 1
            BUGSdeclClassObject <- BUGSdeclClass$new() ## record the line number (a running count of non-`{` lines) for use in naming nodeFunctions later
            BUGSdeclClassObject$setup(code[[i]], contextID, lineNumber)
            declInfo[[iAns]] <<- BUGSdeclClassObject
        }
        if(code[[i]][[1]] == 'for') {        ## e.g. (for i in 1:N).  New context (for-loop info) needed
            indexVarExpr <- code[[i]][[2]]   ## This is the `i`
            indexRangeExpr <- code[[i]][[3]] ## This is the `1:N`
            if(nimbleOptions$prioritizeColonLikeBUGS) indexRangeExpr <- reprioritizeColonOperator(indexRangeExpr)
            nextContextID <- length(contexts) + 1
            forCode <- code[[i]][1:3]        ## This is the (for i in 1:N) without the code block
            
            sinlgeContexts <- c(if(contextID == 1) NULL else contexts[[contextID]]$singleContexts, ## concatenate any current contexts
                                list(BUGSsingleContextClass$new(indexVarExpr = indexVarExpr,       ## Add the new context
                                                                indexRangeExpr = indexRangeExpr,
                                                                forCode = forCode)))
            BUGScontextClassObject <- BUGScontextClass$new()
            BUGScontextClassObject$setup(singleContexts = sinlgeContexts)
            contexts[[nextContextID]] <<- BUGScontextClassObject
            if(length(code[[i]][[4]])==1) {
                stop(paste0('Error, not sure what to do with ', deparse(code[[i]])))
            }
            recurseCode <- if(code[[i]][[4]][[1]] == '{') {
                code[[i]][[4]]
            } else {
                substitute( {ONELINE}, list(ONELINE = code[[i]][[4]]))
            }
            lineNumber <- processBUGScode(recurseCode, nextContextID, lineNumber = lineNumber)  ## Recursive call to process the contents of the for loop
        }
        ## could consider definition-time if-then
    }
    lineNumber
})

modelDefClass$methods(splitConstantsAndData = function() {
    # removes items from constantsNamesList that appear as variables in declInfo
    # also, move detected data to 'data'
    # this deals with case when 'data' are passed in as 'constants'
    if(length(constantsNamesList)) {
        vars <- sapply(declInfo, function(x) x$targetVarName)
        constantsNames <- as.character(constantsNamesList)
        newDataVars <- constantsNames[constantsNames %in% vars]
        if(length(newDataVars)) {
            cat("Detected", paste(newDataVars, collapse = ','), "as data within 'constants'.\n")
            constantsNamesList <<- constantsNamesList[!constantsNames %in% vars]
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
    if(!any(code[[2]] == names(dimensionsList))) {
        ## dimension information was NOT provided for this variable
        ## let's check to make sure all indexes are present
        if(any(unlist(lapply(as.list(code), is.blank)))) {
            stop(paste0('Opps! This part of NIMBLE is still under development.', '\n',
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

modelDefClass$methods(removeTruncationWrapping = function() {
    ## pulls bounds out of T() syntax and puts in lower and upper fields
    for(i in seq_along(declInfo)) {
        
        BUGSdecl <- declInfo[[i]]
        if(BUGSdecl$type != 'stoch' || !(BUGSdecl$valueExpr[[1]] == "T" || BUGSdecl$valueExpr[[1]] == "I")) next

        BUGSdecl$truncation <- list(
            lower = if(BUGSdecl$valueExpr[[3]] == "") -Inf else BUGSdecl$valueExpr[[3]],
            upper = if(BUGSdecl$valueExpr[[4]] == "") Inf else BUGSdecl$valueExpr[[4]]
            # type = as.character(BUGSdecl$valueExpr[[1]]) # T or I? could be used later to check that I() only used for top-level nodes
            )
        if(BUGSdecl$valueExpr[[1]] == "I")
            warning(paste0("Interpreting I(,) as truncation (equivalent to T(,)) in ", deparse(BUGSdecl$code), "; this is only valid when ", deparse(BUGSdecl$targetExpr), " has no unobserved (stochastic) parents."))
                
        newCode <- BUGSdecl$code
        if(BUGSdecl$truncation$lower == -Inf && BUGSdecl$truncation$upper == Inf) {  # user specified no bounds so not really truncated
            BUGSdecl$truncation <- NULL
        }
        newCode[[3]] <- BUGSdecl$valueExpr[[2]]  # insert the core density function call

        BUGSdeclClassObject <- BUGSdeclClass$new()
        BUGSdeclClassObject$setup(newCode, BUGSdecl$contextID, BUGSdecl$sourceLineNumber, BUGSdecl$truncation)
        declInfo[[i]] <<- BUGSdeclClassObject
    }
})




modelDefClass$methods(expandDistributions = function() {
    ## overwrites declInfo for stochastic nodes: calls match.call() on RHS (uses distributions$matchCallEnv)
    for(i in seq_along(declInfo)) {
        
        BUGSdecl <- declInfo[[i]]
        if(BUGSdecl$type != 'stoch') next
        
        newCode <- BUGSdecl$code
        newCode[[3]] <- eval(BUGSdecl$valueExpr, distributions$matchCallEnv)
        
        BUGSdeclClassObject <- BUGSdeclClass$new()
        BUGSdeclClassObject$setup(newCode, BUGSdecl$contextID, BUGSdecl$sourceLineNumber, BUGSdecl$truncation)
        declInfo[[i]] <<- BUGSdeclClassObject
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
            BUGSdeclClassObject$setup(code, BUGSdecl$contextID, BUGSdecl$sourceLineNumber, BUGSdecl$truncation)
            newDeclInfo[[nextNewDeclInfoIndex]]     <- BUGSdeclClassObject
            
            BUGSdeclClassObject <- BUGSdeclClass$new()
            BUGSdeclClassObject$setup(newCode, BUGSdecl$contextID, BUGSdecl$sourceLineNumber, BUGSdecl$truncation)
            newDeclInfo[[nextNewDeclInfoIndex + 1]] <- BUGSdeclClassObject
            
        } else {    # deterministic node
            newRHS <- linkInverses[[linkText]]
            newRHS[[2]] <- BUGSdecl$code[[3]]
            newLHS <- BUGSdecl$targetNodeExpr
            newCode <- substitute(A <- B, list(A = newLHS, B = newRHS))
            
            BUGSdeclClassObject <- BUGSdeclClass$new()
            BUGSdeclClassObject$setup(newCode, BUGSdecl$contextID, BUGSdecl$sourceLineNumber, BUGSdecl$truncation)
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
        if(!(distName %in% distributions$namesVector))    stop('unknown distribution name: ', distName)      ## error if the distribution isn't something we recognize
        distRule <- distributions[[distName]]
        numArgs <- length(distRule$reqdArgs)
        newValueExpr <- quote(dist())       ## set up a parse tree for the new value expression
        newValueExpr[[1]] <- as.name(distName)     ## add in the distribution name
        newValueExpr[1 + (1:numArgs)] <- rep(NA, numArgs)      ## fill in the new parse tree with required arguments
        names(newValueExpr)[1 + (1:numArgs)] <- distRule$reqdArgs    ## add names for the arguments
        params <- as.list(valueExpr[-1])   ## extract the original distribution parameters
        
        if(identical(sort(names(params)), sort(distRule$reqdArgs))) {
            matchedAlt <- 0
        } else {
            matchedAlt <- NULL; count <- 0
            while(is.null(matchedAlt) && count < distRule$numAlts) {
                count <- count + 1
                if(identical(sort(unique(distRule$alts[[count]])), sort(unique(names(params)))))
                    matchedAlt <- count
            }
            if(is.null(matchedAlt)) stop('Error: no available re-parameterization found for: ', deparse(valueExpr), '.')
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
            for(nm in c(nonReqdArgs, distRule$reqdArgs))
                # loop thru possible non-canonical parameters in the expression for the canonical parameter
                transformedParameterPT <- parseTreeSubstitute(pt = transformedParameterPT, pattern = as.name(nm), replacement = params[[nm]])
            newValueExpr[[iArg + 1]] <- transformedParameterPT
        }
        
        ## hold onto the expressions for non-required args
        nonReqdArgExprs <- params[nonReqdArgs]    ## grab the non-required args from the original params list
        names(nonReqdArgExprs) <- if(length(nonReqdArgExprs) > 0) paste0('.', names(nonReqdArgExprs)) else character(0)  ## append '.' to the front of all the old (reparameterized away) param names
        newValueExpr <- as.call(c(as.list(newValueExpr), nonReqdArgExprs))
        
        newCode <- BUGSdecl$code
        newCode[[3]] <- newValueExpr
        
        BUGSdeclClassObject <- BUGSdeclClass$new()
        BUGSdeclClassObject$setup(newCode, BUGSdecl$contextID, BUGSdecl$sourceLineNumber, BUGSdecl$truncation)
        declInfo[[i]] <<- BUGSdeclClassObject
    }  # close loop over declInfo
})
modelDefClass$methods(addRemainingDotParams = function() {
    for(iDecl in seq_along(declInfo)) {
        BUGSdecl <- declInfo[[iDecl]]     ## grab this current BUGS declation info object
        if(BUGSdecl$type == 'determ')  next  ## skip deterministic nodes
        valueExpr <- BUGSdecl$valueExpr   ## grab the RHS (distribution)
        newValueExpr <- valueExpr
        defaultParamExprs <- distributions[[as.character(newValueExpr[[1]])]]$altParams
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
        BUGSdeclClassObject$setup(newCode, BUGSdecl$contextID, BUGSdecl$sourceLineNumber, BUGSdecl$truncation)
        declInfo[[iDecl]] <<- BUGSdeclClassObject
    }
})

modelDefClass$methods(insertDistributionBounds = function() {
    ## pulls bounds out of T() syntax and puts in lower and upper fields
    for(i in seq_along(declInfo)) {
        
        BUGSdecl <- declInfo[[i]]
        if(BUGSdecl$type != 'stoch' || is.null(BUGSdecl$truncation)) next

        newValueExpr <- BUGSdecl$valueExpr   ## grab the RHS (distribution)

        nParams <- length(newValueExpr) 
        newValueExpr[[nParams + 1]] <- BUGSdecl$truncation$lower
        newValueExpr[[nParams + 2]] <- BUGSdecl$truncation$upper
        names(newValueExpr)[(nParams + 1):(nParams + 2)] <- c('lower', 'upper')

        newCode <- BUGSdecl$code
        newCode[[3]] <- newValueExpr

        BUGSdeclClassObject <- BUGSdeclClass$new()
        BUGSdeclClassObject$setup(newCode, BUGSdecl$contextID, BUGSdecl$sourceLineNumber, BUGSdecl$truncation)

        declInfo[[i]] <<- BUGSdeclClassObject
    }
})




modelDefClass$methods(replaceAllConstants = function() {
    ## overwrites declInfo with constants replaced; only replaces scalar constants
    ## does both LHS and RHS of each BUGSdecl code
    for(i in seq_along(declInfo)) {
        newCode <- replaceConstantsRecurse(declInfo[[i]]$code, constantsEnv, constantsNamesList)$code
        
        BUGSdeclClassObject <- BUGSdeclClass$new()
        BUGSdeclClassObject$setup(newCode, declInfo[[i]]$contextID, declInfo[[i]]$sourceLineNumber, declInfo[[i]]$truncation)
        declInfo[[i]] <<- BUGSdeclClassObject
    }
})

neverReplaceable <- list(chol = TRUE, inverse = TRUE) ## only the names matter, any non-null value will do.

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
            if(!any(code[[1]] == distributions$namesVector)) {
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
            
            for(iParam in seq_along(params)) {
                if(grepl('^\\.', names(params)[iParam]))   next        ## skips '.param' names; we do NOT lift these
                paramExpr <- params[[iParam]]
                if(!isExprLiftable(paramExpr))    next     ## if this param isn't an expression, go ahead to next parameter
                newNodeNameExpr <- as.name(paste0('lifted_', nameMashupFromExpr(paramExpr, colonsOK = TRUE)))   ## create the name of the new node
                newNodeNameExprIndexed <- addNecessaryIndexingToNewNode(newNodeNameExpr, paramExpr, contexts[[BUGSdecl$contextID]]$indexVarExprs)  ## add indexing if necessary
                
                newValueExpr[[iParam + 1]] <- newNodeNameExprIndexed  ## update the newValueExpr
                
                newNodeCode <- substitute(LHS <- RHS, list(LHS = newNodeNameExprIndexed, RHS = paramExpr))     ## create code line for declaration of new node
                if(!checkForDuplicateNodeDeclaration(newNodeCode, newNodeNameExprIndexed, newDeclInfo)) {
                    
                    BUGSdeclClassObject <- BUGSdeclClass$new()
                    BUGSdeclClassObject$setup(newNodeCode, BUGSdecl$contextID, BUGSdecl$sourceLineNumber, BUGSdecl$truncation)   ## keep new declaration in the same context, regardless of presence/absence of indexing
                    newDeclInfo[[nextNewDeclInfoIndex]] <- BUGSdeclClassObject
                    
                    nextNewDeclInfoIndex <- nextNewDeclInfoIndex + 1     ## update for lifting other nodes, and re-adding BUGSdecl at the end
                }
            }    # closes loop over params
        }
        
        newCode <- BUGSdecl$code
        newCode[[3]] <- newValueExpr
        
        BUGSdeclClassObject <- BUGSdeclClass$new()
        BUGSdeclClassObject$setup(newCode, BUGSdecl$contextID, BUGSdecl$sourceLineNumber, BUGSdecl$truncation)
        newDeclInfo[[nextNewDeclInfoIndex]] <- BUGSdeclClassObject    ## regardless of anything, add BUGSdecl itself in
    }    # closes loop over declInfo
    declInfo <<- newDeclInfo
})
isExprLiftable <- function(paramExpr) {
    ## determines whether a parameter expression is worthy of lifiting up to a new node
    if(is.name(paramExpr))       return(FALSE)
    if(is.numeric(paramExpr))    return(FALSE)
    if(is.call(paramExpr)) {
        if(paramExpr[[1]] == 'chol')        return(TRUE)    ## do lift calls to chol(...)
        if(paramExpr[[1]] == 'inverse')     return(TRUE)    ## do lift calls to inverse(...)
        if(length(paramExpr) == 1)          return(FALSE)   ## don't generally lift function calls:   fun(...) ## this comment seems incorrect
        if(getCallText(paramExpr) == '[')   return(FALSE)   ## don't lift simply indexed expressions:  x[...]
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
    usedIndexVarsList <- indexVarExprs[indexVarExprs %in% all.vars(paramExpr)]    # this extracts any index variables which appear in 'paramExpr'
    vectorizedIndexExprsList <- extractAnyVectorizedIndexExprs(paramExpr)    # creates a list of any vectorized (:) indexing expressions appearing in 'paramExpr'
    neededIndexExprsList <- c(usedIndexVarsList, vectorizedIndexExprsList)
    if(length(neededIndexExprsList) == 0)  return(newNodeNameExpr)  # no index variables, or vectorized indexing, return the (un-indexed) name expression
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
            ## we've found a node declaration with exactly the same LHS
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
    
    nimFunNames <- distributions$namesExprList
    
    for(i in seq_along(declInfo)){
        declInfo[[i]]$genSymbolicParentNodes(constantsNamesList, contexts[[declInfo[[i]]$contextID]], nimFunNames)
    }
})
modelDefClass$methods(genReplacementsAndCodeReplaced = function() {
    ## sets fields declInfo[[i]]$replacements, $codeReplaced, and $replacementNameExprs
    
    nimFunNames <- distributions$namesExprList
    
    for(i in seq_along(declInfo)) {
        declInfo[[i]]$genReplacementsAndCodeReplaced(constantsNamesList, contexts[[declInfo[[i]]$contextID]], nimFunNames)
    }
})

modelDefClass$methods(genReplacedTargetValueAndParentInfo = function() {
    
    nimFunNames <- distributions$namesExprList
    
    for(i in seq_along(declInfo)) {
        declInfo[[i]]$genReplacedTargetValueAndParentInfo(constantsNamesList, contexts[[declInfo[[i]]$contextID]], nimFunNames)
    }
    NULL
})

modelDefClass$methods(genAltParamsModifyCodeReplaced = function() {
    ## setsfield declInfo[[i]]$altParams, and modifies $codeReplaced to not include .param arguments (if stochastic)
    for(i in seq_along(declInfo)) {
        declInfo[[i]]$genAltParamsModifyCodeReplaced()
    }
})

## genNodeInfo3
## There are no longer nodeInfo objects as we used to have.  This does the related processing in the new system.
modelDefClass$methods(genNodeInfo3 = function(debug = FALSE) {
    if(debug) browser()
    ## Iterate by context (for loops), where 1st context will always be no-for-loop
    for(i in seq_along(contexts)) {
        boolContext <- unlist(lapply(declInfo, function(x) x$contextID == i))                    ## TRUE for BUGS lines (declInfo elements) that use this context
        allReplacements <- do.call('c', lapply(declInfo[boolContext], `[[`, 'replacements'))     ## Get a list of replacement expressions
        if(length(allReplacements) > 0) {
            allReplacementNameExprs <- do.call('c', lapply(declInfo[boolContext], `[[`, 'replacementNameExprs')) ## names of allReplacements as expressions
            boolNotDup <- !duplicated(allReplacements)

            ## This makes an environment with a vector of each replacement and for-loop index from executing the for-loops
            unrolledContextAndReplacementsEnv <- expandContextAndReplacements(allReplacements[boolNotDup], allReplacementNameExprs[boolNotDup], contexts[[i]], constantsEnv)
            contexts[[i]]$replacementsEnv <<- unrolledContextAndReplacementsEnv
        } else {
            contexts[[i]]$replacementsEnv <<- NULL
        }
    }

    for(i in seq_along(declInfo)) {
        ## extract information needed for each declInfo and each parentExpr
        ## note that sometimes lifted nodes don't need some or any of the indices from the context
        ## e.g. there might be a lifted node involving no indices, but it is still embedded in the same context
        BUGSdecl <- declInfo[[i]]
        context <- contexts[[BUGSdecl$contextID]]
        BUGSdecl$replacementsEnv <- context$replacementsEnv
        if(is.null(BUGSdecl$replacementsEnv)) {
            BUGSdecl$unrolledIndicesMatrix <- matrix(nrow = 0, ncol = 0)
            next
        }
        BUGSdecl$outputSize <- BUGSdecl$replacementsEnv$outputSize

        useContext <- unlist(lapply(context$singleContexts, function(x) isNameInExprList(x$indexVarExpr, BUGSdecl$indexExpr)))
        BUGSdecl$replacementsEnv$cbindArgExprs <- unique(c(context$indexVarExprs[useContext], BUGSdecl$replacementNameExprs))
        
        BUGSdecl$unrolledIndicesMatrix <- with(BUGSdecl$replacementsEnv, do.call('cbind', cbindArgExprs))
        rm(list = 'cbindArgExprs', envir = BUGSdecl$replacementsEnv)
        if(!all(useContext)) {
            if(!any(useContext)) {
                BUGSdecl$replacementsEnv <- NULL
                BUGSdecl$unrolledIndicesMatrix <- matrix(nrow = 0, ncol = 0)
                next
            }
            boolUse <- !duplicated(BUGSdecl$unrolledIndicesMatrix[, context$indexVarNames[useContext] ] )
            BUGSdecl$unrolledIndicesMatrix <- BUGSdecl$unrolledIndicesMatrix[boolUse, , drop = FALSE]

            BUGSdecl$replacementsEnv <- list2env(as.data.frame(BUGSdecl$unrolledIndicesMatrix))
            BUGSdecl$replacementsEnv$outputSize <- nrow(BUGSdecl$unrolledIndicesMatrix)
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


expandContextAndReplacements <- function(allReplacements, allReplacementNameExprs, context, constantsEnv) {
 ##   browser()
    numContexts <- length(context$singleContexts)
    if(numContexts == 0) { ## it has no indices or known indices
        replacementsEnv <<- NULL
        return(NULL)
    }

    namesToRemoveAtEnd <- ls(constantsEnv)
    constantsEnvCopy <- list2env(as.list(constantsEnv))
    ## some replacements like min(j:100) should no longer be needed but are still there

    ## If this all works, useContext can be removed
    useContext <- rep(TRUE, numContexts)
    
    valueVarNames <- paste0("INDEXVALUE_", 1:numContexts, "_")
## indexRecordingCode gives lines of code like "INDEXVALUE_1_[iAns] <- i". This will later have its name changed to "i"
    indexRecordingCode <- vector('list', length = numContexts)
    for(i in 1:numContexts) {
        if(useContext[i])
            indexRecordingCode[[i]] <- substitute(V[iAns] <- index, list(V = as.name(valueVarNames[i]), index = context$singleContexts[[i]]$indexVarExpr))
    }

    ##    numReplacements <- length(BUGSdecl$replacements)
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
            replacementRecordingCode[[i]] <- substitute(A[iAns] <- B, list(A = allReplacementNameExprs[[i]], B = allReplacements[[i]])) 
    }
    ## From here through the while loop combines the for loops from the contexts, with the replacementRecordingCode and indexRecordingCode in the innermost
    innerLoopCode <- as.call(c(list(quote(`{`)), replacementRecordingCode, indexRecordingCode, quote(iAns <- iAns + 1)))

    iContext <- numContexts
    while(iContext >= 1) {
        if(useContext[iContext]) {
            newCode <- context$singleContexts[[iContext]]$forCode
            newCode[[4]] <- innerLoopCode
            innerLoopCode <- newCode
        }
        iContext <- iContext - 1
    }
    ## at this point "innerLoopCode" has the full loop  ## determineContextSize does something similar -- creates and executes nested for loops -- only for the purpose of counting how big the result will be
    outputSize <- determineContextSize(context, useContext, constantsEnvCopy)
    for(i in 1:numContexts) {
        if(useContext[i])
            assign(valueVarNames[i], numeric(outputSize), constantsEnvCopy)
    }
    for(i in seq_along(replacementRecordingCode)) {
        if(useReplacement[i])
            assign(names(allReplacements)[i], numeric(outputSize), constantsEnvCopy)
    }
    assign("iAns", 1, constantsEnvCopy)
    eval(innerLoopCode, constantsEnvCopy)
    for(i in 1:numContexts){
        if(useContext[i]) {
            constantsEnvCopy[[ as.character(context$singleContexts[[i]]$indexVarExpr) ]] <- constantsEnvCopy[[ valueVarNames[i] ]]
            rm(list = valueVarNames[i], envir = constantsEnvCopy)
        }
    }
    rm(list = c(namesToRemoveAtEnd, 'iAns'), envir = constantsEnvCopy)
    assign("outputSize", outputSize, constantsEnvCopy)
    return(constantsEnvCopy) ## becomes replacementsEnv
}

determineContextSize <- function(context, useContext = rep(TRUE, length(context$singleContexts)), evalEnv = new.env()) {
    ## could improve this by checking for nested loops that don't use indices from outer loops
    innerLoopCode <- quote(iAns <- iAns + 1)
    iContext <- length(context$singleContexts)
    while(iContext >= 1) {
        if(useContext[iContext]) {
            newCode <- context$singleContexts[[iContext]]$forCode
            newCode[[4]] <- innerLoopCode
            innerLoopCode <- newCode
        }
        iContext <- iContext - 1
    }
    assign("iAns", 0L, evalEnv)
    eval(innerLoopCode, evalEnv)
    ans <- evalEnv$iAns
    rm(list = c('iAns', context$indexVarNames[useContext]), envir = evalEnv)
    return(ans)
}

modelDefClass$methods(genNodeInfo = function() {
    ## sets fields: declInfo[[i]]$indexedNodeInfo, and nodeInfo. these never change.
    
    for(i in seq_along(declInfo)) {
        genNodeInfo_singleDeclaration(declInfo[[i]], contexts[[declInfo[[i]]$contextID]], constantsEnv)
    }
    
    ## gather all the indexedNodeNames from each BUGSdecl object; throw an error if any duplicates exist
    allIndexedNodeNames <- unlist(lapply(declInfo, `[[`, 'indexedNodeNames'))
    duplicateNodes <- unique(allIndexedNodeNames[duplicated(allIndexedNodeNames)])
    if(length(duplicateNodes) > 0)   stop(paste0('duplicate node declarations: ', paste0(duplicateNodes, collapse=', ')))
    
    nodeInfoTemp <- unlist(lapply(declInfo, `[[`, 'indexedNodeInfo'))
    nodeInfo <<- if(is.null(nodeInfoTemp)) list() else nodeInfoTemp
})

#### wrapAsNumeric <- function(code) substitute(as.numeric(X), list(X = code))     ## as.numeric() flattens everything to a vector
wrapAsNumeric <- function(code) substitute({ value <- CODE;   storage.mode(value) <- 'numeric';   value }, list(CODE = code))

genNodeInfo_singleDeclaration <- function(BUGSdecl, context, constantsEnv) {
    ## This function turns out to be a computational bottleneck, so I have added some efficiencies.  There is room for more.
    constantsEnvCopy <- list2env(as.list(constantsEnv))
    indexVariableExprs <- context$indexVarExprs
    indexVarValues <- context$genIndexVarValues(constantsEnvCopy)
    BUGSdecl$indexedNodeInfo <- vector(mode = 'list', length = length(indexVarValues))
    targetNodeExprHasBracket <- hasBracket(BUGSdecl$targetNodeExpr)
    targetNodeExpressionIsVectorized <- unlist(lapply(as.list(BUGSdecl$targetNodeExpr), is.vectorized))
    substituteCodeReadyForValues <- substitute(substitute(CODE, replacementValuesScalar), list(CODE=BUGSdecl$codeReplaced))
    substituteAltParamsReadyForValues <- lapply(BUGSdecl$altParamExprs, function(expr) substitute(substitute(EXPR, replacementValuesScalar), list(EXPR=expr)))
    logProbSubstituteCodeReadyForValues <- substitute(substitute(LOGPROBNODE, replacementValuesScalar), list(LOGPROBNODE=BUGSdecl$logProbNodeExpr))
    doLogProbSubstitute <- !is.null(BUGSdecl$logProbNodeExpr)
    logProbNodeReplacedWithValues <- NULL
    logProbIndexValues <- NULL
    replacementsOneBlock <- as.call(c(list(as.name('list')), lapply(BUGSdecl$replacements, wrapAsNumeric)))

    for(i in seq_along(indexVarValues)) {
        indexValues <- indexVarValues[[i]]
        if(length(indexValues) > 0)     list2env(indexValues, constantsEnvCopy)
        ## evalBracketArgs could operate on args of ':' directly
        evaledTargetNodeExpr <- if(targetNodeExprHasBracket) evalBracketArgsKnownBracket(BUGSdecl$targetNodeExpr, constantsEnvCopy, targetNodeExpressionIsVectorized) else BUGSdecl$targetNodeExpr
        evaledSymbolicParentNodes <- lapply(BUGSdecl$symbolicParentNodes, evalBracketArgs, constantsEnvCopy)
        targetNodeIndexValuesList <- if(targetNodeExprHasBracket) exprAsListDrop2(evaledTargetNodeExpr)
        targetNodeIndexSizes <- unlist(lapply(targetNodeIndexValuesList, function(ind) max(eval(ind)) - min(eval(ind)) + 1))
        replacementValues <- eval(replacementsOneBlock, constantsEnvCopy)
        replacementValuesScalar <- removeNonScalarElementsFromList(replacementValues)
        codeReplacedWithValues <- eval(substituteCodeReadyForValues)
        altParamExprsWithValues <- lapply(substituteAltParamsReadyForValues, function(expr) eval(expr))
        if(doLogProbSubstitute) {
            logProbNodeReplacedWithValues <- eval(logProbSubstituteCodeReadyForValues)
            logProbIndexValues <- if(length(logProbNodeReplacedWithValues)==0)   numeric(0)   else   unlist(as.list(logProbNodeReplacedWithValues)[-c(1,2)])
        }
        targetNodeName <- if(targetNodeExprHasBracket) deparse(evaledTargetNodeExpr) else BUGSdecl$targetNodeName
        
        BUGSdecl$indexedNodeInfo[[i]] <- nodeInfoClass$new(code = BUGSdecl$code,
                                                           type = BUGSdecl$type,
                                                           targetNodeExpr = evaledTargetNodeExpr,
                                                           targetNodeName =  targetNodeName,
                                                           targetNodeIndexValuesList = targetNodeIndexValuesList,
                                                           targetNodeIndexSizes = targetNodeIndexSizes,
                                                           logProbIndexValues = logProbIndexValues,
                                                           targetVarName = BUGSdecl$targetVarName,
                                                           parentNodeExprs = evaledSymbolicParentNodes,
                                                           parentNodeNames = lapply(evaledSymbolicParentNodes, deparse), ## could avoid some deparses here
                                                           nodeFunctionName = targetNodeName,
                                                           indexVariableExprs = indexVariableExprs,
                                                           indexVariableValues = indexValues,
                                                           codeReplaced = BUGSdecl$codeReplaced,
                                                           replacementValues = replacementValues,
                                                           codeReplacedWithValues = codeReplacedWithValues,
                                                           logProbNodeReplacedWithValues = logProbNodeReplacedWithValues,
                                                           altParamExprs = BUGSdecl$altParamExprs,
                                                           altParamExprsWithValues = altParamExprsWithValues
        )
    }   # closes loop over indexVarValues
    indexedNodeNames <- unlist(lapply(BUGSdecl$indexedNodeInfo, `[[`, 'nodeFunctionName'))
    names(BUGSdecl$indexedNodeInfo) <- if(is.null(indexedNodeNames)) character(0) else indexedNodeNames
    
    ## experiment:  this will remove duplicate node names
    ## idea being, this would allow for code like
    ## for(i in 1:10) {
    ##    mu ~ dnorm(0,1)
    ## }
    ## to only produce a single 'mu' variable, for this BUGSdeclClass object
    BUGSdecl$indexedNodeInfo <- BUGSdecl$indexedNodeInfo[unique(names(BUGSdecl$indexedNodeInfo))]
    
    BUGSdecl$indexedNodeNames <- names(BUGSdecl$indexedNodeInfo)
}
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

makeVertexNamesFromIndexArray2 <- function(indArr, minInd = 1, varName) {
    dims <- dim(indArr)
    nDim <- length(dims)
    indArr[indArr < minInd] <- NA
    if(all(is.na(indArr))) return(list(indices = integer(), names = character()))
    arrayWithIndices <- vector('list', length = nDim)
    arrayWithIndices[[1]] <- array( rep(1:dims[1], prod(dims[-1])), dims)
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
    splits <- lapply(arrayWithIndices, split, indArr)

    info <- lapply(splits, lapply, function(x) {r <- range(x); s <- r[1]!=r[2]; c <- if(!s) diff(r)+1 == length(unique(x)) else TRUE; c(r, as.numeric(s), as.numeric(c))})
    ## each entry is (min, max, 0/1 for vector, 0/1 for contiguous)
    dimStrings <- lapply(info, function(x) {
        all <- do.call('rbind', x)
        seps <- rep(':', nrow(all))
        scal <- all[,3]==0
        seps[scal] <- ''
        seps[all[,4]==0] <- '-s-'
        maxStrs <- as.character(all[,2])
        maxStrs[scal] <- ''
        paste0(all[,1], seps, maxStrs)
    })
    dimStrings[['sep']] <- ', '
    newNames <- paste0(varName, '[',  do.call('paste', dimStrings), ']')
    list(indices = as.integer(names(splits[[1]])), names = newNames)
##    contigs <- lapply(splits, lapply, function(x) 
}

makeVertexNamesFromIndexArray <- function(indArr, minInd = 1, varName) {
    ## Obvious potential for speedup via C
    uniqueInds <- unique(as.numeric(indArr))
    if(minInd > 1) uniqueInds <- uniqueInds[uniqueInds >= minInd]
    uniqueInds <- uniqueInds[!is.na(uniqueInds)]
    newNames <- character(length(uniqueInds))
    strfun <- function(x) if(x[1] == x[2]) as.character(x[1]) else paste0(x[1], ':', x[2])
    for(iI in seq_along(uniqueInds)) {
        i <- uniqueInds[iI]
        locationsInArray <- which(indArr == i, TRUE)
        ranges <- apply(locationsInArray, 2, range)
        strings <- apply(ranges, 2, strfun)
        contiguous <- nrow(locationsInArray) == prod(diff(ranges) + 1)
        newNames[iI] <- paste0( varName, if(contiguous) character(0) else "split", '[', paste0(strings, collapse = ", "), ']')
    }
    list(indices = uniqueInds, names = newNames)
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
            var2vertexID[ var2vertexID == ID ] <- nextVertexID
            nextVertexID <- nextVertexID + 1
        }
    }
    list(var2vertexID = var2vertexID, nextVertexID = nextVertexID) 
}

splitVertices <- function(var2vertexID, unrolledBUGSindices, indexExprs = NULL, indexNames = NULL, parentExpr, parentExprReplaced = NULL, parentIndexNamePieces, replacementNameExprs, nextVertexID, maxVertexID, debug = FALSE) {
        if(debug) browser()
    ## 1. Determine which indexExprs are in parentExpr
    useContext <- unlist(lapply(indexExprs, isNameInExprList, parentExpr))
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
    vectorIndices <- lapply(parentIndexNamePieces, function(x) {if(is.list(x)) {allScalar <<- FALSE; return(TRUE)}; FALSE})

    ## step 4 evaporated    
    currentVertexCounts <- tabulate(var2vertexID, maxVertexID)
    ## 5. Set up initial table of vertexIDcounts

    ## 6. All scalar case: iterate or vectorize via cbind and put new vertexIDs over -1s
    if(allScalar) {
        if(anyContext) {
            boolIndexNamePiecesExprs <- !unlist(lapply(parentIndexNamePieces, is.numeric)) 
            if(all(boolIndexNamePiecesExprs)) {
                test <- try(varIndicesToUse <- unrolledBUGSindices[ boolUseUnrolledRow, unlist(parentIndexNamePieces) ])
                if(inherits(test, 'try-error')) browser()
            } else {
                varIndicesToUse <- matrix(nrow = sum(boolUseUnrolledRow), ncol = length(parentIndexNamePieces))
                varIndicesToUse[, boolIndexNamePiecesExprs] <- unrolledBUGSindices[ boolUseUnrolledRow, unlist(parentIndexNamePieces)[boolIndexNamePiecesExprs]]
                indexPieceNumericInds <- which(!boolIndexNamePiecesExprs)
                for(iii in seq_along(indexPieceNumericInds)) varIndicesToUse[, indexPieceNumericInds[iii] ] <- parentIndexNamePieces[[ indexPieceNumericInds[iii] ]]
            }
        }        
        else {
            if(length(parentIndexNamePieces)==1) varIndicesToUse <- 1
            else {
                varIndicesToUse <- matrix(0, nrow = 1, ncol = length(parentIndexNamePieces))
                for(iI in 1:ncol(varIndicesToUse)) varIndicesToUse[1, iI] <- as.numeric(parentExprReplaced[[iI+2]])
            }
        }
        currentVertexIDs <- var2vertexID[ varIndicesToUse ]
        needsVertexID <- is.na(currentVertexIDs)
        numNewVertexIDs <- sum(needsVertexID)
        if(numNewVertexIDs > 0) {
            var2vertexID[ varIndicesToUse ][needsVertexID] <- nextVertexID - 1 + 1:numNewVertexIDs
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
            var2vertexID[ varIndicesToUse][ oldIndices][needsSplit] <- nextVertexID - 1 + 1:numNeedSplit
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
            uniqueCurrentVertexIDs <- unique(currentVertexIDblock)
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
    if(debug) browser()
    anyContext <- ncol(unrolledBUGSindices) > 0 
    if(length(anyContext)==0) browser()
    
    allScalar <- TRUE
    vectorIndices <- lapply(parentIndexNamePieces, function(x) {if(is.list(x)) {allScalar <<- FALSE; return(TRUE)}; FALSE})
        
    if(allScalar) {
        if(is.null(parentIndexNamePieces)) varIndicesToUse <- rep(1, length(targetIDs))
        else {
            if(anyContext) {
                boolIndexNamePiecesExprs <- !unlist(lapply(parentIndexNamePieces, is.numeric)) ##!is.numeric(parentIndexNamePieces)
                if(all(boolIndexNamePiecesExprs)) {
                    test <- try(varIndicesToUse <- unrolledBUGSindices[ , unlist(parentIndexNamePieces) ])
                    if(inherits(test, 'try-error')) browser()
                } else {
                    varIndicesToUse <- matrix(nrow = nrow(unrolledBUGSindices), ncol = length(parentIndexNamePieces))
                    varIndicesToUse[, boolIndexNamePiecesExprs] <- unrolledBUGSindices[ , unlist(parentIndexNamePieces)[boolIndexNamePiecesExprs]]
                    indexPieceNumericInds <- which(!boolIndexNamePiecesExprs)
                    for(iii in seq_along(indexPieceNumericInds)) varIndicesToUse[, indexPieceNumericInds[iii] ] <- parentIndexNamePieces[[ indexPieceNumericInds[iii] ]]
                }
            } else {
                if(length(parentIndexNamePieces)==1) varIndicesToUse <- rep(1, length(targetIDs))
                else {
                    varIndicesToUse <- matrix(0, nrow = 1, ncol = length(parentIndexNamePieces))
                    for(iI in 1:ncol(varIndicesToUse)) varIndicesToUse[1, iI] <- as.numeric(parentExprReplaced[[iI+2]])
                }
            }
        }

        edgesFrom <- var2vertexID[ varIndicesToUse ]
        edgesTo <- targetIDs
        
    } else {
        if(anyContext) {
           if(length(ncol(unrolledBUGSindices))==0) browser()
            colNums <- 1:ncol(unrolledBUGSindices)
            names(colNums) <- dimnames(unrolledBUGSindices)[[2]]
            newIndexExprs <- lapply(replacementNameExprs, function(x) substitute(unrolledBUGSindices[iRow, iCol], list(X = x, iCol = colNums[as.character(x)])))
            accessExpr <- eval( substitute( substitute(AE, newIndexExprs), list(AE = parentExprReplaced) ) )
            iRowRange <- (1:nrow(unrolledBUGSindices))
            maxNewEdges <- nrow(unrolledBUGSindices) * prod(dim(var2vertexID))
        } else {
            accessExpr <- parentExprReplaced
            iRowRange <- 1
            maxNewEdges <- prod(dim(var2vertexID))
        }
       accessExpr[[2]] <- quote(var2vertexID)

       iNextNewEdge <- 1
       edgesFrom <- integer(maxNewEdges)
       edgesTo <- integer(maxNewEdges)
       
        for(iRow in iRowRange) {
            currentVertexIDblock <- eval(accessExpr)
            uniqueCurrentVertexIDs <- unique(currentVertexIDblock)
            numNewEdges <- length(uniqueCurrentVertexIDs)
            edgesFrom[iNextNewEdge - 1 + 1:numNewEdges] <- uniqueCurrentVertexIDs
            edgesTo[iNextNewEdge - 1 + 1:numNewEdges] <- targetIDs[iRow]
            iNextNewEdge <- iNextNewEdge + numNewEdges
        }
       edgesFrom <- edgesFrom[1:(iNextNewEdge-1)]
       edgesTo <- edgesTo[1:(iNextNewEdge-1)]
   }
    list(edgesFrom = edgesFrom, edgesTo = edgesTo)    
}



modelDefClass$methods(genExpandedNodeAndParentNames3 = function(debug = FALSE) {
    if(debug) browser()
    ## 1. initialize origMaps:
    vars_2_nodeOrigID <- new.env()     ## IDs for node function labels, e.g. for x[1:3], x[4]: NA, NA, NA, 2
    vars_2_vertexOrigID <- new.env()  ## IDs for node labels, e.g. x[1:3]: 1, 1, 1, 2
    vars2LogProbName <- new.env()
    vars2LogProbID <- new.env()
    
    for(iV in seq_along(varInfo)) {
        varName <- varInfo[[iV]]$varName
##        logProbVarName <- makeLogProbName(varName)
        
        if(varInfo[[iV]]$nDim > 0) {
            vars_2_nodeOrigID[[varName]] <- array(as.numeric(NA), dim = varInfo[[iV]]$maxs)
            vars2LogProbName[[varName]] <- array(dim = varInfo[[iV]]$maxs)
            vars2LogProbID[[varName]] <- array(dim = varInfo[[iV]]$maxs)
            storage.mode(vars2LogProbName[[varName]]) <- 'character'
        } else {
            vars_2_nodeOrigID[[varName]] <- as.numeric(NA)
            vars2LogProbName[[varName]] <- as.character(NA)
            vars2LogProbID[[varName]] <- as.numeric(NA)
        }
    }

    ## 2. collect names and do eval to create variables in vars_2_nodeOrigID
    next_origID <- 1
    allNodeNames <- character()
    types <- character()
    ## LHS declarations create nodes:
    for(iDI in seq_along(declInfo)) {

        BUGSdecl <- declInfo[[iDI]]
        lhsVar <- BUGSdecl$targetVarName
        nDim <- varInfo[[lhsVar]]$nDim
        if(nDim > 0) {
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
            BUGSdecl$nodeFunctionNames <- paste0(lhsVar, '[', do.call('paste', pieces), ']') 
            allNodeNames <- c(allNodeNames, BUGSdecl$nodeFunctionNames)
            types <- c(types, rep(BUGSdecl$type, length(BUGSdecl$nodeFunctionNames) ) )
            BUGSdecl$origIDs <- next_origID -1 + (1:length(BUGSdecl$nodeFunctionNames))
            next_origID <- next_origID + length(BUGSdecl$nodeFunctionNames)
            
            if(is.environment(BUGSdecl$replacementsEnv)) {
                BUGSdecl$replacementsEnv[['origIDs']] <- BUGSdecl$origIDs
                IDassignCode <- insertSubIndexExpr(BUGSdecl$targetExprReplaced, quote(iAns))
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
            BUGSdecl$nodeFunctionNames <- lhsVar
            allNodeNames <- c(allNodeNames, BUGSdecl$nodeFunctionNames)
            types <- c(types, BUGSdecl$type)
            BUGSdecl$origIDs <- next_origID
            vars_2_nodeOrigID[[lhsVar]] <- next_origID
            next_origID <- next_origID + 1
        }
    }
    nodeNamesLHSall <- allNodeNames
    maxOrigNodeID <- next_origID - 1
    numNodeFunctions <<- maxOrigNodeID
    allNewVertexNames <- character(0)
    allNewVertexIDs <- integer(0)

    ## 2b. Collect logProbNames and do eval to create variables in vars2LogProbName and vars2LogProbID
    if(debug) browser()
    nextLogProbID <- 1
    logProbNames <- character()
    logProbIDs <- integer()
    for(iDI in seq_along(declInfo)) {
        BUGSdecl <- declInfo[[iDI]]
        if(BUGSdecl$type == 'determ') next
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
                BUGSdecl$replacementsEnv[[lhsVar]] <- vars2LogProbID[[lhsVar]]
                eval(forCode, envir = BUGSdecl$replacementsEnv)
                vars2LogProbID[[lhsVar]] <- BUGSdecl$replacementsEnv[[lhsVar]]
                rm(list = c(lhsVar, 'logProbIDs'), envir = BUGSdecl$replacementsEnv)

                BUGSdecl$replacementsEnv[['logProbIDs']] <- newLogProbNames
                BUGSdecl$replacementsEnv[[lhsVar]] <- vars2LogProbName[[lhsVar]]
                eval(forCode, envir = BUGSdecl$replacementsEnv)
                vars2LogProbName[[lhsVar]] <- BUGSdecl$replacementsEnv[[lhsVar]]
                rm(list = c(lhsVar, 'logProbIDs'), envir = BUGSdecl$replacementsEnv)
            } else {
                ## If no replacementsEnv was set up, then there were no index variables (only numerics)
                eval(substitute(A <- B, list(A = targetExprWithMins, B = newLogProbIDs)), envir = vars2LogProbID)
                eval(substitute(A <- B, list(A = targetExprWithMins, B = newLogProbNames)), envir = vars2LogProbName)
            }
        } else { ## nDim == 0
            logProbNames <- c(logProbNames, logProbVarName)
            logProbIDs <- c(logProbIDs, nextLogProbID)
            vars2LogProbID[[lhsVar]] <- nextLogProbID
            vars2LogProbName[[lhsVar]] <- logProbVarName
            nextLogProbID <- nextLogProbID + 1
        }
    }
    
    ## 3. determine total model size from all objects
    ##    copy vars_2_nodeOrigID into vars_2_vertexOrigID
    totModelSize <- 0
    for(iV in seq_along(varInfo)) {
        totModelSize <- totModelSize + if(varInfo[[iV]]$nDim == 0) 1 else prod(varInfo[[iV]]$maxs)
        varName <- varInfo[[iV]]$varName
        vars_2_vertexOrigID[[ varName ]] <- vars_2_nodeOrigID[[ varName ]]
    }

    ## 4. Use RHS pieces to split vertices in vars_2_vertexOrigID
    nextVertexID <- maxOrigNodeID+1
    for(iDI in seq_along(declInfo)) {
        BUGSdecl <- declInfo[[iDI]]
        rhsVars <- BUGSdecl$rhsVars
        for(iV in seq_along(rhsVars)) {
            rhsVar <- rhsVars[iV]
            nDim <- varInfo[[rhsVar]]$nDim
            if(nDim > 0) {
                splitAns <- splitVertices(vars_2_vertexOrigID[[rhsVar]], BUGSdecl$unrolledIndicesMatrix,
                                          contexts[[BUGSdecl$contextID]]$indexVarExprs, contexts[[BUGSdecl$contextID]]$indexVarNames,
                                          BUGSdecl$symbolicParentNodes[[iV]], BUGSdecl$symbolicParentNodesReplaced[[iV]],
                                          BUGSdecl$parentIndexNamePieces[[iV]], BUGSdecl$replacementNameExprs, nextVertexID, totModelSize)
                
                vars_2_vertexOrigID[[rhsVar]] <- splitAns[[1]]
                nextVertexID <- splitAns[[2]]
            } else {
                if(is.na(vars_2_vertexOrigID[[rhsVar]])) {
                    vars_2_vertexOrigID[[rhsVar]] <- nextVertexID
                    nextVertexID <- nextVertexID + 1
                }
            }
        }
    }

    ## 5. In cases where a vertex was split leaving an incomplete part of an original node, complete the split by creating an inferred vertex for the remaining part
    if(debug) browser()
    for(iV in seq_along(varInfo)) {
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
## Need to handle split original nodes before names
##    numVertices <- vars_2_vertexOrigID[['nextVertexID']]-1
    numVertices <- nextVertexID - 1
  ##  vars_2_vertexOrigID[['vertexID_2_nodeID']] <- vars_2_vertexOrigID[['vertexID_2_nodeID']][ 1:numVertices]
    vertexID_2_nodeID <- c(1:maxOrigNodeID, integer(numVertices - maxOrigNodeID) )
    ## for nodeNamesRHSonly and types, we will collect all RHS nodes now as "RHSonly" and later we'll use the vectorID_2_nodeID to relabel LHSinferred
    for(iV in seq_along(varInfo)) {
        varName <- varInfo[[iV]]$varName
        
        if(varInfo[[iV]]$nDim > 0) {
            newVertexNames <- makeVertexNamesFromIndexArray2(vars_2_vertexOrigID[[varName]], maxOrigNodeID + 1, varName = varName)
            ##   test <- makeVertexNamesFromIndexArray2(vars_2_vertexOrigID[[varName]], maxOrigNodeID + 1, varName = varName)
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
    ## ## there can be skips, so we need a map from origGraphIDs to origVertexIDs
    if(debug) browser()

    ## 7. Re-order vertexIDs to make them contiguous
    contigID_2_origVertexID <- c(1:maxOrigNodeID, sort(allNewVertexIDs))
    numContigVertices <- length(contigID_2_origVertexID)
    origVertexID_2_contigID <- numeric(nextVertexID-1)
    origVertexID_2_contigID[contigID_2_origVertexID] <- 1:length(contigID_2_origVertexID)

     if(debug) browser()
    allVertexNames <- allVertexNames[contigID_2_origVertexID]

    if(length(allVertexNames) > maxOrigNodeID) {
        nodeNamesRHSonly <- allVertexNames[ (maxOrigNodeID + 1) : length(allVertexNames) ]
        types <- c(types, rep('RHSonly', length(allVertexNames) - maxOrigNodeID))
    } else
        nodeNamesRHSonly <- character()

    for(iV in seq_along(varInfo)) {
        temp <- vars_2_vertexOrigID[[varInfo[[iV]]$varName]]
        if(!all(is.na(temp))) vars_2_vertexOrigID[[varInfo[[iV]]$varName]][] <- origVertexID_2_contigID[temp] 
    }


    ## 8. Collect edges from RHS vars to LHS nodes
    if(debug) browser()
    edgesFrom <- numeric(0)
    edgesTo <- numeric(0)
    edgesParentExprID <- numeric(0)
    for(iDI in seq_along(declInfo)) {
        BUGSdecl <- declInfo[[iDI]]
        rhsVars <- BUGSdecl$rhsVars
        for(iV in seq_along(rhsVars)) {
            
            rhsVar <- rhsVars[iV]
            nDim <- varInfo[[rhsVar]]$nDim
            newEdges <- collectEdges(vars_2_vertexOrigID[[rhsVar]], BUGSdecl$unrolledIndicesMatrix, BUGSdecl$origIDs, contexts[[BUGSdecl$contextID]]$indexVarExprs,
                                      BUGSdecl$symbolicParentNodesReplaced[[iV]],
                                      BUGSdecl$parentIndexNamePieces[[iV]], BUGSdecl$replacementNameExprs)
            
            edgesFrom <- c(edgesFrom, newEdges[[1]])
            edgesTo <- c(edgesTo, newEdges[[2]])
            edgesParentExprID <- c(edgesParentExprID, rep(iV, length(newEdges[[1]])))
        }
    }

    ## 9. Collect edges from original nodes to inferred vertices
    if(debug) browser()
    for(iV in seq_along(varInfo)) {
        varName <- varInfo[[iV]]$varName
        newEdges <- collectInferredVertexEdges(vars_2_nodeOrigID[[varName]], vars_2_vertexOrigID[[varName]])
        edgesFrom <- c(edgesFrom, newEdges[[1]])
        edgesTo <- c(edgesTo, newEdges[[2]])
        edgesParentExprID <- c(edgesParentExprID, rep(NA, length(newEdges[[1]])))
        vertexID_2_nodeID[newEdges[[2]]] <- newEdges[[1]]
    }
    
   if(debug) browser()
    
    vertexID_2_nodeID <- vertexID_2_nodeID[1:numContigVertices]
    types[ types == 'RHSonly' & vertexID_2_nodeID != 0] <- 'LHSinferred' ## The types == 'RHSonly' could be obtained more easily, since it will be a single set of FALSES followed by a single set of TRUES
    
    ## 10. Build the graph
    if(debug) browser()
    require(igraph)
    graph <<- graph.empty()
    graph <<- add.vertices(graph, length(allVertexNames), name = allVertexNames)
    allEdges <- as.numeric(t(cbind(edgesFrom, edgesTo)))
    graph <<- add.edges(graph, allEdges)
    if(debug) browser()


    ## 11. Topologically sort and re-index all objects with vertex IDs
    newGraphID_2_oldGraphID <- topological.sort(graph, mode = 'out')
    oldGraphID_2_newGraphID <- sort(newGraphID_2_oldGraphID, index = TRUE)$ix
    graph <<- permute.vertices(graph, oldGraphID_2_newGraphID)  #* ## topological sort
    if(debug) browser()
    ## permute all the maps
    vars_2_nodeID <- new.env()
    vars_2_vertexID <- new.env()
    for(iV in seq_along(varInfo)) {
        temp <- vars_2_nodeOrigID[[varInfo[[iV]]$varName]]
        if(!all(is.na(temp))) vars_2_nodeID[[varInfo[[iV]]$varName]] <- oldGraphID_2_newGraphID[temp] else vars_2_nodeID[[varInfo[[iV]]$varName]] <- temp
        temp <- vars_2_vertexOrigID[[varInfo[[iV]]$varName]]
        if(all(is.na(temp))) cat(paste('Something weird: all vertex IDs NA for variable', varInfo[[iV]]$varName))
        temp[temp==-1] <- NA
        vars_2_vertexID[[varInfo[[iV]]$varName]] <- oldGraphID_2_newGraphID[temp]

        if(varInfo[[iV]]$nDim > 0) dim(vars_2_nodeID[[varInfo[[iV]]$varName]]) <- dim(vars_2_vertexID[[varInfo[[iV]]$varName]]) <- dim(vars_2_vertexOrigID[[varInfo[[iV]]$varName]])        
    }

    ## 11b. re-index the graphIDs in the BUGSdecl objects and record graphID_2_declID
    graphID_2_declID <- numeric(numVertices)
    for(iDI in seq_along(declInfo)) {
        BUGSdecl <- declInfo[[iDI]]
        BUGSdecl$graphIDs <- oldGraphID_2_newGraphID[ BUGSdecl$origIDs ]
        graphID_2_declID[ BUGSdecl$graphIDs ] <- iDI
    }

    ## 12. Set up things needed for maps.
    maps <<- mapsClass$new()
    maps$graphID_2_declID <<- graphID_2_declID
    maps$graphID_2_nodeName <<- allVertexNames[newGraphID_2_oldGraphID]
    maps$types <<- types[newGraphID_2_oldGraphID]
    maps$notStoch <<- maps$types != 'stoch'
    maps$nodeNamesLHSall <<- nodeNamesLHSall
    maps$nodeNamesRHSonly <<- nodeNamesRHSonly
    maps$nodeNames <<- maps$graphID_2_nodeName
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
    maps$vars2LogProbID <<- vars2LogProbID
    maps$logProbIDs_2_LogProbName <<- logProbNames

    maps$edgesFrom <<- oldGraphID_2_newGraphID[edgesFrom]
    maps$edgesTo <<- oldGraphID_2_newGraphID[edgesTo]
    maps$edgesParentExprID <<- edgesParentExprID
    fedgesFrom <- factor(maps$edgesFrom, levels = if(length(maps$edgesFrom) > 0) 1:max(maps$edgesFrom) else numeric(0)) ## setting levels ensures blanks inserted into the splits correctly
    maps$edgesFrom2To <<- split(maps$edgesTo, fedgesFrom)
    maps$edgesFrom2ParentExprID <<- split(maps$edgesParentExprID, fedgesFrom)
    maps$graphIDs <<- 1:length(maps$graphID_2_nodeName)
    NULL
})

modelDefClass$methods(genVarInfo3 = function() {
    ## First set up varInfo's for all LHS variables and collect anyStoch.
    ## That allows determination of when logProb information needs to be collected
    for(iDI in seq_along(declInfo)) {
        BUGSdecl <- declInfo[[iDI]]
        ## LHS:
        lhsVar <- BUGSdecl$targetVarName
        if(!(lhsVar %in% names(varInfo))) {
            nDim <- if(length(BUGSdecl$targetNodeExpr)==1) 0 else length(BUGSdecl$targetNodeExpr)-2
            varInfo[[lhsVar]] <<- varInfoClass$new(varName = lhsVar,
                                                   mins = rep(10000000, nDim),
                                                   maxs = rep(0, nDim),
                                                   nDim = nDim,
                                                   anyStoch = FALSE)
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
            if(!(rhsVar %in% names(varInfo))) {
                nDim <- if(length(BUGSdecl$symbolicParentNodes[[iV]])==1) 0 else length(BUGSdecl$symbolicParentNodes[[iV]])-2
                varInfo[[rhsVar]] <<- varInfoClass$new(varName = rhsVar,
                                                       mins = rep(100000, nDim),
                                                       maxs = rep(0, nDim),
                                                       nDim = nDim,
                                                       anyStoch = FALSE)
            }
            if(varInfo[[rhsVar]]$nDim > 0) {
                for(iDim in 1:varInfo[[rhsVar]]$nDim) {
                    indexNamePieces <- BUGSdecl$parentIndexNamePieces[[iV]][[iDim]] 
                    if(is.list(indexNamePieces)) { ## a list would be made if there is a ':' operator in the index expression
                        indsLow <- if(is.numeric(indexNamePieces[[1]])) indexNamePieces[[1]] else BUGSdecl$replacementsEnv[[ indexNamePieces[[1]] ]]
                        indsHigh <- if(is.numeric(indexNamePieces[[2]])) indexNamePieces[[2]] else BUGSdecl$replacementsEnv[[ indexNamePieces[[2]] ]]
                        varInfo[[rhsVar]]$mins[iDim] <<- min(varInfo[[rhsVar]]$mins[iDim], min(indsLow))
                        varInfo[[rhsVar]]$maxs[iDim] <<- max(varInfo[[rhsVar]]$maxs[iDim], max(indsHigh))
                    } else {
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
        if(length(dimensionsList[[dimVarName]]) != varInfo[[dimVarName]]$nDim)   stop('inconsistent dimensions')
        if(any(dimensionsList[[dimVarName]] < varInfo[[dimVarName]]$maxs))  stop(paste0('dimensions specified are smaller than model specification for variable \'', dimVarName, '\''))
        varInfo[[dimVarName]]$maxs <<- dimensionsList[[dimVarName]]
    }
})

modelDefClass$methods(genVarInfo = function() {
    ## uses nodeInfo to set field: varInfo
    for(iDI in seq_along(declInfo)) {
        BUGSdecl <- declInfo[[iDI]]
        ## LHS:
        lhsVar <- BUGSdecl$targetVarName
        if(!(lhsVar %in% names(varInfo))) {
            nDim <- if(length(BUGSdecl$targetNodeExpr)==1) 0 else length(BUGSdecl$targetNodeExpr)-2
            varInfo[[lhsVar]] <<- varInfoClass$new(varName = lhsVar,
                                                   mins = rep(100000, nDim),
                                                   maxs = rep(0, nDim),
                                                   nDim = nDim,
                                                   anyStoch = FALSE)
        }
        varInfo[[lhsVar]]$anyStoch <<- varInfo[[lhsVar]]$anyStoch | (BUGSdecl$type == 'stoch')
        if(varInfo[[lhsVar]]$nDim > 0) {
            for(iDim in 1:varInfo[[lhsVar]]$nDim) {
                indExprsList <- lapply(BUGSdecl$indexedNodeInfo, function(x) x$targetNodeExpr[[iDim + 2]])
                inds <- unlist(lapply(indExprsList, eval))
                rangeInds <- range(inds)
                varInfo[[lhsVar]]$mins[iDim] <<- min(varInfo[[lhsVar]]$mins[iDim], rangeInds[1])
                varInfo[[lhsVar]]$maxs[iDim] <<- max(varInfo[[lhsVar]]$maxs[iDim], rangeInds[2])
            }
        }
        
        ## RHS:
        rhsDeps <- BUGSdecl$indexedNodeInfo[[1]]$parentNodeExprs
        rhsVars <- unlist(lapply(rhsDeps, function(x) if(length(x) == 1) as.character(x) else as.character(x[[2]])))
        
        for(iV in seq_along(rhsVars)) {
            rhsVar <- rhsVars[iV]
            if(!(rhsVar %in% names(varInfo))) {
                
                nDim <- if(length(rhsDeps[[iV]])==1) 0 else length(rhsDeps[[iV]])-2
                varInfo[[rhsVar]] <<- varInfoClass$new(varName = rhsVar,
                                                       mins = rep(100000, nDim),
                                                       maxs = rep(0, nDim),
                                                       nDim = nDim,
                                                       anyStoch = FALSE)
            }
            if(varInfo[[rhsVar]]$nDim > 0) {
                for(iDim in 1:varInfo[[rhsVar]]$nDim) {
                    indExprsList <- lapply(BUGSdecl$indexedNodeInfo, function(x) x$parentNodeExprs[[iV]][[iDim + 2]])
                    inds <- unlist(lapply(indExprsList, eval))
                    rangeInds <- range(inds)
                    varInfo[[rhsVar]]$mins[iDim] <<- min(varInfo[[rhsVar]]$mins[iDim], rangeInds[1])
                    varInfo[[rhsVar]]$maxs[iDim] <<- max(varInfo[[rhsVar]]$maxs[iDim], rangeInds[2])
                }
            }
        }
    }
    
    ## now use dimensionsList, to check / update varInfo
    for(i in seq_along(dimensionsList)) {
        dimVarName <- names(dimensionsList)[i]
        if(!(dimVarName %in% names(varInfo))) next
        if(length(dimensionsList[[dimVarName]]) != varInfo[[dimVarName]]$nDim)   stop('inconsistent dimensions')
        if(any(dimensionsList[[dimVarName]] < varInfo[[dimVarName]]$maxs))  stop(paste0('dimensions specified are smaller than model specification for variable \'', dimVarName, '\''))
        varInfo[[dimVarName]]$maxs <<- dimensionsList[[dimVarName]]
    }
})
modelDefClass$methods(genLogProbVarInfo = function() {
    ## uses varInfo to set field: logProbVarInfo
    anyStoch = unlist(lapply(varInfo, `[[`, 'anyStoch'))
    logProbVarInfo <<- lapply(varInfo[anyStoch], function(x)
        varInfoClass$new(varName = makeLogProbName(x$varName),
                         mins = rep(100000, x$nDim),
                         maxs = rep(0,      x$nDim),
                         nDim = x$nDim,
                         anyStoch = FALSE))
    names(logProbVarInfo) <<- lapply(logProbVarInfo, `[[`, 'varName')
    
    for(iDI in seq_along(declInfo)) {
        BUGSdecl <- declInfo[[iDI]]
        if(BUGSdecl$type != 'stoch')  next
        
        lhsVar <- BUGSdecl$targetVarName
        lhsLogProbVar <- makeLogProbName(lhsVar)
        
        
        if(varInfo[[lhsVar]]$nDim > 0) {
            if(logProbVarInfo[[lhsLogProbVar]]$nDim != varInfo[[lhsVar]]$nDim)   stop('mismatch in nDim, this should never occur')
            
            for(iDim in 1:varInfo[[lhsVar]]$nDim) {
                inds <- unlist(lapply(BUGSdecl$indexedNodeInfo, function(x) x$logProbIndexValues[iDim]))
                rangeInds <- range(inds)
                logProbVarInfo[[lhsLogProbVar]]$mins[iDim] <<- min(logProbVarInfo[[lhsLogProbVar]]$mins[iDim], rangeInds[1])
                logProbVarInfo[[lhsLogProbVar]]$maxs[iDim] <<- max(logProbVarInfo[[lhsLogProbVar]]$maxs[iDim], rangeInds[2])
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
modelDefClass$methods(buildIgraph = function() {
    
    graph <<- graph.empty()
    nodesLHS <- unique(unlist(lapply(declInfo, function(x) x$allTargetNodeNames())))
    nodesLHSVec <- nodesLHS[grepl(':', nodesLHS)]
    nodesLHSInferred <- nl_vectorizedExpandNodeIndex(nodesLHSVec)
    if(any(nodesLHSInferred %in% nodesLHS))    stop('duplicate declaration of some nodes')
    nodesRHSExprs <- unique(unlist(lapply(declInfo, function(x) x$allParentNodeExprs())))
    nodesRHSAll <- nl_vectorizedExpandNodeIndexExprs(nodesRHSExprs) #*
    
    nodesAll <- unique(c(nodesLHS, nodesLHSInferred, nodesRHSAll))
    graph <<- add.vertices(graph, length(nodesAll), name = nodesAll)
    
    allEdges <- unlist(lapply(declInfo, function(x) x$allEdges())) #**
    graph <<- add.edges(graph, allEdges)
    
    graph <<- permute.vertices(graph, sort(topological.sort(graph, mode = 'out'), index = TRUE)$ix)  #* ## topological sort
})

modelDefClass$methods(buildMaps2 = function() {
	
    ## ditto - and this is repeated work from buildIgraph
    nodesLHS <- unique(unlist(lapply(declInfo, function(x) x$allTargetNodeNames())))
    nodesLHSVec <- nodesLHS[grepl(':', nodesLHS)]
    nodesLHSInferred <- nl_vectorizedExpandNodeIndex(nodesLHSVec)
    if(any(nodesLHSInferred %in% nodesLHS))    stop('duplicate declaration of some nodes')
    nodesLHSAll <- c(nodesLHS, nodesLHSInferred)
    
    nodesLHFInferredLookupList <- list()
    for(nnVec in nodesLHSVec) nodesLHFInferredLookupList[nl_expandNodeIndex(nnVec)] <- nnVec
    
    ##ditto    - and ditto on repeated work
    nodesRHSExprs <- unique(unlist(lapply(declInfo, function(x) x$allParentNodeExprs())))
    nodesRHSAll <- nl_vectorizedExpandNodeIndexExprs(nodesRHSExprs)
    
    nodeNamesAll <- unique(c(nodesLHSAll, nodesRHSAll))
    nodeNamesRHSOnly <- setdiff(nodesRHSAll, nodesLHSAll)
    

    ## this new version cuts out the "middle person" of graphNodesList
    nodeNames <- V(graph)$name ## formerly graphNames
    graphIDs <- seq_along(nodeNames)
    numCase2 <- length(nodesLHSInferred)
    numCase3 <- length(nodeNamesRHSOnly)
    
    case2Names <- unlist(lapply(nodesLHSInferred, function(x) nodeInfo[[ nodesLHFInferredLookupList[[x]] ]]$nodeFunctionName))
    nodeFunctionNamesRaw <- c(unlist(lapply(nodesLHS, function(x) nodeInfo[[x]]$nodeFunctionName)),
                              case2Names,
                              rep('DOES_NOT_EXIST', numCase3))
    names(nodeFunctionNamesRaw) <- c(nodesLHS, nodesLHSInferred, nodeNamesRHSOnly)
    nodeFunctionNamesRaw <- nodeFunctionNamesRaw[ nodeNames ]
    
    originNodeNamesRaw <- nodeNames
    names(originNodeNamesRaw) <- nodeNames
    originNodeNamesRaw[nodesLHSInferred] <- case2Names

    types <- c(unlist(lapply(nodesLHS, function(x) nodeInfo[[x]]$type)),
               rep('LHSinferred', numCase2),
               rep('RHSonly', numCase3))
    names(types) <-  c(nodesLHS, nodesLHSInferred, nodeNamesRHSOnly)
    types <- types[ nodeNames ]
##    nodeFunctionNames  not needed
    
    maps <<- mapsClass$new()
    maps$setup2(nodeNames, graphIDs, nodeFunctionNamesRaw, originNodeNamesRaw, types, graph, varInfo, nodeInfo)

})

modelDefClass$methods(newModel = function(data = list(), inits = list(), where = globalenv(), modelName = character()) {
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
    model$buildNodeFunctions(where = where)
    model$buildNodesList() ## This step makes RStudio choke, we think from circular reference classes -- fixed, by not displaying Global Environment in RStudio
    model$setData(data)
    # prevent overwriting of data values by inits
    for(varName in intersect(names(inits), model$getVarNames())) {
        dataVars <- model$isData(varName)
        if(sum(dataVars) && !identical(data[[varName]][dataVars],
                                      inits[[varName]][dataVars])) {
            inits[[varName]][dataVars] <- data[[varName]][dataVars]
            nonNAinits <- !is.na(inits[[varName]][dataVars])
            # only warn if user passed conflicting actual values
            if(!identical(data[[varName]][dataVars][nonNAinits],
                                      inits[[varName]][dataVars][nonNAinits]))
                warning("newModel: Conflict between 'data' and 'inits' for ", varName, "; using values from 'data'.\n")
        }
    }
    nonVarIndices <- !names(inits) %in% model$getVarNames()
    if(sum(nonVarIndices))
        warning("newModel: ", paste(names(inits)[nonVarIndices], collapse = ','),
                " ", ifelse(sum(nonVarIndices) > 1, "are", "is"), " not ", ifelse(sum(nonVarIndices) > 1, "variables", "a variable"), " in the model; initial value ignored.")
    model$setInits(inits[!nonVarIndices])

	# Below is the code that checks if an index is missing    
#    allVarNames <- model$getVarNames()
#    for(var in allVarNames){
#    	varExpandedNodeNameLength = length(model$expandNodeNames(var, returnScalarComponents = TRUE))
#       	varValuesLength = length(model[[var]])
#        if(varExpandedNodeNameLength != varValuesLength){
#        	warningText <- paste('missing node detected, i.e. something like x[2] declared but not x[1]')
#        	warningText <- paste(warningText, '\n variable name = ', var)
#         	warningText <- paste(warningText, '\n Can be remedied by adding dummy nodes, i.e. x[1] <- 0')
#           	warning(warningText)
#            }
#         }
    ## fixing the problem with RStudio hanging: over-writing the str() method for this model class
    ## added by DT, April 2015
    thisClassName <- as.character(class(model))
    eval(substitute(METHOD <- function(object, ...) str(NULL), list(METHOD = as.name(paste0('str.', thisClassName)))), envir = globalenv())
    return(model)
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
        cat(paste0(deparse(BUGSdecl$code, width.cutoff=500L)))
        cat(paste0(rep('   }}}', length(contexts[[BUGSdecl$contextID]]$singleContexts)), collapse=''))
        cat('\n')
    }
})



modelDefClass$methods(nodeName2GraphIDs = function(nodeName, nodeFunctionID = TRUE){
	if(length(nodeName) == 0)
		return(NULL)	
	if(nodeFunctionID)
		output <- unique(unlist(sapply(nodeName, parseEvalNumeric, env = maps$vars2GraphID_functions, USE.NAMES = FALSE)))
	else
		output <- unlist(sapply(nodeName, parseEvalNumeric, env = maps$vars2GraphID_values, USE.NAMES = FALSE))	
	return(output[!is.na(output)])
})

## next two functions work for properly formed nodeNames.
modelDefClass$methods(nodeName2LogProbName = function(nodeName){ ## used in 2 places: MCMC_build and cppInterfaces_models
	if(length(nodeName) == 0)
		return(NULL)
	output <- unique(unlist(sapply(nodeName, parseEvalCharacter, env = maps$vars2LogProbName, USE.NAMES = FALSE)))
	return(output[!is.na(output)])
})

modelDefClass$methods(nodeName2LogProbID = function(nodeName){ ## used only in cppInterfaces_models
	if(length(nodeName) == 0)
		return(NULL)
	output <- unique(unlist(sapply(nodeName, parseEvalNumeric, env = maps$vars2LogProbID, USE.NAMES = FALSE) ) ) 
	return(output[!is.na(output)])
})

parseEvalNumeric <- function(x, env){
	ans <- eval(parse(text = x, keep.source = FALSE)[[1]], envir = env)
	as.numeric(ans)
}
parseEvalCharacter <- function(x, env){
	ans <- eval(parse(text = x, keep.source = FALSE)[[1]], envir = env)
	as.character(ans)
}



#The class we use keep track of graphIDs
nodeVector <- setRefClass(Class = "nodeVector",
						fields = 
							list(origNodeNames 				=	'ANY',
								expandedNodeFunctionNames 	=	'ANY',
								expandedNodeValuesNames 	=	'ANY',
								origGraphIDs_values 		=	'ANY',
								origGraphIDs_functions		=	'ANY',
								sortedGraphIDs_values 		=	'ANY',
								sortedGraphIDs_functions	=	'ANY',
								sortOrder_values 			= 	'ANY',
								sortOrder_functions 		= 	'ANY',
								model 						=	'ANY'),
						methods = 	#These methods need to check if objects are initiated. If not, it needs
									#to initiate. If so, it just returns object
							list(getOrigNames = function(){ 
									if(!inherits(origNodeNames, 'uninitializedField') ) 
										origNodeNames 
									else stop('origNodeNames never initialized!')
									},
								getExpandedFunctionNames = function(){ 
									if(!inherits(expandedNodeFunctionNames, 'uninitializedField') ) 
										expandedNodeFunctionNames 
									else {
										expandedNodeFunctionNames <<- unique(model$modelDef$maps$graphID_2_nodeFunctionName[origGraphIDs_functions])
										expandedNodeFunctionNames <<- expandedNodeFunctionNames[expandedNodeFunctionNames != 'DOES_NOT_EXIST']
										expandedNodeFunctionNames
										}
									},
								getExpandedValuesNames = function(){ 
									if(!inherits(expandedNodeValuesNames, 'uninitializedField') ) 
										expandedNodeValuesNames 
									else {
										expandedNodeValuesNames <<- unique(model$modelDef$maps$graphID_2_nodeName[origGraphIDs_values])
										expandedNodeFunctionNames <<- expandedNodeFunctionNames[expandedNodeFunctionNames != 'DOES_NOT_EXIST']
										expandedNodeValuesNames
										}
									},
								getOrigIDs_values = function() { 
									if(!inherits(origGraphIDs_values, 'uninitializedField') ) 
										origGraphIDs_values 
									else stop('origGraphIDs_values never initialized! (this indicates a bug in the nimble source code: should not be able to build a nodeVector without initializing OrigGraphIDs_values)')
									},
								getOrigIDs_functions = function() { 
									if(!inherits(origGraphIDs_functions, 'uninitializedField') ) 
										origGraphIDs_functions 
									else stop('origGraphIDs_functions never initialized! (this indicates a bug in the nimble source code: should not be able to build a nodeVector without initializing OrigGraphIDs_values)')
									},
								getSortedIDs_values = function() { 
									if(!inherits(sortedGraphIDs_values, 'uninitializedField') ) 
										sortedGraphIDs_values 
									else{ 
										if(inherits(sortOrder_values, 'uninitializedField') ) 
											sortOrder_values <<- order(origGraphIDs_values)
 										sortedGraphIDs_values <<- origGraphIDs_values[sortOrder_values]
										sortedGraphIDs_values
									}
								},
								getSortedIDs_functions = function() { 
									if(!inherits(sortedGraphIDs_functions, 'uninitializedField') ) 
										sortedGraphIDs_functions 
									else{
										if(inherits(sortOrder_functions, 'uninitializedField') ) 
											sortOrder_functions <<- order(origGraphIDs_functions)
										sortedGraphIDs_functions <<- origGraphIDs_functions[sortOrder_functions]
										sortedGraphIDs_functions
									}
								},
								getSortOrder_values = function() { 
									if(!inherits(sortOrder_values, 'uninitializedField') ) 
										sortOrder_values 
									else{
										sortOrder_values <<- order(origGraphIDs_values)
										sortOrder_values	
									}
								},	
								getSortOrder_functions = function() { 
									if(!inherits(sortOrder_functions, 'uninitializedField') ) 
										sortOrder_functions 
									else {
										sortOrder_functions <<- order(origGraphIDs_functions)
										sortOrder_functions
									}
								},	
								initialize = function(...){
									callSuper(...)
									if(inherits(model, 'uninitializedField'))	
										stop('model missing when creating nodeVector')
									if(inherits(origGraphIDs_values, 'uninitializedField') )	{
										if(!inherits(origNodeNames, 'uninitializedField') ){
											origGraphIDs_values <<- model$modelDef$nodeName2GraphIDs(origNodeNames, FALSE)
											origGraphIDs_functions <<- model$modelDef$nodeName2GraphIDs(origNodeNames, TRUE)
										}
										else{
											origNodeNames <<- model$modelDef$maps$graphID_2_nodeName[origGraphIDs_functions]
											origGraphIDs_values <<- model$modelDef$nodeName2GraphIDs(origNodeNames, FALSE)
										}
									}										
									else if(inherits(origGraphIDs_functions, 'uninitializedField')	){
											origNodeNames <<- model$modelDef$maps$graphID_2_nodeName[origGraphIDs_values]
											origGraphIDs_functions <<- model$modelDef$nodeName2GraphIDs(origNodeNames, TRUE)
										}

								}))

