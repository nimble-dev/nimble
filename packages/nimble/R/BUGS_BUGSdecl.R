## The BUGSdeclClass contains the pulled-apart content of a BUGS declaration line

## nimbleOrRfunctionNames is used to determine what can be evaluated in R if every argument is known OR in C++ (nimble) if arguments are other nodes
nimbleOrRfunctionNames <- c('+','-','/','*','(','exp','log','pow','^','%%','%*%',
                            'equals','inprod','nimbleEquals',
                            'sqrt', 'logit', 'expit', 'ilogit', 'probit', 'iprobit', 'phi', 'cloglog', 'icloglog', 'chol', 'step', 'nimbleStep', 'inverse',
                            'sin','cos','tan','asin','acos','atan','cosh','sinh','tanh', 'asinh', 'acosh', 'atanh',
                            'cube', 'abs', 'lgamma', 'loggam', 'log1p', 'lfactorial', ##'factorial', 'gamma',
                            'ceiling', 'floor', 'round', 'trunc',
                            'mean','sum','max','min','prod',
                            'asRow', 'asCol',
                            '>', '<', '>=', '<=', '==', '!=', '&', '|')

#' BUGSdeclClass contains the information extracted from one BUGS declaration
BUGSdeclClass <- setRefClass('BUGSdeclClass',
                             
                             fields = list(
                                 ###### the following are set in setup(), and never change.
                                 contextID = 'ANY',
                                 sourceLineNumber = 'ANY',
                                 code = 'ANY',        ## original BUGS code line
                                 type = 'ANY',
                                 targetExpr = 'ANY',   ## LHS of code
                                 valueExpr = 'ANY',    ## RHS of code
                                 transExpr = 'ANY',
                                 indexExpr = 'ANY',
                                 targetVarExpr = 'ANY',
                                 targetNodeExpr = 'ANY',
                                 targetVarName = 'ANY',
                                 targetNodeName = 'ANY',
                                 
                                 ## truncation information
                                 truncation = 'ANY',
                                 
                                 ## set in setIndexVariableExprs(), and never changes.
                                 indexVariableExprs = 'ANY',
                                 
                                 ## set in genSymbolicParentNodes(), and never changes.
                                 symbolicParentNodes = 'ANY',
                                 
                                 ##### the following are set in genReplacementsAndCodeReplaced(), and never change, with one exception
                                 replacements = 'ANY',
                                 codeReplaced = 'ANY',     ## *** MODIFIED **** in genAltParamsModifyCodeReplaced() to remove the .params
                                 replacementNameExprs = 'ANY',
                                 logProbNodeExpr = 'ANY',
                                 
                                 ##### the following is set in genAltParamsModifyCodeReplaced(), and never changes
                                 altParamExprs = 'ANY',   ## contains the "smart" expression from each altParam: from original parameterization, whenever possible
                                 
                                 ## the following are set in modelDefClass$genNodeInfo(), and never change:
                                 indexedNodeInfo = 'ANY',
                                 indexedNodeNames = 'ANY',

                                 ## new for version 3
                                 targetExprReplaced = 'ANY',
                                 valueExprReplaced = 'ANY',
                                 symbolicParentNodesReplaced = 'ANY',
                                 rhsVars = 'ANY',
                                 targetIndexNamePieces = 'ANY',
                                 parentIndexNamePieces = 'ANY',
                                 replacementsEnv = 'ANY',
                                 nodeFunctionNames = 'ANY',

                                 outputSize = 'ANY', ## should match nrow(unrolledIndicesMatrix)
                                 origIDs = 'ANY',
                                 graphIDs = 'ANY',
                                 unrolledIndicesMatrix = 'ANY',
                                 numUnrolledNodes = 'ANY' ## differs from outputSize ONLY for a no-context singleton, a ~ dnorm(b, c), so numUnrolledNodes is 1 but outputSize is 0
                             ),   
                             
                             methods = list(
                                 setup                          = function() {},
                                 setIndexVariableExprs          = function() {},
                                 genSymbolicParentNodes         = function() {},
                                 genReplacementsAndCodeReplaced = function() {},
                                 genAltParamsModifyCodeReplaced = function() {},

                                 genReplacedTargetValueAndParentInfo = function() {},

                                 allParentVarNames = function() {
                                     unlist(lapply(symbolicParentNodes, function(x) if(length(x) == 1) as.character(x) else if(x[[1]] == '[') as.character(x[[2]]) else stop('Error in allParentVarNames')))
                                 },
                                 
                                 allTargetNodeNames = function() {
                                     unlist(lapply(indexedNodeInfo, `[[`, 'targetNodeName'))
                                 },
                                 
                                 allParentNodeNames = function() {
                                     unique(unlist(lapply(indexedNodeInfo, `[[`, 'parentNodeNames')))
                                 },

                                 allParentNodeExprs = function() {
                                     unique(unlist(lapply(indexedNodeInfo, `[[`, 'parentNodeExprs')))
                                 },
                                 allEdges = function() {
                                     edgesIn  <- unlist(lapply(indexedNodeInfo, 
                                                               function(x) {
                                                                   expandedNodeIndices <- nl_vectorizedExpandNodeIndexExprs(x$parentNodeExprs)
                                                                   if(length(expandedNodeIndices) > 0) 
                                                                       strsplit(paste(expandedNodeIndices, 
                                                                                      x$targetNodeName, 
                                                                                      sep = '&', collapse = '&'), '&') else NULL}))
                                     edgesOut <- unlist(lapply(indexedNodeInfo, 
                                                               function(x) if(is.vectorized(x$targetNodeName)) 
                                                               strsplit(paste(x$targetNodeName, 
                                                                                  nl_expandNodeIndexExpr(x$targetNodeExpr), 
                                                                                  sep = '&', collapse = '&'), '&') else NULL))
                                     return(c(edgesIn, edgesOut))
                                 },
                                 getDistribution = function() {
                                     if(type != 'stoch')  stop('getting distribution of non-stochastic node')
                                     return(as.character(valueExprReplaced[[1]]))
                                 }
                             )
)


BUGSdeclClass$methods(setup = function(code, contextID, sourceLineNum, truncation = NULL) {
    ## master entry function.
    ## uses 'contextID' to set the field: contextID.
    ## uses 'code' argument, to set the fields:
    ## code
    ## targetExpr, valueExpr
    ## targetVarExpr, targetNodeExpr
    ## targetVarName, targetNodeName
    
    contextID <<- contextID
    sourceLineNumber <<- sourceLineNum
    code <<- code
    truncation <<- truncation
    
    if(code[[1]] == '~') {
        type <<- 'stoch'
        if(!is.call(code[[3]]) || (!any(code[[3]][[1]] == getDistributionsInfo('namesVector')) && code[[3]][[1]] != "T" && code[[3]][[1]] != "I"))
            stop(paste0('Improper syntax for stochastic declaration: ', deparse(code)))
    } else if(code[[1]] == '<-') {
        type <<- 'determ'
        if( is.call(code[[3]]) &&  any(code[[3]][[1]] == getDistributionsInfo('namesVector')))
            stop(paste0('Improper syntax for determistic declaration: ', deparse(code)))
    } else {
        stop(paste0('Improper syntax for declaration: ', deparse(code)))
    }
    
    targetExpr <<- code[[2]]
    valueExpr <<- code[[3]]
    
    transExpr <<- NULL
    indexExpr <<- NULL
    
    if(length(targetExpr) > 1) {
        ## There is a tranformation and/or a subscript
        if(targetExpr[[1]] == '[') {
            ## It is a subscript only
            indexExpr <<- as.list(targetExpr[-c(1,2)]) 
            targetVarExpr <<- targetExpr[[2]]
            targetNodeExpr <<- targetExpr
        } else {
            ## There is a transformation, possibly with a subscript
            transExpr <<- targetExpr[[1]]
            targetNodeExpr <<- targetExpr[[2]]
            if(length(targetNodeExpr)>1) {
                ## There are subscripts inside the transformation
                if(targetNodeExpr[[1]] != '[') {
                    print(paste("Invalid subscripting for", deparse(targetExpr)))
                }
                indexExpr <<- as.list(targetNodeExpr[-c(1,2)])
                targetVarExpr <<- targetNodeExpr[[2]]
            } else {
                targetVarExpr <<- targetNodeExpr
            }
        }
    } else {
        ## no tranformation or subscript
        targetVarExpr <<- targetExpr
        targetNodeExpr <<- targetVarExpr
    }
    
    targetVarName <<- deparse(targetVarExpr)
    targetNodeName <<- deparse(targetNodeExpr)
})


BUGSdeclClass$methods(setIndexVariableExprs = function(exprs) {
    indexVariableExprs <<- exprs
})


BUGSdeclClass$methods(genSymbolicParentNodes = function(constantsNamesList, context, nimFunNames) {
    ## sets the field symbolicparentNodes
    symbolicParentNodes <<- unique(getSymbolicParentNodes(valueExpr, constantsNamesList, context$indexVarExprs, nimFunNames)) 
})

## move this to a util file when everything is working.  It is convenient here for now
makeIndexNamePieces <- function(indexCode) {
    if(length(indexCode) == 1) return(if(is.numeric(indexCode)) indexCode else as.character(indexCode))
    ## Diagnostic for messed up indexing here
    if(as.character(indexCode[[1]] != ':')) stop(paste0("Error processing model: something is wrong with the index ", deparse(indexCode),". Note that any variables in index expressions must be provided as constants.  NIMBLE does not yet allow indices that are model nodes."), call. = FALSE)
    p1 <- indexCode[[2]]
    p2 <- indexCode[[3]]
    list( if(is.numeric(p1)) p1 else as.character(p1),
      if(is.numeric(p2)) p2 else as.character(p2))
    ## e.g. makeIndexNamePieces(quote(i))
    ##      makeIndexNamePieces(quote(i:100))
}

BUGSdeclClass$methods(genReplacedTargetValueAndParentInfo = function(constantsNamesList, context, nimFunNames) { ## assuming codeReplaced is there
    ## generate hasBracket info
    
    targetExprReplaced <<- codeReplaced[[2]] ## shouldn't have any link functions at this point
    valueExprReplaced <<- codeReplaced[[3]]

    symbolicParentNodesReplaced <<- unique(getSymbolicParentNodes(valueExprReplaced, constantsNamesList, c(context$indexVarExprs, replacementNameExprs), nimFunNames))
    rhsVars <<- unlist(lapply(symbolicParentNodesReplaced,  function(x) if(length(x) == 1) as.character(x) else as.character(x[[2]])))

    ## note that makeIndexNamePieces is designed only for indices that are a single name or number or a `:` operator with single name or number for each argument
    ## This relies on the fact that any expression will have been lifted by this point and what it has been replaced with is simply a name
    ## This means makeIndexNamePieces can include a diagnostic
    targetIndexNamePieces <<- try(if(length(targetExprReplaced) > 1) lapply(targetExprReplaced[-c(1,2)], makeIndexNamePieces) else NULL)
    if(inherits(targetIndexNamePieces, 'try-error')) stop(paste('Error occurred defining ', deparse(targetExprReplaced)), call. = FALSE)
    parentIndexNamePieces <<- lapply(symbolicParentNodesReplaced, function(x) if(length(x) > 1) lapply(x[-c(1,2)], makeIndexNamePieces) else NULL)
    NULL
})
                      
BUGSdeclClass$methods(genReplacementsAndCodeReplaced = function(constantsNamesList, context, nimFunNames) {
    replacementsAndCode <- genReplacementsAndCodeRecurse(code, c(constantsNamesList,context$indexVarExprs), nimFunNames)
    replacements <<- replacementsAndCode$replacements
    codeReplaced <<- replacementsAndCode$codeReplaced
    
    if(type == 'determ') { logProbNodeExpr <<- NULL }
    if(type == 'stoch') {
        logProbNodeExprAndReplacements <- genLogProbNodeExprAndReplacements(code, codeReplaced, context$indexVarExprs)
        logProbNodeExpr <<- logProbNodeExprAndReplacements$logProbNodeExpr
        replacements <<- c(replacements, logProbNodeExprAndReplacements$replacements)
    }
    
    replacementNameExprs <<- lapply(as.list(names(replacements)), as.name)
    names(replacementNameExprs) <<- names(replacements)
})

## only affects stochastic nodes
## removes any params in codeReplaced which begin with '.'
## generates the altParamExprs list, which contains the expression for each alternate parameter,
## which is taken from the .param expression (no longer taken from getDistribution(distName)$altParams, ever)
BUGSdeclClass$methods(genAltParamsModifyCodeReplaced = function() {
    
    altParamExprs <<- list()
    
    if(type == 'stoch') {
        RHSreplaced <- codeReplaced[[3]]
        paramNamesAll <- names(RHSreplaced)
        paramNamesDotLogicalVector <- grepl('^\\.', paramNamesAll)
        RHSreplacedWithoutDotParams <- RHSreplaced[!paramNamesDotLogicalVector]    ## removes all parameters whose name begins with '.' from distribution
        codeReplaced[[3]] <<- RHSreplacedWithoutDotParams
        
        altParamExprs <<- if(any(paramNamesDotLogicalVector)) as.list(RHSreplaced[paramNamesDotLogicalVector]) else list()
        names(altParamExprs) <<- gsub('^\\.', '', names(altParamExprs))    ## removes the '.' from each name
#         dotParamNames <- names(dotParamExprs)
#         distRuleAltParamExprs <- getDistribution(as.character(RHSreplaced[[1]]))$altParams
#         for(altParam in names(distRuleAltParamExprs)) {
#             if(altParam %in% dotParamNames) {
#                 altParamExprs[[altParam]] <<- dotParamExprs[[altParam]]
#             } else {
#                 defaultParamExpr <- getDistributions(as.character(RHSreplaced[[1]]))$altParams[[altParam]]
#                 subParamExpr <- eval(substitute(substitute(EXPR, as.list(RHSreplaced)[-1]), list(EXPR=defaultParamExpr)))
#                 altParamExprs[[altParam]] <<- subParamExpr
#             }
#         }
    }
})

getSymbolicParentNodes <- function(code, constNames = list(), indexNames = list(), nimbleFunctionNames = list(), addDistNames = FALSE) {
    ## replaceConstants looks to see if name of a function exists in R
    ## getSymbolicVariables requires a list of nimbleFunctionNames.
    ## The latter could take the former approach
    if(addDistNames) nimbleFunctionNames <- c(nimbleFunctionNames, getDistributionsInfo('namesExprList'))
    ans <- getSymbolicParentNodesRecurse(code, constNames, indexNames, nimbleFunctionNames)
    return(ans$code)
}
getSymbolicParentNodesRecurse <- function(code, constNames = list(), indexNames = list(), nimbleFunctionNames = list()) {
    ## Takes as input some code and returns the variables in it
    ## Expects one line of code, no '{'s
    ## However, indexNames and constNames are not identified as separate variables
    ## e.g. x[i] is returned as 'x' and 'i' or as 'x[i]' if i is an indexName
    ## indexNames and constNames can be substituted at compile time, such as a block index variable
    ## every function EXCEPT those in nimbleFunctionNames can be evaluated at compile time
    ## constNames, indexNames and nimbleFunctionNames should be lists of names
    if(is.numeric(code)) {
        return(list(code = NULL, replaceable = TRUE, hasIndex = FALSE))
    }
    cLength <- length(code)
    if(cLength == 1) {
        if(is.name(code)) {
            if(any(code == indexNames)) {
                return(list(code = NULL, replaceable = TRUE, hasIndex = TRUE))
            }
            if(any(code == constNames)) {
                return(list(code = NULL, replaceable = TRUE, hasIndex = FALSE))
            }
            ## just something regular: not constant or index
            return(list(code = list(code), replaceable = FALSE, hasIndex = FALSE))
        }
    }

    if(is.call(code)) {
        if(code[[1]] == '[') {
            contents <- lapply(code[-c(1,2)], function(x) getSymbolicParentNodesRecurse(x, constNames, indexNames, nimbleFunctionNames))
            contentsCode <- unlist(lapply(contents, function(x) x$code), recursive = FALSE)
            contentsHasIndex <- unlist(lapply(contents, function(x) x$hasIndex))
            contentsReplaceable <- unlist(lapply(contents, function(x) x$replaceable))
            variable <- getSymbolicParentNodesRecurse(code[[2]], constNames, indexNames, nimbleFunctionNames)
            
            if(variable$hasIndex) stop('Error: Variable', deparse(code[[2]]), 'on outside of [ contains a BUGS code index.')
            if(variable$replaceable) {
                return(list(code = contentsCode,
                            replaceable = all(contentsReplaceable),
                            hasIndex = any(contentsHasIndex)))
            } else {
                if(all(contentsReplaceable)) {
                    return(list(code = c(contentsCode, list(code)),
                                replaceable = FALSE,
                                hasIndex = any(contentsHasIndex)))
                } else {
                    return(list(code = c(contentsCode, list(code[[2]])),
                                replaceable = FALSE,
                                hasIndex = any(contentsHasIndex)))
                }
            }
        } else {
            if(cLength > 1) {
                contents <- lapply(code[-1], function(x) getSymbolicParentNodesRecurse(x, constNames, indexNames, nimbleFunctionNames))
                contentsCode <- unlist(lapply(contents, function(x) x$code), recursive = FALSE)
                contentsHasIndex <- unlist(lapply(contents, function(x) x$hasIndex))
                contentsReplaceable <- unlist(lapply(contents, function(x) x$replaceable))
                allContentsReplaceable <- all(contentsReplaceable)
            } else {
                contentsCode <- NULL
                contentsHasIndex <- FALSE
                allContentsReplaceable <- TRUE
            }
            isRfunction <- !any(code[[1]] == nimbleFunctionNames)
            isRonly <- isRfunction &
                !checkNimbleOrRfunctionNames(deparse(code[[1]]))
#                !any(deparse(code[[1]]) == nimbleOrRfunctionNames)
            if(isRonly & !allContentsReplaceable) stop(paste('Error, R function', deparse(code[[1]]),' has non-replaceable node values as arguments.  Must be a nimble function.'))
            
            return(list(code = contentsCode,
                        replaceable = allContentsReplaceable & isRfunction,
                        hasIndex = any(contentsHasIndex)))
        }
    }
    stop(paste('Something went wrong in getSymbolicVariablesRecurse with', deparse(code)))
}

checkNimbleOrRfunctionNames <- function(functionName) {
    if(any(functionName == nimbleOrRfunctionNames)) return(TRUE)
    if(exists(functionName) && is.rcf(get(functionName))) return(TRUE)  # FIXME: deal with finding by R's scoping rules here and in genCpp_sizeProcessing (currently line 139)
    return(FALSE)
}


## The replaceVariableLHS arg avoids replacement of x[i] if x[i] is on the LHS of <- or ~
genReplacementsAndCodeRecurse <- function(code, constAndIndexNames, nimbleFunctionNames, replaceVariableLHS = TRUE, debug = FALSE) {
    if(debug) browser()
    if(is.numeric(code))  return(list(codeReplaced = code, replacements = list(), replaceable = TRUE))
    cLength <- length(code)
    if(cLength == 1) {
        if(is.name(code)) {
            if(any(code == constAndIndexNames) & replaceVariableLHS) return(replaceAllCodeSuccessfully(code))
            else  return(list(codeReplaced = code, replacements = list(), replaceable = FALSE))
        }
    }
    if(is.call(code)) {
        if(code[[1]] == '[') {
            contents <- lapply(code[-c(1,2)], function(x) genReplacementsAndCodeRecurse(x, constAndIndexNames, nimbleFunctionNames, debug = debug))
            contentsCodeReplaced <- lapply(contents, function(x) x$codeReplaced)
            contentsReplacements <- lapply(contents, function(x) x$replacements)
            contentsReplaceable  <- unlist(lapply(contents, function(x) x$replaceable))
            if(replaceVariableLHS) {
                variable <- genReplacementsAndCodeRecurse(code[[2]], constAndIndexNames, nimbleFunctionNames, debug = debug)
                if(variable$replaceable && all(contentsReplaceable))  return(replaceAllCodeSuccessfully(code))
            }
            return(replaceWhatWeCan(code, contentsCodeReplaced, contentsReplacements, contentsReplaceable, startingAt=3))
        }
        assignment <- any(code[[1]] == c('<-', '~'))
        if(cLength > 1) {
            if(assignment) {
                ## In an assignment, prevent the outermost variable on the LHS from being replaced 
                contents <- c(list(genReplacementsAndCodeRecurse(code[[2]], constAndIndexNames, nimbleFunctionNames, replaceVariableLHS = FALSE, debug)),
                              lapply(code[-c(1,2)], function(x) genReplacementsAndCodeRecurse(x, constAndIndexNames, nimbleFunctionNames, debug = debug)))
            } else {
                contents <- lapply(code[-1], function(x) genReplacementsAndCodeRecurse(x, constAndIndexNames, nimbleFunctionNames, debug = debug))
            }
            contentsCodeReplaced <- lapply(contents, function(x) x$codeReplaced)
            contentsReplacements <- lapply(contents, function(x) x$replacements)
            contentsReplaceable  <- unlist(lapply(contents, function(x) x$replaceable))
            allContentsReplaceable <- all(contentsReplaceable)
        } else {
            contentsCodeReplaced <- list()
            contentsReplacements <- list()
            contentsReplaceable  <- list()
            allContentsReplaceable <- TRUE
        }
        if(code[[1]] == ':')   return(replaceWhatWeCan(code, contentsCodeReplaced, contentsReplacements, contentsReplaceable, startingAt=2, replaceable=allContentsReplaceable))
        if(assignment)         return(replaceWhatWeCan(code, contentsCodeReplaced, contentsReplacements, contentsReplaceable, startingAt=2))
        isRfunction <- !any(code[[1]] == nimbleFunctionNames)
#        isRonly <- isRfunction & !any(deparse(code[[1]]) == nimbleOrRfunctionNames)
        isRonly <- isRfunction & !checkNimbleOrRfunctionNames(deparse(code[[1]]))
        if(isRonly & !allContentsReplaceable) stop(paste0('Error, R function \"', deparse(code[[1]]),'\" has non-replaceable node values as arguments.  Must be a nimble function.'))
        if(isRfunction & allContentsReplaceable)   return(replaceAllCodeSuccessfully(code))
        return(replaceWhatWeCan(code, contentsCodeReplaced, contentsReplacements, contentsReplaceable, startingAt=2))
    }
    stop(paste('Something went wrong in genReplacementsAndCodeRecurse with', deparse(code)))
}
replaceAllCodeSuccessfully <- function(code) {
    deparsedCode <- Rname2CppName(code, colonsOK = TRUE) ## nameMashup
    replacements <- list()
    replacements[[deparsedCode]] <- code
    return(list(codeReplaced = as.name(deparsedCode), replacements = replacements, replaceable = TRUE))
}
replaceWhatWeCan <- function(code, contentsCodeReplaced, contentsReplacements, contentsReplaceable, startingAt, replaceable=FALSE) {
    replacements <- list()
    codeReplaced <- code
    if(length(code) >= startingAt) for(i in seq_along(contentsReplaceable)) {
        replacements <- c(replacements, contentsReplacements[[i]])
        codeReplaced[[i+startingAt-1]] <- contentsCodeReplaced[[i]]
    }
    replacements <- replacements[unique(names(replacements))]
    list(codeReplaced = codeReplaced, replacements = replacements, replaceable = replaceable)
}
genLogProbNodeExprAndReplacements <- function(code, codeReplaced, indexVarExprs) {
    logProbNodeExpr <- codeReplaced[[2]]   ## initially, we'll use the replaced version
    replacements <- list()
    
    if(length(logProbNodeExpr) == 1) {
        ## no indexing present
        logProbNodeExpr <- as.name(makeLogProbName(logProbNodeExpr))
        
    } else {
        ## indexing on the LHS node
        if(logProbNodeExpr[[1]] != '[')    stop('something wrong')
        logProbNodeExpr[[2]] <- as.name(makeLogProbName(logProbNodeExpr[[2]]))
        
        origLHS <- code[[2]]
        for(i in seq_along(origLHS)[-c(1,2)]) {
            origIndex <- origLHS[[i]]
            if(is.vectorized(origIndex)) {
                if(any(indexVarExprs %in% all.vars(origIndex))) {
                    ## the vectorized index includes a loop-indexing variable; we will create a replacement, for a memberData, for each nodeFunction
                    replacementExpr <- substitute(min(EXPR), list(EXPR=origIndex))
                    replacementName <- Rname2CppName(replacementExpr, colonsOK = TRUE)##nameMashup
                    logProbNodeExpr[[i]] <- as.name(replacementName)
                    replacements[[replacementName]] <- replacementExpr
                } else {
                    ## no loop-indexing variables present in the vectorized index.  this index should be constant for all instances of this nodeFunction
                    logProbIndexValue <- as.numeric(min(eval(origIndex)))   ## if this eval() causes an error... then is some case I haven't thought of.
                    logProbNodeExpr[[i]] <- logProbIndexValue
                }
            }
        }
    }
    
    list(logProbNodeExpr = logProbNodeExpr, replacements = replacements)
}
